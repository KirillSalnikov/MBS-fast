// Two stages, two existing queues, no third polygon buffer. After clipping the
// next queue temporarily contains clipped apertures with parent optical state
// and a target facet. Optics transforms each slot in place. Rejected/depth-limit
// reflections become inactive slots and are skipped on the next level.
struct ClippedAperture { Vec vertices[MaxVertices];int count; };
// Both children have the same aperture. Optics borrows its immutable vertices
// until the outgoing beam is saved and the reflection replaces the staged item.
struct OpticalView {
    const Vec* vertices;
    Vec direction,basis;
    Cx j[4];
    double path,front;
    int count,facet,depth;
    unsigned char track[20];
};
__device__ __forceinline__ void storeBeam(Beam& dest,const OpticalView& source) {
    dest.count=source.count;dest.facet=source.facet;dest.depth=source.depth;
    dest.direction=source.direction;dest.basis=source.basis;dest.path=source.path;dest.front=source.front;
    for(int j=0;j<4;++j)dest.j[j]=source.j[j];
    if(dest.vertices!=source.vertices)
        for(int v=0;v<source.count;++v)dest.vertices[v]=source.vertices[v];
    for(int d=0;d<=source.depth;++d)dest.track[d]=source.track[d];
}
__device__ __forceinline__ void addNonzeroStatus(Status* states,int orientation,const Status& s) {
    if(s.failure)atomicExch(&states[orientation].failure,s.failure);
    if(s.count)atomicAdd(&states[orientation].count,s.count);
    if(s.evaluated)atomicAdd(&states[orientation].evaluated,s.evaluated);
    if(s.rejected)atomicAdd(&states[orientation].rejected,s.rejected);
    if(s.rejectedJones)atomicAdd(&states[orientation].rejectedJones,s.rejectedJones);
    if(s.rejectedArea)atomicAdd(&states[orientation].rejectedArea,s.rejectedArea);
    if(s.rejectedImportance)atomicAdd(&states[orientation].rejectedImportance,s.rejectedImportance);
}
template<int lanes> __global__ void clipLevel(const Input* inputs,const QueueItem* current,const unsigned* currentSize,
        QueueItem* next,unsigned* nextSize,unsigned queueCapacity,Status* states) {
    const unsigned lane=threadIdx.x&(lanes-1),base=(threadIdx.x&31u)-lane;
    const unsigned mask=(0xffffffffu>>(32-lanes))<<base;
    const unsigned group=(blockIdx.x*blockDim.x+threadIdx.x)/lanes;
    const unsigned stride=blockDim.x*gridDim.x/lanes,size=min(*currentSize,queueCapacity);
    for(unsigned index=group;index<size;index+=stride){
        const QueueItem& item=current[index];const int orientation=item.orientation;
        const Input& in=inputs[orientation];const Beam& beam=item.beam;
        Status status={};int active=0;
        if(lane==0)active=!atomicAdd(&states[orientation].failure,0)
            &&beam.depth<in.maxDepth&&!pruned(in,beam,status);
        active=__shfl_sync(mask,active,0,lanes);
        if(!active)continue;
        for(int f=lane;f<in.faceCount&&!status.failure;f+=lanes){
            if(f==beam.facet)continue;
            const double cosine=maximum(-1.0,minimum(1.0,dot(in.faces[f].outward,beam.direction)));
            if(cosine<=Epsilon)continue;
            ClippedAperture result;intersect(in,beam,f,result,status);
            if(status.failure||result.count<3)continue;
            const unsigned slot=atomicAdd(nextSize,1u);
            if(slot>=queueCapacity){status.failure=3;break;}
            QueueItem& target=next[slot];target.orientation=orientation;Beam& dest=target.beam;
            dest.count=result.count;dest.facet=f;dest.depth=beam.depth;
            dest.direction=beam.direction;dest.basis=beam.basis;dest.path=beam.path;dest.front=beam.front;
            for(int j=0;j<4;++j)dest.j[j]=beam.j[j];
            for(int d=0;d<=beam.depth;++d)dest.track[d]=beam.track[d];
            for(int v=0;v<result.count;++v)dest.vertices[v]=result.vertices[v];
        }
        for(int step=lanes/2;step;step/=2)status.failure=max(status.failure,__shfl_down_sync(mask,status.failure,step,lanes));
        status.evaluated=warpSum<lanes>(status.evaluated,mask);
        if(lane==0)addNonzeroStatus(states,orientation,status);
    }
}
__global__ void splitLevel(const Input* inputs,QueueItem* staged,const unsigned* stagedSize,unsigned queueCapacity,
        Output* output,unsigned* total,unsigned capacity,Status* states) {
    const unsigned size=min(*stagedSize,queueCapacity);
    for(unsigned index=blockIdx.x*blockDim.x+threadIdx.x;index<size;index+=blockDim.x*gridDim.x){
        QueueItem& item=staged[index];const int orientation=item.orientation;
        if(atomicAdd(&states[orientation].failure,0))continue;
        const Input& in=inputs[orientation];const Beam& beam=item.beam;const int f=beam.facet;
        const double cosine=maximum(-1.0,minimum(1.0,dot(in.faces[f].outward,beam.direction)));
        Status status={};OpticalView reflected={},transmitted;
        reflected.count=beam.count;reflected.vertices=beam.vertices;
        // split uses the parent's directions, Jones matrix, front and path;
        // its original polygon/facet is not needed after successful clipping.
        const bool outgoing=split(in,beam,in.faces[f].outward,cosine,reflected,transmitted);
        reflected.facet=transmitted.facet=f;reflected.depth=transmitted.depth=beam.depth+1;
        for(int p=0;p<=beam.depth;++p)reflected.track[p]=transmitted.track[p]=beam.track[p];
        reflected.track[beam.depth+1]=transmitted.track[beam.depth+1]=f;
        if(!finiteBeam(reflected)||(outgoing&&!finiteBeam(transmitted))){status.failure=4;addNonzeroStatus(states,orientation,status);continue;}
        if(outgoing){
            const unsigned slot=atomicAdd(total,1u);
            if(slot>=capacity){status.failure=2;addNonzeroStatus(states,orientation,status);continue;}
            storeBeam(output[slot].beam,transmitted);output[slot].orientation=orientation;output[slot].order=0;++status.count;
        }
        const bool rejected=pruned(in,reflected,status);
        if(!rejected&&reflected.depth<in.maxDepth)storeBeam(item.beam,reflected);
        else item.beam.depth=in.maxDepth; // inactive: do not process this slot again
        addNonzeroStatus(states,orientation,status);
    }
}
