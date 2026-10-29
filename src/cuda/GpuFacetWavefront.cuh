// Included in GpuConvexQueue's anonymous namespace. One subgroup visits one
// aperture; its 4..32 lanes test candidate facets independently. Only integer queue
// counters are reduced here: optical amplitudes never cross orientations.
template<int lanes,class T> __device__ T warpSum(T value,unsigned mask) {
    for(int step=lanes/2;step;step/=2)value+=__shfl_down_sync(mask,value,step,lanes);
    return value;
}
template<int lanes> __global__ void traceFacets(const Input* inputs,const QueueItem* current,const unsigned* currentSize,
        QueueItem* next,unsigned* nextSize,unsigned queueCapacity,
        Output* output,unsigned* total,unsigned capacity,Status* states) {
    const unsigned lane=threadIdx.x&(lanes-1),base=(threadIdx.x&31u)-lane;
    const unsigned mask=(0xffffffffu>>(32-lanes))<<base;
    const unsigned warp=(blockIdx.x*blockDim.x+threadIdx.x)/lanes;
    const unsigned stride=blockDim.x*gridDim.x/lanes,size=min(*currentSize,queueCapacity);
    for(unsigned index=warp;index<size;index+=stride){
        const QueueItem& item=current[index];const int orientation=item.orientation;
        const Input& in=inputs[orientation];const Beam& beam=item.beam;
        Status status={};int active=0;
        if(lane==0)active=!atomicAdd(&states[orientation].failure,0)
            &&beam.depth<in.maxDepth&&!pruned(in,beam,status);
        active=__shfl_sync(mask,active,0,lanes);
        if(!active)continue;
        for(int f=lane;f<in.faceCount&&!status.failure;f+=lanes)
            traceFacet(in,beam,orientation,f,next,nextSize,queueCapacity,output,total,capacity,status);
        for(int step=lanes/2;step;step/=2)status.failure=max(status.failure,__shfl_down_sync(mask,status.failure,step,lanes));
        status.count=warpSum<lanes>(status.count,mask);
        status.evaluated=warpSum<lanes>(status.evaluated,mask);status.rejected=warpSum<lanes>(status.rejected,mask);
        status.rejectedJones=warpSum<lanes>(status.rejectedJones,mask);status.rejectedArea=warpSum<lanes>(status.rejectedArea,mask);
        status.rejectedImportance=warpSum<lanes>(status.rejectedImportance,mask);
        if(lane==0)addStatus(states,orientation,status);
    }
}
