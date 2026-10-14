#include "GpuConvexQueueMath.cuh"
#include "GpuProcessLock.h"
#include <cuda_runtime.h>
#include <algorithm>
#include <chrono>
#include <cstdio>
#include <climits>
#include <cstdlib>

namespace GpuConvexQueue {
namespace {
__global__ void residentTrace(const Input* inputs,int count,Beam* stacks,
                             Output* output,unsigned* total,unsigned capacity,Status* states) {
    int orientation=blockIdx.x*blockDim.x+threadIdx.x;if(orientation>=count)return;
    const Input& in=inputs[orientation];Status status={};
    int size=in.seedCount;
    // Stride the per-orientation stacks to expose independent orientations.
    for(int i=0;i<size;++i)stacks[i*count+orientation]=in.seeds[i];
    while(size&&!status.failure){
        Beam beam=stacks[(--size)*count+orientation];
        if(beam.depth>=in.maxDepth||pruned(in,beam,status))continue;
        for(int f=0;f<in.faceCount&&!status.failure;++f){
            if(f==beam.facet)continue;
            Vec normal=in.faces[f].outward;
            double cosine=maximum(-1.0,minimum(1.0,dot(normal,beam.direction)));
            if(cosine<=Epsilon)continue;
            Beam reflected,transmitted;
            intersect(in,beam,f,reflected,status);
            if(status.failure||reflected.count<3)continue;
            bool outgoing=split(in,beam,normal,cosine,reflected,transmitted);
            reflected.facet=transmitted.facet=f;
            reflected.depth=transmitted.depth=beam.depth+1;
            for(int p=0;p<=beam.depth;++p)reflected.track[p]=transmitted.track[p]=beam.track[p];
            reflected.track[beam.depth+1]=transmitted.track[beam.depth+1]=f;
            if(!finiteBeam(reflected)||(outgoing&&!finiteBeam(transmitted))){status.failure=4;break;}
            if(outgoing){
                unsigned slot=atomicAdd(total,1u);
                if(slot>=capacity){status.failure=2;break;}
                output[slot].beam=transmitted;output[slot].orientation=orientation;output[slot].order=status.count++;
            }
            if(!pruned(in,reflected,status)){
                if(size>=StackCapacity){status.failure=3;break;}
                stacks[size++*count+orientation]=reflected;
            }
        }
    }
    states[orientation]=status;
}
struct QueueItem { Beam beam;int orientation; };
__global__ void seedQueue(const Input* inputs,int count,QueueItem* queue,unsigned* size) {
    int i=blockIdx.x*blockDim.x+threadIdx.x;if(i>=count)return;
    unsigned start=atomicAdd(size,static_cast<unsigned>(inputs[i].seedCount));
    for(int s=0;s<inputs[i].seedCount;++s){queue[start+s].beam=inputs[i].seeds[s];queue[start+s].orientation=i;}
}
__global__ void traceLevel(const Input* inputs,const QueueItem* current,const unsigned* currentSize,
                          QueueItem* next,unsigned* nextSize,unsigned queueCapacity,
                          Output* output,unsigned* total,unsigned capacity,Status* states) {
    unsigned index=blockIdx.x*blockDim.x+threadIdx.x;
    if(index>=*currentSize||index>=queueCapacity)return;
    const QueueItem& item=current[index];int orientation=item.orientation;
    if(atomicAdd(&states[orientation].failure,0))return;
    const Input& in=inputs[orientation];const Beam& beam=item.beam;Status status={};
    if(beam.depth>=in.maxDepth||pruned(in,beam,status))return;
    for(int f=0;f<in.faceCount&&!status.failure;++f){
        if(f==beam.facet)continue;
        Vec normal=in.faces[f].outward;
        double cosine=maximum(-1.0,minimum(1.0,dot(normal,beam.direction)));
        if(cosine<=Epsilon)continue;
        Beam reflected,transmitted;intersect(in,beam,f,reflected,status);
        if(status.failure||reflected.count<3)continue;
        bool outgoing=split(in,beam,normal,cosine,reflected,transmitted);
        reflected.facet=transmitted.facet=f;reflected.depth=transmitted.depth=beam.depth+1;
        for(int p=0;p<=beam.depth;++p)reflected.track[p]=transmitted.track[p]=beam.track[p];
        reflected.track[beam.depth+1]=transmitted.track[beam.depth+1]=f;
        if(!finiteBeam(reflected)||(outgoing&&!finiteBeam(transmitted))){status.failure=4;break;}
        if(outgoing){
            unsigned slot=atomicAdd(total,1u);
            if(slot>=capacity){status.failure=2;break;}
            output[slot].beam=transmitted;output[slot].orientation=orientation;output[slot].order=0;++status.count;
        }
        if(!pruned(in,reflected,status)&&reflected.depth<in.maxDepth){
            unsigned slot=atomicAdd(nextSize,1u);
            if(slot>=queueCapacity){status.failure=3;break;}
            next[slot].beam=reflected;next[slot].orientation=orientation;
        }
    }
    if(status.failure)atomicExch(&states[orientation].failure,status.failure);
    atomicAdd(&states[orientation].count,status.count);
    atomicAdd(&states[orientation].evaluated,status.evaluated);atomicAdd(&states[orientation].rejected,status.rejected);
    atomicAdd(&states[orientation].rejectedJones,status.rejectedJones);atomicAdd(&states[orientation].rejectedArea,status.rejectedArea);
    atomicAdd(&states[orientation].rejectedImportance,status.rejectedImportance);
}
struct Buffers {
    Input* inputs=nullptr;Beam* stacks=nullptr;Output* output=nullptr;Status* status=nullptr;unsigned* total=nullptr;
    cudaStream_t stream=nullptr;cudaEvent_t start=nullptr,stop=nullptr;
    QueueItem* queues[2]={nullptr,nullptr};unsigned* queueSizes=nullptr;
    ~Buffers(){if(stream)cudaStreamSynchronize(stream);cudaFree(inputs);cudaFree(stacks);cudaFree(output);cudaFree(status);cudaFree(total);
        cudaFree(queues[0]);cudaFree(queues[1]);cudaFree(queueSizes);
        if(start)cudaEventDestroy(start);if(stop)cudaEventDestroy(stop);if(stream)cudaStreamDestroy(stream);}
};
}
bool Trace(const std::vector<Input>& inputs,std::vector<Output>& output,std::vector<size_t>& order,std::vector<Status>& status,Metrics& metrics,int threads) {
    output.clear();order.clear();status.clear();metrics=Metrics();if(inputs.empty())return true;
    const auto started=std::chrono::steady_clock::now();
    GpuProcessLock lock;Buffers b;
    size_t freeBytes=0,totalBytes=0;
    const int count=static_cast<int>(inputs.size());
    const char* levelsEnv=std::getenv("MBS_GPU_TRACE_CONVEX_LEVELS");
    bool levels=levelsEnv&&levelsEnv[0]=='1';
    unsigned queueCapacity=static_cast<unsigned>(count)*512;
    size_t stackBytes=levels?2*size_t(queueCapacity)*sizeof(QueueItem):size_t(count)*StackCapacity*sizeof(Beam);
    size_t inputBytes=inputs.size()*sizeof(Input);
    if(cudaMemGetInfo(&freeBytes,&totalBytes)!=cudaSuccess)return false;
    // At most one quarter of currently free memory; leave diffraction room.
    size_t budget=freeBytes/4;
    if(budget<stackBytes+inputBytes+sizeof(Output)*size_t(count)*64)return false;
    size_t cap=std::min(size_t(count)*2048,(budget-stackBytes-inputBytes)/sizeof(Output));
    if(cap>UINT_MAX)return false;
    metrics.bytes=stackBytes+inputBytes+cap*sizeof(Output)+inputs.size()*sizeof(Status)+sizeof(unsigned);
    const auto check=[](cudaError_t e){if(e!=cudaSuccess){std::fprintf(stderr,"GPU convex queue: %s; complete CPU fallback\n",cudaGetErrorString(e));return false;}return true;};
    if(!check(cudaStreamCreateWithFlags(&b.stream,cudaStreamNonBlocking))
       ||!check(cudaEventCreate(&b.start))||!check(cudaEventCreate(&b.stop))
       ||!check(cudaMalloc(&b.inputs,inputBytes))
       ||!check(cudaMalloc(&b.output,cap*sizeof(Output)))||!check(cudaMalloc(&b.status,inputs.size()*sizeof(Status)))
       ||!check(cudaMalloc(&b.total,sizeof(unsigned))))return false;
    if(levels){
        if(!check(cudaMalloc(&b.queues[0],queueCapacity*sizeof(QueueItem)))
           ||!check(cudaMalloc(&b.queues[1],queueCapacity*sizeof(QueueItem)))||!check(cudaMalloc(&b.queueSizes,2*sizeof(unsigned)))
           ||!check(cudaMemsetAsync(b.queueSizes,0,2*sizeof(unsigned),b.stream))
           ||!check(cudaMemsetAsync(b.status,0,inputs.size()*sizeof(Status),b.stream)))return false;
    }else if(!check(cudaMalloc(&b.stacks,stackBytes)))return false;
    if(!check(cudaMemcpyAsync(b.inputs,inputs.data(),inputBytes,cudaMemcpyHostToDevice,b.stream))
       ||!check(cudaMemsetAsync(b.total,0,sizeof(unsigned),b.stream))||!check(cudaEventRecord(b.start,b.stream)))return false;
    if(levels){
        seedQueue<<<(count+63)/64,64,0,b.stream>>>(b.inputs,count,b.queues[0],b.queueSizes);
        int depth=0;for(const auto& in:inputs)depth=std::max(depth,in.maxDepth);
        for(int d=0;d<depth;++d){int current=d%2,next=1-current;
            if(!check(cudaMemsetAsync(b.queueSizes+next,0,sizeof(unsigned),b.stream)))return false;
            traceLevel<<<(queueCapacity+63)/64,64,0,b.stream>>>(b.inputs,b.queues[current],b.queueSizes+current,
                b.queues[next],b.queueSizes+next,queueCapacity,b.output,b.total,static_cast<unsigned>(cap),b.status);
            if(!check(cudaGetLastError()))return false;
        }
    }else residentTrace<<<(count+63)/64,64,0,b.stream>>>(b.inputs,count,b.stacks,b.output,b.total,static_cast<unsigned>(cap),b.status);
    if(!check(cudaGetLastError())||!check(cudaEventRecord(b.stop,b.stream)))return false;
    unsigned used=0;
    status.resize(inputs.size());
    if(!check(cudaMemcpyAsync(&used,b.total,sizeof(unsigned),cudaMemcpyDeviceToHost,b.stream))
       ||!check(cudaMemcpyAsync(status.data(),b.status,inputs.size()*sizeof(Status),cudaMemcpyDeviceToHost,b.stream))
       ||!check(cudaStreamSynchronize(b.stream)))return false;
    float ms=0;if(!check(cudaEventElapsedTime(&ms,b.start,b.stop)))return false;metrics.kernelSeconds=ms/1000.0;
    // An overflow invalidates only orientations flagged by the device. Entries
    // written earlier by other orientations remain complete and independently valid.
    const auto downloadStarted=std::chrono::steady_clock::now();
    output.resize(std::min(size_t(used),cap));
    if(!check(cudaMemcpyAsync(output.data(),b.output,output.size()*sizeof(Output),cudaMemcpyDeviceToHost,b.stream))
       ||!check(cudaStreamSynchronize(b.stream)))return false;
    metrics.downloadSeconds=std::chrono::duration<double>(std::chrono::steady_clock::now()-downloadStarted).count();
    const auto sortStarted=std::chrono::steady_clock::now();
    std::vector<size_t> offsets(count+1,0);
    for(const auto& item:output)++offsets[item.orientation+1];
    for(int i=0;i<count;++i)offsets[i+1]+=offsets[i];
    std::vector<size_t> cursor=offsets;order.resize(output.size());
    for(size_t i=0;i<output.size();++i)order[cursor[output[i].orientation]++]=i;
    #pragma omp parallel for num_threads(threads) schedule(dynamic,4)
    for(int orientation=0;orientation<count;++orientation)
    std::sort(order.begin()+offsets[orientation],order.begin()+offsets[orientation+1],[levels,&output](size_t ai,size_t bi){
        const Output& a=output[ai];const Output& b=output[bi];
        if(a.orientation!=b.orientation)return a.orientation<b.orientation;
        if(!levels)return a.order<b.order;
        // Restore the original CPU DFS order: parent first, descending child
        // facets, and ascending outgoing facets within each visited parent.
        int length=std::min(a.beam.depth,b.beam.depth);
        for(int i=0;i<length;++i)if(a.beam.track[i]!=b.beam.track[i])return a.beam.track[i]>b.beam.track[i];
        if(a.beam.depth!=b.beam.depth)return a.beam.depth<b.beam.depth;
        return a.beam.facet<b.beam.facet;
    });
    metrics.sortSeconds=std::chrono::duration<double>(std::chrono::steady_clock::now()-sortStarted).count();
    metrics.wallSeconds=std::chrono::duration<double>(std::chrono::steady_clock::now()-started).count();
    return true;
}
}
