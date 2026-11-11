#include "GpuConvexQueueMath.cuh"
#include "GpuTraceCapacity.h"
#include "GpuProcessLock.h"
#include <cuda_runtime.h>
#include <algorithm>
#include <atomic>
#include <chrono>
#include <cstdio>
#include <climits>
#include <cstdlib>
#include <list>
#include <mutex>
#include <memory>
#include <string>
#include <stdexcept>

namespace GpuConvexQueue {
#include "GpuTraceReuse.cuh"
namespace {
__global__ void residentTrace(const Input* inputs,int count,Beam* stacks,
                             Output* output,unsigned* total,unsigned capacity,Status* states,int stackCapacity) {
    int orientation=blockIdx.x*blockDim.x+threadIdx.x;if(orientation>=count)return;
    const Input& in=inputs[orientation];Status status={};
    int size=in.seedCount;
    // Stride the per-orientation stacks to expose independent orientations.
    for(int i=0;i<size;++i)stacks[size_t(i)*count+orientation]=in.seeds[i];
    while(size&&!status.failure){
        Beam beam=stacks[size_t(--size)*count+orientation];
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
                if(size>=stackCapacity){status.failure=3;break;}
                stacks[size_t(size++)*count+orientation]=reflected;
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
__device__ __forceinline__ void traceFacet(const Input& in,const Beam& beam,int orientation,int f,
        QueueItem* next,unsigned* nextSize,unsigned queueCapacity,
        Output* output,unsigned* total,unsigned capacity,Status& status) {
        if(f==beam.facet)return;
        Vec normal=in.faces[f].outward;
        double cosine=maximum(-1.0,minimum(1.0,dot(normal,beam.direction)));
        if(cosine<=Epsilon)return;
        Beam reflected,transmitted;intersect(in,beam,f,reflected,status);
        if(status.failure||reflected.count<3)return;
        bool outgoing=split(in,beam,normal,cosine,reflected,transmitted);
        reflected.facet=transmitted.facet=f;reflected.depth=transmitted.depth=beam.depth+1;
        for(int p=0;p<=beam.depth;++p)reflected.track[p]=transmitted.track[p]=beam.track[p];
        reflected.track[beam.depth+1]=transmitted.track[beam.depth+1]=f;
        if(!finiteBeam(reflected)||(outgoing&&!finiteBeam(transmitted))){status.failure=4;return;}
        if(outgoing){
            unsigned slot=atomicAdd(total,1u);
            if(slot>=capacity){status.failure=2;return;}
            output[slot].beam=transmitted;output[slot].orientation=orientation;output[slot].order=0;++status.count;
        }
        if(!pruned(in,reflected,status)&&reflected.depth<in.maxDepth){
            unsigned slot=atomicAdd(nextSize,1u);
            if(slot>=queueCapacity){status.failure=3;return;}
            next[slot].beam=reflected;next[slot].orientation=orientation;
        }
}
__device__ __forceinline__ void addFailure(int* destination,int failure){
    int previous=atomicAdd(destination,0);
    for(;;){const int combined=CombineFailure(previous,failure);if(combined==previous)return;
        const int observed=atomicCAS(destination,previous,combined);if(observed==previous)return;previous=observed;}
}
__device__ __forceinline__ void addStatus(Status* states,int orientation,const Status& status) {
    if(status.failure)addFailure(&states[orientation].failure,status.failure);
    atomicAdd(&states[orientation].count,status.count);
    atomicAdd(&states[orientation].evaluated,status.evaluated);atomicAdd(&states[orientation].rejected,status.rejected);
    atomicAdd(&states[orientation].rejectedJones,status.rejectedJones);atomicAdd(&states[orientation].rejectedArea,status.rejectedArea);
    atomicAdd(&states[orientation].rejectedImportance,status.rejectedImportance);
}
__device__ __forceinline__ void traceQueuedBeam(unsigned index,const Input* inputs,const QueueItem* current,
        QueueItem* next,unsigned* nextSize,unsigned queueCapacity,
        Output* output,unsigned* total,unsigned capacity,Status* states) {
    const QueueItem& item=current[index];int orientation=item.orientation;
    if(atomicAdd(&states[orientation].failure,0))return;
    const Input& in=inputs[orientation];const Beam& beam=item.beam;Status status={};
    if(beam.depth>=in.maxDepth||pruned(in,beam,status))return;
    for(int f=0;f<in.faceCount&&!status.failure;++f)
        traceFacet(in,beam,orientation,f,next,nextSize,queueCapacity,output,total,capacity,status);
    addStatus(states,orientation,status);
}
template<bool compact> __global__ void traceLevel(const Input* inputs,const QueueItem* current,const unsigned* currentSize,
        QueueItem* next,unsigned* nextSize,unsigned queueCapacity,
        Output* output,unsigned* total,unsigned capacity,Status* states) {
    const unsigned size=min(*currentSize,queueCapacity);
    for(unsigned index=blockIdx.x*blockDim.x+threadIdx.x;index<size;index+=blockDim.x*gridDim.x){
        traceQueuedBeam(index,inputs,current,next,nextSize,queueCapacity,output,total,capacity,states);
        if(!compact)break;
    }
}
#include "GpuFacetWavefront.cuh"
#include "GpuSplitWavefront.cuh"
struct Buffers {
    Input* inputs=nullptr;Beam* stacks=nullptr;Output* output=nullptr;Status* status=nullptr;unsigned* total=nullptr;
    cudaStream_t stream=nullptr;cudaEvent_t start=nullptr,stop=nullptr;
    QueueItem* queues[2]={nullptr,nullptr};unsigned* queueSizes=nullptr;
    enum Slot {Inputs,OutputData,States,Total,Stack,Queue0,Queue1,QueueSizes,Slots};
    size_t sizes[Slots]={};int device=0,reused=0,allocated=0;
    bool reuseLarger=false,healthy=true;
    explicit Buffers(int id):device(id){}
    template<class T> cudaError_t Allocate(T** address,size_t bytes,Slot slot){
        if(*address&&(sizes[slot]==bytes||(reuseLarger&&sizes[slot]>=bytes))){++reused;return cudaSuccess;}
        if(*address){const auto error=cudaFree(*address);if(error!=cudaSuccess){healthy=false;return error;}*address=nullptr;sizes[slot]=0;}
        const auto error=cudaMalloc(address,bytes);if(error==cudaSuccess){sizes[slot]=bytes;++allocated;}return error;
    }
    template<class T> void Release(T*& address,Slot slot){
        if(address){const auto error=cudaFree(address);if(error!=cudaSuccess){healthy=false;return;}}
        address=nullptr;sizes[slot]=0;
    }
    template<class T> cudaError_t DiscardMismatch(T*& address,Slot slot,size_t desired){
        if(sizes[slot]==desired)return cudaSuccess;
        if(address){const auto error=cudaFree(address);if(error!=cudaSuccess){healthy=false;return error;}}
        address=nullptr;sizes[slot]=0;return cudaSuccess;
    }
    size_t Bytes()const{size_t total=0;for(size_t n:sizes)total+=n;return total;}
    void ReleaseScratch(){Release(stacks,Stack);Release(queues[0],Queue0);Release(queues[1],Queue1);}
    void Trim(size_t limit,bool keepScratch=false){
        if(!keepScratch||Bytes()>limit)ReleaseScratch();
        if(Bytes()>limit)Release(output,OutputData);
        if(Bytes()>limit)Release(inputs,Inputs);
        if(Bytes()>limit){Release(status,States);Release(total,Total);Release(queueSizes,QueueSizes);}
    }
    bool ReleaseAllOwned() noexcept {
        if(!Bytes()&&!stream&&!start&&!stop)return true;
        int before=0;
        if(cudaGetDevice(&before)!=cudaSuccess)return false;
        if(before!=device&&cudaSetDevice(device)!=cudaSuccess)return false;
        bool ok=!stream||cudaStreamSynchronize(stream)==cudaSuccess;
        if(ok){
            healthy=true;
            Release(inputs,Inputs);Release(output,OutputData);Release(status,States);Release(total,Total);
            ReleaseScratch();Release(queueSizes,QueueSizes);
            ok=healthy;
        }
        if(ok&&start){ok=cudaEventDestroy(start)==cudaSuccess;if(ok)start=nullptr;}
        if(ok&&stop){ok=cudaEventDestroy(stop)==cudaSuccess;if(ok)stop=nullptr;}
        if(ok&&stream){ok=cudaStreamDestroy(stream)==cudaSuccess;if(ok)stream=nullptr;}
        if(before!=device&&cudaSetDevice(before)!=cudaSuccess)ok=false;
        return ok;
    }
    // Normal runtime destruction follows a checked, complete retirement.
    // A quarantined owner reaches this best-effort path only at process exit.
    ~Buffers(){ReleaseAllOwned();}
};
#include "GpuTraceWorkspace.cuh"
}
bool ReleaseWorkspace(){
    auto& cache=RawCache();std::unique_lock<std::mutex> lock(cache.mutex);
    if(WorkspaceQuarantined().load(std::memory_order_acquire))return false;
    auto& saved=SavedWorkspace();const size_t bytes=saved?saved->Bytes():0;
    if(!RetireWorkspace(saved))return false;
    if(bytes)std::fprintf(stderr,"GPU convex buffer pool: idle_release_MiB=%.3f\n",bytes/1048576.);
    return true;
}
ChunkLimit ResidentChunkLimit(int requested,int faces,int maximumDepth,bool allowWide){
    ChunkLimit result;int device=0;
    if(cudaGetDevice(&device)!=cudaSuccess){cudaGetLastError();return result;}
    auto& cache=RawCache();std::unique_lock<std::mutex> lock(cache.mutex);
    size_t free=0,total=0,available=0;
    if(cudaMemGetInfo(&free,&total)!=cudaSuccess){cudaGetLastError();return result;}
    RequireHealthyWorkspace();
    const auto& workspace=SavedWorkspace();
    const size_t retained=workspace&&workspace->device==device?workspace->Bytes():0;
    if(!Capacity::Add(free,retained,available))return result;
    result.queried=true;
    // Trace creates its stream/events before querying free memory. Leave a
    // small additional margin here; concurrent allocations still require its
    // existing checked allocations and complete CPU fallback.
    const size_t quarter=available/4,margin=4*1024*1024;
    result.budgetBytes=quarter>margin?quarter-margin:0;
    const char* levelsValue=std::getenv("MBS_GPU_TRACE_CONVEX_LEVELS");
    const bool levels=levelsValue&&levelsValue[0]=='1';
    const int fitting=Capacity::FittingChunk(requested,faces,maximumDepth,levels,
        result.budgetBytes,sizeof(QueueItem),result.reservedBytes,allowWide);
    result.fits=fitting>0;result.count=std::max(1,fitting);return result;
}
static bool TraceImpl(const std::vector<Input>& inputs,std::vector<Output>& output,std::vector<size_t>& order,std::vector<Status>& status,Metrics& metrics,int threads,
           const std::vector<Output>* external,
           const std::function<bool(const Input*,const Output*,size_t)>& consumeDevice,const std::string& capacityIdentity) {
    output.clear();order.clear();status.clear();metrics=Metrics();if(inputs.empty())return true;
    const auto started=std::chrono::steady_clock::now();
    const bool deepCapacity=Capacity::DeepEnabled();
    if(inputs.size()>size_t(INT_MAX)-63)return false;
    for(const auto& in:inputs)if(!Capacity::ValidInput(in,deepCapacity)){
        std::fprintf(stderr,"GPU convex queue fallback: invalid input/depth bounds; accepted_GPU=0\n");return false;}
    if(external){
        if(external->size()>size_t(UINT_MAX))return false;
        for(const auto& item:*external){const auto& beam=item.beam;
            if(item.orientation<0||size_t(item.orientation)>=inputs.size()||beam.count<3||beam.count>MaxVertices
               ||beam.depth<0||beam.depth>inputs[item.orientation].maxDepth
               ||(beam.facet!=INT_MAX&&(beam.facet<0||beam.facet>=inputs[item.orientation].faceCount)))return false;
        }
    }
    const char* scheduleEnv=std::getenv("MBS_GPU_TRACE_CONVEX_SCHEDULE");
    const std::string schedule=scheduleEnv?scheduleEnv:"legacy";
    const bool splitSchedule=schedule=="split"||schedule=="split8";
    if(schedule!="legacy"&&schedule!="compact"&&schedule!="facets"&&schedule!="facets4"&&schedule!="facets8"&&schedule!="facets16"&&!splitSchedule)
        throw std::runtime_error("MBS_GPU_TRACE_CONVEX_SCHEDULE must be legacy, compact, facets4, facets8, facets16, facets, split or split8");
    GpuProcessLock lock;
    auto& cache=RawCache();std::unique_lock<std::mutex> cacheLock(cache.mutex);
    RequireHealthyWorkspace();
    std::string cacheKey;int device=0;cudaGetDevice(&device);
    const char* levelsValue=std::getenv("MBS_GPU_TRACE_CONVEX_LEVELS");
    const bool hintActive=Capacity::HintEnabled()&&cache.capacityHint.Select(capacityIdentity,device,schedule,levelsValue&&levelsValue[0]=='1');
    if(!hintActive)cache.capacityHint.Reset();
    Capacity::HintAttempt hintAttempt(hintActive?&cache.capacityHint:nullptr);
    if(consumeDevice&&cache.limit){
        const size_t inputCount=inputs.size(),externalCount=external?external->size():0;
        cacheKey.append(reinterpret_cast<const char*>(&inputCount),sizeof(inputCount));
        cacheKey.append(reinterpret_cast<const char*>(inputs.data()),inputs.size()*sizeof(Input));
        cacheKey.append(reinterpret_cast<const char*>(&externalCount),sizeof(externalCount));
        if(externalCount)cacheKey.append(reinterpret_cast<const char*>(external->data()),externalCount*sizeof(Output));
        const char* levels=std::getenv("MBS_GPU_TRACE_CONVEX_LEVELS");cacheKey+=levels?levels:"0";cacheKey+=" schedule="+schedule;
        for(auto i=cache.entries.begin();i!=cache.entries.end();++i)if((*i)->device==device&&(*i)->key==cacheKey){
            auto& saved=**i;status=saved.status;metrics.bytes=saved.bytes;
            const bool ok=consumeDevice(saved.inputs,saved.output,saved.count);
            metrics.wallSeconds=std::chrono::duration<double>(std::chrono::steady_clock::now()-started).count();
            std::fprintf(stderr,"GPU raw trace reuse: hit orientations=%zu raw_beams=%zu cached_MiB=%.3f\n",inputs.size(),saved.count,cache.bytes/1048576.);
            if(ok)cache.entries.splice(cache.entries.begin(),cache.entries,i);
            else{cache.bytes-=saved.bytes;cache.entries.erase(i);}
            hintAttempt.successful=ok;
            return ok;
        }
    }
    WorkspaceLease workspace(device);Buffers& b=*workspace.buffer;
    if(!b.healthy)return false;
    size_t freeBytes=0,totalBytes=0;
    const int count=static_cast<int>(inputs.size());
    const char* levelsEnv=std::getenv("MBS_GPU_TRACE_CONVEX_LEVELS");
    bool levels=levelsEnv&&levelsEnv[0]=='1';
    if(splitSchedule&&!levels)throw std::runtime_error("split scheduling requires MBS_GPU_TRACE_CONVEX_LEVELS=1");
    const unsigned initial=external?static_cast<unsigned>(external->size()):0;
    const auto check=[](cudaError_t e){if(e!=cudaSuccess){std::fprintf(stderr,"GPU convex queue: %s; complete CPU fallback\n",cudaGetErrorString(e));return false;}return true;};
    if((!b.stream&&!check(cudaStreamCreateWithFlags(&b.stream,cudaStreamNonBlocking)))
       ||(!b.start&&!check(cudaEventCreate(&b.start)))||(!b.stop&&!check(cudaEventCreate(&b.stop))))return false;
    // Query after stream/event creation. Retained typed buffers are included
    // exactly once; context/driver allocations remain outside available RAM.
    size_t beforeReservation=0;
    if(!check(cudaMemGetInfo(&freeBytes,&totalBytes))||!Capacity::Add(freeBytes,b.Bytes(),beforeReservation))return false;
    const size_t budget=beforeReservation/4;
    size_t queuePerOrientation=0,beamsPerOrientation=0,queueWanted=0,outputWanted=0,minimumOutput=0;
    Capacity::PerOrientation(inputs,deepCapacity,queuePerOrientation,beamsPerOrientation);
    if(!Capacity::Multiply(inputs.size(),queuePerOrientation,queueWanted)
       ||!Capacity::Multiply(inputs.size(),beamsPerOrientation,outputWanted)
       ||!Capacity::Multiply(inputs.size(),64,minimumOutput))return false;
    Capacity::Plan plan;
    if(!Capacity::MakePlan(inputs.size(),levels,queueWanted,outputWanted,std::max(minimumOutput,size_t(initial)),
                           budget,sizeof(QueueItem),plan,deepCapacity?StackCapacity:LegacyStackCapacity)){
        std::fprintf(stderr,"GPU convex queue fallback: initial capacity exceeds budget orientations=%d budget_MiB=%.3f accepted_GPU=0\n",count,budget/1048576.);return false;}
    const Capacity::Plan baselinePlan=plan;
    size_t maximumQueue=0,maximumOutput=0;
    if(!Capacity::Multiply(plan.queue,4,maximumQueue)||!Capacity::Multiply(plan.output,4,maximumOutput))return false;
    if(hintActive){Capacity::Plan hinted;
        const bool usedHint=cache.capacityHint.Apply(plan,inputs.size(),levels,budget,sizeof(QueueItem),hinted);
        if(usedHint)plan=hinted;
        std::fprintf(stderr,"GPU convex capacity hint: used=%d queue_per_orientation=%zu output_per_orientation=%zu baseline_queue=%zu baseline_output=%zu selected_queue=%zu selected_output=%zu budget_MiB=%.3f\n",
            usedHint?1:0,cache.capacityHint.queuePerOrientation,cache.capacityHint.outputPerOrientation,
            baselinePlan.queue,baselinePlan.output,plan.queue,plan.output,budget/1048576.);}
    const size_t inputBytes=plan.inputBytes;
    unsigned used=0,retries=0;size_t cap=0;
    for(;;){
    const unsigned queueCapacity=static_cast<unsigned>(plan.queue);cap=plan.output;
    const size_t stackBytes=plan.scratchBytes;
    if(deepCapacity)std::fprintf(stderr,"GPU convex capacity: attempt=%u orientations=%d queue_capacity=%u output_capacity=%zu planned_MiB=%.3f budget_MiB=%.3f\n",
        retries,count,queueCapacity,cap,plan.bytes/1048576.,budget/1048576.);
    // Drop every mismatched retained allocation before growing any slot. This
    // prevents old output/scratch plus its replacement exceeding the plan.
    const size_t queueBytes=levels?size_t(queueCapacity)*sizeof(QueueItem):0;
    static_assert(Buffers::Slots==8,"physical budget must cover every trace allocation");
    const Capacity::SlotBytes required{{inputBytes,cap*sizeof(Output),plan.statusBytes,sizeof(unsigned),
        levels?0:stackBytes,queueBytes,queueBytes,levels?2*sizeof(unsigned):0}};
    Capacity::SlotBytes existing{},physical{};std::copy(b.sizes,b.sizes+Buffers::Slots,existing.begin());
    size_t physicalBytes=0;
    if(!Capacity::PhysicalSlots(required,existing,budget,b.reuseLarger,physical,physicalBytes))return false;
    metrics.bytes=std::max(metrics.bytes,physicalBytes);
    if(workspace.pooled)std::fprintf(stderr,"GPU convex buffer pool plan: logical_MiB=%.3f physical_MiB=%.3f budget_MiB=%.3f\n",
        plan.bytes/1048576.,physicalBytes/1048576.,budget/1048576.);
    if(!check(b.DiscardMismatch(b.inputs,Buffers::Inputs,physical[Buffers::Inputs]))
       ||!check(b.DiscardMismatch(b.output,Buffers::OutputData,physical[Buffers::OutputData]))
       ||!check(b.DiscardMismatch(b.status,Buffers::States,physical[Buffers::States]))
       ||!check(b.DiscardMismatch(b.total,Buffers::Total,physical[Buffers::Total]))
       ||!check(b.DiscardMismatch(b.stacks,Buffers::Stack,physical[Buffers::Stack]))
       ||!check(b.DiscardMismatch(b.queues[0],Buffers::Queue0,physical[Buffers::Queue0]))
       ||!check(b.DiscardMismatch(b.queues[1],Buffers::Queue1,physical[Buffers::Queue1]))
       ||!check(b.DiscardMismatch(b.queueSizes,Buffers::QueueSizes,physical[Buffers::QueueSizes])))return false;
    if(!check(b.Allocate(&b.inputs,inputBytes,Buffers::Inputs))
       ||!check(b.Allocate(&b.output,cap*sizeof(Output),Buffers::OutputData))||!check(b.Allocate(&b.status,inputs.size()*sizeof(Status),Buffers::States))
       ||!check(b.Allocate(&b.total,sizeof(unsigned),Buffers::Total)))return false;
    if(levels){
        if(!check(b.Allocate(&b.queues[0],queueCapacity*sizeof(QueueItem),Buffers::Queue0))
           ||!check(b.Allocate(&b.queues[1],queueCapacity*sizeof(QueueItem),Buffers::Queue1))||!check(b.Allocate(&b.queueSizes,2*sizeof(unsigned),Buffers::QueueSizes))
           ||!check(cudaMemsetAsync(b.queueSizes,0,2*sizeof(unsigned),b.stream))
           ||!check(cudaMemsetAsync(b.status,0,inputs.size()*sizeof(Status),b.stream)))return false;
    }else if(!check(b.Allocate(&b.stacks,stackBytes,Buffers::Stack)))return false;
    if(initial&&!check(cudaMemcpyAsync(b.output,external->data(),initial*sizeof(Output),cudaMemcpyHostToDevice,b.stream)))return false;
    if(!check(cudaMemcpyAsync(b.inputs,inputs.data(),inputBytes,cudaMemcpyHostToDevice,b.stream))
       ||!check(cudaMemcpyAsync(b.total,&initial,sizeof(unsigned),cudaMemcpyHostToDevice,b.stream))||!check(cudaEventRecord(b.start,b.stream)))return false;
    if(levels){
        seedQueue<<<(count+63)/64,64,0,b.stream>>>(b.inputs,count,b.queues[0],b.queueSizes);
        int processors=0;if(!check(cudaDeviceGetAttribute(&processors,cudaDevAttrMultiProcessorCount,device)))return false;
        const unsigned compactBlocks=std::max(1,std::min(int((queueCapacity+63)/64),processors*4));
        const unsigned facetBlocks=std::max(1,std::min(int((queueCapacity+1)/2),processors*4));
        const unsigned splitBlocks=std::max(1,std::min(int((queueCapacity+63)/64),processors*16));
        int depth=0;for(const auto& in:inputs)depth=std::max(depth,in.maxDepth);
        for(int d=0;d<depth;++d){int current=d%2,next=1-current;
            if(!check(cudaMemsetAsync(b.queueSizes+next,0,sizeof(unsigned),b.stream)))return false;
            if(splitSchedule){
                if(schedule=="split8")clipLevel<8><<<splitBlocks,64,0,b.stream>>>(b.inputs,b.queues[current],b.queueSizes+current,
                    b.queues[next],b.queueSizes+next,queueCapacity,b.status);
                else clipLevel<1><<<splitBlocks,64,0,b.stream>>>(b.inputs,b.queues[current],b.queueSizes+current,
                    b.queues[next],b.queueSizes+next,queueCapacity,b.status);
                if(!check(cudaGetLastError()))return false;
                splitLevel<<<splitBlocks,64,0,b.stream>>>(b.inputs,b.queues[next],b.queueSizes+next,queueCapacity,
                    b.output,b.total,static_cast<unsigned>(cap),b.status);
            }
            else if(schedule=="facets4")traceFacets<4><<<facetBlocks,64,0,b.stream>>>(b.inputs,b.queues[current],b.queueSizes+current,
                b.queues[next],b.queueSizes+next,queueCapacity,b.output,b.total,static_cast<unsigned>(cap),b.status);
            else if(schedule=="facets8")traceFacets<8><<<facetBlocks,64,0,b.stream>>>(b.inputs,b.queues[current],b.queueSizes+current,
                b.queues[next],b.queueSizes+next,queueCapacity,b.output,b.total,static_cast<unsigned>(cap),b.status);
            else if(schedule=="facets16")traceFacets<16><<<facetBlocks,64,0,b.stream>>>(b.inputs,b.queues[current],b.queueSizes+current,
                b.queues[next],b.queueSizes+next,queueCapacity,b.output,b.total,static_cast<unsigned>(cap),b.status);
            else if(schedule=="facets")traceFacets<32><<<facetBlocks,64,0,b.stream>>>(b.inputs,b.queues[current],b.queueSizes+current,
                b.queues[next],b.queueSizes+next,queueCapacity,b.output,b.total,static_cast<unsigned>(cap),b.status);
            else if(schedule=="compact")traceLevel<true><<<compactBlocks,64,0,b.stream>>>(b.inputs,b.queues[current],b.queueSizes+current,
                b.queues[next],b.queueSizes+next,queueCapacity,b.output,b.total,static_cast<unsigned>(cap),b.status);
            else traceLevel<false><<<(queueCapacity+63)/64,64,0,b.stream>>>(b.inputs,b.queues[current],b.queueSizes+current,
                b.queues[next],b.queueSizes+next,queueCapacity,b.output,b.total,static_cast<unsigned>(cap),b.status);
            if(!check(cudaGetLastError()))return false;
        }
    }else residentTrace<<<(count+63)/64,64,0,b.stream>>>(b.inputs,count,b.stacks,b.output,b.total,static_cast<unsigned>(cap),b.status,plan.stackCapacity);
    if(!check(cudaGetLastError())||!check(cudaEventRecord(b.stop,b.stream)))return false;
    used=0;
    status.resize(inputs.size());
    if(!check(cudaMemcpyAsync(&used,b.total,sizeof(unsigned),cudaMemcpyDeviceToHost,b.stream))
       ||!check(cudaMemcpyAsync(status.data(),b.status,inputs.size()*sizeof(Status),cudaMemcpyDeviceToHost,b.stream))
       ||!check(cudaStreamSynchronize(b.stream)))return false;
    float ms=0;if(!check(cudaEventElapsedTime(&ms,b.start,b.stop)))return false;metrics.kernelSeconds+=ms/1000.0;
    int failures[5]={};bool failed=false,otherFailure=false;
    for(const auto& state:status)if(state.failure){failed=true;
        if(state.failure>0&&state.failure<=4)++failures[state.failure];
        if(state.failure!=2&&state.failure!=3)otherFailure=true;}
    if((consumeDevice||deepCapacity)&&(used>cap||failed)){
        Capacity::Plan next;
        const bool retry=deepCapacity&&Capacity::Retry(plan,inputs.size(),levels,failures[3]>0,
            used>cap||failures[2]>0,otherFailure,retries,budget,sizeof(QueueItem),next)
            &&(!hintActive||(next.queue<=maximumQueue&&next.output<=maximumOutput));
        std::fprintf(stderr,"GPU convex capacity %s: attempt=%u raw_beams=%u output_capacity=%zu queue_capacity=%u vertex_limit=%d output_limit=%d work_queue_limit=%d nonfinite=%d accepted_GPU=0\n",
            retry?"retry":"fallback",retries,used,cap,queueCapacity,failures[1],failures[2],failures[3],failures[4]);
        if(!retry)return false;
        plan=next;++retries;continue;
    }
    break;
    }
    if(consumeDevice){
        // Tracing scratch is no longer needed. Keep only inputs and raw output
        // until the consumer completes; no raw beam crosses PCIe.
        const bool keptScratch=workspace.KeepScratchForConsumer();
        if(!keptScratch)b.ReleaseScratch();
        if(!b.healthy)return false;
        if(workspace.pooled)std::fprintf(stderr,"GPU convex buffer pool consumer: scratch_retained=%d physical_MiB=%.3f\n",
            keptScratch?1:0,b.Bytes()/1048576.);
        const bool ok=consumeDevice(b.inputs,b.output,used);
        const size_t cacheBytes=b.sizes[Buffers::Inputs]+size_t(used)*sizeof(Output)+cacheKey.size()+status.size()*sizeof(Status)+sizeof(TraceCacheEntry);
        size_t compactReservation=0;
        const bool compactWithinBudget=Capacity::Add(b.Bytes(),std::max(size_t(1),size_t(used))*sizeof(Output),compactReservation)&&compactReservation<=budget;
        if(ok&&!cacheKey.empty()&&cacheBytes<=cache.limit&&((!deepCapacity&&!workspace.pooled)||compactWithinBudget)){
            while(cache.bytes+cacheBytes>cache.limit&&!cache.entries.empty()){cache.bytes-=cache.entries.back()->bytes;cache.entries.pop_back();}
            Output* compact=nullptr;
            // Cache allocation is optional and must never invalidate a valid result.
            if(cudaMalloc(&compact,std::max(size_t(1),size_t(used))*sizeof(Output))==cudaSuccess){
                if(cudaMemcpy(compact,b.output,size_t(used)*sizeof(Output),cudaMemcpyDeviceToDevice)==cudaSuccess){
                    std::unique_ptr<TraceCacheEntry> saved(new TraceCacheEntry);
                    saved->key=std::move(cacheKey);saved->device=device;saved->bytes=cacheBytes;saved->count=used;saved->status=status;
                    saved->inputs=b.inputs;b.inputs=nullptr;b.sizes[Buffers::Inputs]=0;saved->output=compact;cache.bytes+=cacheBytes;
                    cache.entries.push_front(std::move(saved));
                }else{cudaFree(compact);cudaGetLastError();}
            }else cudaGetLastError();
        }
        metrics.wallSeconds=std::chrono::duration<double>(std::chrono::steady_clock::now()-started).count();
        if(deepCapacity)std::fprintf(stderr,"GPU convex capacity result: accepted_GPU=%d retries=%u raw_beams=%u peak_planned_MiB=%.3f budget_MiB=%.3f\n",
            ok?count:0,retries,used,metrics.bytes/1048576.,budget/1048576.);
        if(ok&&hintActive)cache.capacityHint.Learn(plan,inputs.size());
        hintAttempt.successful=ok;
        workspace.success=ok;return ok;
    }
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
    if(deepCapacity)std::fprintf(stderr,"GPU convex capacity result: accepted_GPU=%d retries=%u raw_beams=%u peak_planned_MiB=%.3f budget_MiB=%.3f\n",
        count,retries,used,metrics.bytes/1048576.,budget/1048576.);
    workspace.success=std::all_of(status.begin(),status.end(),[](const Status& s){return s.failure==0;});
    if(workspace.success&&hintActive)cache.capacityHint.Learn(plan,inputs.size());
    hintAttempt.successful=workspace.success;return true;
}
bool Trace(const std::vector<Input>& inputs,std::vector<Output>& output,std::vector<size_t>& order,
           std::vector<Status>& status,Metrics& metrics,int threads,const std::vector<Output>* external,
           const std::function<bool(const Input*,const Output*,size_t)>& consumeDevice,const std::string& capacityIdentity){
    RequireHealthyWorkspace();
    const bool ok=TraceImpl(inputs,output,order,status,metrics,threads,external,consumeDevice,capacityIdentity);
    // TraceImpl's lease has retired before accepting any result or CPU retry.
    // Cleanup failure must not masquerade as a valid result or an empty cache.
    RequireHealthyWorkspace();
    return ok;
}
}
