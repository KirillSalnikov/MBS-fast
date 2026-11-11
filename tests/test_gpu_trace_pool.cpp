#include "cuda/GpuTraceCapacity.h"
#include <cassert>
#include <atomic>
#include <cstdio>
#include <memory>
#include <stdexcept>

using cudaError_t=int;
static const int cudaSuccess=0;
static size_t freeBytes=8ULL*1024*1024*1024;
static int currentDevice=0,queryFailure=0,alive=0,retireFailure=0,retireCalls=0,trimFailure=0;
static int cudaGetDevice(int* device){*device=currentDevice;return queryFailure;}
static int cudaMemGetInfo(size_t* free,size_t* total){*free=freeBytes;*total=16ULL*1024*1024*1024;return queryFailure;}
namespace GpuConvexQueue {
// Test the actual lease's ownership/zero-override/device/query decisions.
// The mock records trimming requests; no CUDA call or allocation is made.
struct Buffers {
    enum Slot {Inputs,OutputData,States,Total,Stack,Queue0,Queue1,QueueSizes,Slots};
    size_t sizes[Slots]={};int device,reused=0,allocated=0,trimCalls=0;
    bool healthy=true,reuseLarger=false,lastKeepScratch=false;
    explicit Buffers(int id):device(id){++alive;}
    ~Buffers(){--alive;}
    size_t Bytes()const{size_t n=0;for(size_t s:sizes)n+=s;return n;}
    bool ReleaseAllOwned() noexcept {++retireCalls;if(retireFailure)return false;
        for(auto& n:sizes)n=0;
        return true;}
    void Trim(size_t limit,bool keep=false){++trimCalls;lastKeepScratch=keep;
        if(trimFailure){healthy=false;return;}
        for(auto& size:sizes){const size_t total=Bytes();if(total<=limit)break;size-=std::min(size,total-limit);}}
};
#include "cuda/GpuTraceWorkspace.cuh"
}
int main(){
    using namespace GpuConvexQueue;
    const size_t mib=1024*1024,gib=1024*mib,maximum=std::numeric_limits<size_t>::max();
    Capacity::SlotBytes required{{10,20,3,4,0,50,60,8}},existing{{100,200,30,40,900,500,600,80}},physical{};
    size_t bytes=0;
    assert(Capacity::PhysicalSlots(required,existing,2000,true,physical,bytes));
    assert(physical[4]==0&&bytes==1550); // unused previous mode never retained
    assert(Capacity::PhysicalSlots(required,existing,200,true,physical,bytes));
    assert(physical==required&&bytes==155); // shrink every excess before growth
    assert(Capacity::PhysicalSlots(required,existing,2000,false,physical,bytes));
    assert(physical==required&&bytes==155);
    assert(!Capacity::PhysicalSlots(required,existing,154,true,physical,bytes));
    auto larger=required;larger[0]=500;
    assert(Capacity::PhysicalSlots(larger,existing,2000,true,physical,bytes)&&physical[0]==500&&bytes==1950);
    auto overflow=required;overflow[0]=maximum;
    assert(!Capacity::PhysicalSlots(overflow,existing,maximum,true,physical,bytes));
    overflow=existing;overflow[0]=maximum;
    assert(!Capacity::PhysicalSlots(required,overflow,maximum,true,physical,bytes));
    assert(Capacity::RetentionLimit(true,8*gib,gib,4*gib)==2*gib);
    assert(Capacity::RetentionLimit(true,2*gib,gib,2*gib)==3*gib/4);
    assert(Capacity::RetentionLimit(true,8*gib,gib,64*mib)==64*mib);
    assert(!Capacity::RetentionLimit(false,8*gib,gib,2*gib));
    assert(!Capacity::RetentionLimit(true,511*mib,gib,2*gib));
    assert(!Capacity::RetentionLimit(true,maximum,1,2*gib));

    unsetenv("MBS_GPU_TRACE_CONVEX_WORKSPACE_MIB");setenv("MBS_GPU_TRACE_BUFFER_POOL","1",1);
    {WorkspaceLease w(0);assert(w.pooled&&w.buffer->reuseLarger&&w.limit==2*gib);
        w.buffer->sizes[Buffers::Queue0]=256*mib;w.buffer->sizes[Buffers::OutputData]=512*mib;
        assert(w.KeepScratchForConsumer());w.success=true;}
    assert(SavedWorkspace()&&alive==1&&SavedWorkspace()->lastKeepScratch);
    {WorkspaceLease w(0);assert(!SavedWorkspace()&&alive==1&&w.buffer->Bytes()==768*mib);
        freeBytes=900*mib;assert(!w.KeepScratchForConsumer());freeBytes=8*gib;w.success=true;}
    // Explicit zero must drop the existing reservation and defeat opt-in.
    setenv("MBS_GPU_TRACE_CONVEX_WORKSPACE_MIB","0",1);
    {WorkspaceLease w(0);assert(!w.pooled&&!w.buffer->reuseLarger&&!w.limit&&!SavedWorkspace());w.success=true;}
    assert(!SavedWorkspace()&&!alive);
    setenv("MBS_GPU_TRACE_CONVEX_WORKSPACE_MIB","64",1);
    {WorkspaceLease w(0);assert(w.pooled&&w.limit==64*mib);w.buffer->sizes[Buffers::OutputData]=32*mib;w.success=true;}
    assert(SavedWorkspace()&&SavedWorkspace()->Bytes()==32*mib);
    // Device changes release the prior slot; they cannot add its bytes to a
    // different GPU's free-memory budget.
    currentDevice=1;
    {WorkspaceLease w(1);assert(w.buffer->device==1&&!w.buffer->Bytes()&&alive==1);w.success=true;}
    queryFailure=1;
    {WorkspaceLease w(1);assert(!w.limit&&!w.KeepScratchForConsumer());w.success=true;}
    assert(!SavedWorkspace()&&!alive);queryFailure=0;
    unsetenv("MBS_GPU_TRACE_CONVEX_WORKSPACE_MIB");
    {WorkspaceLease w(1);w.buffer->sizes[Buffers::Queue0]=64*mib;w.success=true;freeBytes=511*mib;}
    assert(!SavedWorkspace()&&!alive);freeBytes=8*gib;
    {WorkspaceLease w(1);w.buffer->sizes[Buffers::OutputData]=gib;/* failed computation */}
    assert(!SavedWorkspace()&&!alive);
    {WorkspaceLease w(1);w.success=true;currentDevice=2;}
    assert(!SavedWorkspace()&&!alive);
    // Default off preserves the old small workspace semantics.
    setenv("MBS_GPU_TRACE_BUFFER_POOL","0",1);setenv("MBS_GPU_TRACE_CONVEX_WORKSPACE_MIB","64",1);
    {WorkspaceLease w(2);assert(!w.pooled&&!w.buffer->reuseLarger&&w.limit==64*mib);w.success=true;}
    assert(SavedWorkspace()&&!SavedWorkspace()->lastKeepScratch);SavedWorkspace().reset();assert(!alive);
    setenv("MBS_GPU_TRACE_CONVEX_WORKSPACE_MIB","2048",1);bool threw=false;
    try{WorkspaceLease w(2);}catch(const std::runtime_error&){threw=true;}assert(threw);
    setenv("MBS_GPU_TRACE_BUFFER_POOL","1",1);
    {WorkspaceLease w(2);assert(w.limit==2*gib);}
    setenv("MBS_GPU_TRACE_CONVEX_WORKSPACE_MIB","2049",1);threw=false;
    try{WorkspaceLease w(2);}catch(const std::runtime_error&){threw=true;}assert(threw&&!alive);
    // Failed free has one fixed owner, never retries in a destructor, and
    // rejects every next lease before allocating anything new.
    setenv("MBS_GPU_TRACE_CONVEX_WORKSPACE_MIB","0",1);retireFailure=1;
    {WorkspaceLease w(2);w.buffer->sizes[Buffers::OutputData]=gib;w.success=true;}
    assert(WorkspaceQuarantined()&&SavedWorkspace()&&alive==1&&SavedWorkspace()->Bytes()==gib);
    const int beforeRetire=retireCalls;
    for(int i=0;i<5;++i){threw=false;try{WorkspaceLease w(2);}catch(const std::runtime_error&){threw=true;}
        assert(threw&&alive==1&&retireCalls==beforeRetire&&SavedWorkspace()->Bytes()==gib);}
    // Test-only teardown reset: production intentionally requires a new process.
    retireFailure=0;assert(RetireWorkspace(SavedWorkspace()));WorkspaceQuarantined()=false;
    assert(!alive);
    setenv("MBS_GPU_TRACE_CONVEX_WORKSPACE_MIB","64",1);
    {WorkspaceLease w(2);w.buffer->sizes[Buffers::Queue0]=32*mib;w.success=true;}
    trimFailure=1;threw=false;
    try{WorkspaceLease w(2);}catch(const std::runtime_error&){threw=true;}
    assert(threw&&WorkspaceQuarantined()&&alive==1&&SavedWorkspace()->Bytes()==32*mib);
    trimFailure=0;SavedWorkspace()->healthy=true;
    assert(RetireWorkspace(SavedWorkspace()));WorkspaceQuarantined()=false;assert(!alive);
    // A failure before retention must quarantine even when total bytes are zero.
    {WorkspaceLease w(2);w.success=true;w.buffer->healthy=false;}
    assert(WorkspaceQuarantined()&&SavedWorkspace()&&alive==1);
    SavedWorkspace()->healthy=true;assert(RetireWorkspace(SavedWorkspace()));WorkspaceQuarantined()=false;
    assert(!alive);
    setenv("MBS_GPU_TRACE_CONVEX_WORKSPACE_MIB","00",1);assert(Capacity::WorkspaceExplicitlyDisabled());
    setenv("MBS_GPU_TRACE_CONVEX_WORKSPACE_MIB","64",1);assert(!Capacity::WorkspaceExplicitlyDisabled());
    setenv("MBS_GPU_TRACE_CONVEX_WORKSPACE_MIB","0x",1);assert(!Capacity::WorkspaceExplicitlyDisabled());
    std::puts("PASS bounded quarantine/fail-closed lease, physical eight-slot budgets/grow/shrink/overflow, scratch lease, explicit0, device/query/pressure/failure release, default-off policy");
}
