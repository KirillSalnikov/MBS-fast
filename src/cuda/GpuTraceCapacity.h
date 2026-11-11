#pragma once
#include "GpuConvexQueue.h"
#include <algorithm>
#include <climits>
#include <cmath>
#include <cstdlib>
#include <limits>
#include <array>

// Host-only, allocation-free policy. Capacity failures never authorize using
// partial optical results. The caller retains one fixed byte budget on retry.
namespace GpuConvexQueue { namespace Capacity {
inline bool WorkspaceExplicitlyDisabled(){
    const char* value=std::getenv("MBS_GPU_TRACE_CONVEX_WORKSPACE_MIB");
    if(!value)return false;
    char* end=nullptr;const long mib=std::strtol(value,&end,10);
    return end!=value&&!*end&&mib==0;
}
inline bool DeepEnabled(){
    const char* value=std::getenv("MBS_GPU_TRACE_DEEP_CAPACITY");
    return value&&value[0]=='1'&&value[1]=='\0';
}
inline bool WideChunksEnabled(){
    const char* value=std::getenv("MBS_GPU_TRACE_WIDE_CHUNKS");
    return DeepEnabled()&&value&&value[0]=='1'&&value[1]=='\0';
}
inline bool HintEnabled(){
    const char* value=std::getenv("MBS_GPU_TRACE_CAPACITY_HINT");
    return DeepEnabled()&&value&&value[0]=='1'&&value[1]=='\0';
}
inline bool Add(size_t a,size_t b,size_t& out){
    if(b>std::numeric_limits<size_t>::max()-a)return false;
    out=a+b;return true;
}
inline bool Multiply(size_t a,size_t b,size_t& out){
    if(a&&b>std::numeric_limits<size_t>::max()/a)return false;
    out=a*b;return true;
}
inline bool ValidInput(const Input& in,bool deep){
    const int maximum=deep?MaxDepth:LegacyMaxDepth;
    if(in.maxDepth<0||in.maxDepth>maximum||in.faceCount<0||in.faceCount>MaxFacets
       ||in.seedCount<0||in.seedCount>in.faceCount)return false;
    // A non-direct batch may contain an unsupported orientation marked empty.
    if(!in.faceCount)return !in.seedCount;
    if(!(in.ri>0)||!std::isfinite(in.ri))return false;
    for(int f=0;f<in.faceCount;++f)
        if(in.faces[f].count<3||in.faces[f].count>MaxVertices)return false;
    for(int s=0;s<in.seedCount;++s){const auto& beam=in.seeds[s];
        if(beam.count<3||beam.count>MaxVertices||beam.depth<0||beam.depth>in.maxDepth
           ||beam.facet<0||beam.facet>=in.faceCount)return false;
        for(int d=0;d<=beam.depth;++d)if(beam.track[d]>=in.faceCount)return false;
    }
    return true;
}
inline void ExpandPerOrientation(int faces,int depth,bool deep,size_t& queue,size_t& output){
    if(deep||faces>16){
        const size_t complexity=size_t(faces)*size_t(std::max(1,depth));
        queue=std::max(queue,std::min(size_t(8192),complexity*16));
        output=std::max(output,std::min(size_t(32768),complexity*32));
    }
}
inline void PerOrientation(const std::vector<Input>& inputs,bool deep,size_t& queue,size_t& output){
    queue=512;output=2048;
    for(const auto& in:inputs)ExpandPerOrientation(in.faceCount,in.maxDepth,deep,queue,output);
}
struct Plan {
    size_t queue=0,output=0,inputBytes=0,statusBytes=0,scratchBytes=0,bytes=0;
    int stackCapacity=StackCapacity;
};
inline bool MakePlan(size_t count,bool levels,size_t queue,size_t outputWanted,
                     size_t outputMinimum,size_t budget,size_t queueItemBytes,Plan& result,int stackCapacity=StackCapacity){
    if(!count||count>size_t(INT_MAX)-63||queue>size_t(UINT_MAX)-63||outputWanted>size_t(UINT_MAX)
       ||outputMinimum>outputWanted||stackCapacity<=0||stackCapacity>StackCapacity)return false;
    Plan p;p.queue=queue;p.stackCapacity=stackCapacity;
    size_t fixed=0,temporary=0;
    if(!Multiply(count,sizeof(Input),p.inputBytes)||!Multiply(count,sizeof(Status),p.statusBytes)
       ||!Add(p.inputBytes,p.statusBytes,fixed)||!Add(fixed,sizeof(unsigned)*(levels?3:1),fixed))return false;
    if(levels){
        if(!Multiply(queue,queueItemBytes,temporary)||!Multiply(temporary,2,p.scratchBytes))return false;
    }else if(!Multiply(count,size_t(stackCapacity),temporary)||!Multiply(temporary,sizeof(Beam),p.scratchBytes))return false;
    if(!Add(fixed,p.scratchBytes,fixed)||fixed>budget)return false;
    p.output=std::min(outputWanted,(budget-fixed)/sizeof(Output));
    if(p.output<outputMinimum||!Multiply(p.output,sizeof(Output),temporary)||!Add(fixed,temporary,p.bytes))return false;
    result=p;return true;
}
inline bool Retry(const Plan& previous,size_t count,bool levels,bool queueOverflow,
                  bool outputOverflow,bool otherFailure,unsigned attempt,size_t budget,
                  size_t queueItemBytes,Plan& next){
    if(attempt>=2||otherFailure||(!queueOverflow&&!outputOverflow)||(!levels&&queueOverflow))return false;
    size_t queue=previous.queue,output=previous.output;
    if(queueOverflow&&!Multiply(queue,2,queue))return false;
    if(outputOverflow&&!Multiply(output,2,output))return false;
    // Keep the other capacity: trading it away would merely move the failure.
    return MakePlan(count,levels,queue,output,output,budget,queueItemBytes,next,previous.stackCapacity);
}
inline int FittingChunk(int requested,int faces,int depth,bool levels,size_t budget,
                        size_t queueItemBytes,size_t& reservedBytes,bool wide=false){
    reservedBytes=0;
    if(requested<1||faces<1||faces>MaxFacets||depth<0||depth>MaxDepth)return 0;
    size_t queue=512,output=2048;ExpandPerOrientation(faces,depth,true,queue,output);
    // Two retries can double both capacities twice. No clipping of output is
    // permitted in this reservation; every intermediate retry must fit too.
    if((levels&&!Multiply(queue,4,queue))||!Multiply(output,4,output))return 0;
    for(int count=std::min(wide?512:128,requested);count>0;--count){
        size_t totalQueue=0,totalOutput=0;Plan plan;
        if(Multiply(queue,size_t(count),totalQueue)&&Multiply(output,size_t(count),totalOutput)
           &&MakePlan(size_t(count),levels,totalQueue,totalOutput,totalOutput,budget,queueItemBytes,plan)){
            reservedBytes=plan.bytes;return count;}
    }
    return 0;
}
struct Hint {
    std::string key;
    size_t queuePerOrientation=0,outputPerOrientation=0;
    static const size_t MaximumKeyBytes=16384;
    void Reset(){std::string().swap(key);queuePerOrientation=outputPerOrientation=0;}
    bool Select(const std::string& recipe,int device,const std::string& schedule,bool levels){
        if(recipe.empty()||recipe.size()>MaximumKeyBytes-128||schedule.size()>32){Reset();return false;}
        const std::string next=std::to_string(device)+" "+schedule+" "+(levels?"1":"0")+"\n"+recipe;
        if(next!=key){std::string(next).swap(key);queuePerOrientation=outputPerOrientation=0;}
        return true;
    }
    void Learn(const Plan& success,size_t count){
        if(key.empty()||!count)return;
        queuePerOrientation=std::max(queuePerOrientation,success.queue/count+(success.queue%count!=0));
        outputPerOrientation=std::max(outputPerOrientation,success.output/count+(success.output%count!=0));
    }
    bool Apply(const Plan& baseline,size_t count,bool levels,size_t budget,size_t queueItemBytes,Plan& chosen)const{
        if(key.empty()||!count)return false;
        size_t q=0,o=0,qMaximum=0,oMaximum=0;
        if(!Multiply(queuePerOrientation,count,q)||!Multiply(outputPerOrientation,count,o)
           ||!Multiply(baseline.queue,4,qMaximum)||!Multiply(baseline.output,4,oMaximum))return false;
        q=std::max(baseline.queue,std::min(q,qMaximum));o=std::max(baseline.output,std::min(o,oMaximum));
        if(q==baseline.queue&&o==baseline.output)return false;
        return MakePlan(count,levels,q,o,o,budget,queueItemBytes,chosen,baseline.stackCapacity);
    }
};
struct HintAttempt {
    Hint* hint;bool successful=false;
    explicit HintAttempt(Hint* value):hint(value){}
    ~HintAttempt(){if(hint&&!successful)hint->Reset();}
};
inline bool BufferPoolEnabled(){
    const char* value=std::getenv("MBS_GPU_TRACE_BUFFER_POOL");
    return value&&value[0]=='1'&&value[1]=='\0';
}
inline size_t RetentionLimit(bool queried,size_t free,size_t owned,size_t requested){
    size_t available=0;
    if(!queried||free<512ULL*1024*1024||!Add(free,owned,available))return 0;
    return std::min(std::min(requested,size_t(2ULL*1024*1024*1024)),available/4);
}
using SlotBytes=std::array<size_t,8>;
inline bool PhysicalSlots(const SlotBytes& required,const SlotBytes& existing,size_t budget,
                          bool reuseLarger,SlotBytes& physical,size_t& bytes){
    size_t exact=0,retained=0;
    SlotBytes choice{};
    for(size_t i=0;i<required.size();++i){
        if(!Add(exact,required[i],exact))return false;
        choice[i]=required[i]?(reuseLarger?std::max(required[i],existing[i]):required[i]):0;
        if(!Add(retained,choice[i],retained))return false;
    }
    if(exact>budget)return false;
    // Free excess capacities before any new allocation. Never discard an
    // active slot merely to claim a smaller physical-memory reservation.
    if(retained>budget){physical=required;bytes=exact;}
    else{physical=choice;bytes=retained;}
    return true;
}
static_assert(TrackCapacity>MaxDepth,"track storage must include the final child");
static_assert(StackCapacity>=MaxFacets+MaxDepth*(MaxFacets-2),"DFS stack bound");
}}
