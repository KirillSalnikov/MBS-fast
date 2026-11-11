#include "cuda/GpuTraceCapacity.h"
#include <cassert>
#include <cstdio>
#include <cstdint>

using namespace GpuConvexQueue;
struct QueueElement {Beam beam;int orientation;};
static Input input(int faces,int depth){
    Input in={};in.faceCount=faces;in.seedCount=faces;in.maxDepth=depth;in.ri=1.3116;
    for(int f=0;f<faces;++f){in.faces[f].count=4;in.seeds[f].count=4;
        in.seeds[f].facet=f;in.seeds[f].track[0]=f;}
    return in;
}
int main(){
    const size_t mib=1024*1024,budget=2048*mib,qBytes=sizeof(QueueElement);
    assert(!setenv("MBS_GPU_TRACE_DEEP_CAPACITY","0",1));assert(!Capacity::DeepEnabled());
    assert(!setenv("MBS_GPU_TRACE_DEEP_CAPACITY","1",1));assert(Capacity::DeepEnabled());
    assert(!setenv("MBS_GPU_TRACE_DEEP_CAPACITY","10",1));assert(!Capacity::DeepEnabled());
    assert(!setenv("MBS_GPU_TRACE_WIDE_CHUNKS","1",1));assert(!Capacity::WideChunksEnabled());
    assert(!setenv("MBS_GPU_TRACE_DEEP_CAPACITY","1",1));assert(Capacity::WideChunksEnabled());
    for(const char* value: {"0","true","10",""}){
        assert(!setenv("MBS_GPU_TRACE_WIDE_CHUNKS",value,1));assert(!Capacity::WideChunksEnabled());
    }
    assert(!unsetenv("MBS_GPU_TRACE_WIDE_CHUNKS"));assert(!Capacity::WideChunksEnabled());
    auto prism=input(8,18);assert(Capacity::ValidInput(prism,false));
    prism.maxDepth=22;assert(!Capacity::ValidInput(prism,false)&&Capacity::ValidInput(prism,true));
    prism.maxDepth=23;assert(!Capacity::ValidInput(prism,true));prism.maxDepth=22;
    auto invalid=prism;invalid.seedCount=33;assert(!Capacity::ValidInput(invalid,true));
    invalid=prism;invalid.seedCount=9;assert(!Capacity::ValidInput(invalid,true));
    invalid=prism;invalid.faceCount=33;assert(!Capacity::ValidInput(invalid,true));
    invalid=prism;invalid.seeds[0].depth=24;assert(!Capacity::ValidInput(invalid,true));
    invalid=prism;invalid.seeds[0].depth=-1;assert(!Capacity::ValidInput(invalid,true));
    invalid=prism;invalid.seeds[0].facet=8;assert(!Capacity::ValidInput(invalid,true));
    invalid=prism;invalid.seeds[0].depth=22;invalid.seeds[0].track[22]=7;
    assert(Capacity::ValidInput(invalid,true));invalid.seeds[0].track[22]=8;assert(!Capacity::ValidInput(invalid,true));
    invalid=prism;invalid.faces[0].count=25;assert(!Capacity::ValidInput(invalid,true));
    invalid=prism;invalid.ri=std::numeric_limits<double>::infinity();assert(!Capacity::ValidInput(invalid,true));
    Input empty={};assert(Capacity::ValidInput(empty,true));empty.seedCount=1;assert(!Capacity::ValidInput(empty,true));

    size_t queue=0,output=0;
    std::vector<Input> inputs(128,input(8,18));
    Capacity::PerOrientation(inputs,false,queue,output);assert(queue==512&&output==2048);
    Capacity::PerOrientation(inputs,true,queue,output);assert(queue==2304&&output==4608);
    Capacity::Plan first,retry,second;
    assert(Capacity::MakePlan(inputs.size(),true,queue*128,output*128,128*64,budget,qBytes,first));
    assert(first.bytes<=budget&&first.output==output*128);
    assert(Capacity::Retry(first,128,true,true,false,false,0,budget,qBytes,retry));
    assert(retry.queue==2*first.queue&&retry.output==first.output&&retry.bytes<=budget);
    assert(!Capacity::Retry(first,128,true,true,false,true,0,budget,qBytes,retry));
    assert(!Capacity::Retry(first,128,true,false,false,false,0,budget,qBytes,retry));
    assert(!Capacity::Retry(first,128,false,true,false,false,0,budget,qBytes,retry));
    assert(!Capacity::Retry(first,128,true,true,true,false,2,budget,qBytes,retry));
    assert(Capacity::Retry(first,128,true,false,true,false,0,budget,qBytes,retry));
    assert(retry.queue==first.queue&&retry.output==2*first.output);
    // One fixed tight budget: a retry must not silently expand it or borrow
    // output space to conceal a queue overflow.
    assert(!Capacity::Retry(first,128,true,true,false,false,0,first.bytes,qBytes,retry));
    assert(Capacity::MakePlan(1,true,64,64,64,8*mib,qBytes,first));
    assert(Capacity::Retry(first,1,true,true,true,false,0,8*mib,qBytes,retry));
    assert(Capacity::Retry(retry,1,true,true,true,false,1,8*mib,qBytes,second));
    assert(!Capacity::Retry(second,1,true,true,true,false,2,8*mib,qBytes,retry));

    size_t result=123;const size_t maximum=std::numeric_limits<size_t>::max();
    assert(!Capacity::Multiply(maximum,2,result)&&result==123);
    assert(!Capacity::Add(maximum,1,result)&&result==123);
    assert(!Capacity::MakePlan(maximum,true,64,64,1,budget,qBytes,first));
    assert(!Capacity::MakePlan(1,true,size_t(UINT_MAX)+1,64,1,budget,qBytes,first));
    assert(!Capacity::MakePlan(1,true,size_t(UINT_MAX),64,1,budget,qBytes,first));
    assert(!Capacity::MakePlan(1,true,64,size_t(UINT_MAX)+1,1,budget,qBytes,first));
    assert(!Capacity::MakePlan(1,true,64,64,65,budget,qBytes,first));
    assert(!Capacity::MakePlan(1,true,64,64,1,budget,maximum,first));
    assert(!Capacity::MakePlan(1,true,64,64,1,1,qBytes,first));
    assert(Capacity::MakePlan(1,false,512,2048,64,budget,qBytes,first,LegacyStackCapacity));
    assert(first.scratchBytes==LegacyStackCapacity*sizeof(Beam));
    assert(Capacity::MakePlan(1,false,512,2048,64,budget,qBytes,first));
    assert(first.scratchBytes==StackCapacity*sizeof(Beam));
    size_t reserved=999;
    assert(Capacity::FittingChunk(1024,8,22,true,budget,qBytes,reserved)>0&&reserved<=budget);
    assert(Capacity::FittingChunk(7,8,22,true,budget,qBytes,reserved)==7);
    assert(Capacity::FittingChunk(INT_MAX,8,6,true,16*budget,qBytes,reserved)==128);
    assert(Capacity::FittingChunk(INT_MAX,8,6,true,16*budget,qBytes,reserved,true)==512);
    assert(Capacity::FittingChunk(73,8,6,true,16*budget,qBytes,reserved,true)==73);
    assert(!Capacity::FittingChunk(512,8,22,true,1,qBytes,reserved,true)&&reserved==0);
    assert(!Capacity::FittingChunk(512,8,22,true,budget,maximum,reserved,true));
    assert(!Capacity::FittingChunk(128,8,22,true,1,qBytes,reserved)&&reserved==0);
    assert(!Capacity::FittingChunk(0,8,22,true,budget,qBytes,reserved));
    assert(!Capacity::FittingChunk(128,33,22,true,budget,qBytes,reserved));
    assert(!Capacity::FittingChunk(128,8,23,true,budget,qBytes,reserved));
    assert(!Capacity::FittingChunk(128,8,22,true,budget,maximum,reserved));

    Capacity::Hint hint;
    assert(!hint.Select("",0,"legacy",true));
    assert(hint.Select("shape-A ri1.31 depth22 cutoffoff",0,"legacy",true));
    assert(Capacity::MakePlan(16,true,2816*16,5632*16,5632*16,budget,qBytes,first));
    assert(Capacity::Retry(first,16,true,false,true,false,0,budget,qBytes,retry));
    {Capacity::HintAttempt attempt(&hint);hint.Learn(retry,16);attempt.successful=true;}
    assert(hint.outputPerOrientation==11264&&hint.queuePerOrientation==2816);
    assert(Capacity::MakePlan(115,true,2816*115,5632*115,5632*115,budget,qBytes,first));
    assert(hint.Apply(first,115,true,budget,qBytes,retry));
    assert(retry.output==11264*115&&retry.queue==first.queue&&retry.bytes<=budget);
    assert(!hint.Apply(first,115,true,first.bytes,qBytes,retry)); // baseline still fits
    assert(hint.Select("shape-A ri1.31 depth22 cutoffoff",0,"legacy",true)&&hint.outputPerOrientation==11264);
    {Capacity::HintAttempt failed(&hint);}
    assert(hint.key.empty()&&!hint.outputPerOrientation); // exceptions/fatal errors cannot teach
    for(const auto& key: {"new shape","new material","new depth","new cutoff"}){
        assert(hint.Select(key,0,"legacy",true));hint.Learn(first,115);
        assert(hint.Select(std::string(key)+" changed",0,"legacy",true)&&!hint.outputPerOrientation);
    }
    hint.Select("key",0,"legacy",true);hint.Learn(first,115);
    assert(hint.Select("key",1,"legacy",true)&&!hint.outputPerOrientation);
    hint.Learn(first,115);assert(hint.Select("key",1,"split8",true)&&!hint.outputPerOrientation);
    hint.Learn(first,115);assert(hint.Select("key",1,"split8",false)&&!hint.outputPerOrientation);
    assert(!hint.Select(std::string(20000,'x'),0,"legacy",true)&&hint.key.empty());
    hint.Select("bounded",0,"legacy",true);hint.queuePerOrientation=maximum;hint.outputPerOrientation=maximum;
    assert(!hint.Apply(first,115,true,budget,qBytes,retry));
    hint.queuePerOrientation=2816*64;hint.outputPerOrientation=5632*64;
    assert(hint.Apply(first,115,true,16*budget,qBytes,retry));
    assert(retry.queue==first.queue*4&&retry.output==first.output*4); // original absolute ceiling

    // Failure merge must be commutative/associative regardless of warp and
    // atomic arrival order. Capacity 2/3 can never hide geometry 1 or NaN 4.
    for(int a=0;a<=4;++a)for(int b=0;b<=4;++b){
        assert(CombineFailure(a,b)==CombineFailure(b,a));
        for(int c=0;c<=4;++c)
            assert(CombineFailure(CombineFailure(a,b),c)==CombineFailure(a,CombineFailure(b,c)));
    }
    assert(CombineFailure(1,3)==1&&CombineFailure(3,4)==4);

    // Many small and large plans, including insufficient budgets. Account
    // every allocated byte, not merely the output/queue estimates.
    for(size_t count: {size_t(1),size_t(64),size_t(128),size_t(129),size_t(512),size_t(1024)})
    for(int faces: {4,8,20,32})for(int depth: {6,18,20,22})for(bool levels: {false,true}){
        std::vector<Input> batch(1,input(faces,depth));Capacity::PerOrientation(batch,true,queue,output);
        for(size_t limit: {16*mib,128*mib,1024*mib,4096*mib})for(bool wide: {false,true}){
            size_t reserved=0;const int fitting=Capacity::FittingChunk(int(count),faces,depth,levels,limit,qBytes,reserved,wide);
            assert(fitting>=0&&fitting<=int(count)&&fitting<=(wide?512:128)&&reserved<=limit);
            if(fitting){
                Capacity::Plan base,one,two;
                assert(Capacity::MakePlan(size_t(fitting),levels,queue*fitting,output*fitting,output*fitting,limit,qBytes,base));
                assert(Capacity::Retry(base,size_t(fitting),levels,levels,true,false,0,limit,qBytes,one));
                assert(Capacity::Retry(one,size_t(fitting),levels,levels,true,false,1,limit,qBytes,two));
                assert(two.bytes==reserved&&two.bytes<=limit);
                if(fitting<std::min(wide?512:128,int(count)))
                    assert(!Capacity::MakePlan(size_t(fitting+1),levels,queue*(fitting+1)*(levels?4:1),
                        output*(fitting+1)*4,output*(fitting+1)*4,limit,qBytes,base));
            }
            Capacity::Plan plan;
            if(!Capacity::MakePlan(count,levels,queue*count,output*count,64*count,limit,qBytes,plan))continue;
            const size_t exact=count*sizeof(Input)+count*sizeof(Status)+sizeof(unsigned)*(levels?3:1)
                +plan.scratchBytes+plan.output*sizeof(Output);
            assert(exact==plan.bytes&&exact<=limit&&plan.output>=64*count);
            if(Capacity::Retry(plan,count,levels,levels,true,false,0,limit,qBytes,retry))assert(retry.bytes<=limit);
        }
    }
    std::printf("PASS deep-capacity opt-in, depth22/track bounds, seed guards, byte budgets, bounded retries, sticky failures; sizeof Beam=%zu Output=%zu\n",sizeof(Beam),sizeof(Output));
}
