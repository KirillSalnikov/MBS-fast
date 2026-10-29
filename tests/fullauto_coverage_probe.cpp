// Forecast-only regression: a zero-byte index exercises the original exact
// candidate search. No optical work or production checkpoint mutation.
#include "../src/FullAuto.cpp"
using namespace FullAuto;
static void Same(const WorkForecast& a,const WorkForecast& b){
    Require(a.seconds==b.seconds&&a.jobs==b.jobs&&a.profiled==b.profiled,"coverage index changed forecast");
}
static void Init(Engine& e){
    e.controls=false;e.mirror=false;e.phi=2;e.angularCache=true;e.workQueuePolicy="on";
    e.incremental=true;e.orientationBlock=64;e.seeds.assign(production,production+8);
    e.theta={170,172,174,176,178,180};
}
static void Add(Engine& e,unsigned seed,int depth,int offset,const std::vector<double>& angles,
                const std::string& id,int paired=0){
    Job j;e.SetSampling(j,64,seed,depth,angles,false);
    j.identity=id;j.folder=id;j.offset=offset;j.pairedDepth=paired;
    j.traceSeconds=1;j.diffractionSeconds=2;j.seconds=3;e.CacheJob(j);
}
int main(){try{
    ArgPP args;RunConfig config;Engine e(args,config);Init(e);
    Add(e,e.seeds[0],6,0,{170,172,174,176},"first");
    Add(e,e.seeds[0],6,0,{174,176,178},"overlap");
    Add(e,e.seeds[0],6,64,{172,174},"offset");
    Add(e,e.seeds[1],6,0,{172,176},"second_seed");
    Add(e,e.seeds[0],8,0,{170,180},"paired",6);
    e.phi=4;Add(e,e.seeds[0],6,0,{170,180},"other_phi");e.phi=2;
    e.mirror=true;Add(e,e.seeds[0],6,0,{172,180},"other_mirror");e.mirror=false;
    e.controls=true;Add(e,e.seeds[0],6,0,{172,180},"other_controls");e.controls=false;
    e.base={"different_recipe"};Add(e,e.seeds[0],6,0,{172,180},"other_recipe");e.base.clear();
    const std::vector<std::vector<Engine::WorkRequest>> requests={
        {{64,0,6,0,{170,180}}},{{64,0,6,0,{174,176}}},
        {{128,0,6,0,{170,172,174,180}}},{{64,0,8,6,{170,180}}},
        {{64,0,6,0,{170,174}},{128,0,6,0,{174,180}}}};
    for(size_t limit:{size_t(0),size_t(7*sizeof(double)-1),size_t(7*sizeof(double)),
                     size_t(7*sizeof(double)+1),size_t(8*1024*1024)}){
        AdditionalWork original(e,e.seeds.data(),false,false,false,0);
        AdditionalWork indexed(e,e.seeds.data(),false,false,false,limit);
        for(const auto& r:requests){Same(original.Forecast(r),indexed.Forecast(r));Require(indexed.CoverageIndexBytes()<=limit,"coverage byte cap exceeded");}
        for(auto r=requests.rbegin();r!=requests.rend();++r)Same(original.Forecast(*r),indexed.Forecast(*r));
        if(limit==7*sizeof(double))Require(indexed.CoverageIndexBytes()==limit,"exact coverage boundary was not used");
    }
    Engine empty(args,config);Init(empty);
    AdditionalWork emptyIndex(empty,empty.seeds.data()),emptyOriginal(empty,empty.seeds.data(),false,false,false,0);
    for(const auto& r:requests)Same(emptyOriginal.Forecast(r),emptyIndex.Forecast(r));
    Require(emptyIndex.CoverageIndexBytes()==0,"empty coverage allocated angles");
    AdditionalWork ignored(e,e.seeds.data(),true),ignoredOriginal(e,e.seeds.data(),true,false,false,0);
    for(const auto& r:requests)Same(ignoredOriginal.Forecast(r),ignored.Forecast(r));
    Require(ignored.CoverageIndexBytes()==0,"ignored coverage allocated angles");
    // Long completed prefixes, a hole, a partial tail and overlapping groups.
    // The zero-budget path cannot skip prefixes and remains an independent
    // exact leaf-by-leaf reference for the scheduling accelerator.
    Engine longRun(args,config);Init(longRun);
    for(unsigned seed:longRun.seeds)for(int offset=0;offset<1024;offset+=64){
        const auto angles=offset==384&&seed==longRun.seeds.back()?std::vector<double>{170,180}:std::vector<double>{170,172,174,180};
        Add(longRun,seed,6,offset,angles,Number(seed)+"_"+Number(offset));
        Add(longRun,seed,8,offset,angles,"paired_"+Number(seed)+"_"+Number(offset),6);
    }
    std::vector<std::vector<Engine::WorkRequest>> longRequests;
    for(int n:{64,1024,256,512,96,2048,384}){
        longRequests.push_back({{n,0,6,0,{170,180}}});
        longRequests.push_back({{n,0,6,0,{172,174}}});
        longRequests.push_back({{n,0,8,6,{170,174}},{256,0,6,0,{170,180}}});
    }
    for(size_t limit:{size_t(0),size_t(127),size_t(128),size_t(256),size_t(8192),size_t(8*1024*1024)}){
        AdditionalWork reference(longRun,longRun.seeds.data(),false,false,false,0),fast(longRun,longRun.seeds.data(),false,false,false,limit);
        for(const auto& request:longRequests){Same(reference.Forecast(request),fast.Forecast(request));Require(fast.CoverageIndexBytes()<=limit,"leading coverage exceeded shared budget");}
    }
    // Each phi child must spend from the same bound, not receive a new budget.
    e.thetaPhi={{0,2},{1,4},{2,2},{3,4},{4,2},{5,4}};
    for(size_t limit:{size_t(1),size_t(56),size_t(80),size_t(1024)}){
        AdditionalWork original(e,e.seeds.data(),false,false,false,0),indexed(e,e.seeds.data(),false,false,false,limit);
        for(const auto& r:requests){Same(original.Forecast(r),indexed.Forecast(r));Require(indexed.CoverageIndexBytes()<=limit,"phi children multiplied coverage cap");}
    }
    std::cout<<"PASS: bounded exact coverage union; byte boundaries, overlap, empty, ignored, incompatible sampling, paired depths and shared phi budget\n";
    return 0;
}catch(const std::exception& error){std::cerr<<error.what()<<'\n';return 1;}}
