#define FULLAUTO_STALL_HELPERS_ONLY
#include "fullauto_stall_probe.cpp"
int main(int argc,char** argv){try{
    Require(argc==2,"OUTPUT");ArgPP args;RunConfig config;Engine e(args,config);SetupMemory(e,argv[1],"/bin/true",false);
    e.cutoffPolicy="staged";e.theta={170,175,180};e.orientationBlock=4096;e.balancedCorrection=true;
    auto timing=[&](const std::string& suffix,int count,double seconds){Job j;e.SetSampling(j,count,9999,14,e.theta,false);
        j.executionHash=Hash(e.ExecutionHash()+suffix);j.identity=Hash(suffix);j.folder=e.root+"/"+j.identity;j.traceSeconds=seconds;j.seconds=seconds;e.CacheJob(j);};
    timing(" staged-fast-child-v1",4096,4096);timing(" off-child-v1",64,256);
    AdditionalWork incomplete(e,e.seeds.data());Require(!incomplete.Forecast({{4096,0,14,0,e.theta}}).profiled,"missing balanced stage treated as free");
    timing(" staged-balanced-child-v1",512,1024);AdditionalWork complete(e,e.seeds.data());
    const auto a=complete.Forecast({{4096,0,14,0,e.theta}});Require(a.profiled&&a.jobs==40&&std::fabs(a.seconds-8*6016)<1e-9,"staged cost omitted a term");
    const auto b=complete.Forecast({{512,0,14,0,e.theta}});Require(b.profiled&&b.jobs==32&&std::fabs(b.seconds-8*960)<1e-9,"duplicate balanced minimum cost");
    std::cout<<"PASS staged timing: all levels, unknown profiles, minimum-count deduplication\n";return 0;
}catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}}
