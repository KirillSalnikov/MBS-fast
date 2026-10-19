#include "../src/FullAuto.cpp"
using namespace FullAuto;
void Setup(Engine& e,const std::string& root){
    e.root=root;e.signature="fixes-probe";Directory(root);e.exe="/bin/true";e.binaryHash=FileHash(e.exe);
    e.controls=false;e.depthSharing=false;e.gpus={""};e.seeds.assign(production,production+8);
    e.phi=1;e.minPhi=1;e.maxPhi=2;e.pilot=64;e.maxCount=1024;e.target=2;e.initial=256;e.theta={170,180};
}
void Calibration(const std::string& root,bool allUnresolved){
    ArgPP args;RunConfig config;Engine e(args,config);Setup(e,root);
    const double signs[8]={-1,-1,-1,-1,1,1,1,1};
    for(int phi:{1,2})for(int n:{64,128,256,512})for(int seed=0;seed<8;++seed){
        e.phi=phi;Job j;e.SetSampling(j,n,validation[seed],1,e.theta,false);j.angles=e.theta;j.weightHash=Hash("");j.matrix=Matrix(2);
        for(auto& row:j.matrix){row[0]=100+(allUnresolved?(phi==1?240:130):(phi==1?130:10))*signs[seed];row[5]=20;}
        j.identity=Number(phi)+"_"+Number(n)+"_"+Number(seed);j.folder=root+"/fake_"+j.identity;j.seconds=n*.001*(phi==1?1:2);e.CacheJob(j);
    }
    e.Calibrate(e.theta);Require(e.phi==2,"phi selection used cheapest unresolved ratio instead of informative score");
    const auto report=Read(root+"/calibration.json");
    Require(report.find(allUnresolved?"\"selection_resolved\":false":"\"selection_resolved\":true")!=std::string::npos,"phi resolution receipt wrong");
}
void BoundedWork(const std::string& root){
    ArgPP args;RunConfig config;Engine e(args,config);Setup(e,root);e.maxCount=65536;e.local=true;e.incremental=true;e.orientationBlock=256;
    e.theta.clear();std::vector<size_t> ids;for(size_t i=0;i<32;++i){e.theta.push_back(170+i*.25);ids.push_back(i);}
    Job j;e.SetSampling(j,256,9999,2,e.theta,false);j.identity="timing-profile";j.folder=root+"/profile";j.angles=e.theta;
    j.seconds=100;j.traceSeconds=1;j.diffractionSeconds=99;e.CacheJob(j);
    Budgets current(1,std::vector<int>(32,256)),proposed(1,std::vector<int>(32,16384));
    Statistics stats;stats.interval.resize(32);stats.ratioInterval.resize(32);for(size_t i=0;i<32;++i)stats.interval[i][0]=.1+i*.01;
    AdditionalWork work(e,e.seeds.data());const auto original=proposed;
    BoundRoundWork(e,work,current,proposed,{2},ids,stats);
    Require(proposed!=original&&proposed!=current,"long round was not split with progress");
    for(size_t i=0;i<32;++i)Require(proposed[0][i]>=current[0][i]&&proposed[0][i]<=original[0][i]&&(proposed[0][i]&(proposed[0][i]-1))==0,"split damaged prefix counts");
    Require(Read(root+"/round_work_limit.json").find("\"applied\":true")!=std::string::npos,"round split receipt missing");
}
void SharedGpuResources(const std::string& root){
    ArgPP args;RunConfig config;config.useGpu=true;Engine e(args,config);
    e.gpus={"0"};e.hostBudgetMb=4096;e.dispatchThreads=6;
    const char* names[]={"MBS_GPU_MEM_FRACTION","MBS_HOST_MEM_BUDGET_MB"};
    std::vector<std::pair<bool,std::string>> saved;for(auto name:names){const char* value=getenv(name);saved.push_back({value!=nullptr,value?value:""});unsetenv(name);}
    e.dispatchGpus={"0","0"};e.ValidateDispatchResources();
    Require(e.WorkersOnGpu("0")==2&&std::abs(e.WorkerGpuFraction("0")-.44)<1e-12,"shared VRAM budget incorrect");
    const std::string shared=e.WorkerIdentity("0");e.dispatchGpus={"0"};Require(shared!=e.WorkerIdentity("0"),"resource changes retained old worker");
    e.dispatchGpus={"0","0"};
    auto rejected=[&](){try{e.ValidateDispatchResources();return false;}catch(const std::exception&){return true;}};
    e.dispatchGpus={"0","0","0"};Require(rejected(),"three workers per GPU accepted");e.dispatchGpus={"0","0"};
    e.dispatchThreads=0;Require(rejected(),"shared GPU accepted without CPU pool bound");e.dispatchThreads=6;
    for(auto value:{"0.88","nan","0.1","-1","bad"}){setenv(names[0],value,1);Require(rejected(),"unsafe explicit VRAM budget accepted");}
    setenv(names[0],"0.3",1);e.ValidateDispatchResources();e.ConfigureWorkerEnvironment("0");
    Require(std::string(getenv(names[0]))=="0.3","explicit VRAM override changed");
    unsetenv(names[0]);unsetenv(names[1]);e.ConfigureWorkerEnvironment("0");
    Require(std::abs(std::stod(getenv(names[0]))-.44)<1e-12&&std::string(getenv(names[1]))=="4096","shared child environment incorrect");
    setenv(names[1],"8192",1);Require(rejected(),"unsafe explicit RAM budget accepted");setenv(names[1],"2048",1);e.ValidateDispatchResources();
    for(size_t i=0;i<2;++i){if(saved[i].first)setenv(names[i],saved[i].second.c_str(),1);else unsetenv(names[i]);}
    unsetenv("CUDA_VISIBLE_DEVICES");
    unsetenv("MBS_HOST_MEM_BUDGET_MB");unsetenv("MBS_GPU_MEM_FRACTION");
    e.root=root;Directory(root);e.signature="dispatch-budget-test";e.gpuCandidates={"0"};e.hostBudgetMb=40000;e.SaveTuning();
    const std::string tuning=Read(root+"/execution_tuning.tsv");e.hostBudgetMb=512;e.dispatchGpus={"0","0"};
    Require(e.RestoreTuning()&&e.hostBudgetMb==512,"checkpoint restored stale per-worker RAM budget");
    Require(Read(root+"/execution_tuning.tsv")==tuning,"dispatch budget changed saved sampling recipe");
    for(size_t i=0;i<2;++i){if(saved[i].first)setenv(names[i],saved[i].second.c_str(),1);else unsetenv(names[i]);}

}
int main(int argc,char** argv){try{
    Require(argc==2,"output directory required");
    using namespace FullAutoScheduling;
    Require(!Resolved(1e100)&&Resolved(.2),"unbounded interval ranked as resolved");
    Require(PilotLimit(512,1<<20,false)==4096&&PilotLimit(512,1<<20,true)==512,"pilot bounds/override changed");
    Require(ShortPacket(9,32,4096,256,false)&&ShortPacket(147,1,4096,0,false),"real short packets excluded");
    Require(!ShortPacket(9,32,4096,256,true)&&!ShortPacket(6318,64,8192,0,false),"explicit or wide packet modified");
    Require(!ShortPacket(4,32,4096,512,false,512)&&ShortPacket(4,32,4096,0,false,512),"float cap or existing cap ignored");
    Require(!RetryPhi(1024,2048,1)&&RetryPhi(1024,4096,1)&&!RetryPhi(1024,4096,3),"phi retries not bounded");
    {ArgPP args;RunConfig config;Engine cpu(args,config);Require(cpu.MergeInterpolationWork(),"CPU grouping optimization disabled");
     config.useGpu=true;Engine gpu(args,config);
#ifdef MBS_GPU_FP32
     Require(!gpu.MergeInterpolationWork(),"unvalidated FP32 grouping enabled");
#else
     Require(gpu.MergeInterpolationWork(),"double grouping optimization disabled");
#endif
     gpu.go=true;Require(gpu.MergeInterpolationWork(),"GO grouping optimization disabled");}
    SharedGpuResources(std::string(argv[1])+"/dispatch_resources");
    BoundedWork(std::string(argv[1])+"/round_work");
    Calibration(std::string(argv[1])+"/resolved",false);Calibration(std::string(argv[1])+"/unresolved",true);
    std::cout<<"fullauto scheduling regression checks passed\n";return 0;
}catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}}
