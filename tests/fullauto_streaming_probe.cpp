// Frozen optical workloads exercise queue overlap without changing sampling.
#include "../src/FullAuto.cpp"
using namespace FullAuto;
void SetupStream(Engine& e,const std::string& root,const std::string& binary,const std::string& mode){
    e.root=root;e.signature="streaming-regression-v1";Directory(root);Directory(root+"/jobs");e.exe=binary;e.binaryHash=FileHash(binary);
    e.controls=false;e.gpus={"",""};e.seeds.assign(production,production+8);e.workerThreads=1;e.hostBudgetMb=512;
    e.phi=1;e.target=6;e.theta={170,175,180};e.initial=9216;e.orientationBlock=512;e.workerPolicy=mode;e.resourcePolicy="off";e.workQueuePolicy="on";
    e.base={"--method","po","--backend","cpu","--particle","1",".2",".1","--ri","1.3116","0","--wavelength-um",".532","--symmetry","1","1","--cutoff-profile","off"};
    e.balancedCorrection=true;e.selectedCutoff="off";WriteGrid(root+"/requested_theta.csv",e.theta);
}
void Preparation(const std::string& root){
    Require(FullAutoPreparation::ThreadCap(12,12,1)==12,"12 cores truncated");
    Require(FullAutoPreparation::ThreadCap(64,16,4)==16,"EPYC partition divided twice");
    Require(FullAutoPreparation::ThreadCap(64,8,4)==8,"affinity ignored");
    Require(FullAutoPreparation::ThreadCap(12,12,2)==6,"two cards oversubscribed");
    ArgPP args;RegisterCliOptions(args);RunConfig config;config.useGpu=true;Engine e(args,config);e.root=root;Directory(root);
    e.workerThreads=4;e.prepareThreadCap=12;e.controls=false;e.phi=1;e.gpus={"0"};
    Job j;j.identity="balanced-small";j.folder=root+"/sample";j.count=512;j.depth=8;j.phi=1;j.inputAngles={170,180};
    j.executionHash=Hash(e.ExecutionHash()+" balanced-child-v1");j.preparationThreads=4;j.traceSeconds=1;j.diffractionSeconds=.01;j.processingSeconds=1;
    for(int i=0;i<4;++i){j.identity="baseline"+Number(i);j.folder=root+"/sample"+Number(i);e.CacheJob(j);}
    Require(e.AutoPrepareThreads(512,8,0,2,j.executionHash)==12,"child profile ignored");
    Require(e.AutoPrepareThreads(512,8,0,2)==0,"balanced profile leaked into off parent");
    std::vector<std::string> command={"--threads","4","--po-prepare-threads","4"};e.ApplyAutomaticPreparation(j,command);
    Require(j.preparationThreads==12&&command.back()=="12","actual child not tuned at dispatch");
    j.preparationThreads=12;j.processingSeconds=1.2;
    for(int i=0;i<4;++i)e.ObservePreparation(j);
    Require(e.AutoPrepareThreads(512,8,0,2,j.executionHash)==0,"slower 12-thread candidate retained");
    j.processingSeconds=.6;for(int i=0;i<8;++i)e.ObservePreparation(j);
    Require(e.AutoPrepareThreads(512,8,0,2,j.executionHash)==12,"faster measured candidate ignored");
    j.traceSeconds=.01;j.diffractionSeconds=2;
    for(int i=0;i<32;++i)e.ObservePreparation(j);
    e.ApplyAutomaticPreparation(j,command);Require(j.preparationThreads==0&&command.size()==2,"fresh child timings ignored");
    for(const char* flag:{"--threads","--po-prepare-threads"}){
        ArgPP manual;RegisterCliOptions(manual);const char* av[]={"probe",flag,"3"};manual.Parse(3,av);Engine fixed(manual,config);fixed.root=root;fixed.prepareThreadCap=12;
        const auto original=command;fixed.ApplyAutomaticPreparation(j,command);Require(command==original,"manual flag was overridden");
    }
    std::cout<<"PASS affinity cap, child-profile isolation/live refresh, manual flags\n";
}
int main(int argc,char** argv){try{
    Require(argc>=2,"OUTPUT [BINARY MODE]");Preparation(std::string(argv[1])+"/preparation");if(argc==2)return 0;
    Require(argc==4,"OUTPUT BINARY MODE");ArgPP args;RunConfig config;Engine e(args,config);SetupStream(e,argv[1],argv[2],argv[3]);
    const bool restore=getenv("TASK_STREAM_RESTORE");if(restore)e.WarmCache();
    const auto before=e.cache.size();const double begin=Seconds();
    e.Prefetch({Engine::WorkRequest{e.initial,0,6,4,e.theta}},e.seeds.data());
    auto result=e.PrefixStage(e.initial,e.seeds.data(),8,6,e.theta);
    WriteMatrix(e.root+"/frozen_matrix.dat",e.theta,Mean(e.Data(result)));
    IntegralReader reader;std::ostringstream ints;ints<<std::setprecision(17);for(const auto& j:result){auto i=reader.Resolve(j.folder,6);ints<<i.area<<' '<<i.ext<<'\n';}Atomic(e.root+"/frozen_integrals.txt",ints.str());
    if(restore)Require(e.cache.size()==before,"resume issued new work");
    Atomic(e.root+"/timing.json","{\"seconds\":"+Number(Seconds()-begin)+",\"cache_jobs\":"+Number(e.cache.size())+"}");
    std::cout<<"PASS frozen streamed allocation and integral provenance\n";return 0;
}catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}}
