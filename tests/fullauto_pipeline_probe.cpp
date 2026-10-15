// Real CUDA/CPU workers on a frozen allocation. The historical profile fixture
// supplies timing evidence only; its optical matrices are never reused.
#include "../src/FullAuto.cpp"
int main(int argc,const char* argv[]){try{
    using namespace FullAuto;
    ArgPP args;RegisterCliOptions(args);args.Parse(argc,argv);
#ifdef USE_CUDA
    const bool cuda=true;
#else
    const bool cuda=false;
#endif
    RunConfig config=RunConfig::FromCommandLine(args,false,cuda);Engine e(args,config);e.Init(argc,argv);
    const char* binary=getenv("TASK_MBS_WORKER");Require(binary,"TASK_MBS_WORKER required");e.exe=binary;e.binaryHash=FileHash(binary);
    e.workerThreads=config.threads>0?config.threads:1;e.prepareThreadCap=8;e.phi=args.GetIntValue("nphi");
    e.kernel.pipeline=args.IsCatched("orientation_pipeline");e.SetMirror(false);
    WriteGrid(e.root+"/requested_theta.csv",e.theta);Directory(e.root+"/jobs");
    signal(SIGTERM,StopWorkers);signal(SIGINT,StopWorkers);
    if(getenv("TASK_MBS_WARM_CACHE"))e.WarmCache();
    if(const char* fixture=getenv("TASK_MBS_PROFILES")){
        std::ifstream input(fixture);Require(bool(input),"cannot open saved timing profiles");int n,rows,threads,records=0;double trace,diffraction;
        while(input>>n>>rows>>threads>>trace>>diffraction){Job j;j.identity="timing_only_"+Number(records++);j.folder=j.identity;
            j.executionHash=e.ExecutionHash();j.count=n;j.depth=4;j.phi=64;j.seed=0;j.inputAngles.resize(rows);
            j.preparationThreads=threads;j.traceSeconds=trace;j.diffractionSeconds=diffraction;e.CacheJob(j);}
        Require(input.eof()&&records>0,"invalid saved timing profiles");
    }
    const bool legacy=getenv("TASK_MBS_LEGACY_PREPARE")!=nullptr;
    if(legacy)e.resourcePolicy="off"; // Reproduce the old, measured one-thread decision.
    const int n=e.initial,depth=config.maxReflections;const double started=Seconds();
    e.Prefetch({Engine::WorkRequest{n,0,depth,0,e.theta}},e.seeds.data());
    const auto result=e.PrefixStage(n,e.seeds.data(),8,depth,e.theta);
    WriteMatrix(e.root+"/frozen_matrix.dat",e.theta,Mean(e.Data(result)));
    Atomic(e.root+"/frozen_profile.json","{\"seconds\":"+Number(Seconds()-started)+",\"orientations_per_seed\":"+Number(n)
        +",\"gpus\":"+Number(e.gpus.size())+",\"theta_rows\":"+Number(e.theta.size())+",\"legacy_preparation_decision\":"+(legacy?"true":"false")+"}\n");
    std::cout<<Read(e.root+"/frozen_profile.json");return 0;
}catch(const std::exception& error){std::cerr<<error.what()<<std::endl;return 1;}}
