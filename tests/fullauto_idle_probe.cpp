// Real saved-block replay: assemble another angular subset without new optics.
#include "../src/FullAuto.cpp"
int main(int argc,const char* argv[]){try{
    using namespace FullAuto;
    if(argc==2&&std::string(argv[1])=="--queue-policy"){
        ArgPP args;RunConfig config;Engine cpu(args,config);cpu.gpus={""};Require(!cpu.QueueEnabled(),"single CPU unexpectedly queued");
        config.useGpu=true;Engine gpu(args,config);gpu.gpus={"0"};Require(gpu.QueueEnabled(),"single GPU auto queue disabled");
        gpu.workQueuePolicy="off";Require(!gpu.QueueEnabled(),"explicit queue off ignored");
        gpu.workQueuePolicy="on";gpu.angularCache=false;Require(!gpu.QueueEnabled(),"queue enabled without angular cache");
        std::cout<<"PASS: single GPU automatic queue, explicit off, CPU fallback\n";return 0;
    }
    ArgPP args;RegisterCliOptions(args);args.Parse(argc,argv);
    RunConfig config=RunConfig::FromCommandLine(args,false,false);Engine e(args,config);e.Init(argc,argv);
    const char* worker=getenv("TASK_MBS_WORKER");Require(worker,"TASK_MBS_WORKER required");e.exe=worker;e.binaryHash=FileHash(e.exe);
    Directory(e.root+"/jobs");WriteGrid(e.root+"/requested_theta.csv",e.theta);
    e.phi=args.GetIntValue("nphi");e.SetMirror(false);e.sharedMeans=e.controls;e.workerPolicy="persistent";e.workQueuePolicy="on";
    if(e.controls){e.weightTheta={e.theta.front(),e.theta.back()};e.weights={{{.7,1.1}},{{1.2,.9}}};}
    if(getenv("TASK_MBS_WARM_CACHE"))e.WarmCache();
    const int n=e.initial,depth=config.maxReflections;
    e.Prefetch({Engine::WorkRequest{n,0,depth,0,e.theta}},e.seeds.data());
    // Use the original full grid only as a numerical reference, never count it
    // as an independent statistical sample.
    auto full=e.PrefixStage(n,e.seeds.data(),8,depth,e.theta);
    std::vector<double> subset;std::vector<size_t> ids;
    for(size_t i=0;i<e.theta.size();++i)if(i%3==0||i+1==e.theta.size()){subset.push_back(e.theta[i]);ids.push_back(i);}
    const size_t before=e.cache.size();const double started=Seconds();
    auto rows=e.PrefixStage(n,e.seeds.data(),8,depth,subset);const double elapsed=Seconds()-started;
    double error=0;for(int s=0;s<8;++s)for(size_t r=0;r<ids.size();++r)for(int j=0;j<16;++j)
        error=std::max(error,std::fabs(rows[s].matrix[r][j]-full[s].matrix[ids[r]][j])/full[s].matrix[ids[r]][0]);
    Require(error<1e-12,"saved prefix selection changed optical matrices");
    // Integral source traversal must remain complete even when intermediate
    // angular views are omitted. All native sources have the same windows.
    IntegralReader reader;for(int s=0;s<8;++s){auto ref=reader.Resolve(full[s].folder,depth),value=reader.Resolve(rows[s].folder,depth);
        Require(std::fabs(value.area-ref.area)<1e-10*ref.area&&std::fabs(value.ext-ref.ext)<1e-10*std::max(1.,std::fabs(ref.ext)),"prefix integral provenance changed");}
    WriteMatrix(e.root+"/frozen_matrix.dat",subset,Mean(e.Data(rows)));
    Atomic(e.root+"/replay_profile.json","{\"seconds\":"+Number(elapsed)+",\"new_derived_jobs\":"+Number(e.cache.size()-before)+",\"all16_difference\":"+Number(error)+",\"blocks_per_seed\":"+Number(n/e.orientationBlock)+",\"theta\":"+Number(subset.size())+",\"integral_sources_complete\":true}\n");
    std::cout<<Read(e.root+"/replay_profile.json");return 0;
}catch(const std::exception& error){std::cerr<<error.what()<<std::endl;return 1;}}
