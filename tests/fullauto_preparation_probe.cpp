// Regression for the preparation scheduler: historical GPU stalls must not
// pin CPU tracing to one thread, and fresh evidence must replace old profiles.
#include "../src/FullAuto.cpp"
int main(int argc,const char* argv[]){try{
    using namespace FullAuto;
    Require(argc==2||argc==3,"usage: fullauto_preparation_probe OUTPUT [PRODUCTION_PROFILES_TSV]");Directory(argv[1]);
    ArgPP args;RegisterCliOptions(args);RunConfig config;config.useGpu=true;
    Engine e(args,config);e.root=argv[1];e.workerThreads=1;e.prepareThreadCap=8;e.controls=false;e.phi=64;
    if(argc==3){
        std::ifstream input(argv[2]);Require(bool(input),"cannot read production profiles");int n,rows,threads;double trace,diffraction;int records=0;
        while(input>>n>>rows>>threads>>trace>>diffraction){Job j;j.identity="saved"+Number(records++);j.folder=j.identity;
            j.count=n;j.depth=4;j.phi=64;j.executionHash=e.ExecutionHash();j.inputAngles.resize(rows);
            j.preparationThreads=threads;j.traceSeconds=trace;j.diffractionSeconds=diffraction;e.CacheJob(j);}
        Require(input.eof()&&records>0,"invalid/empty production profiles");
        Require(e.AutoPrepareThreads(8192,4,0,26)==8,"real 26-row production profile still limits CPU tracing to one thread");
        Require(e.AutoPrepareThreads(8192,4,0,44)==8,"real 44-row production profile still limits CPU tracing to one thread");
        std::cout<<"PASS: replayed "<<records<<" production profiles; 26-row and 44-row jobs select 8 preparation threads\n";return 0;
    }
    auto sample=[&](const std::string& name,int rows,double trace,double diffraction){
        Job j;j.identity=name;j.folder=name;j.count=8192;j.depth=4;j.phi=64;j.block=64;
        j.executionHash=e.ExecutionHash();j.inputAngles.resize(rows);j.preparationThreads=1;
        j.traceSeconds=trace;j.diffractionSeconds=diffraction;return j;
    };
    // Measured epyc1 production values, seconds per complete native block.
    e.CacheJob(sample("comparable",26,1.890,.693));
    for(int i=0;i<20;++i)e.CacheJob(sample("old_two_row_stall"+Number(i),2,1.743,7.8));
    Require(e.AutoPrepareThreads(8192,4,0,26)==8,"unrelated two-row GPU stalls disabled preparation pool");
    // A transient slow GPU in the same bucket cannot dominate its median.
    for(int i=0;i<5;++i)e.CacheJob(sample("steady"+Number(i),26,1.890,.693));
    e.CacheJob(sample("stall",26,1.890,20));e.preparationProfiles.clear();
    Require(e.AutoPrepareThreads(8192,4,0,26)==8,"one stalled job disabled preparation pool");
    for(int i=0;i<32;++i)e.ObservePreparation(sample("new_gpu_bound",26,1.890,10));
    Require(e.AutoPrepareThreads(8192,4,0,26)==0,"fresh GPU-bound work did not replace restored profile");
    for(int i=0;i<32;++i)e.ObservePreparation(sample("new_cpu_bound",26,1.890,.693));
    Require(e.AutoPrepareThreads(8192,4,0,26)==8,"fresh CPU-bound work did not replace stale profile");
    Require(e.AutoPrepareThreads(512,4,0,26)==0,"tiny pilots changed preparation resources");
    e.resourcePolicy="off";Require(e.AutoPrepareThreads(8192,4,0,26)==0,"explicit resource policy ignored");e.resourcePolicy="auto";
    e.go=true;Require(e.AutoPrepareThreads(8192,4,0,26)==0,"PO preparation applied to GO");e.go=false;
    config.useGpu=false;Require(e.AutoPrepareThreads(8192,4,0,26)==0,"CPU worker resource override");config.useGpu=true;
    for(const char* flag:{"--threads","--po-prepare-threads"}){
        ArgPP fixed;RegisterCliOptions(fixed);const char* av[]={"probe",flag,"2"};fixed.Parse(3,av);
        Engine explicitEngine(fixed,config);explicitEngine.workerThreads=1;explicitEngine.prepareThreadCap=8;
        Require(explicitEngine.AutoPrepareThreads(8192,4,0,26)==0,"explicit thread override ignored");
    }
    for(size_t n:{size_t(1),size_t(2),size_t(3),size_t(26),size_t(8192)}){
        size_t b=FullAutoPreparation::Bucket(n);Require(b<=n&&b>n/2,"invalid comparable-work bucket");
    }
    std::cout<<"PASS: comparable production profiles, outlier resistance, live refresh, explicit overrides, CPU/GO guards\n";
    return 0;
}catch(const std::exception& error){std::cerr<<error.what()<<std::endl;return 1;}}
