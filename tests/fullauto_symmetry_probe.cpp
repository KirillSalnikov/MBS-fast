// Exercise the actual controller's geometry/numerical rotation selection without
// running the unrelated production convergence loop. Uses a real MBS worker.
#include "../src/FullAuto.cpp"
#include "../src/ArgPP.h"
#include "../src/CliOptions.h"
#include "../src/RunConfig.h"
#include <iostream>
int main(int argc,const char*argv[]){
    try{
        ArgPP args;RegisterCliOptions(args);args.Parse(argc,argv);
#ifdef USE_CUDA
        const bool cuda=true;
#else
        const bool cuda=false;
#endif
        RunConfig config=RunConfig::FromCommandLine(args,false,cuda);
        FullAuto::Engine e(args,config);e.Init(argc,argv);
        const char* worker=getenv("TASK_MBS_WORKER");FullAuto::Require(worker,"TASK_MBS_WORKER required");
        e.exe=worker;e.binaryHash=FullAuto::FileHash(e.exe);
        FullAuto::WriteGrid(e.root+"/requested_theta.csv",e.theta);FullAuto::Directory(e.root+"/jobs");
        signal(SIGTERM,FullAuto::StopWorkers);signal(SIGINT,FullAuto::StopWorkers);
        e.controls=false;e.phi=args.GetIntValue("nphi");e.TuneRotations(e.theta);
        const int meanThreads=e.analyticMeanThreads;e.SaveTuning();e.analyticMeanThreads=0;
        FullAuto::Require(e.RestoreTuning()&&e.analyticMeanThreads==meanThreads,"mean thread policy was not restored");
        e.State("symmetry_test_complete",0,{128},e.theta.size(),0,0,0);
        std::cout<<FullAuto::Read(e.root+"/rotation_verification.json")<<std::endl;return 0;
    }catch(const std::exception& error){std::cerr<<error.what()<<std::endl;return 1;}
}
