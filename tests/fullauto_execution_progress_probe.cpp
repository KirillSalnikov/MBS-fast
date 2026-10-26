#define main native_publication_original_main
#include "fullauto_native_publication_probe.cpp"
#undef main
int main(int argc,char** argv){try{
    Require(argc==4,"OUTPUT WORKER POLICY");ArgPP args;RunConfig config;Engine e(args,config);
    SetupPublication(e,argv[1],argv[2],false);e.workerPolicy=argv[3];e.balancedCorrection=false;
    e.SetExecutionPhase("native_execution",4);
    const auto result=e.DirectStage(64,e.seeds.data(),8,8,e.theta,true);
    Require(result.size()==8&&e.nativeCommitted==8&&e.nativeReturned==8,"durable work counter changed");
    Require(e.nativeSubmitted==8+e.nativeRetries&&e.nativeOrientations==512,"retry counted as new samples");
    Require(e.liveNative.empty()&&e.dispatchProgress.empty(),"stale active workers after completion");
    e.SetExecutionPhase("test_complete",4);
    std::ostringstream proof;proof<<std::setprecision(17);
    for(const auto& job:result){std::set<std::string> verified;ValidateStoredJob(job.folder,verified);
        for(const auto& row:job.matrix)for(double v:row)proof<<v<<' ';}
    Atomic(e.root+"/proof.txt",proof.str());
    std::cout<<"PASS completed workers, durable counters, no stale active dispatch\n";return 0;
}catch(const std::exception& error){std::cerr<<error.what()<<'\n';return 1;}}
