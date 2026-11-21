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
    Require(e.phaseNativeCommitted==8&&e.phaseAngleOrientationSamples==8*64*4*4,"phase work lost or double counted");
    Require(e.recentNativeBlocks.size()==8,"missing durable angular progress");
    for(const auto& block:e.recentNativeBlocks)Require(Exists(block.folder+"/complete.txt")&&block.rows==4&&block.lo==170&&block.hi==180,"reported block was not committed");
    e.SetExecutionPhase("test_complete",4);
    Require(e.phaseNativeCommitted==0&&e.phaseAngleOrientationSamples==0&&e.recentNativeBlocks.empty(),"new phase retained old completion counters");
    {
        Engine broken(args,config);
        SetupPublication(broken,e.root+"/failure",argv[2],false);
        broken.workerPolicy=argv[3];broken.balancedCorrection=false;
        broken.base.push_back("--invalid-test-worker-option");
        bool failed=false;
        try{broken.DirectStage(64,broken.seeds.data(),8,8,broken.theta,true);}
        catch(const std::exception& error){failed=std::string(error.what()).find("worker failed")!=std::string::npos;}
        Require(failed&&broken.nativeAbandoned>0&&broken.nativeAbandoned==broken.nativeSubmitted,
            "failed native batch still active or counted as returned");
        const auto snapshot=Read(broken.root+"/execution_progress.json");
        Require(broken.nativeCommitted==0&&broken.liveNative.empty()&&
            snapshot.find("\"active_workers\":[]")!=std::string::npos&&
            snapshot.find("\"batch_active\":false")!=std::string::npos,
            "failure left stale published worker or dispatch");
    }
    if(e.workerPolicy=="persistent"){
        std::vector<Job> reference;
        for(bool preflight:{false,true}){
            Engine rosette(args,config);
            SetupPublication(rosette,e.root+(preflight?"/preflight":"/fallback"),argv[2],false);
            rosette.base={"--method","po","--backend","cpu","--particle","3","1","1","--rosette-joined",
                "--ri","1.3116","0","--wavelength-um",".532","--symmetry","2","2","--cutoff-profile","off"};
            rosette.nonconvexGeometry=preflight;
            rosette.depthSharing=rosette.DepthSharingRequested();
            rosette.target=4;rosette.pilot=8;rosette.balancedCorrection=false;
            rosette.VerifyDepthSharing(rosette.theta);
            Require(!rosette.depthSharing&&rosette.nativeSubmitted==(preflight?0:1),
                "nonconvex depth capability did not skip exactly the unsupported worker request");
            const auto jobs=rosette.DirectStage(8,rosette.seeds.data(),2,4,rosette.theta);
            if(!preflight){reference=jobs;continue;}
            for(size_t index=0;index<jobs.size();++index)
                Require(jobs[index].matrix==reference.at(index).matrix,"preflight changed rosette matrix");
        }
    }
    std::ostringstream proof;proof<<std::setprecision(17);
    for(const auto& job:result){std::set<std::string> verified;ValidateStoredJob(job.folder,verified);
        for(const auto& row:job.matrix)for(double v:row)proof<<v<<' ';}
    Atomic(e.root+"/proof.txt",proof.str());
    std::cout<<"PASS completed workers, durable counters, no stale active dispatch\n";return 0;
}catch(const std::exception& error){std::cerr<<error.what()<<'\n';return 1;}}
