// Test harness for a frozen allocation, independent of timing-based adaptation.
// TASK_MBS_WORKER names the real worker executable; this probe never runs optics.
#include "../src/FullAuto.cpp"
int main(int argc,const char* argv[]){
    try{
        if(argc==3&&std::string(argv[1])=="--test-file-hash"){
            const std::string path=argv[2];FullAuto::Atomic(path,std::string(4*1024*1024,'a'));
            const double begin=FullAuto::Seconds();std::string digest;for(int i=0;i<128;++i)digest=FullAuto::Hash(FullAuto::Read(path));const double cold=FullAuto::Seconds()-begin;
            const double cachedBegin=FullAuto::Seconds();for(int i=0;i<128;++i)FullAuto::Require(FullAuto::FileHash(path)==digest,"digest mismatch");const double cached=FullAuto::Seconds()-cachedBegin;
            struct stat stamp;stat(path.c_str(),&stamp);std::fstream file(path.c_str(),std::ios::in|std::ios::out);file.put('b');file.close();
            struct timespec times[2]={stamp.st_atim,stamp.st_mtim};utimensat(AT_FDCWD,path.c_str(),times,0);
            FullAuto::Require(FullAuto::FileHash(path)!=digest,"same-size mutation with restored mtime escaped validation");
            std::remove(path.c_str());std::cout<<"{\"uncached_seconds\":"<<cold<<",\"cached_seconds\":"<<cached<<",\"mutation_rejected\":true}\n";return 0;
        }
        if(argc==3&&std::string(argv[1])=="--test-depth-promotion"){
            ArgPP args;RunConfig config;FullAuto::Engine e(args,config);
            e.root=argv[2];e.signature="test_physics_checkpoint";e.autoDepth=true;e.target=6;e.depthStart=6;e.depthMaximum=18;e.depthBudget=.0025;
            FullAuto::Directory(e.root);FullAuto::Directory(e.root+"/jobs");
            FullAuto::Atomic(e.root+"/physics_validation.json","{\"depth_pass\":false,\"cutoff_pass\":true}");
            FullAuto::Atomic(e.root+"/adaptive_checkpoint.tsv","allocation_at_n6");
            FullAuto::Atomic(e.root+"/execution_tuning.tsv","tuning_at_n6");
            FullAuto::Require(FullAuto::PromoteAutomaticDepth(e,4)&&e.target==8,"failed to promote automatic depth");
            FullAuto::Require(!FullAuto::Exists(e.root+"/adaptive_checkpoint.tsv")&&FullAuto::Read(e.root+"/depth_attempts/n6/adaptive_checkpoint.tsv")=="allocation_at_n6","old estimator was reused or lost");
            FullAuto::Require(!FullAuto::Exists(e.root+"/execution_tuning.tsv"),"old depth tuning remained live");
            FullAuto::Atomic(e.root+"/physics_validation.json","{\"depth_pass\":false,\"cutoff_pass\":false}");
            FullAuto::Require(!FullAuto::PromoteAutomaticDepth(e,4)&&e.target==8,"unrelated cutoff failure changed depth");
            e.target=18;FullAuto::Require(!FullAuto::PromoteAutomaticDepth(e,4),"automatic depth exceeded ceiling");
            e.target=6;e.autoDepth=false;FullAuto::Require(!FullAuto::PromoteAutomaticDepth(e,4),"explicit target changed");
            std::cout<<"Automatic depth promotion archives old allocations and respects explicit physics and ceiling\n";return 0;
        }
        ArgPP args;RegisterCliOptions(args);args.Parse(argc,argv);
#ifdef USE_CUDA
        const bool cuda=true;
#else
        const bool cuda=false;
#endif
        RunConfig config=RunConfig::FromCommandLine(args,false,cuda);
        FullAuto::Engine e(args,config);e.Init(argc,argv);
        const char* worker=getenv("TASK_MBS_WORKER");FullAuto::Require(worker,"TASK_MBS_WORKER is required");e.exe=worker;e.binaryHash=FullAuto::FileHash(e.exe);
        FullAuto::WriteGrid(e.root+"/requested_theta.csv",e.theta);FullAuto::Directory(e.root+"/jobs");
        signal(SIGTERM,FullAuto::StopWorkers);signal(SIGINT,FullAuto::StopWorkers);
        e.phi=args.GetIntValue("nphi");e.depthSharing=e.coherent&&e.On("fullauto_depth_share",false);e.kernel.pipeline=args.IsCatched("orientation_pipeline");
        if(getenv("TASK_MBS_WARM_CACHE"))e.WarmCache();
        const std::vector<int> depths={1,2,config.maxReflections};FullAuto::Budgets counts(3,std::vector<int>(e.theta.size()));std::vector<size_t> ids;
        for(size_t i=0;i<e.theta.size();++i){ids.push_back(i);counts[0][i]=e.initial;counts[1][i]=i%2?e.initial/2:e.initial/4;counts[2][i]=e.initial/4;}
        const double begin=FullAuto::Seconds();auto result=FullAuto::EvaluateBudget(e,depths,counts,ids,e.seeds.data());
        FullAuto::WriteMatrix(e.root+"/frozen_matrix.dat",e.theta,FullAuto::Mean(FullAuto::Combine(result.levels)));
        FullAuto::Atomic(e.root+"/frozen_seconds.txt",FullAuto::Number(FullAuto::Seconds()-begin)+"\n");return 0;
    }catch(const std::exception& error){std::cerr<<error.what()<<std::endl;return 1;}
}
