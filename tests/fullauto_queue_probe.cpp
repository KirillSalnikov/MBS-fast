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
        FullAuto::AdditionalWork forecast(e,e.seeds.data());const auto planned=forecast.Forecast(FullAuto::BudgetRequests(e,depths,counts,ids));
        const double begin=FullAuto::Seconds();auto result=FullAuto::EvaluateBudget(e,depths,counts,ids,e.seeds.data());
        size_t native=0;for(const auto& item:e.cache)if(!item.second.assembled&&std::find(e.seeds.begin(),e.seeds.end(),item.second.seed)!=e.seeds.end())++native;
        FullAuto::AdditionalWork cached(e,e.seeds.data());const auto remaining=cached.Forecast(FullAuto::BudgetRequests(e,depths,counts,ids));
        if(e.QueueEnabled()){
            FullAuto::Require(remaining.jobs==0,"completed frozen allocation forecast new jobs");
            if(!getenv("TASK_MBS_WARM_CACHE"))FullAuto::Require(planned.jobs==native,"forecast differs from native frozen job count");
        }
        FullAuto::Atomic(e.root+"/work_forecast.json","{\"planned_native_jobs\":"+FullAuto::Number(planned.jobs)+",\"saved_native_jobs\":"+FullAuto::Number(native)+",\"remaining_native_jobs\":"+FullAuto::Number(remaining.jobs)+"}\n");
        FullAuto::WriteMatrix(e.root+"/frozen_matrix.dat",e.theta,FullAuto::Mean(FullAuto::Combine(result.levels)));
        FullAuto::Atomic(e.root+"/frozen_seconds.txt",FullAuto::Number(FullAuto::Seconds()-begin)+"\n");
        if(getenv("TASK_MBS_AUDIT_FAILURE")){
            FullAuto::AdaptiveState state;state.counts=counts;state.knots=ids;state.valid.assign(ids.size(),false);state.coarse.assign(ids.size(),false);state.previous.resize(ids.size());
            e.statistical=1e-10;e.eps=1e-9;e.local=true;e.localPolicy="off";
            FullAuto::Require(!FullAuto::IndependentAudit(e,state,depths,result,ids),"deliberately tight real audit unexpectedly passed");
            FullAuto::Require(state.counts!=counts,"real failed audit did not schedule work");FullAuto::SaveAdaptive(e,state);
            FullAuto::AdaptiveState restored=state;FullAuto::Require(FullAuto::RestoreAdaptive(e,restored),"cannot restore failed real audit");
            const unsigned previousSerial=restored.auditSerial;
            FullAuto::Require(!FullAuto::IndependentAudit(e,restored,depths,result,ids)&&restored.auditSerial>previousSerial,"restored audit reused seeds or passed without precision");
        }
        return 0;
    }catch(const std::exception& error){std::cerr<<error.what()<<std::endl;return 1;}
}
