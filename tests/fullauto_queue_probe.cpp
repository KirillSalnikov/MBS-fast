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
        const std::vector<int> depths={1,2,config.maxReflections};FullAuto::Budgets counts(3,std::vector<int>(e.theta.size()));std::vector<size_t> ids;
        for(size_t i=0;i<e.theta.size();++i){ids.push_back(i);counts[0][i]=e.initial;counts[1][i]=i%2?e.initial/2:e.initial/4;counts[2][i]=e.initial/4;}
        const double begin=FullAuto::Seconds();auto result=FullAuto::EvaluateBudget(e,depths,counts,ids,e.seeds.data());
        FullAuto::WriteMatrix(e.root+"/frozen_matrix.dat",e.theta,FullAuto::Mean(FullAuto::Combine(result.levels)));
        FullAuto::Atomic(e.root+"/frozen_seconds.txt",FullAuto::Number(FullAuto::Seconds()-begin)+"\n");return 0;
    }catch(const std::exception& error){std::cerr<<error.what()<<std::endl;return 1;}
}
