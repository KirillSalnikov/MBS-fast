#define main old_lazy_main
#include "fullauto_lazy_restore_probe.cpp"
#undef main
int main(int argc,char** argv){try{
    Require(argc==3,"new output root and native CPU worker binary required");Require(!Exists(argv[1]),"use a new output directory");
    for(bool memory:{false,true}){
        ArgPP args;RegisterCliOptions(args);const auto traceFlag=CliCanonicalFlag("trace_cache_mb");const char* command[]={"probe",traceFlag.c_str(),"0"};args.Parse(3,command);
        RunConfig config;Engine e(args,config);const std::string root=std::string(argv[1])+(memory?"/memory":"/legacy");SetupMemory(e,root,argv[2],false);
        e.base.insert(e.base.end(),{traceFlag,"0"});e.theta={150,170};e.kernelCount=16;e.initial=65536;e.kernel.pipeline=false;e.fixedPipeline=true;
        e.resourcePolicy="auto";e.checkpointInterval=300;e.shortPacketVerifying=false;e.memoryMode=memory;
        if(memory){e.memoryStore.reset(new FullAutoMemoryStore::Store(root,8*1024*1024,root));e.memoryScope.reset(new FullAutoMemoryStore::Scope(*e.memoryStore));}
        // A real scientific job already has exactly the calibration seed/N.
        // Fresh timing must neither reuse its old seconds nor replace it.
        const unsigned seed=4001;const Job scientific=e.DirectStage(16,&seed,1,6,e.theta,false)[0];
        const auto bytes=Read(scientific.data),marker=Read(scientific.folder+"/complete.txt");const auto originalBase=e.base;
        const size_t cacheSize=e.cache.size();const auto recipe=e.ExecutionHash();uint64_t committed=e.nativeCommitted;
        for(int repeat=0;repeat<2;++repeat){
            e.CalibrateKernel(e.theta);Require(e.nativeCommitted>committed,"calibration replayed a cached wall time");committed=e.nativeCommitted;
            Require(e.base==originalBase&&e.resourcePolicy=="auto"&&e.checkpointInterval==300&&!e.shortPacketVerifying&&e.kernelTrialNonce.empty()
                &&e.ExecutionHash()==recipe,"calibration changed the native production recipe");
            Require(e.cache.size()==cacheSize&&Read(scientific.data)==bytes&&Read(scientific.folder+"/complete.txt")==marker,"calibration changed scientific estimator sources");
            if(memory)Require(e.memoryStore->List(root+"/kernel_trials",true).empty(),"memory timing outputs retained");
            else {DIR* dir=opendir((root+"/kernel_trials").c_str());Require(dir,"missing private trial parent");size_t count=0;while(auto* item=readdir(dir))if(std::string(item->d_name)!="."&&std::string(item->d_name)!="..")++count;closedir(dir);Require(count==0,"legacy timing outputs retained");}
            Require(Read(root+"/kernel_calibration.json").find("\"trial_matrices_retained\":false")!=std::string::npos,"kernel audit missing");
        }
        std::cout<<"PASS native CPU storage="<<(memory?"memory":"legacy")<<": two fresh calibrations, scientific all16/integrity untouched, no trial matrices, recipe restored\n";
    }return 0;
}catch(const std::exception& error){std::cerr<<error.what()<<'\n';return 1;}}
