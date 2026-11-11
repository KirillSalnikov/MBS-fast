#include "../src/FullAuto.cpp"
using namespace FullAuto;
namespace {
Job Sample(double wall,double processing=.01){Job j;j.seconds=wall;j.processingSeconds=processing;j.matrix.resize(2);
    for(auto& row:j.matrix){row.fill(0.);row[0]=2;row[1]=row[4]=-.1;row[5]=1.5;row[10]=1.;row[15]=.5;}return j;}
void PurePolicy(){
    Kernel baseline;const auto candidates=KernelCalibration::Candidates(baseline,false,false,false);
    Require(candidates.size()==5,"default layouts not deduplicated");
    Kernel fixed;fixed.pipeline=true;fixed.block=128;fixed.warp="thread";
    Require(KernelCalibration::Candidates(fixed,true,true,true).size()==1,"explicit layouts were overridden or duplicated");
    for(const auto& k:KernelCalibration::Candidates(fixed,false,true,true))Require(k.block==128&&k.warp=="thread","explicit environment settings changed");
    auto ignore=[](const Kernel&,int,bool,const Job&,double){};
    double clock=0;size_t requests=0;
    auto slow=[&](const Kernel&,int,bool){++requests;clock+=8;return Sample(8,8);};
    const auto stopped=KernelCalibration::Select(baseline,candidates,256,16,8192,1e-5,true,false,slow,ignore,[&]{return clock;});
    Require(requests==1&&KernelCalibration::Key(stopped.selected)==KernelCalibration::Key(baseline)
        &&stopped.elapsed>stopped.budget,"slow nonpreemptible bootstrap did not stop further launches");
    clock=0;requests=0;std::map<std::string,int> measured;size_t warmups=0;
    auto cold=[&](const Kernel& k,int,bool warm){++requests;
        if(warm){++warmups;clock+=.3;return Sample(.3,.01);}
        ++measured[KernelCalibration::Key(k)];const double seconds=!k.pipeline?.16:k.block==128?.04:k.block==256?.08:k.warp=="thread"?.09:.07;
        clock+=seconds;return Sample(seconds,seconds*.9);};
    const auto tuned=KernelCalibration::Select(baseline,candidates,256,16,65536,1e-5,true,false,cold,ignore,[&]{return clock;});
    Require(tuned.selected.pipeline&&tuned.selected.block==128,"cold warmup was included in layout ranking");
    Require(tuned.requests<=12&&warmups==5&&measured[KernelCalibration::Key(tuned.selected)]==2,"unbounded or missing top-two repeat");
    clock=0;requests=0;
    auto wrong=[&](const Kernel& k,int,bool warm){++requests;double seconds=warm?.002:!k.pipeline?.1:k.block==128?.0001:.05;
        clock+=seconds;auto j=Sample(seconds,seconds);if(!warm&&k.block==128)j.matrix[0][15]+=.1;return j;};
    const auto rejected=KernelCalibration::Select(baseline,candidates,256,16,65536,1e-5,true,false,wrong,ignore,[&]{return clock;});
    Require(rejected.selected.block!=128,"fast numerically invalid layout won");
    for(const auto& trial:rejected.trials)if(trial.kernel.block==128)Require(!trial.accepted&&trial.seconds.size()==1,"invalid layout was repeated or accepted");
    auto a=Sample(.01),b=a;a.matrix[0][1]=std::numeric_limits<double>::quiet_NaN();
    Require(KernelCalibration::Error(a,b,true)>1,"nonfinite matrix passed numerical check");
    a=b;a.matrix[0][0]*=1.01;Require(KernelCalibration::Error(a,b,true)>.005,"ratio normalization change not detected");
    clock=0;requests=0;auto noWork=[&](const Kernel&,int,bool)->Job{++requests;throw std::runtime_error("unexpected fixed-layout work");};
    const auto allFixed=KernelCalibration::Select(fixed,{fixed},256,16,8192,1e-5,true,false,noWork,ignore,[&]{return clock;});
    Require(requests==0&&KernelCalibration::Key(allFixed.selected)==KernelCalibration::Key(fixed),"fixed layout launched tuning");
    // Tiny production jobs must not fund a full five-layout search.
    clock=0;requests=0;auto tiny=[&](const Kernel&,int,bool){++requests;clock+=.02;return Sample(.02,.019);};
    const auto limited=KernelCalibration::Select(baseline,candidates,256,16,32,1e-5,true,false,tiny,ignore,[&]{return clock;});
    Require(requests==1&&KernelCalibration::Key(limited.selected)==KernelCalibration::Key(baseline),"small production budget did not retain default");
    std::cout<<"PASS bounded calibration: cold warmup excluded, complete request costs, top-two repeats, explicit settings, numerical reject, slow/tiny fallback\n";
}
void PrivateCleanup(const std::string& path,bool memory){
    ArgPP args;RunConfig config;Engine e(args,config);e.root=path;Directory(path);e.memoryMode=memory;
    if(memory){e.memoryStore.reset(new FullAutoMemoryStore::Store(path,1024*1024,path));e.memoryScope.reset(new FullAutoMemoryStore::Scope(*e.memoryStore));}
    const std::string scientific=path+"/jobs/science",trial=path+"/kernel_trials/t012345";
    Directory(scientific+"/result");Directory(trial+"/result");
    Atomic(scientific+"/result/result.dat","science\n");Atomic(scientific+"/complete.txt","untouched scientific marker\n");
    Atomic(trial+"/complete.txt","temporary\n");Atomic(trial+"/result/result.dat","trial\n");
    e.kernelTrialNonce="fresh-timing";Job job;job.identity="trial";job.folder=trial;job.count=64;job.seed=4001;job.depth=6;job.matrix=Sample(.1).matrix;
    e.CacheJob(job);Require(e.cache.empty()&&e.folderIndex.empty()&&e.samplingIndex.empty()&&e.payloads.Bytes()==0,"timing job became a scientific source");
    e.RemoveKernelTrial(trial);Require(!Exists(trial)&&Read(scientific+"/result/result.dat")=="science\n","trial cleanup changed science");
    for(const auto& bad:{scientific,path+"/kernel_trials/t../jobs/science",path+"/kernel_trials/"}){bool rejected=false;try{e.RemoveKernelTrial(bad);}catch(const std::exception&){rejected=true;}Require(rejected,"cleanup accepted nonprivate namespace");}
    e.kernelTrialNonce.clear();job.identity="scientific";job.folder=scientific;e.CacheJob(job);Require(e.cache.size()==1,"production cache did not resume after timing");
    std::cout<<"PASS private trial cleanup storage="<<(memory?"memory":"legacy")<<"; scientific files and indices preserved\n";
}
void RestorePolicy(){
    for(bool explicitCache:{false,true})for(bool failed:{false,true}){
        ArgPP args;RegisterCliOptions(args);const auto traceFlag=CliCanonicalFlag("trace_cache_mb");if(explicitCache){const char* command[]={"probe",traceFlag.c_str(),"32"};args.Parse(3,command);}
        RunConfig config;Engine e(args,config);e.base={"--threads","7","--orientation-chunk","128"};
        if(explicitCache){e.base.push_back(traceFlag);e.base.push_back("32");}
        const auto before=e.base;e.resourcePolicy="auto";e.checkpointInterval=300;e.packJobs=true;e.shortPacketVerifying=false;
        e.kernel.pipeline=true;e.kernel.block=128;e.kernel.warp="thread";const auto initial=KernelCalibration::Key(e.kernel);
        try{Engine::KernelCalibrationScope scope(e);Require(e.resourcePolicy=="off"&&e.checkpointInterval==0&&!e.packJobs&&e.shortPacketVerifying,"timing resources were not frozen");
            if(explicitCache)Require(e.base==before,"explicit prepared cache overridden");
            else Require(e.base.size()==before.size()+2&&e.base[e.base.size()-2]==traceFlag&&e.base.back()=="0","automatic cache not disabled for fresh orientation comparison");
            e.kernelTrialNonce="fixture";e.kernel.block=256;if(failed)throw std::runtime_error("injected timing failure");scope.keepKernel=true;
        }catch(const std::runtime_error&){}
        Require(e.base==before&&e.resourcePolicy=="auto"&&e.checkpointInterval==300&&e.packJobs&&!e.shortPacketVerifying&&e.kernelTrialNonce.empty(),"kernel calibration changed production options");
        Require(failed?KernelCalibration::Key(e.kernel)==initial:e.kernel.block==256,"failed calibration did not restore default or committed selection was lost");
    }
    std::cout<<"PASS scope restores base/resources/checkpoint/pack/cache on success and failure; explicit cache retained\n";
}
}
int main(int argc,char** argv){try{Require(argc==2,"new output path required");Require(!Exists(argv[1]),"use a new output path");PurePolicy();RestorePolicy();PrivateCleanup(std::string(argv[1])+"/legacy",false);PrivateCleanup(std::string(argv[1])+"/memory",true);return 0;}
catch(const std::exception& error){std::cerr<<error.what()<<'\n';return 1;}}
