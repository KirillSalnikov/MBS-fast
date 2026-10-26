// Real immutable native jobs, paired depths and signed cutoff correction.
// Baseline/candidate builds use the same worker and frozen sampling.
#ifndef MBS_TEST_FULLAUTO_PATH
#define MBS_TEST_FULLAUTO_PATH "../src/FullAuto.cpp"
#endif
#include MBS_TEST_FULLAUTO_PATH
using namespace FullAuto;

void SetupPublication(Engine& e,const std::string& root,const std::string& binary,bool controls){
    e.root=root;e.signature="native-publication-v1";Directory(root);Directory(root+"/jobs");
    e.exe=binary;e.binaryHash=FileHash(binary);e.compactIndex=true;e.controls=controls;e.sharedMeans=controls;
    e.depthSharing=true;e.gpus={""};e.seeds.assign(production,production+8);e.workerThreads=1;e.analyticMeanThreads=1;e.hostBudgetMb=1024;
    e.phi=4;e.target=8;e.theta={170,173,176,180};e.initial=512;e.orientationBlock=512;
    e.workerPolicy="persistent";e.resourcePolicy="off";e.balancedCorrection=true;e.selectedCutoff="off";
    e.base={"--method","po","--backend","cpu","--particle","1","2","1","--ri","1.3116",controls?".001":"0",
        "--wavelength-um",".532","--symmetry","1","1","--cutoff-profile","off"};
    WriteGrid(root+"/requested_theta.csv",e.theta);
}
int main(int argc,char** argv){try{
    Require(argc==3,"OUTPUT BINARY");const std::string root=argv[1];Directory(root);
    ArgPP args;RunConfig config;double elapsed=0;size_t flushes=0,jobs=0;
    std::ostringstream proof;proof<<std::setprecision(17);
    for(bool controls:{false,true})for(int paired:{0,6}){
        Engine e(args,config);const auto caseRoot=root+"/controls"+Number(controls)+"_paired"+Number(paired);
        SetupPublication(e,caseRoot,argv[2],controls);
        const double began=Seconds();
        const auto result=e.DirectStage(512,e.seeds.data(),8,8,e.theta,true,false,0,paired);
        elapsed+=Seconds()-began;flushes+=e.indexJournalFlushes;jobs+=e.cache.size();
#ifndef MBS_BASELINE_PROBE
        Require(e.indexJournalFlushes==1,"native correction journal still flushed per child");
#endif
        for(size_t s=0;s<result.size();++s){
            const auto& parent=result[s];std::set<std::string> verified;ValidateStoredJob(parent.folder,verified);
            for(int depth:{8,paired})if(depth){const Job output=depth==8?parent:e.LowerView(parent,depth);
                proof<<controls<<' '<<paired<<' '<<depth<<' '<<output.seed<<' '<<output.count<<'\n';
                for(const auto& row:output.matrix){for(double v:row)proof<<v<<' ';proof<<'\n';}
                for(const auto& row:output.controls){for(double v:row)proof<<v<<' ';proof<<'\n';}
                IntegralReader reader;const auto integral=reader.Resolve(output.folder,depth);proof<<integral.area<<' '<<integral.ext<<'\n';
                // The signed source sampling and order are also part of the proof.
                std::istringstream sources(Read(parent.folder+"/cutoff_sources.tsv"));int d,sign;std::string folder;
                while(sources>>d>>sign>>folder)if(d==depth){const auto& j=*e.folderIndex.at(folder);
                    proof<<d<<' '<<sign<<' '<<j.seed<<' '<<j.count<<' '<<j.offset<<' '<<j.phi<<'\n';}
            }
        }
#ifndef MBS_BASELINE_PROBE
        Require(e.liveNative.empty()&&e.dispatchProgress.empty(),"completed dispatch still shown active");
        Require(e.nativeCommitted==24&&e.nativeSubmitted==24&&e.nativeReturned==24,"native telemetry counts wrong");
        Require(e.nativeOrientations==8*(512+64+64),"orientation counter changed native sampling");
        e.WriteExecutionProgress(true);
        const auto progress=Read(caseRoot+"/execution_progress.json");
        Require(progress.find("\"batch_active\":false")!=std::string::npos&&progress.find("\"accuracy_recomputed\":false")!=std::string::npos,
            "progress invents a fresh accuracy or leaves a stale batch");
        e.executionProgressWritten=Seconds()-11;e.WriteExecutionProgress();
        Require(Read(caseRoot+"/execution_progress.json")!=progress,"heartbeat does not update");
#endif
        // Simulate loss of the entire acceleration journal. Every native and
        // corrected marker must restore without rerunning a single worker.
        for(auto& service:e.services)service.reset();
        Atomic(caseRoot+"/completed_jobs.journal","");
        Engine resumed(args,config);SetupPublication(resumed,caseRoot,argv[2],controls);
        resumed.WarmCache();resumed.forbidNativeWork=true;
        const auto restored=resumed.DirectStage(512,resumed.seeds.data(),8,8,resumed.theta,true,false,0,paired);
        for(size_t s=0;s<result.size();++s)Require(restored[s].matrix==result[s].matrix&&restored[s].controls==result[s].controls,
            "journal fallback changed matrices or controls");
    }
    Atomic(root+"/proof.txt",proof.str());
    Atomic(root+"/timing.json","{\"worker_seconds\":"+Number(elapsed)+",\"index_flushes\":"+Number(flushes)
        +",\"saved_jobs\":"+Number(jobs)+",\"proof_fnv\":"+Quoted(Hash(proof.str()))+"}");
    std::cout<<Read(root+"/timing.json")<<'\n';
    std::cout<<"PASS real all16/controls/signed integrals; native sampling; durable-marker journal recovery\n";return 0;
}catch(const std::exception& error){std::cerr<<error.what()<<'\n';return 1;}}
