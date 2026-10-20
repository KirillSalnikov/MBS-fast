#define main legacy_lazy_restore_test_main
#include "fullauto_lazy_restore_probe.cpp"
#undef main
int main(int argc,char** argv){try{
    Require(argc==3,"OUTPUT BINARY");
    for(bool controls:{false,true}){
        ArgPP args;RunConfig config;Engine e(args,config);
        SetupMemory(e,std::string(argv[1])+"/controls"+Number(controls),argv[2],controls);
        e.ConfigureIndex(true);e.workerPolicy="process";e.workQueuePolicy="off";e.initial=128;e.orientationBlock=64;
        e.payloads.SetLimit(1);
        const auto first=e.SharedStage(128,e.seeds.data(),8,6,4,e.theta);
        const auto upper=e.Data(first.first),lower=e.Data(first.second);
        for(const auto& item:e.cache){
            const auto& j=item.second;
            Require(!Exists(j.folder+"/sampling.index")&&!Exists(j.folder+"/sampling.index.fnv"),"per-job index files remain");
            const auto path=SamplingIndexPath(j.folder);
            Require(DecodeSamplingIndex(path,Read(path))==e.SamplingText(j),"embedded sampling metadata changed");
        }
        const auto journalPath=e.root+"/completed_jobs.journal",journal=Read(journalPath);
        size_t offset=0,last=0;
        while(offset<journal.size()){
            last=offset;const auto end=journal.find('\n',offset);Require(end!=std::string::npos,"bad fixture journal");
            std::istringstream header(journal.substr(offset,end-offset));size_t length;header>>length;offset=end+1+length+1;
        }
        // Crash after a completed job was committed but before its journal
        // append finished. Recover that job from its directory, without optics.
        Atomic(journalPath,journal.substr(0,last)+"120 ");
        Engine resumed(args,config);SetupMemory(resumed,e.root,argv[2],controls);resumed.ConfigureIndex(false);
        resumed.orientationBlock=64;resumed.initial=128;resumed.payloads.SetLimit(1);resumed.forbidNativeWork=true;
        resumed.WarmCache();
        Require(resumed.cache.size()==e.cache.size(),"torn tail lost completed jobs");
        const auto restored=resumed.SharedStage(128,e.seeds.data(),8,6,4,e.theta);
        Require(resumed.Data(restored.first)==upper&&resumed.Data(restored.second)==lower,"compact restore changed all16");
        for(size_t seed=0;seed<8;++seed){
            IntegralReader a,b;
            const auto x=a.Resolve(first.first[seed].folder,6),y=b.Resolve(restored.first[seed].folder,6);
            Require(x.area==y.area&&x.ext==y.ext,"compact integral provenance changed");
        }
        for(const auto& item:e.cache){
            const auto& a=item.second;const auto& b=resumed.cache.at(item.first);
            Require(a.seconds==b.seconds&&a.controlRefinement==b.controlRefinement&&a.controlSetup==b.controlSetup
                &&a.traceSeconds==b.traceSeconds&&a.processingSeconds==b.processingSeconds
                &&a.preparationThreads==b.preparationThreads,"journal changed scheduling diagnostics");
        }
        Engine quick(args,config);SetupMemory(quick,e.root,argv[2],controls);quick.ConfigureIndex(false);quick.payloads.SetLimit(1);
        quick.WarmCache();Require(quick.cache.size()==e.cache.size()&&quick.payloads.Bytes()==0,"journal eagerly restored matrices");
        const auto victim=quick.cache.begin()->second;const auto markerPath=victim.folder+"/complete.txt";
        const auto marker=Read(markerPath);Atomic(markerPath,marker+"\n");
        bool rejected=false;try{quick.Materialize(victim);}catch(const std::exception&){rejected=true;}
        Atomic(markerPath,marker);Require(rejected,"changed journal source marker accepted");
        auto damaged=Read(journalPath);damaged[damaged.find('\n')+1]^=1;Atomic(journalPath,damaged);
        Engine broken(args,config);SetupMemory(broken,e.root,argv[2],controls);broken.ConfigureIndex(false);
        rejected=false;try{broken.WarmCache();}catch(const std::exception&){rejected=true;}
        Require(rejected,"corrupt committed journal record accepted");
        Atomic(journalPath,journal);
        std::cout<<"PASS compact controls="<<controls<<": paired all16/integrals exact; no per-job index files; lazy journal restore; torn append recovery; corrupt marker/record rejected\n";
    }
    return 0;
}catch(const std::exception& error){std::cerr<<error.what()<<'\n';return 1;}}
