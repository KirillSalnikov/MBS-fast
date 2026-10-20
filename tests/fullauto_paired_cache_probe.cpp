#define main lazy_restore_test_main
#include "fullauto_lazy_restore_probe.cpp"
#undef main
int main(int argc,char** argv){try{
    Require(argc==3,"OUTPUT BINARY");
    for(bool controls:{false,true}){
        ArgPP args;RunConfig config;Engine e(args,config);
        SetupMemory(e,std::string(argv[1])+"/controls"+Number(controls),argv[2],controls);
        e.workerPolicy="process";e.workQueuePolicy="off";e.theta={170,175,180};e.initial=256;e.orientationBlock=64;
        e.dispatchThreads=2;e.payloads.SetLimit(1);
        // First assembly exercises the fallback and establishes an exact
        // reference before the reusable prefix exists.
        auto first=e.SharedStage(256,e.seeds.data(),8,6,4,e.theta);
        e.forbidNativeWork=true;
        std::pair<std::vector<Job>,std::vector<Job>> reused;
        Require(e.CachedPairPrefix(256,e.seeds.data(),8,6,4,e.theta,reused),"exact prefix not reused");
        for(size_t s=0;s<8;++s){
            Require(first.first[s].identity==reused.first[s].identity&&first.second[s].identity==reused.second[s].identity,"prefix identity changed");
            Require(first.first[s].matrix==reused.first[s].matrix&&first.second[s].matrix==reused.second[s].matrix,"all16 changed");
            Require(first.first[s].controls==reused.first[s].controls&&first.second[s].controls==reused.second[s].controls,"controls changed");
            IntegralReader reader;
            for(const auto& pair:{std::make_pair(first.first[s],reused.first[s]),std::make_pair(first.second[s],reused.second[s])}){
                const auto a=reader.Resolve(pair.first.folder,pair.first.depth),b=reader.Resolve(pair.second.folder,pair.second.depth);
                Require(a.area==b.area&&a.ext==b.ext,"integrals changed");
            }
        }
        auto view=std::find_if(e.cache.begin(),e.cache.end(),[](const std::pair<const std::string,Job>& entry){return entry.second.pairedDepth==4;});
        Require(view!=e.cache.end(),"missing paired source");
        const auto opens=ValidationFiles().opens;auto lazy=e.LowerView(view->second,4,false);
        Require(lazy.matrix.empty()&&ValidationFiles().opens==opens,"metadata registration read a payload");
        std::pair<std::vector<Job>,std::vector<Job>> absent;
        Require(!e.CachedPairPrefix(512,e.seeds.data(),8,6,4,e.theta,absent)&&absent.first.empty(),"uncomputed prefix reused");
        Require(!e.CachedPairPrefix(256,e.seeds.data(),8,6,4,{170,180},absent),"partial grid incorrectly reused");
        ++e.phi;Require(!e.CachedPairPrefix(256,e.seeds.data(),8,6,4,e.theta,absent),"different phi reused");--e.phi;
        // An absent lower view must fall back before returning any output.
        const auto saved=e.cache.at(lazy.identity);e.cache.erase(lazy.identity);
        Require(!e.CachedPairPrefix(256,e.seeds.data(),8,6,4,e.theta,absent),"new view bypassed row selection");
        e.CacheJob(saved);
        // With payload caching disabled, a corrupt selected leaf must still
        // be rejected through the reused prefix's dependency graph.
        const auto victim=view->second;const auto original=Read(victim.data);Atomic(victim.data,original+"corrupt\n");
        bool rejected=false;try{std::pair<std::vector<Job>,std::vector<Job>> bad;e.CachedPairPrefix(256,e.seeds.data(),8,6,4,e.theta,bad);}catch(const std::exception&){rejected=true;}
        Atomic(victim.data,original);Require(rejected,"corrupt descendant accepted");
        // A new last-exact candidate must change the prefix identity, even
        // if all its values happen to coincide with an earlier candidate.
        Job later=victim;later.identity="zz-new-exact-candidate";later.pairedDepth=0;e.CacheJob(later);
        Require(!e.CachedPairPrefix(256,e.seeds.data(),8,6,4,e.theta,absent),"changed exact selection reused stale prefix");
        e.samplingIndex[Engine::SamplingSlot(later.count,later.seed,later.depth,later.offset)].erase(later.identity);
        e.cache.erase(later.identity);
        Engine restored(args,config);SetupMemory(restored,e.root,argv[2],controls);
        restored.theta=e.theta;restored.initial=256;restored.orientationBlock=64;restored.workerPolicy="process";restored.workQueuePolicy="off";
        restored.payloads.SetLimit(1);restored.WarmCache();restored.forbidNativeWork=true;
        auto resumed=restored.SharedStage(256,e.seeds.data(),8,6,4,e.theta);
        Require(restored.Data(resumed.first)==e.Data(first.first)&&restored.Data(resumed.second)==e.Data(first.second),"resume changed estimator");
        std::cout<<"PASS controls="<<controls<<": eight replicas, all16/controls/integrals exact, lazy registration, fallback, corrupt-source rejection, resume\n";
    }
    return 0;
}catch(const std::exception& error){std::cerr<<error.what()<<'\n';return 1;}}
