#define main original_lazy_main
#include "fullauto_lazy_restore_probe.cpp"
#undef main
int main(int argc,char** argv){try{
    Require(argc==3,"OUTPUT BINARY");
    for(bool controls:{false,true}){
        ArgPP args;RunConfig config;Engine legacy(args,config),batch(args,config);
        SetupMemory(legacy,std::string(argv[1])+"/legacy"+Number(controls),argv[2],controls);
        SetupMemory(batch,std::string(argv[1])+"/batch"+Number(controls),argv[2],controls);
        legacy.theta=batch.theta={0,5,10,15,170,175,180};legacy.orientationBlock=batch.orientationBlock=32;
        WriteGrid(legacy.root+"/requested_theta.csv",legacy.theta);WriteGrid(batch.root+"/requested_theta.csv",batch.theta);
        if(controls)for(auto* e:{&legacy,&batch}){e->weightTheta=e->theta;e->weights.resize(e->theta.size());
            for(size_t i=0;i<e->theta.size();++i)e->weights[i]={{.3+e->theta[i]*.001,.2-e->theta[i]*.0005}};}
        const std::vector<std::vector<double>> grids={{0,5},{10,15},{170,175},{175,180}};
        for(int offset=0;offset<256;offset+=32){auto jobs=legacy.DirectStage(32,legacy.seeds.data(),8,6,legacy.theta,true,false,offset);
            for(const auto& j:jobs)batch.CacheJob(j);}
        legacy.forbidNativeWork=batch.forbidNativeWork=true;
        legacy.payloads.Clear();batch.payloads.Clear();legacy.payloads.SetLimit(48*1024);batch.payloads.SetLimit(48*1024);
        const auto legacyMisses=legacy.payloads.misses,batchMisses=batch.payloads.misses;
        std::vector<std::vector<Job>> references;
        for(const auto& grid:grids)references.push_back(legacy.PrefixStage(256,legacy.seeds.data(),8,6,grid));
        std::vector<Engine::WorkRequest> requests;for(const auto& grid:grids)requests.push_back({256,0,6,0,grid});
        batch.SetPrefixAssemblyPlan(requests,batch.seeds.data());
        for(size_t k=0;k<grids.size();++k){const auto result=batch.PrefixStage(256,batch.seeds.data(),8,6,grids[k]);
            for(size_t s=0;s<8;++s){const auto& a=result[s];const auto& b=references[k][s];
                Require(a.matrix==b.matrix&&a.controls==b.controls&&a.controlRefinement==b.controlRefinement&&a.controlSetup==b.controlSetup,"streamed prefix changed matrix or controls");
                Require(Read(a.folder+"/angular_sources.tsv")==Read(b.folder+"/angular_sources.tsv"),"streamed prefix changed selected source order");}
            if(k==0)Require(batch.preparedPrefixes.size()==3,"remaining prefixes were not streamed together");
        }
        Require(batch.preparedPrefixes.empty()&&batch.preparedPrefixBytes==0,"consumed prefixes retained memory");
        const auto oldReads=legacy.payloads.misses-legacyMisses,newReads=batch.payloads.misses-batchMisses;
        Require(newReads*2<oldReads,"streaming did not reduce repeated payload loads");
        Require(Read(batch.root+"/fullauto_history.jsonl").find("\"completed_prefixes\":4")!=std::string::npos,"batch path not exercised");
        batch.ClearPrefixAssemblyPlan();
        // A changed native candidate set invalidates a previously prepared
        // prefix before it can replace the authoritative source selection.
        batch.SetPrefixAssemblyPlan(requests,batch.seeds.data());
        Engine fresh(args,config);SetupMemory(fresh,std::string(argv[1])+"/generation"+Number(controls),argv[2],controls);
        fresh.orientationBlock=32;fresh.theta=legacy.theta;
        WriteGrid(fresh.root+"/requested_theta.csv",fresh.theta);
        for(const auto& item:legacy.cache)if(item.second.count==32)fresh.CacheJob(item.second);
        fresh.forbidNativeWork=true;fresh.SetPrefixAssemblyPlan(requests,fresh.seeds.data());
        auto first=fresh.TakePreparedPrefix(256,fresh.seeds.data(),8,6,grids[0]);Require(bool(first)&&fresh.preparedPrefixes.size()==3,"generation fixture not prepared");
        Job candidate=legacy.cache.begin()->second;candidate.count=32;candidate.identity="new_native_candidate";candidate.folder=fresh.root+"/unselected_candidate";
        fresh.CacheJob(candidate);Require(fresh.preparedPrefixGeneration!=fresh.prefixLeafGeneration,"native candidate did not invalidate prepared sampling");
        // Disable preparation of replacements: the stale item must be dropped.
        fresh.prefixAssemblyPlans.clear();Require(!fresh.TakePreparedPrefix(256,fresh.seeds.data(),8,6,grids[1]),"stale source selection reused");
        Require(fresh.preparedPrefixes.empty()&&fresh.preparedPrefixBytes==0,"invalidated prefixes retained memory");
        std::cout<<"PASS controls="<<controls<<": four prefixes, eight seeds, exact all16/controls/source order, no native work, candidate invalidation and bounded lifetime; payload reads "<<oldReads<<" -> "<<newReads<<'\n';
    }
    return 0;
}catch(const std::exception& error){std::cerr<<error.what()<<'\n';return 1;}}
