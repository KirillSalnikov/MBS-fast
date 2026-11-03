// A signed multilevel/control-variate estimate can have nonpositive M11
// before its corrections converge. Preserve the estimate and all source
// blocks; refine its sampling instead of clipping it or rejecting controls.
// Return 0 for an ordinary estimate, 1 for scheduled work, 2 at the limit.
int RecoverNonpositiveProduction(Engine& e,AdaptiveState& state,
        const std::vector<int>& depths,const AdaptiveData& stage,
        const std::vector<size_t>& ids,const Samples& samples,const char* phase){
    const auto mean=Mean(samples);Require(mean.size()==ids.size(),"positive sampling row mismatch");
    std::vector<size_t> bad;
    for(size_t r=0;r<ids.size();++r)if(!(mean[r][0]>0))bad.push_back(r);
    if(bad.empty())return 0;
    const auto current=state.counts;auto proposed=current;
    Statistics priority;priority.interval.resize(ids.size());priority.ratioInterval.resize(ids.size());
    for(size_t r:bad){priority.interval[r][0]=1e100;
        for(auto& level:proposed)level[ids[r]]=Doubled(e,level[ids[r]]);}
    if(!e.local)for(auto& level:proposed){const int peak=*std::max_element(level.begin(),level.end());
        std::fill(level.begin(),level.end(),peak);}
    AdditionalWork work(e,e.seeds.data());
    BoundRoundWork(e,work,current,proposed,depths,ids,priority);
    const bool grows=proposed!=current;
    const auto folder=e.root+"/nonpositive_sampling/round"+Number(state.round)+"_"+phase;
    Directory(folder);WriteMatrix(folder+"/mueller_evaluated.dat",Angles(e.theta,ids),mean);
    Atomic(folder+"/sources.tsv",stage.sources);
    std::ostringstream details;details<<"theta_deg,M11";
    for(size_t l=0;l<current.size();++l)details<<",N_before"<<l<<",N_next"<<l;
    details<<'\n';for(size_t r:bad){details<<Number(e.theta[ids[r]])<<','<<Number(mean[r][0]);
        for(size_t l=0;l<current.size();++l)details<<','<<current[l][ids[r]]<<','<<proposed[l][ids[r]];
        details<<'\n';state.valid[ids[r]]=false;state.coarse[ids[r]]=false;}
    Atomic(folder+"/sampling.csv",details.str());state.counts=proposed;state.streak=0;
    const std::string status=grows?"refining_nonpositive_M11":"sampling_limit_nonpositive_M11";
    const std::string report="{\"converged\":false,\"phase\":"+Quoted(status)+",\"round\":"+Number(state.round)
        +",\"snapshot_folder\":"+Quoted(folder)+",\"matrix\":"+Quoted(folder+"/mueller_evaluated.dat")
        +",\"nonpositive_rows\":"+Number(bad.size())+",\"combined_error95\":null,\"accuracy_resolved\":false"
        +",\"matrix_clipped\":false,\"controls_unchanged\":true,\"saved_blocks_retained\":true,\"eight_seeds_retained\":true}";
    e.Record(e.root+"/intermediate_current.json",report);e.Log("nonpositive_production_sampling",report);
    SaveAdaptive(e,state,status);e.State(status,state.round,PeakBudgets(state.counts),ids.size(),1e100,1e100,0);
    return grows?1:2;
}
