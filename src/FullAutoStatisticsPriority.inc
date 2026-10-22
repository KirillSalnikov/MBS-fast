// Interpolation and pointwise sampling have separate failures. In particular,
// a noisy support has zero interpolation residual and must still grow while a
// different guard is waiting for angular refinement. Forecasts only schedule
// work; every acceptance gate continues to use the measured eight replicas.
Budgets PlanPointwiseStatistics(Engine& e,const AdaptiveState& state,const AdaptiveData& stage,
        const std::vector<size_t>& ids,const Statistics& stats,const InterpolationResult& interp){
    Budgets proposed=state.counts;size_t changed=0;
    std::ostringstream report;
    report<<"theta_deg,level,current_count,proposed_count,statistical95,scheduling_tolerance\n";
    for(size_t r=0;r<ids.size();++r){
        const size_t id=ids[r];const double tolerance=RemainingStatistics(e,ResidualAt(e,interp,ids,id));
        const double error=WorstMetric(stats,r,e.ratios).error;
        if(error<=tolerance)continue;
        const auto selected=AllocateRow(e,stage,state.counts,ids,r,stats.mean,tolerance);
        for(size_t l=0;l<selected.size();++l){
            const int current=state.counts[l][id];
            const int cap=current>e.maxCount/4?e.maxCount:4*current;
            proposed[l][id]=std::max(current,std::min(cap,selected[l]));
            if(proposed[l][id]>current){++changed;
                report<<Number(e.theta[id])<<','<<l<<','<<current<<','<<proposed[l][id]
                    <<','<<Number(error)<<','<<Number(tolerance)<<'\n';}
        }
    }
    if(!e.local)for(auto& level:proposed){int peak=0;for(size_t id:ids)peak=std::max(peak,level[id]);
        for(auto& n:level)n=std::max(n,peak);}
    const auto name="pointwise_sampling_priority_round"+Number(state.round)+"_"+Hash(report.str())+".csv";
    Atomic(e.root+"/"+name,report.str());
    e.Log("pointwise_sampling_priority","{\"round\":"+Number(state.round)+",\"growing_level_rows\":"+Number(changed)
        +",\"details\":"+Quoted(name)+",\"runs_during_interpolation_failure\":true,\"accuracy_gates_unchanged\":true}");
    return proposed;
}
void MergeGrowingBudgets(Budgets& target,const Budgets& proposal){
    Require(target.size()==proposal.size(),"sampling proposal level mismatch");
    for(size_t l=0;l<target.size();++l){Require(target[l].size()==proposal[l].size(),"sampling proposal grid mismatch");
        for(size_t id=0;id<target[l].size();++id)target[l][id]=std::max(target[l][id],proposal[l][id]);}
}
Budgets InterpolationSamplingProposal(const AdaptiveState& state,const InterpolationResult& interp){
    Budgets result=state.counts;
    for(const auto& plan:interp.growth)if(!plan.biasDominated&&plan.level>=0)
        for(const auto& next:plan.counts)result.at(size_t(plan.level)).at(next.first)=
            std::max(result.at(size_t(plan.level)).at(next.first),next.second);
    return result;
}
