// A linear interpolation residual estimates the same expectation with its
// own matched prefixes. It need not inherit a noisy production row's N.
// Counts are frozen before reading the residual, saved with AdaptiveState,
// and grow only after a failed measured residual. No forecast accepts data.
void PrepareResidualBudget(Engine& e,AdaptiveState& state){
    if(!e.local)return;
    const int phi=ResidualPhi(e);
    if(!state.interpolationCounts.empty()){
        const bool changed=state.interpolationPhi?state.interpolationPhi!=phi:MixedResidualPhi(e);
        if(changed){
            // Old mixed-phi noise says nothing about the cost of the newly
            // coupled residual. Freeze a small new pilot before reading it;
            // old native blocks and all production budgets remain reusable.
            const auto old=PeakBudgets(state.interpolationCounts);
            const int start=std::max(e.minimum,512);
            for(auto& level:state.interpolationCounts)for(int& n:level)n=std::min(n,start);
            e.Log("interpolation_phi_budget_restart","{\"old_phi\":"+Number(state.interpolationPhi)+",\"new_phi\":"+Number(phi)
                +",\"old_peak_base_count\":"+Number(old[0])+",\"new_peak_base_count\":"+Number(PeakBudgets(state.interpolationCounts)[0])
                +",\"production_unchanged\":true,\"cached_blocks_retained\":true,\"accuracy_not_accepted\":true}");
        }
        state.interpolationPhi=phi;return;
    }
    state.interpolationCounts=state.counts;
    const auto ids=AdaptiveGuards(state.knots,e.theta.size());
    for(size_t l=0;l<state.counts.size();++l){
        int smallest=e.maxCount;for(size_t id:ids)smallest=std::min(smallest,state.counts[l][id]);
        const int start=std::max(e.minimum,std::min(smallest,l?e.minimum:e.initial)/8);
        for(size_t id=0;id<e.theta.size();++id)
            state.interpolationCounts[l][id]=std::min(state.counts[l][id],start);
    }
    if(MixedResidualPhi(e))for(auto& level:state.interpolationCounts)for(int& n:level)n=std::min(n,std::max(e.minimum,512));
    state.interpolationPhi=phi;
    e.Record(e.root+"/interpolation_sampling_policy.json",
        "{\"policy\":\"separate matched residual prefixes v1\",\"replicas\":8,"
        "\"initial_base_fraction\":0.125,\"sampling_choices_frozen_before_residual\":true,"
        "\"growth_requires_measured_failure\":true,\"all16_and_normalized_ratios\":true,"
        "\"tolerance\":"+Number(e.interpolation)+",\"production_counts_preserved\":true,"
        "\"native_blocks_deleted\":0,\"accuracy_gates_unchanged\":true}");
    e.Log("interpolation_budget_initialized","{\"policy\":\"separate matched residual prefixes v1\","
        "\"base_count_per_seed\":"+Number(PeakBudgets(state.interpolationCounts)[0])
        +",\"production_peak_base_count\":"+Number(PeakBudgets(state.counts)[0])
        +",\"accuracy_accepted\":false}");
}
void ResidualBudgetFile(Engine& e,const AdaptiveState& state,const std::vector<size_t>& ids,const std::string& folder){
    if(state.interpolationCounts.empty())return;
    std::ostringstream out;out<<"theta_deg";for(size_t l=0;l<state.interpolationCounts.size();++l)out<<",N_residual_level"<<l;
    out<<",phi_points\n";for(size_t id:ids){out<<Number(e.theta[id]);
        for(const auto& level:state.interpolationCounts)out<<','<<level[id];out<<','<<ResidualPhi(e)<<'\n';}
    Atomic(folder+"/interpolation_orientation_budget.csv",out.str());
}
