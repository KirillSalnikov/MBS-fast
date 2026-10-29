// Forecast production AND matched-prefix interpolation. Only scheduling changes:
// every matrix/ratio/interpolation gate still runs on the full grid.
std::vector<Engine::WorkRequest> AuditRequests(const Engine&,const std::vector<int>&,
    const Budgets&,const std::vector<size_t>&,const std::vector<size_t>&);
size_t RoundPhysicalGpus(const Engine& e){
    const auto slots=e.WorkerGpus();return std::max(size_t(1),std::set<std::string>(slots.begin(),slots.end()).size());
}
struct RoundOutputPriority {double error;std::vector<size_t> nodes;};
void BoundRoundWork(Engine& e,AdditionalWork& work,const Budgets& current,Budgets& proposed,
        const std::vector<int>& depths,const std::vector<size_t>& ids,const Statistics& stats,
        const std::vector<size_t>& knots={},const std::vector<std::vector<size_t>>& cohorts={},const Budgets* residualCounts=nullptr,
        const std::vector<RoundOutputPriority>& outputPriorities={}){
    const bool separatePhi=residualCounts&&!knots.empty()&&MixedResidualPhi(e);
    WorkForecast residualForecast;
    if(separatePhi){ResidualPhiScope scope(e,true);AdditionalWork residualWork(e,e.seeds.data());
        residualForecast=residualWork.Forecast(AuditRequests(e,depths,*residualCounts,ids,knots));}
    auto estimate=[&](const Budgets& counts){
        if(separatePhi){auto result=work.Forecast(BudgetRequests(e,depths,counts,ids));
            result.seconds+=residualForecast.seconds;result.jobs+=residualForecast.jobs;
            result.profiled=result.profiled&&residualForecast.profiled;return result;}
        if(residualCounts&&!knots.empty()){
            auto requests=BudgetRequests(e,depths,counts,ids);
            const auto residual=AuditRequests(e,depths,*residualCounts,ids,knots);
            requests.insert(requests.end(),residual.begin(),residual.end());return work.Forecast(requests);}
        return work.Forecast(knots.empty()?BudgetRequests(e,depths,counts,ids):AuditRequests(e,depths,counts,ids,knots));};
    const auto original=estimate(proposed);const double budget=1800.*RoundPhysicalGpus(e);
    if(!original.profiled||original.seconds<=budget)return;
    const Budgets wanted=proposed;Budgets step=current;
    std::set<size_t> coupled;for(const auto& group:cohorts)coupled.insert(group.begin(),group.end());
    for(const auto& group:outputPriorities)coupled.insert(group.nodes.begin(),group.nodes.end());
    std::vector<std::pair<double,std::vector<size_t>>> rank;
    auto score=[&](const std::vector<size_t>& group,double outputError){double error=outputError;bool grows=false;
        for(size_t id:group){const auto at=std::lower_bound(ids.begin(),ids.end(),id);Require(at!=ids.end()&&*at==id,"round cohort outside evaluated grid");
            error=std::max(error,WorstMetric(stats,at-ids.begin(),e.ratios).error);
            for(size_t l=0;l<current.size();++l)grows=grows||step[l][id]>current[l][id];}
        if(grows)rank.push_back({error,group});};
    for(size_t id:ids)for(size_t l=0;l<current.size();++l){
        Require(wanted[l][id]>=current[l][id],"round work limiter cannot shrink completed prefixes");
        // Preserve the full signed-phi/stencil growth prescription.
        step[l][id]=coupled.count(id)?wanted[l][id]:std::min(wanted[l][id],Doubled(e,current[l][id]));}
    for(const auto& group:cohorts)score(group,0);
    // Rank the error of the actual interpolated output, not just the smaller
    // errors at its individual support nodes. Keep its whole stencil together
    // so latency limiting cannot undo RefineError's coupled growth proposal.
    for(const auto& group:outputPriorities)score(group.nodes,group.error);
    for(size_t id:ids)if(!coupled.count(id))score({id},0);
    std::sort(rank.rbegin(),rank.rend());if(rank.empty())return;
    auto selected=step;auto forecast=estimate(selected);size_t active=rank.size();
    const double fullSeconds=forecast.seconds;bool preservedSharing=false;double partitionOverhead=0;
    // Several hours without feedback do not justify a small tracing saving.
    // Still reject partitions which duplicate most of the tracing work.
    const double overheadLimit=fullSeconds>4*budget?.35:.10;
    for(int attempt=0;e.local&&forecast.profiled&&forecast.seconds>budget&&active>1&&attempt<8;++attempt){
        const size_t nextActive=std::max(size_t(1),size_t(std::ceil(active*std::max(.05,std::min(.5,budget/forecast.seconds)))));Budgets candidate=current;
        std::set<size_t> selectedIds;
        for(size_t k=0;k<nextActive;++k)for(size_t id:rank[k].second){selectedIds.insert(id);
            for(size_t l=0;l<current.size();++l)candidate[l][id]=step[l][id];}
        const auto next=estimate(candidate);
        // Compare the cost of the entire partition, not seconds per selected
        // angle-orientation. Difficult rows can be more expensive and less
        // covered by cached blocks; their cost per unit is not repeated work.
        // The complementary request measures the actual lost trace sharing.
        Budgets remainder=step;
        for(size_t id:selectedIds)for(size_t l=0;l<current.size();++l)remainder[l][id]=current[l][id];
        const auto rest=estimate(remainder);
        const double splitSeconds=next.seconds+rest.seconds;
        partitionOverhead=fullSeconds>0?std::max(0.,splitSeconds/fullSeconds-1):0;
        const bool latencyPriority=fullSeconds>4*budget&&next.profiled&&next.seconds<.9*forecast.seconds;
        if(!next.profiled||!rest.profiled||(splitSeconds>(1+overheadLimit)*fullSeconds&&!latencyPriority)){preservedSharing=true;break;}
        selected=candidate;active=nextActive;forecast=next;
    }
    const bool accepted=forecast.profiled&&forecast.seconds<original.seconds&&selected!=current;
    if(accepted)proposed=selected;
    size_t growing=0;for(size_t id:ids){bool grows=false;for(size_t l=0;l<current.size();++l)grows=grows||proposed[l][id]>current[l][id];growing+=grows;}
    e.Record(e.root+"/round_work_limit.json","{\"policy\":\"complete eight-seed prefix and output/residual cohort scheduling v4\",\"applied\":"+std::string(accepted?"true":"false")
        +",\"output_priority_cohorts\":"+Number(outputPriorities.size())
        +",\"forecast_includes_common_prefix_interpolation\":"+(knots.empty()?"false":"true")+",\"physical_gpus\":"+Number(RoundPhysicalGpus(e))
        +",\"original_worker_seconds\":"+Number(original.seconds)+",\"selected_worker_seconds\":"+Number(accepted?forecast.seconds:original.seconds)
        +",\"advisory_worker_seconds\":"+Number(budget)+",\"selected_growing_rows\":"+Number(growing)+",\"selected_cohorts\":"+Number(active)
        +",\"shared_trace_efficiency_preserved\":"+(preservedSharing?"true":"false")
        +",\"partition_cost_model\":\"selected plus complementary cached-work forecast\",\"partition_overhead_fraction\":"+Number(partitionOverhead)
        +",\"allowed_partition_overhead_fraction\":"+Number(overheadLimit)
        +",\"within_forecast_budget\":"+((accepted?forecast.seconds:original.seconds)<=budget?"true":"false")
        +",\"forecast_is_elapsed_time_guarantee\":false,\"accuracy_gates_unchanged\":true}");
}
