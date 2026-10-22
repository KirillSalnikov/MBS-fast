// Split expensive proposals at complete eight-seed prefixes. Only scheduling
// changes: every matrix/ratio/interpolation gate still runs on the full grid.
void BoundRoundWork(Engine& e,AdditionalWork& work,const Budgets& current,Budgets& proposed,
        const std::vector<int>& depths,const std::vector<size_t>& ids,const Statistics& stats){
    const auto original=work.Forecast(BudgetRequests(e,depths,proposed,ids));
    const double budget=1800.*std::max(size_t(1),e.gpus.size());
    if(!original.profiled||original.seconds<=budget||original.jobs<256)return;
    const Budgets wanted=proposed;Budgets step=current;
    std::vector<std::pair<double,size_t>> rank;
    for(size_t r=0;r<ids.size();++r){const size_t id=ids[r];bool grows=false;
        for(size_t l=0;l<current.size();++l){step[l][id]=std::min(wanted[l][id],Doubled(e,current[l][id]));grows=grows||step[l][id]>current[l][id];}
        if(grows)rank.push_back({WorstMetric(stats,r,e.ratios).error,id});}
    std::sort(rank.rbegin(),rank.rend());if(rank.empty())return;
    auto selected=step;auto forecast=work.Forecast(BudgetRequests(e,depths,selected,ids));size_t active=rank.size();
    const double fullSeconds=forecast.seconds;bool preservedSharing=false;double partitionOverhead=0;
    // At most three additional cost probes; never create incomplete seed groups. A
    // single indivisible angular prefix may exceed this advisory budget.
    for(int attempt=0;e.local&&forecast.profiled&&forecast.seconds>budget&&active>1&&attempt<3;++attempt){
        const size_t nextActive=std::max(size_t(1),size_t(std::ceil(active*std::max(.05,std::min(.5,budget/forecast.seconds)))));Budgets candidate=current;
        for(size_t k=0;k<nextActive;++k)for(size_t l=0;l<current.size();++l)candidate[l][rank[k].second]=step[l][rank[k].second];
        const auto next=work.Forecast(BudgetRequests(e,depths,candidate,ids));
        // Compare the cost of the entire partition, not seconds per selected
        // angle-orientation. Difficult rows can be more expensive and less
        // covered by cached blocks; their cost per unit is not repeated work.
        // The complementary request measures the actual lost trace sharing.
        Budgets remainder=step;
        for(size_t k=0;k<nextActive;++k)for(size_t l=0;l<current.size();++l)
            remainder[l][rank[k].second]=current[l][rank[k].second];
        const auto rest=work.Forecast(BudgetRequests(e,depths,remainder,ids));
        const double splitSeconds=next.seconds+rest.seconds;
        partitionOverhead=fullSeconds>0?std::max(0.,splitSeconds/fullSeconds-1):0;
        if(!next.profiled||!rest.profiled||splitSeconds>1.10*fullSeconds){preservedSharing=true;break;}
        selected=candidate;active=nextActive;forecast=next;
    }
    const bool accepted=forecast.profiled&&forecast.seconds<original.seconds&&selected!=current;
    if(accepted)proposed=selected;
    e.Record(e.root+"/round_work_limit.json","{\"policy\":\"complete eight-seed prefix scheduling\",\"applied\":"+std::string(accepted?"true":"false")
        +",\"original_worker_seconds\":"+Number(original.seconds)+",\"selected_worker_seconds\":"+Number(forecast.seconds)
        +",\"advisory_worker_seconds\":"+Number(budget)+",\"selected_growing_rows\":"+Number(active)
        +",\"shared_trace_efficiency_preserved\":"+(preservedSharing?"true":"false")
        +",\"partition_cost_model\":\"selected plus complementary cached-work forecast\",\"partition_overhead_fraction\":"+Number(partitionOverhead)
        +",\"within_forecast_budget\":"+(forecast.profiled&&forecast.seconds<=budget?"true":"false")
        +",\"forecast_is_elapsed_time_guarantee\":false,\"accuracy_gates_unchanged\":true}");
}
