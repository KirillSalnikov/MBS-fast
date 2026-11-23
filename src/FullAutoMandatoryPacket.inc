// Schedule only work already requested by the level-stability check. A packet
// changes when a level is evaluated; the next round still checks every level.
bool SelectMandatoryPacket(Engine& e,AdditionalWork& work,const Budgets& current,
        const Budgets& wanted,const Budgets& required,const std::vector<int>& depths,
        const std::vector<size_t>& ids,const Statistics& stats,const std::vector<size_t>& knots,
        const Budgets* residual,Budgets& selected,double budget){
    std::vector<int> floors;
    for(size_t l=0;l<depths.size();++l){
        const auto minimum=PrefixMinimumCount(e,l,depths.size());
        if(minimum>e.maxCount)return false;
        floors.push_back(int(minimum));
    }
    long long baselineGain=0;
    for(size_t l=0;l<depths.size();++l)for(size_t id:ids){
        const int gain=selected[l][id]-current[l][id];
        if(gain<=0)continue;
        if(selected[l][id]>floors[l]||selected[l][id]>required[l][id])return false;
        baselineGain+=gain;
    }
    if(!baselineGain)return false;
    auto price=[&](const Budgets& counts){
        auto requests=BudgetRequests(e,depths,counts,ids);
        if(residual&&!knots.empty()){
            const auto extra=AuditRequests(e,depths,*residual,ids,knots);
            requests.insert(requests.end(),extra.begin(),extra.end());
        }
        for(size_t l=0;l<depths.size();++l){
            std::vector<size_t> checked;
            Budgets half(l?2:1,std::vector<int>(e.theta.size(),1));
            for(size_t id:ids)if(counts[l][id]>=floors[l]){
                checked.push_back(id);
                for(auto& level:half)level[id]=counts[l][id]/2;
            }
            if(checked.empty())continue;
            const std::vector<int> checkDepths=l?std::vector<int>{depths[l-1],depths[l]}:std::vector<int>{depths[l]};
            const auto extra=BudgetRequests(e,checkDepths,half,checked);
            requests.insert(requests.end(),extra.begin(),extra.end());
        }
        return work.Forecast(requests);
    };
    const auto baseline=price(selected);
    if(!baseline.profiled||baseline.seconds<=0)return false;
    const double ceiling=std::max(budget,baseline.seconds);
    double best=baselineGain/baseline.seconds;
    struct Packet {
        size_t level;int increase;std::vector<size_t> rows;
        size_t Gain()const{return std::min(rows.size(),size_t(128))*size_t(increase);}
    };
    std::vector<Packet> packets;
    std::vector<std::pair<double,size_t>> ranked;
    for(size_t r=0;r<ids.size();++r)ranked.push_back({WorstMetric(stats,r,e.ratios).error,ids[r]});
    std::stable_sort(ranked.begin(),ranked.end(),[](const std::pair<double,size_t>& a,const std::pair<double,size_t>& b){return a.first>b.first;});
    for(size_t l=0;l<depths.size();++l){
        std::map<std::array<int,6>,std::vector<size_t>> groups;
        for(const auto& row:ranked){
            const size_t id=row.second;
            const int count=std::min({wanted[l][id],required[l][id],floors[l]});
            if(count<=current[l][id])continue;
            const std::array<int,6> key={{e.PhiForTheta(id),current[l][id],count,residual?(*residual)[l][id]:0,
                l?current[l-1][id]:0,l&&residual?(*residual)[l-1][id]:0}};
            groups[key].push_back(id);
        }
        for(auto& entry:groups)packets.push_back({l,entry.first[2]-entry.first[1],std::move(entry.second)});
    }
    std::stable_sort(packets.begin(),packets.end(),[](const Packet& a,const Packet& b){return a.Gain()>b.Gain();});
    size_t selectedRows=0,selectedLevel=0,firstId=0;long long selectedGain=0;WorkForecast chosen;
    for(size_t g=0;g<std::min(size_t(6),packets.size());++g){
        const auto& packet=packets[g];
        // Price one larger group, then compare the original bounded candidates.
        for(size_t width=std::min(size_t(256),packet.rows.size());width;width=width>128?128:width/2){
            auto trial=current;
            for(size_t i=0;i<width;++i)trial[packet.level][packet.rows[i]]+=packet.increase;
            const auto cost=price(trial);
            if(!cost.profiled||cost.seconds>ceiling)continue;
            const long long gain=static_cast<long long>(width)*packet.increase;
            const double efficiency=gain/std::max(1.e-6,cost.seconds);
            if(efficiency>best){
                best=efficiency;chosen=cost;selected=std::move(trial);selectedRows=width;
                selectedLevel=packet.level;selectedGain=gain;firstId=packet.rows.front();
            }
            if(width<=128)break;
        }
    }
    if(!selectedRows)return false;
    e.Record(e.root+"/round_work_limit.json","{\"policy\":\"mandatory level prefix packet v1\",\"applied\":true"
        ",\"reason\":\"more required level samples per forecast worker second\",\"physical_gpus\":1"
        ",\"decision_after_round\":"+Number(e.executionRound)+",\"selected_forecast_profiled\":true"
        ",\"selected_growing_rows\":"+Number(selectedRows)+",\"selected_depth\":"+Number(depths[selectedLevel])
        +",\"count_before_per_seed\":"+Number(current[selectedLevel][firstId])+",\"count_after_per_seed\":"+Number(selected[selectedLevel][firstId])
        +",\"phi_points\":"+Number(e.PhiForTheta(firstId))+",\"selected_mandatory_samples_per_seed\":"+Number(selectedGain)+",\"baseline_mandatory_samples_per_seed\":"+Number(baselineGain)
        +",\"selected_worker_seconds\":"+Number(chosen.seconds)+",\"baseline_worker_seconds\":"+Number(baseline.seconds)
        +",\"advisory_worker_seconds\":"+Number(budget)+",\"within_forecast_budget\":"+(chosen.seconds<=budget?"true":"false")
        +",\"forecast_includes_level_prefix_checks\":true,\"forecast_is_elapsed_time_guarantee\":false,\"accuracy_gates_unchanged\":true}");
    return true;
}
