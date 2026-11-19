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
        const std::vector<RoundOutputPriority>& outputPriorities={},double costFactor=1.,double* selectedRawForecast=nullptr){
    const bool separatePhi=residualCounts&&!knots.empty()&&MixedResidualPhi(e);
    WorkForecast residualForecast;
    if(separatePhi){ResidualPhiScope scope(e,true);AdditionalWork residualWork(e,e.seeds.data());
        residualForecast=residualWork.Forecast(AuditRequests(e,depths,*residualCounts,ids,knots));}
    auto rawEstimate=[&](const Budgets& counts){
        if(separatePhi){auto result=work.Forecast(BudgetRequests(e,depths,counts,ids));
            result.seconds+=residualForecast.seconds;result.jobs+=residualForecast.jobs;
            result.profiled=result.profiled&&residualForecast.profiled;return result;}
        if(residualCounts&&!knots.empty()){
            auto requests=BudgetRequests(e,depths,counts,ids);
            const auto residual=AuditRequests(e,depths,*residualCounts,ids,knots);
            requests.insert(requests.end(),residual.begin(),residual.end());return work.Forecast(requests);}
        return work.Forecast(knots.empty()?BudgetRequests(e,depths,counts,ids):AuditRequests(e,depths,counts,ids,knots));};
    Require(costFactor>=.25&&costFactor<=1.,"invalid round work cost factor");
    auto estimate=[&](const Budgets& counts){auto r=rawEstimate(counts);r.seconds*=costFactor;return r;};
    if(selectedRawForecast)*selectedRawForecast=0;
    const auto original=estimate(proposed);const double budget=1800.*RoundPhysicalGpus(e);
    if(original.profiled&&original.seconds<=budget){
        if(selectedRawForecast)*selectedRawForecast=original.seconds/costFactor;
        return;
    }
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
    std::sort(rank.rbegin(),rank.rend());
    auto selected=step;size_t active=rank.size();
    // The 2x step can use measured native recipes even when the larger
    // proposal requires a new cutoff-correction profile. Test the work we
    // will actually dispatch before reducing it to a single cohort.
    auto forecast=estimate(selected);
    const bool initialStepProfiled=forecast.profiled;
    const bool unknownCostFallback=!forecast.profiled&&e.local&&!rank.empty();
    if(unknownCostFallback){
        // Missing measurements are not zero work. Advance only the highest
        // priority complete cohort until the real worker timings are known.
        // Keep its prescribed signed stencil and all eight seed prefixes.
        selected=current;active=1;
        for(size_t id:rank.front().second)for(size_t l=0;l<current.size();++l)selected[l][id]=step[l][id];
        forecast=estimate(selected);
    }
    size_t knownRecipeForecasts=0,knownRecipeCohorts=0;
    // A missing recipe for one future correction must not serialize all the
    // already measured 2x work. Group complete cohorts by their sampling
    // transitions, then admit a group only after an exact union forecast.
    // Cache coverage may differ inside a group: the representative is never
    // used to assert that the whole group is measured or cheap.
    if(unknownCostFallback&&forecast.profiled&&forecast.seconds<=budget){
        std::map<std::vector<int>,std::vector<size_t>> groups;
        std::vector<std::vector<int>> order;
        for(size_t k=0;k<rank.size();++k){std::vector<int> key;bool grows=false;
            for(size_t id:rank[k].second){key.push_back(e.PhiForTheta(id));
                for(size_t l=0;l<current.size();++l){key.push_back(selected[l][id]);key.push_back(step[l][id]);
                    grows=grows||step[l][id]>selected[l][id];}}
            if(!grows)continue;
            if(!groups.count(key))order.push_back(key);groups[key].push_back(k);
        }
        for(const auto& key:order){if(knownRecipeForecasts>=8)break;
            auto candidate=selected;
            for(size_t k:groups[key])for(size_t id:rank[k].second)for(size_t l=0;l<current.size();++l)
                candidate[l][id]=std::max(candidate[l][id],step[l][id]);
            const auto next=estimate(candidate);++knownRecipeForecasts;
            if(next.profiled&&next.seconds<=budget){selected=std::move(candidate);forecast=next;
                active+=groups[key].size();knownRecipeCohorts+=groups[key].size();}
        }
    }
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
    Budgets packingStep=step;
    bool boundedPrefixStep=false;size_t prefixStepForecasts=0,sharingForecasts=0,sharingCohorts=0;
    std::ostringstream sharingDecisions;sharingDecisions<<'[';
    // A whole output/stencil cohort can prescribe 4x (or more) samples. Merely
    // reducing the number of cohorts leaves even a single angle above the
    // latency budget. First shorten its growth to one complete dyadic prefix
    // step, then spend any remaining budget on other complete ranked cohorts.
    // All eight replicas and every accuracy test are unchanged. Do not touch
    // the unknown-cost fallback or split a 2x step into non-dyadic samples.
    if(e.local&&original.profiled&&forecast.profiled&&forecast.seconds>budget&&!rank.empty()){
        Budgets shorter=step;
        for(size_t id:ids)for(size_t l=0;l<current.size();++l)
            shorter[l][id]=std::min(step[l][id],Doubled(e,current[l][id]));
        if(shorter!=step){
            auto prefix=[&](size_t count){Budgets candidate=current;
                for(size_t k=0;k<count;++k)for(size_t id:rank[k].second)
                    for(size_t l=0;l<current.size();++l)candidate[l][id]=shorter[l][id];
                return candidate;};
            auto candidate=prefix(active);auto next=estimate(candidate);++prefixStepForecasts;
            if(next.profiled&&next.seconds<forecast.seconds){
                selected=candidate;forecast=next;boundedPrefixStep=true;
                // Bounded binary search only changes scheduling. Recheck each
                // proposed width with the real cache/corrected-native forecast;
                // never infer GPU seconds from angle-orientation counts.
                size_t low=active,high=rank.size()+1;
                while(forecast.seconds<=budget&&low+1<high&&prefixStepForecasts<8){
                    const size_t width=low+(high-low)/2;
                    candidate=prefix(width);next=estimate(candidate);++prefixStepForecasts;
                    if(next.profiled&&next.seconds<=budget){
                        low=width;selected=candidate;forecast=next;
                    }else high=width;
                }
                active=low;
                packingStep=shorter;
            }
        }
    }
    // Fill spare capacity for ordinary 2x steps too, not only shortened 4x steps.
    if(e.local&&original.profiled&&forecast.profiled&&forecast.seconds<=budget&&!rank.empty()){
        // A costly next-ranked cohort must not hide a compatible
        // later one. Prefer replacing singleton endpoint padding by
        // useful rows, then sharing already scheduled Sobol offsets.
        // The cheap event sweep ranks candidates only; the exact
        // native/cached-work forecast must still admit every cohort.
        using Events=std::map<int,std::pair<int,int>>;
        using SamplingGroup=std::pair<size_t,int>; // level and phi
        std::set<size_t> tried;
        for(;forecast.seconds<=budget&&sharingForecasts<8;){
            std::map<SamplingGroup,Events> baseEvents;
            for(size_t l=0;l<current.size();++l)for(size_t id:ids)if(selected[l][id]>current[l][id]){
                auto& events=baseEvents[{l,e.PhiForTheta(id)}];
                ++events[current[l][id]].first;--events[selected[l][id]].first;}
            size_t best=rank.size();double bestCost=0,bestShared=0,bestTotal=0,bestFraction=0;
            for(size_t k=0;k<rank.size();++k){if(tried.count(k))continue;
                std::map<SamplingGroup,Events> candidateEvents;
                for(size_t id:rank[k].second)for(size_t l=0;l<current.size();++l){
                    const int begin=selected[l][id],end=packingStep[l][id];if(end<=begin)continue;
                    const SamplingGroup key{l,e.PhiForTheta(id)};auto at=candidateEvents.find(key);
                    if(at==candidateEvents.end())at=candidateEvents.emplace(key,baseEvents[key]).first;
                    ++at->second[begin].second;--at->second[end].second;
                }
                double total=0,shared=0,extra=0;
                for(const auto& group:candidateEvents){int present=0,added=0,last=0;
                    for(const auto& point:group.second){const double span=double(point.first)-last;
                        if(added){total+=span*added;if(present)shared+=span*added;
                            extra+=span*(std::max(2,present+added)-(present?std::max(2,present):0));}
                        present+=point.second.first;added+=point.second.second;last=point.first;}
                }
                if(!total||!shared)continue;
                const double cost=extra/total,fraction=shared/total;
                // Iteration is in decreasing error order: equal
                // sharing retains the original statistical priority.
                if(best==rank.size()||cost<bestCost||(cost==bestCost&&fraction>bestFraction)){
                    best=k;bestCost=cost;bestShared=shared;bestTotal=total;bestFraction=fraction;}
            }
            if(best==rank.size())break;tried.insert(best);auto candidate=selected;
            for(size_t id:rank[best].second)for(size_t l=0;l<current.size();++l)
                candidate[l][id]=std::max(candidate[l][id],packingStep[l][id]);
            const auto next=estimate(candidate);const bool admit=next.profiled&&next.seconds<=budget;
            if(sharingForecasts++)sharingDecisions<<',';
            sharingDecisions<<"{\"cohort_rank\":"<<best<<",\"error_priority\":"<<Number(rank[best].first)
                <<",\"shared_orientation_fraction\":"<<Number(bestFraction)
                <<",\"shared_orientations\":"<<Number(bestShared)<<",\"additional_orientations\":"<<Number(bestTotal)
                <<",\"added_padded_row_work_ratio\":"<<Number(bestCost)
                <<",\"forecast_worker_seconds\":"<<(next.profiled?Number(next.seconds):"null")
                <<",\"accepted\":"<<(admit?"true":"false")<<'}';
            if(admit){selected=candidate;forecast=next;++active;++sharingCohorts;}
        }
    }
    sharingDecisions<<']';
    const bool accepted=selected!=current&&(original.profiled
        ?forecast.profiled&&forecast.seconds<original.seconds:selected!=wanted);
    if(accepted)proposed=selected;
    const auto& chosen=accepted?forecast:original;
    if(selectedRawForecast&&chosen.profiled)*selectedRawForecast=chosen.seconds/costFactor;
    size_t growing=0;for(size_t id:ids){bool grows=false;for(size_t l=0;l<current.size();++l)grows=grows||proposed[l][id]>current[l][id];growing+=grows;}
    e.Record(e.root+"/round_work_limit.json","{\"policy\":\"complete eight-seed prefix and output/residual cohort scheduling v7\",\"applied\":"+std::string(accepted?"true":"false")
        +",\"unknown_cost_fallback\":"+(unknownCostFallback?"true":"false")
        +",\"worker_cost_factor\":"+Number(costFactor)+",\"original_forecast_profiled\":"+(original.profiled?"true":"false")
        +",\"shorter_step_forecast_profiled\":"+(initialStepProfiled?"true":"false")
        +",\"known_recipe_group_forecasts\":"+Number(knownRecipeForecasts)+",\"known_recipe_group_cohorts\":"+Number(knownRecipeCohorts)
        +",\"selected_forecast_profiled\":"+((accepted?forecast.profiled:original.profiled)?"true":"false")
        +",\"no_growing_cohorts\":"+(rank.empty()?"true":"false")
        +",\"cost_model\":\"native physical or balanced children; conservative three-leg correction (five legacy paired legs); no assembly views\""
        +",\"output_priority_cohorts\":"+Number(outputPriorities.size())
        +",\"forecast_includes_common_prefix_interpolation\":"+(knots.empty()?"false":"true")+",\"physical_gpus\":"+Number(RoundPhysicalGpus(e))
        +",\"original_worker_seconds\":"+(original.profiled?Number(original.seconds):"null")+",\"selected_worker_seconds\":"+((accepted?forecast.profiled:original.profiled)?Number(accepted?forecast.seconds:original.seconds):"null")
        +",\"advisory_worker_seconds\":"+Number(budget)+",\"selected_growing_rows\":"+Number(growing)+",\"selected_cohorts\":"+Number(active)
        +",\"bounded_prefix_step_applied\":"+(boundedPrefixStep?"true":"false")+",\"prefix_step_forecasts\":"+Number(prefixStepForecasts)
        +",\"sharing_backfill_forecasts\":"+Number(sharingForecasts)+",\"sharing_backfill_cohorts\":"+Number(sharingCohorts)
        +",\"sharing_backfill_admissions\":"+sharingDecisions.str()
        +",\"shared_trace_efficiency_preserved\":"+(boundedPrefixStep?"null":(preservedSharing?"true":"false"))
        +",\"partition_overhead_scope\":"+Quoted(boundedPrefixStep?"before_prefix_step":"selected_partition")
        +",\"final_partition_overhead_fraction\":"+(boundedPrefixStep?"null":Number(partitionOverhead))
        +",\"partition_cost_model\":\"selected plus complementary cached-work forecast\",\"partition_overhead_fraction\":"+Number(partitionOverhead)
        +",\"allowed_partition_overhead_fraction\":"+Number(overheadLimit)
        +",\"within_forecast_budget\":"+((accepted?forecast.profiled:original.profiled)&&(accepted?forecast.seconds:original.seconds)<=budget?"true":"false")
        +",\"forecast_is_elapsed_time_guarantee\":false,\"accuracy_gates_unchanged\":true}");
}
