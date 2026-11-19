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
    const bool observe=e.unknownRoundFeedback&&!e.go&&e.executionRound>=0&&RoundPhysicalGpus(e)==1;
    if(observe)e.FinishTimedCohort(Seconds());
    else{e.unknownRound={};e.recentCohortTimes.clear();}
    const auto original=estimate(proposed);const double budget=1800.*RoundPhysicalGpus(e);
    if(original.profiled&&original.seconds<=budget){e.unknownRound={};return;}
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
    double measuredSecondsPerRow=0;
    size_t empiricalLimit=1,backfilled=0;
    std::string timedRecipe;
    std::set<size_t> timedRows;
    if(!original.profiled&&e.local&&!rank.empty()){
        // Start with one complete cohort. Repeated short rounds on the same
        // optical recipe provide a bounded way out of an incomplete cost model.
        selected=current;active=1;
        if(observe&&rank.front().second.size()==1){
            auto rowKey=[&](size_t id){
                std::ostringstream key;
                key<<e.PhiForTheta(id)<<' '<<(e.theta[id]<=1?-1:e.theta[id]>=179?1:0);
                for(size_t l=0;l<current.size();++l){
                    key<<' '<<depths[l]<<' '<<current[l][id]<<' '<<step[l][id];
                    if(residualCounts)key<<' '<<residualCounts->at(l).at(id);
                }
                return key.str();
            };
            const auto first=rowKey(rank.front().second.front());
            std::ostringstream signature;signature<<e.ExecutionHash()<<' '<<e.target<<' '<<separatePhi<<' '<<first;
            for(size_t id:ids)signature<<' '<<id;
            signature<<" /";for(size_t id:knots)signature<<' '<<id;
            signature<<" resources "<<e.kernel.block<<' '<<e.kernel.pipeline<<' '<<e.kernel.warp<<' '<<e.controls
                     <<' '<<e.sharedMeans<<' '<<e.mirror<<' '<<e.dispatchThreads<<' '<<e.orientationBlock;
            if(residualCounts)for(const auto& level:*residualCounts){
                signature<<" residual";for(size_t id:ids)signature<<' '<<level.at(id);
            }
            timedRecipe=Hash(signature.str());
            const auto measured=e.recentCohortTimes.find(timedRecipe);
            if(measured!=e.recentCohortTimes.end()){
                const auto& t=measured->second;measuredSecondsPerRow=t.latest;
                if(t.previous>0){
                    const double unit=std::max(t.previous,t.latest);
                    empiricalLimit=size_t(std::max(1.,std::min(double(std::min(size_t(32),2*t.rows)),budget/(8*unit))));
                }
            }
            if(empiricalLimit==1){
                timedRows.insert(rank.front().second.front());
            }else{
                decltype(rank) admitted,deferred;
                admitted.reserve(empiricalLimit);deferred.reserve(rank.size());
                for(auto& item:rank){
                    const auto& group=item.second;
                    if(admitted.size()<empiricalLimit&&group.size()==1&&rowKey(group.front())==first
                            &&timedRows.insert(group.front()).second){
                        if(!deferred.empty())++backfilled;
                        admitted.push_back(std::move(item));
                    }else deferred.push_back(std::move(item));
                }
                active=admitted.size();
                admitted.insert(admitted.end(),std::make_move_iterator(deferred.begin()),std::make_move_iterator(deferred.end()));
                rank.swap(admitted);
            }
        }
        for(size_t k=0;k<active;++k)for(size_t id:rank[k].second)
            for(size_t l=0;l<current.size();++l)selected[l][id]=step[l][id];
    }else e.unknownRound={};
    // An expensive unpriced future upgrade need not serialize already priced
    // dyadic steps. Try the full nearest step first, then complete cohorts in
    // priority order. Admit only fully priced work, with time headroom.
    size_t pricedProbes=0,pricedAdmissions=0;
    if(observe&&e.pricedSubround&&!original.profiled&&active<rank.size()&&!rank.empty()){
        const double limit=budget*.5;
        auto test=[&](const Budgets& proposal){
            ++pricedProbes;
            const auto value=estimate(proposal);
            return value.profiled&&value.seconds<=limit;
        };
        if(test(step)){
            selected=step;pricedAdmissions=rank.size()-active;active=rank.size();
        }else{
            // Non-monotonic costs are possible when two sets share traces.
            // Failure only stops speculation; it never changes a known estimate.
            const size_t ceiling=std::min(rank.size(),size_t(32));
            size_t upper=ceiling;
            bool failed=false;
            while(active<upper){
                const size_t next=failed?active+(upper-active+1)/2:std::min(upper,2*active);
                auto trial=selected;
                for(size_t k=active;k<next;++k)for(size_t id:rank[k].second)
                    for(size_t l=0;l<current.size();++l)trial[l][id]=step[l][id];
                if(!test(trial)){failed=true;upper=next-1;continue;}
                selected=std::move(trial);pricedAdmissions+=next-active;active=next;
            }
        }
        // These rounds have measured component profiles. They must not be
        // recorded as observations of the empirically timed singleton recipe.
        if(pricedAdmissions)timedRows.clear();
    }
    auto forecast=estimate(selected);
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
    size_t growing=0;for(size_t id:ids){bool grows=false;for(size_t l=0;l<current.size();++l)grows=grows||proposed[l][id]>current[l][id];growing+=grows;}
    if(!timedRows.empty()&&growing==timedRows.size()){
        bool same=true;
        for(size_t id:timedRows)for(size_t l=0;l<current.size();++l)same=same&&proposed[l][id]==step[l][id];
        if(same){
            e.unknownRound.recipe=timedRecipe;e.unknownRound.round=e.executionRound;
            e.unknownRound.rows=growing;e.unknownRound.started=Seconds();
            e.unknownRound.workerSeconds=e.nativeWorkerSeconds;e.unknownRound.nativeBlocks=e.nativeCommitted;
        }
    }
    e.Record(e.root+"/round_work_limit.json","{\"policy\":\"complete eight-seed prefix and output/residual cohort scheduling v7\",\"applied\":"+std::string(accepted?"true":"false")
        +",\"observed_previous_seconds_per_row\":"+Number(measuredSecondsPerRow)+",\"priced_subround_probes\":"+Number(pricedProbes)+",\"priced_cohort_admissions\":"+Number(pricedAdmissions)+",\"empirical_cohort_limit\":"+Number(empiricalLimit)+",\"timed_cohorts_backfilled\":"+Number(backfilled)
        +",\"unknown_cost_fallback\":"+(original.profiled?"false":"true")
        +",\"original_forecast_profiled\":"+(original.profiled?"true":"false")
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
