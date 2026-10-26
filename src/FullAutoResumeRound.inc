// Published diagnostics only propose an allocation. They never certify a
// matrix: the ordinary source validation and all accuracy gates run again.
bool PublishedRoundBudget(Engine& e,const AdaptiveState& state,Budgets& counts,
        Statistics& priorities,std::string& folder,bool& passed){
    const std::string path=e.root+"/intermediate_current.json";if(!Exists(path))return false;
    const auto text=CompactJson(Read(path));const std::string key="\"snapshot_folder\":\"";
    const auto at=text.find(key);if(at==std::string::npos)return false;
    const auto begin=at+key.size(),end=text.find('"',begin);if(end==std::string::npos)return false;
    folder=text.substr(begin,end-begin);
    // Escaped/nonlocal paths are optional scheduling hints we decline to use.
    if(folder.find('\\')!=std::string::npos||folder.find("..")!=std::string::npos
            ||folder.compare(0,e.root.size()+17,e.root+"/round_snapshots/")!=0)return false;
    if(!Exists(folder+"/parameters_accuracy.csv")||!Exists(folder+"/theta_support.csv"))return false;
    if(ReadGrid(folder+"/theta_support.csv")!=Angles(e.theta,state.knots))return false;
    const auto ids=AdaptiveGuards(state.knots,e.theta.size());counts=state.counts;
    priorities.interval.resize(ids.size());priorities.ratioInterval.resize(ids.size());
    std::istringstream file(Read(folder+"/parameters_accuracy.csv"));std::string line;
    if(!std::getline(file,line))return false;
    std::ostringstream expected;expected<<"theta_deg,is_support,interpolation_residual95,previous_available,half_prefix_stable";
    for(size_t l=0;l<counts.size();++l)expected<<",N_level"<<l;
    for(int j=0;j<16;++j)expected<<",estimated95_M"<<j/4+1<<j%4+1<<"_scaled_M11";
    for(int j=1;j<16;++j)expected<<",estimated95_M"<<j/4+1<<j%4+1<<"_over_M11";
    expected<<",phi_points";
    if(line!=expected.str())return false;
    size_t row=0;while(std::getline(file,line)){if(row>=ids.size())return false;
        std::replace(line.begin(),line.end(),',',' ');std::istringstream values(line);
        double theta,residual;int support,valid,coarse,phi;std::string extra;const size_t id=ids[row];
        if(!(values>>theta>>support>>residual>>valid>>coarse)||theta!=e.theta[id]
                ||support!=int(std::binary_search(state.knots.begin(),state.knots.end(),id)))return false;
        for(auto& level:counts){int n;if(!(values>>n)||n<=0||n>level[id]||(n&(n-1)))return false;level[id]=n;}
        for(double& error:priorities.interval[row])if(!(values>>error)||!std::isfinite(error)||error<0)return false;
        for(int j=1;j<16;++j){double& error=priorities.ratioInterval[row][j];if(!(values>>error)||!std::isfinite(error)||error<0)return false;}
        if(!(values>>phi)||phi!=e.PhiForTheta(id)||(values>>extra))return false;
        ++row;
    }
    passed=text.find("\"numerical_round_passed\":true")!=std::string::npos;
    return row==ids.size();
}
void RevisitPublishedPass(Engine& e,AdaptiveState& state){
    if(e.finalChecks||e.stablePasses!=1||state.streak<1)return;
    Budgets counts;Statistics priorities;std::string folder;bool passed=false;
    if(!PublishedRoundBudget(e,state,counts,priorities,folder,passed)||!passed||counts==state.counts)return;
    // This is the last already-passing production estimator, not a search for
    // an older favourable prefix. Its matrices and errors are not loaded here.
    const auto old=state.counts;state.counts=counts;
    e.Record(e.root+"/resume_round_work.json","{\"policy\":\"revisit last passing production allocation\",\"snapshot\":"+Quoted(folder)
        +",\"required_production_passes\":1,\"native_blocks_deleted\":0,\"accuracy_accepted_from_diagnostics\":false,\"all_accuracy_gates_recomputed\":true}");
    e.Log("resume_published_allocation","{\"previous_peak\":"+Number(PeakBudgets(old)[0])+",\"revisited_peak\":"+Number(PeakBudgets(counts)[0])
        +",\"saved_blocks_retained\":true,\"convergence_decision_pending\":true}");
}
// A scheduling hint may span changed support. Values are never accepted from
// diagnostics; new/changed-phi rows start at the ordinary initial allocation.
bool PublishedSchedulingBudget(Engine& e,const AdaptiveState& state,Budgets& counts,
        Statistics& priorities,std::string& folder,bool& passed,bool& residualPending){
    const auto path=e.root+"/intermediate_current.json";if(!Exists(path))return false;
    const auto text=CompactJson(Read(path));const std::string key="\"snapshot_folder\":\"";const auto at=text.find(key);
    if(at==std::string::npos)return false;
    const auto begin=at+key.size(),end=text.find('"',begin);
    if(end==std::string::npos)return false;
    folder=text.substr(begin,end-begin);
    if(folder.find('\\')!=std::string::npos||folder.find("..")!=std::string::npos
            ||folder.compare(0,e.root.size()+17,e.root+"/round_snapshots/")!=0)return false;
    if(!Exists(folder+"/parameters_accuracy.csv"))return false;
    struct Row {std::vector<int> counts;Mueller interval{},ratio{};int phi=0;};std::map<size_t,Row> rows;
    std::istringstream file(Read(folder+"/parameters_accuracy.csv"));std::string line;
    if(!std::getline(file,line))return false;
    std::ostringstream expected;expected<<"theta_deg,is_support,interpolation_residual95,previous_available,half_prefix_stable";
    for(size_t l=0;l<state.counts.size();++l)expected<<",N_level"<<l;
    for(int j=0;j<16;++j)expected<<",estimated95_M"<<j/4+1<<j%4+1<<"_scaled_M11";
    for(int j=1;j<16;++j)expected<<",estimated95_M"<<j/4+1<<j%4+1<<"_over_M11";
    expected<<",phi_points";
    if(line!=expected.str())return false;
    while(std::getline(file,line)){
        std::replace(line.begin(),line.end(),',',' ');std::istringstream values(line);double theta,residual;int support,valid,coarse;std::string extra;
        if(!(values>>theta>>support>>residual>>valid>>coarse)||!std::isfinite(theta)||!std::isfinite(residual)||residual<0)return false;
        const auto q=std::lower_bound(e.theta.begin(),e.theta.end(),theta);if(q==e.theta.end()||*q!=theta)return false;
        const size_t id=q-e.theta.begin();Row row;row.counts.resize(state.counts.size());
        for(int& n:row.counts)if(!(values>>n)||n<=0||n>e.maxCount||(n&(n-1)))return false;
        for(double& v:row.interval)if(!(values>>v)||!std::isfinite(v)||v<0)return false;
        for(int j=1;j<16;++j)if(!(values>>row.ratio[j])||!std::isfinite(row.ratio[j])||row.ratio[j]<0)return false;
        if(!(values>>row.phi)||row.phi<=0||row.phi>e.maxPhi||(values>>extra)||!rows.emplace(id,row).second)return false;
    }
    if(rows.empty())return false;
    const auto ids=AdaptiveGuards(state.knots,e.theta.size());
    counts=state.counts;
    priorities.interval.resize(ids.size());priorities.ratioInterval.resize(ids.size());size_t newRows=0;
    for(size_t r=0;r<ids.size();++r){const size_t id=ids[r];auto q=rows.lower_bound(id);if(q==rows.end())--q;
        const bool same=q->first==id&&q->second.phi==e.PhiForTheta(id);newRows+=!same;
        priorities.interval[r]=q->second.interval;priorities.ratioInterval[r]=q->second.ratio;
        for(size_t l=0;l<counts.size();++l)counts[l][id]=std::min(state.counts[l][id],same?q->second.counts[l]:(l?e.minimum:e.initial));
    }
    passed=text.find("\"numerical_round_passed\":true")!=std::string::npos;residualPending=false;
    const std::string residualKey="\"interpolation_residual95\":";const auto rp=text.find(residualKey);
    if(rp!=std::string::npos){std::istringstream value(text.substr(rp+residualKey.size()));double error;
        if(value>>error)residualPending=std::isfinite(error)&&error>e.interpolation;}
    e.Log("resume_scheduling_snapshot","{\"new_or_changed_phi_rows\":"+Number(newRows)+",\"acceptance_from_diagnostics\":false,\"changed_support_allowed_for_scheduling_only\":true}");
    return true;
}
void BoundResumedRound(Engine& e,AdaptiveState& state,const std::vector<int>& depths){
    const auto correction=ReadCorrectionProbe(e,state);
    if(correction.attempt&&!correction.complete){
        // This small in-flight probe is keyed to the exact production sources.
        // Finish it before rescheduling, rather than invalidate its receipt.
        e.Log("resume_round_bound_deferred","{\"reason\":\"finish pending correction subround\",\"round\":"+Number(state.round)+"}");return;}
    Budgets baseline;Statistics priorities;std::string folder;bool passed=false;
    bool residualPending=false;
    if(!PublishedSchedulingBudget(e,state,baseline,priorities,folder,passed,residualPending)||passed)return;
    const auto ids=AdaptiveGuards(state.knots,e.theta.size());
    // Keep every complete eight-seed prefix already produced in the interrupted
    // round. Incomplete blocks also remain in the cache for subsequent rounds.
    const auto completed=CachedCompletedBudgets(e,depths,baseline,ids);
    for(size_t l=0;l<baseline.size();++l)for(size_t id:ids)
        baseline[l][id]=std::max(baseline[l][id],std::min(state.counts[l][id],completed[l][id]));
    for(const auto& row:correction.rows)baseline[row.level][row.id]=std::max(baseline[row.level][row.id],row.after);
    auto selected=state.counts;
    const bool residualFirst=residualPending&&!state.interpolationCounts.empty();
    if(residualFirst)selected=baseline; // Pending residual growth remains in its own saved allocation.
    else {AdditionalWork work(e,e.seeds.data());
        BoundRoundWork(e,work,baseline,selected,depths,ids,priorities,state.knots,{},state.interpolationCounts.empty()?nullptr:&state.interpolationCounts);}
    if(selected==state.counts)return;
    state.counts=selected;
    e.Record(e.root+"/resume_round_work.json","{\"policy\":\"bound unfinished round using last published allocation and complete cached prefixes\",\"snapshot\":"+Quoted(folder)
        +",\"residual_first\":"+(residualFirst?"true":"false")+",\"native_blocks_deleted\":0,\"complete_prefixes_retained\":true,\"accuracy_accepted_from_diagnostics\":false,\"all_accuracy_gates_recomputed\":true}");
    e.Log("resume_round_bounded","{\"snapshot\":"+Quoted(folder)+",\"saved_blocks_retained\":true,\"convergence_decision_pending\":true}");
}

// The empirical N is a useful convergence forecast, but it can make the first
// production estimate take days after the angular pilot expands the grid.
// Freeze a cost-bounded first allocation before looking at any production
// values. Preserve explicit starts and every reusable complete seed prefix.
void BoundInitialProduction(Engine& e,AdaptiveState& state,const std::vector<int>& depths){
    if(state.round!=0||!e.local||e.args.IsCatched("fullauto_initial")||e.config.adaptive.loadedFromFile
            ||std::any_of(state.valid.begin(),state.valid.end(),[](bool value){return value;}))return;
    const auto ids=AdaptiveGuards(state.knots,e.theta.size());
    if(ids.empty()||state.counts.empty())return;
    const Budgets wanted=state.counts;
    AdditionalWork work(e,e.seeds.data());
    auto estimate=[&](const Budgets& counts){return work.Forecast(BudgetRequests(e,depths,counts,ids));};
    const double limit=1800.*RoundPhysicalGpus(e);const auto original=estimate(wanted);
    if(original.profiled&&original.seconds<=limit)return;
    Budgets selected=wanted;
    for(size_t id=0;id<e.theta.size();++id)selected[0][id]=std::min(wanted[0][id],512);
    // Materialized prefixes can vary by angle. No matrix values or confidence
    // intervals enter this metadata-only selection.
    const auto cached=CachedCompletedBudgets(e,depths,selected,ids);
    for(size_t l=0;l<selected.size();++l)for(size_t id:ids)
        selected[l][id]=std::max(selected[l][id],std::min(wanted[l][id],cached[l][id]));
    // During the very first prefetch, long prefixes may exist only as native
    // leaves. Recognize a complete eight-seed prefix before it is assembled.
    // At and above the saved block size, the coverage test is monotone.
    int commonMaximum=e.maxCount;for(size_t id:ids)commonMaximum=std::min(commonMaximum,wanted[0][id]);
    int completeBase=0;
    for(int n=Power2(std::max(512,e.orientationBlock));n<=commonMaximum;){
        const Budgets candidate(1,std::vector<int>(e.theta.size(),n));
        const auto remaining=work.Forecast(BudgetRequests(e,{depths[0]},candidate,ids));
        if(remaining.jobs)break;
        completeBase=n;
        for(size_t id:ids)selected[0][id]=std::max(selected[0][id],n);
        if(n>commonMaximum/2)break;
        n*=2;
    }
    auto chosen=estimate(selected);
    // Grow by full powers of two only when existing timings cover every leg.
    // Missing timings retain the small pilot; they never mean zero cost.
    while(chosen.profiled&&chosen.seconds<=limit){
        Budgets candidate=selected;bool changed=false;
        for(size_t id:ids){const int n=std::min(wanted[0][id],Doubled(e,selected[0][id]));
            changed=changed||n!=selected[0][id];candidate[0][id]=n;}
        if(!changed)break;
        const auto next=estimate(candidate);if(!next.profiled||next.seconds>limit)break;
        selected=std::move(candidate);chosen=next;
    }
    if(selected==wanted)return;
    state.counts=selected;
    e.Record(e.root+"/initial_production_work.json",std::string("{\"policy\":\"cost-bounded first production estimate v1\",\"round\":0")
        +",\"requested_base_peak\":"+Number(PeakBudgets(wanted)[0])+",\"selected_base_peak\":"+Number(PeakBudgets(selected)[0])
        +",\"complete_native_leaf_prefix\":"+Number(completeBase)+",\"forecast_profiled\":"+(chosen.profiled?"true":"false")
        +",\"forecast_worker_seconds\":"+Number(chosen.seconds)+",\"advisory_worker_seconds\":"+Number(limit)
        +",\"explicit_initial_unchanged\":true,\"correction_minima_preserved\":true,\"cached_blocks_retained\":true"
        +",\"allocation_reads_matrix_values\":false,\"accuracy_gates_unchanged\":true}");
}
