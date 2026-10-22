// Every depth component must be stable in the complete physical estimator.
// These are pointwise empirical prefix checks, not rigorous error bounds.
const char* OutputCachePolicy(){return "complete-compatible-prefixes-v1";}
const char* PrefixStabilityPolicy(){return "all-level-paired-prefix-v1";}
double LevelPrefixChange(const Mueller& total,const Mueller& current,const Mueller& alternate,bool ratios){
    Mueller other{};for(int j=0;j<16;++j){
        if(!std::isfinite(total[j])||!std::isfinite(current[j])||!std::isfinite(alternate[j]))return 1e100;
        other[j]=total[j]-current[j]+alternate[j];if(!std::isfinite(other[j]))return 1e100;
    }
    if(!(total[0]>0&&other[0]>0))return 1e100;
    double change=0;for(int j=0;j<16;++j){
        change=std::max(change,std::fabs(other[j]-total[j])/total[0]);
        if(ratios&&j)change=std::max(change,std::fabs(other[j]/other[0]-total[j]/total[0]));
    }return change;
}
long long PrefixMinimumCount(const Engine& e,size_t level){
    if(!level)return 2LL*e.initial;
    // Do not accept a cheap correction at its first apparently quiet pilot.
    // Four doublings of the correction pilot are capped by the base startup
    // check; this avoids imposing a huge base count on every correction.
    return std::max(2LL*e.minimum,std::min(2LL*e.initial,16LL*e.minimum));
}
struct LevelPrefixDecision {
    Budgets required;std::vector<std::vector<double>> errors;bool passed=true;
    size_t failed=0,failedJoint=0;std::string sources,table,jointTable;
};
double PairedLevelPrefixChange(const Engine& e,const Samples& total,const Samples& current,const Samples& alternate,size_t row){
    Require(total.size()==8&&current.size()==8&&alternate.size()==8,"prefix check needs eight complete replicas");
    Samples reference(8,Matrix(1)),candidate(8,Matrix(1));
    for(size_t s=0;s<8;++s){Require(row<total[s].size()&&row<current[s].size()&&row<alternate[s].size(),"prefix check row missing");
        reference[s][0]=total[s][row];for(int j=0;j<16;++j){
            const double a=current[s][row][j],b=alternate[s][row][j],t=reference[s][0][j];
            if(!std::isfinite(a)||!std::isfinite(b)||!std::isfinite(t))return 1e100;
            candidate[s][0][j]=t-a+b;if(!std::isfinite(candidate[s][0][j]))return 1e100;}}
    const Mueller mean=Mean(reference)[0];
    const double shift=LevelPrefixChange(mean,mean,Mean(candidate)[0],e.ratios);
    if(shift>=1e99)return shift;
    return std::max(shift,e.PairedBias(candidate,reference));
}
LevelPrefixDecision CheckLevelPrefixes(Engine& e,const std::vector<int>& depths,AdaptiveState& state,
        const AdaptiveData& stage,const std::vector<size_t>& ids){
    LevelPrefixDecision result;result.required=state.counts;
    result.errors.assign(depths.size(),std::vector<double>(ids.size(),1e100));
    const Samples total=Combine(stage.levels);Samples joint=total;const double tolerance=NumericalBudget(e);
    std::fill(state.coarse.begin(),state.coarse.end(),false);
    std::vector<bool> passed(ids.size(),true),jointAvailable(ids.size(),true);std::ostringstream table;
    table<<"theta_deg,level,depth,count_per_seed,reference_count_per_seed,minimum_stable_count,paired_change95,budget,passed,required_count,phi_points\n";
    for(size_t level=0;level<depths.size();++level){
        std::vector<size_t> checked;Budgets half(level?2:1,std::vector<int>(e.theta.size(),1));
        for(size_t id:ids)if(state.counts[level][id]>=PrefixMinimumCount(e,level)){
            checked.push_back(id);for(auto& n:half)n[id]=state.counts[level][id]/2;}
        if(!checked.empty()){
            const std::vector<int> checkDepths=level?std::vector<int>{depths[level-1],depths[level]}:std::vector<int>{depths[0]};
            // Persist complete replicas, never acceptance. The key includes
            // counts, phi, coefficients, execution and seed identities; restore
            // also validates every source graph before any value is reused.
            const auto key=InterpolationGroupKey(e,checkDepths,half,checked,e.seeds.data());AdaptiveData previous;
            if(!RestoreInterpolationGroup(e,key,checkDepths.size(),checked.size(),previous)){
                previous=EvaluateBudget(e,checkDepths,half,checked,e.seeds.data());SaveInterpolationGroup(e,key,previous);}
            result.sources+=previous.sources;Samples alternate=stage.levels[level];
            const auto& old=previous.levels.back();
            for(size_t k=0;k<checked.size();++k){const size_t row=std::lower_bound(ids.begin(),ids.end(),checked[k])-ids.begin();
                for(size_t s=0;s<8;++s){alternate[s][row]=old[s][k];
                    for(int j=0;j<16;++j)joint[s][row][j]+=old[s][k][j]-stage.levels[level][s][row][j];}
                result.errors[level][row]=PairedLevelPrefixChange(e,total,stage.levels[level],alternate,row);
            }
        }
        for(size_t row=0;row<ids.size();++row){const size_t id=ids[row];const int n=state.counts[level][id];
            jointAvailable[row]=jointAvailable[row]&&n>=PrefixMinimumCount(e,level);
            const bool ok=result.errors[level][row]<=tolerance;passed[row]=passed[row]&&ok;
            if(!ok){++result.failed;result.required[level][id]=Doubled(e,n);}
            table<<Number(e.theta[id])<<','<<level<<','<<depths[level]<<','<<n<<','<<(n>=PrefixMinimumCount(e,level)?n/2:0)
                <<','<<PrefixMinimumCount(e,level)<<','<<Number(result.errors[level][row])<<','<<Number(tolerance)<<','<<ok
                <<','<<result.required[level][id]<<','<<e.PhiForTheta(id)<<'\n';
        }
    }
    std::ostringstream jointRows;jointRows<<"theta_deg,all_half_prefixes_available,paired_change95,budget,passed,growth_level\n";
    for(size_t row=0;row<ids.size();++row){
        const double error=jointAvailable[row]?PairedLevelPrefixChange(e,total,total,joint,row):1e100;
        const bool ok=jointAvailable[row]&&error<=tolerance;int growth=-1;
        if(!ok){++result.failedJoint;
            // Individual failures already schedule their own levels. If only
            // their sum fails, grow the largest measured contribution first.
            if(passed[row])for(size_t l=0;l<depths.size();++l)if(state.counts[l][ids[row]]<e.maxCount
                    &&(growth<0||result.errors[l][row]>result.errors[size_t(growth)][row]))growth=int(l);
            if(growth>=0)result.required[size_t(growth)][ids[row]]=Doubled(e,state.counts[size_t(growth)][ids[row]]);
        }
        passed[row]=passed[row]&&ok;
        jointRows<<Number(e.theta[ids[row]])<<','<<jointAvailable[row]<<','<<Number(error)<<','<<Number(tolerance)<<','<<ok<<','<<growth<<'\n';
    }
    for(size_t row=0;row<ids.size();++row){state.coarse[ids[row]]=passed[row];result.passed=result.passed&&passed[row];}
    result.table=table.str();Atomic(e.root+"/level_prefix_stability.csv",result.table);
    result.jointTable=jointRows.str();Atomic(e.root+"/joint_prefix_stability.csv",result.jointTable);
    Atomic(e.root+"/level_prefix_sources.tsv",result.sources);
    e.Log("level_prefix_stability","{\"policy\":"+Quoted(PrefixStabilityPolicy())+",\"round\":"+Number(state.round)
        +",\"passed\":"+(result.passed?"true":"false")+",\"failed_level_rows\":"+Number(result.failed)
        +",\"failed_joint_rows\":"+Number(result.failedJoint)
        +",\"all16_and_full_matrix_ratios\":true,\"paired_eight_seed_noise\":true,\"recomputed_for_current_estimator\":true,\"tolerance\":"+Number(tolerance)+"}");
    return result;
}
void PreparePrefixStabilityPolicy(Engine& e,AdaptiveState& state,bool resumed){
    const auto path=e.root+"/prefix_stability_policy.json";
    const bool changed=!Exists(path)||Read(path).find(PrefixStabilityPolicy())==std::string::npos;
    if(changed&&resumed){state.streak=0;std::fill(state.coarse.begin(),state.coarse.end(),false);e.independentPassed=false;
        SaveAdaptive(e,state);e.Log("prefix_stability_policy_migration","{\"production_blocks_retained\":true,\"old_stability_acceptance_invalidated\":true}");}
    e.Record(path,"{\"policy\":"+Quoted(PrefixStabilityPolicy())+",\"all_depth_levels\":true,\"joint_level_sum\":true,\"recompute_each_round\":true,\"comparison\":\"full estimator with each level and then all levels replaced by half prefixes; paired eight-seed mean shift plus noise95\",\"correction_minimum\":\"max(2*minimum,min(2*initial,16*minimum))\",\"statistical_not_rigorous_bound\":true}");
}
void RequirePrefixGrowth(const LevelPrefixDecision& check,const std::vector<size_t>& ids,Budgets& proposed){
    for(size_t level=0;level<proposed.size();++level)for(size_t id:ids)proposed[level][id]=std::max(proposed[level][id],check.required[level][id]);
}
bool CheckCachedLevelPrefixes(Engine& e,const std::vector<int>& depths,AdaptiveState& state,
        AdaptiveData& stage,const std::vector<size_t>& ids,bool adoptComplete){
    // Select candidate counts solely from compatible complete coverage, before
    // inspecting any matrix. Production adopts complete larger prefixes in
    // this round; the legacy diagnostic mode grows only inconsistent levels.
    if(adoptComplete)e.Record(e.root+"/output_cache_policy.json",
        "{\"policy\":"+Quoted(OutputCachePolicy())+",\"round\":"+Number(state.round)
        +",\"selection\":\"complete compatible metadata before values\",\"replicas\":8,"
        "\"same_phi_and_paired_depths\":true,\"native_work_forbidden_during_promotion\":true,"
        "\"output_statistics_recomputed\":true,\"geometry_counts_frozen_before_output_promotion\":true}");
    const auto counts=CachedCompletedBudgets(e,depths,state.counts,ids);
    const std::string header="theta_deg,level,current_count,cached_count,paired_change95,joint_change95,budget,promoted,reason,phi_points\n";
    if(counts==state.counts){
        Atomic(e.root+"/cached_prefix_consistency.csv",header);Atomic(e.root+"/cached_prefix_consistency_sources.tsv","");
        e.Log("cached_prefix_consistency","{\"round\":"+Number(state.round)+",\"larger_complete_candidates\":false,\"promoted_level_rows\":0}");return false;}
    const auto key=InterpolationGroupKey(e,depths,counts,ids,e.seeds.data());AdaptiveData cached;
    if(!RestoreInterpolationGroup(e,key,depths.size(),ids.size(),cached)){
        struct NoNative {Engine& e;bool saved;NoNative(Engine& engine):e(engine),saved(e.forbidNativeWork){e.forbidNativeWork=true;}
            ~NoNative(){e.forbidNativeWork=saved;}} guard(e);
        cached=EvaluateBudget(e,depths,counts,ids,e.seeds.data(),false);SaveInterpolationGroup(e,key,cached);
    }
    const Samples total=Combine(stage.levels),joint=Combine(cached.levels);Budgets promoted=state.counts;size_t changed=0;
    std::ostringstream rows;rows<<header;
    for(size_t row=0;row<ids.size();++row){const size_t id=ids[row];
        const double jointError=PairedLevelPrefixChange(e,total,total,joint,row);
        std::vector<double> errors(depths.size(),0);std::vector<std::string> reasons(depths.size());bool any=false;int largest=-1;
        for(size_t level=0;level<depths.size();++level){if(counts[level][id]<=state.counts[level][id])continue;
            errors[level]=PairedLevelPrefixChange(e,total,stage.levels[level],cached.levels[level],row);
            if(adoptComplete)reasons[level]="complete_cached_output";
            else if(state.counts[level][id]<PrefixMinimumCount(e,level))reasons[level]="minimum_prefix";
            if(errors[level]>NumericalBudget(e))reasons[level]="individual_discrepancy";
            any=any||!reasons[level].empty();
            if(largest<0||errors[level]>errors[size_t(largest)])largest=int(level);
        }
        if(!any&&jointError>NumericalBudget(e)&&largest>=0)reasons[size_t(largest)]="joint_discrepancy";
        for(size_t level=0;level<depths.size();++level){if(counts[level][id]<=state.counts[level][id])continue;
            const bool grow=!reasons[level].empty();if(grow){promoted[level][id]=counts[level][id];state.coarse[id]=false;++changed;}
            rows<<Number(e.theta[id])<<','<<level<<','<<state.counts[level][id]<<','<<counts[level][id]<<','<<Number(errors[level])
                <<','<<Number(jointError)<<','<<Number(NumericalBudget(e))<<','<<grow<<','<<(grow?reasons[level]:"consistent")<<','<<e.PhiForTheta(id)<<'\n';
        }
    }
    Atomic(e.root+"/cached_prefix_consistency.csv",rows.str());Atomic(e.root+"/cached_prefix_consistency_sources.tsv",cached.sources);
    // A promotion reassembles this round before its ordinary snapshot. Keep
    // the triggering evidence even when the next consistency check is quiet.
    std::string decisionFile,sourcesFile;
    if(changed){const auto stem="cached_prefix_consistency_round"+Number(state.round)+"_"+Hash(rows.str()+cached.sources);
        decisionFile=stem+".csv";sourcesFile=stem+"_sources.tsv";
        Atomic(e.root+"/"+decisionFile,rows.str());Atomic(e.root+"/"+sourcesFile,cached.sources);}
    e.Log("cached_prefix_consistency","{\"round\":"+Number(state.round)+",\"promoted_level_rows\":"+Number(changed)
        +",\"details_file\":"+Quoted(decisionFile)+",\"sources_file\":"+Quoted(sourcesFile)
        +",\"candidate_counts_selected_before_values\":true,\"native_work_forbidden\":true,\"joint_level_sum_checked\":true,\"adopt_all_complete_prefixes\":"+(adoptComplete?"true":"false")
        +",\"growth_reasons\":\"complete cached output, individual or joint discrepancy, or minimum prefix\"}");
    if(!changed)return false;
    state.counts=promoted;
    if(adoptComplete){
        // Adopt once in this round, in place. Restarting interpolation here
        // recursively spreads the largest stencil count into its neighbours.
        // The caller recomputes all output statistics and half-prefix checks,
        // retaining the frozen common-prefix geometry check for this round.
        stage=std::move(cached);
    }else state.streak=0;
    SaveAdaptive(e,state);return true;
}
