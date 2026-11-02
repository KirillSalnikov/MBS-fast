// Adjacent support estimates must resolve their angular DIFFERENCE as well as
// their individual values. This checks sampling noise, not physical slope or
// curvature: deterministic sharp structure is never flattened or penalized.
// No matrices are modified. Failed pairs become bounded scheduling cohorts.
const char* AngularSamplingPolicy(){return "incoherent-paired-support-noise-v1";}
struct AngularSamplingDecision {
    bool enabled=false,passed=true;double maximum=0,tolerance=0;
    size_t pairs=0,failed=0;Budgets proposed;
    std::vector<RoundOutputPriority> priorities;
    std::string table;
};
double AngularNoise(const std::array<double,8>& values){
    double mean=0,variance=0;for(double v:values){if(!std::isfinite(v))return 1e100;mean+=v/8;}
    for(double v:values)variance+=(v-mean)*(v-mean)/7;
    return t95*std::sqrt(variance/8);
}
AngularSamplingDecision CheckAngularSampling(const Engine& e,const AdaptiveState& state,
        const AdaptiveData& stage,const std::vector<size_t>& ids,const Samples& samples,const Statistics& stats){
    AngularSamplingDecision result;result.proposed=state.counts;
    result.enabled=!e.go&&!e.coherent&&e.local;result.tolerance=e.interpolation;
    std::ostringstream table;table<<"theta_left_deg,theta_right_deg,noise95,allowance,passed,element,normalized_ratio";
    for(size_t l=0;l<state.counts.size();++l)table<<",N_left_level"<<l<<",N_right_level"<<l<<",next_left_level"<<l<<",next_right_level"<<l;
    table<<'\n';
    if(!result.enabled){result.table=table.str();return result;}
    Require(result.tolerance>0&&samples.size()==8&&stats.mean.size()==ids.size(),"invalid angular sampling input");
    Require(!stage.levels.empty()&&stage.levels.size()==state.counts.size(),"angular sampling level mismatch");
    for(const auto& sample:samples)Require(sample.size()==ids.size(),"angular sampling row mismatch");
    if(e.ratios)Require(stats.ratioInterval.size()==ids.size(),"missing angular ratio intervals");
    for(size_t l=0;l<stage.levels.size();++l){
        Require(state.counts[l].size()==e.theta.size()&&stage.levels[l].size()==8,"angular sampling budget/replica mismatch");
        for(const auto& sample:stage.levels[l])Require(sample.size()==ids.size(),"angular level row mismatch");
    }
    struct Pair {size_t left,right;double noise;bool passed;int element;bool ratio;};
    std::vector<Pair> pairs;
    for(size_t k=1;k<state.knots.size();++k){
        const size_t left=state.knots[k-1],right=state.knots[k];
        const auto la=std::lower_bound(ids.begin(),ids.end(),left),ra=std::lower_bound(ids.begin(),ids.end(),right);
        Require(la!=ids.end()&&*la==left&&ra!=ids.end()&&*ra==right,"angular support absent from production");
        const size_t a=la-ids.begin(),b=ra-ids.begin();
        const auto& ma=stats.mean[a];const auto& mb=stats.mean[b];
        const bool positive=std::isfinite(ma[0])&&std::isfinite(mb[0])&&ma[0]>0&&mb[0]>0;
        // Nonpositive intermediate estimates already fail the ordinary gate.
        // Keep them refinable instead of terminating the controller here.
        const double scale=positive?std::min(ma[0],mb[0]):1.;
        double worst=-1;int element=0;bool ratio=false;
        std::array<double,8> contrast{};
        auto value=[&](const Mueller& va,const Mueller& vb,int j,bool normalized){
            if(!normalized)return (vb[j]-va[j])/scale;
            return (vb[j]-(mb[j]/mb[0])*vb[0])/mb[0]
                -(va[j]-(ma[j]/ma[0])*va[0])/ma[0];
        };
        for(int j=0;j<16;++j)for(int normalized=0;normalized<=(positive&&e.ratios&&j?1:0);++normalized){
            std::array<double,8> v{};for(int s=0;s<8;++s)v[s]=value(samples[s][a],samples[s][b],j,normalized);
            double noise=AngularNoise(v);
            // A linearized ratio contrast cannot rescue an unresolved
            // denominator. The ordinary Fieller gate remains authoritative.
            if(normalized&&(!std::isfinite(stats.ratioInterval[a][j])||!std::isfinite(stats.ratioInterval[b][j])
                    ||stats.ratioInterval[a][j]>=1e50||stats.ratioInterval[b][j]>=1e50))noise=1e100;
            if(noise>worst){worst=noise;element=j;ratio=normalized;contrast=v;}
        }
        if(!positive)worst=1e100;
        ++result.pairs;result.maximum=std::max(result.maximum,worst);
        const bool passed=worst<=result.tolerance;
        if(!passed){
            result.passed=false;++result.failed;
            // Positive covariance identifies levels that contribute noise to
            // the actual combined contrast, including cutoff corrections.
            // Freeze all proposals from the current counts; no recursive
            // propagation of a remote maximum along the support grid.
            double mean=0;for(double v:contrast)mean+=v/8;
            std::vector<double> covariance(stage.levels.size());
            for(size_t l=0;l<stage.levels.size();++l){
                Require(stage.levels[l].size()==8,"angular sampling requires eight replicas");
                std::array<double,8> v{};double average=0;
                for(int s=0;s<8;++s){v[s]=value(stage.levels[l][s][a],stage.levels[l][s][b],element,ratio);average+=v[s]/8;}
                for(int s=0;s<8;++s)covariance[l]+=(v[s]-average)*(contrast[s]-mean)/7;
            }
            size_t dominant=0;for(size_t l=1;l<covariance.size();++l)if(covariance[l]>covariance[dominant])dominant=l;
            bool synchronized=false;
            for(size_t l=0;l<covariance.size();++l){
                if(covariance[dominant]>0&&covariance[l]<covariance[dominant]/8)continue;
                const int na=state.counts[l][left],nb=state.counts[l][right];if(na==nb)continue;
                const size_t lower=na<nb?left:right;const int n=std::min(na,nb);
                const int cap=n>e.maxCount/4?e.maxCount:4*n;
                const int next=std::min(std::max(na,nb),cap);
                result.proposed[l][lower]=std::max(result.proposed[l][lower],next);synchronized=synchronized||next>n;
            }
            if(!synchronized){
                // Matching prefixes already cancel most common fluctuations.
                // If noise remains, grow both ends of the dominant level.
                for(size_t id:{left,right})result.proposed[dominant][id]=
                    std::max(result.proposed[dominant][id],Doubled(e,state.counts[dominant][id]));
            }
            result.priorities.push_back({worst,{left,right}});
        }
        pairs.push_back({left,right,worst,passed,element,ratio});
    }
    // Report the final combined proposal after overlapping pairs have voted.
    for(const auto& pair:pairs){
        const size_t left=pair.left,right=pair.right;const double worst=pair.noise;
        const bool passed=pair.passed,ratio=pair.ratio;const int element=pair.element;
        table<<Number(e.theta[left])<<','<<Number(e.theta[right])<<','<<Number(worst)<<','<<Number(result.tolerance)
            <<','<<passed<<",M"<<element/4+1<<element%4+1<<','<<ratio;
        for(size_t l=0;l<state.counts.size();++l)table<<','<<state.counts[l][left]<<','<<state.counts[l][right]
            <<','<<result.proposed[l][left]<<','<<result.proposed[l][right];
        table<<'\n';
    }
    result.table=table.str();return result;
}
std::string AngularSamplingFields(const AngularSamplingDecision& result){
    return ",\"angular_sampling_policy\":"+Quoted(AngularSamplingPolicy())
        +",\"angular_sampling_enabled\":"+(result.enabled?"true":"false")
        +",\"angular_sampling_pass\":"+(result.passed?"true":"false")
        +",\"angular_sampling_noise95\":"+Number(result.maximum)
        +",\"angular_sampling_tolerance\":"+Number(result.tolerance)
        +",\"angular_sampling_failed_pairs\":"+Number(result.failed);
}
void SaveAngularSampling(Engine& e,const AngularSamplingDecision& result,const std::string& folder){
    Atomic(folder+"/angular_sampling.csv",result.table);
    const auto report="{\"schema\":1"+AngularSamplingFields(result)
        +",\"checked_pairs\":"+Number(result.pairs)
        +",\"replicas\":8,\"physical_slope_constrained\":false,\"matrices_modified\":false"
         ",\"scope\":\"paired support-increment sampling noise; ratios use linearized influence with an unresolved-Fieller guard; not a simultaneous confidence bound\"}";
    Atomic(folder+"/angular_sampling.json",report);e.Log("angular_sampling",report);
}

void PrepareAngularSamplingPolicy(Engine& e,AdaptiveState& state,bool resumed){
    const auto path=e.root+"/angular_sampling_policy.json";
    const bool enabled=!e.go&&!e.coherent&&e.local;
    const bool changed=!Exists(path)||Read(path).find(AngularSamplingPolicy())==std::string::npos;
    if(enabled&&changed&&resumed){
        state.streak=0;e.independentPassed=false;SaveAdaptive(e,state);
        e.Log("angular_sampling_policy_migration","{\"production_blocks_retained\":true,\"old_pass_streak_reset\":true}");
    }
    e.Record(path,"{\"policy\":"+Quoted(AngularSamplingPolicy())
        +",\"scope\":\"incoherent PO with local orientation budgets\",\"tolerance\":"+Number(e.interpolation)
        +",\"native_blocks_deleted\":0,\"sampling_noise_only\":true,\"all16\":true,\"ratios_follow_fullauto_ratios\":true}");
}
// Saved diagnostics are scheduling hints only. Never restore acceptance from
// them. The current samples, means and accuracy gates are recomputed on resume.
std::vector<RoundOutputPriority> SavedAngularPriorities(const Engine& e,const AdaptiveState& state,
        const std::string& folder){
    if(e.go||e.coherent||!e.local||!Exists(folder+"/angular_sampling.json")
            ||Read(folder+"/angular_sampling.json").find(AngularSamplingPolicy())==std::string::npos
            ||!Exists(folder+"/angular_sampling.csv"))return {};
    std::vector<RoundOutputPriority> priorities;
    std::istringstream file(Read(folder+"/angular_sampling.csv"));std::string line;
    if(!std::getline(file,line)||line.find("theta_left_deg,theta_right_deg,noise95,allowance,passed,element,")!=0)return {};
    while(std::getline(file,line)){
        std::replace(line.begin(),line.end(),',',' ');std::istringstream values(line);
        double a,b,noise,tolerance;int passed;
        if(!(values>>a>>b>>noise>>tolerance>>passed)||!std::isfinite(a)||!std::isfinite(b)
                ||!std::isfinite(noise)||noise<0||!(tolerance>0)||!(a<b)||(passed!=0&&passed!=1))return {};
        if(passed)continue;
        const auto left=std::lower_bound(e.theta.begin(),e.theta.end(),a),right=std::lower_bound(e.theta.begin(),e.theta.end(),b);
        if(left==e.theta.end()||*left!=a||right==e.theta.end()||*right!=b)continue;
        const size_t i=left-e.theta.begin(),j=right-e.theta.begin();
        const auto knot=std::lower_bound(state.knots.begin(),state.knots.end(),i);
        if(knot==state.knots.end()||*knot!=i||knot+1==state.knots.end()||*(knot+1)!=j)continue;
        priorities.push_back({noise,{i,j}});
    }
    return priorities;
}
