// Adjacent support estimates must resolve their angular DIFFERENCE as well as
// their individual values. This checks sampling noise, not physical slope or
// curvature: deterministic sharp structure is never flattened or penalized.
// No matrices are modified. Failed pairs/triples become bounded cohorts.
const char* AngularSamplingPolicy(){return "incoherent-component-shape-noise-v3";}
struct AngularSamplingDecision {
    bool enabled=false,passed=true;double maximum=0,tolerance=0;
    double maximumRatio=0;Mueller envelope{},allowance{};
    size_t pairs=0,failed=0;Budgets proposed;
    size_t shapes=0,failedShapes=0;double maximumShapeRatio=0;
    std::vector<RoundOutputPriority> priorities;
    std::string table,shapeTable;
};
double AngularNoise(const std::array<double,8>& values){
    double mean=0,variance=0;for(double v:values){if(!std::isfinite(v))return 1e100;mean+=v/8;}
    for(double v:values)variance+=(v-mean)*(v-mean)/7;
    const double noise=t95*std::sqrt(variance/8);
    return std::isfinite(noise)?noise:1e100;
}
AngularSamplingDecision CheckAngularSampling(const Engine& e,const AdaptiveState& state,
        const AdaptiveData& stage,const std::vector<size_t>& ids,const Samples& samples,const Statistics& stats){
    AngularSamplingDecision result;result.proposed=state.counts;
    result.enabled=!e.go&&!e.coherent&&e.local;result.tolerance=e.interpolation;
    std::ostringstream table;table<<"theta_left_deg,theta_right_deg,noise95,allowance,passed,element,normalized_ratio";
    for(size_t l=0;l<state.counts.size();++l)table<<",N_left_level"<<l<<",N_right_level"<<l<<",next_left_level"<<l<<",next_right_level"<<l;
    table<<'\n';
    std::ostringstream shapeTable;shapeTable<<"theta_left_deg,theta_center_deg,theta_right_deg,span_support_steps,noise95,allowance,passed,element,normalized_ratio,mean_contrast";
    for(size_t l=0;l<state.counts.size();++l)shapeTable<<",N_left_level"<<l<<",N_center_level"<<l<<",N_right_level"<<l
        <<",next_left_level"<<l<<",next_center_level"<<l<<",next_right_level"<<l;
    shapeTable<<'\n';
    if(!result.enabled){result.table=table.str();result.shapeTable=shapeTable.str();return result;}
    Require(result.tolerance>0&&samples.size()==8&&stats.mean.size()==ids.size(),"invalid angular sampling input");
    Require(!stage.levels.empty()&&stage.levels.size()==state.counts.size(),"angular sampling level mismatch");
    for(const auto& sample:samples)Require(sample.size()==ids.size(),"angular sampling row mismatch");
    if(e.ratios)Require(stats.ratioInterval.size()==ids.size(),"missing angular ratio intervals");
    // The M11-scaled contract alone can accept visibly noisy weak elements.
    // Reserve at most epsilon/4 for angular noise relative to each component's
    // amplitude. v2 used the whole epsilon and still accepted visible noise.
    // A stricter explicit interpolation tolerance is respected as well.
    // epsilon is the amplitude floor; no division by a local zero is needed.
    // These are stopping diagnostics, not weights applied to the estimator.
    Require(e.eps>0&&e.eps<1,"invalid angular component tolerance");
    result.envelope.fill(e.eps);
    for(const auto& mean:stats.mean)if(std::isfinite(mean[0])&&mean[0]>0)
        for(int j=0;j<16;++j)if(std::isfinite(mean[j]))
            result.envelope[j]=std::max(result.envelope[j],std::fabs(mean[j]/mean[0]));
    const double angularFraction=std::min(result.tolerance,e.eps/4);
    for(int j=0;j<16;++j)result.allowance[j]=std::min(result.tolerance,angularFraction*result.envelope[j]);
    for(size_t l=0;l<stage.levels.size();++l){
        Require(state.counts[l].size()==e.theta.size()&&stage.levels[l].size()==8,"angular sampling budget/replica mismatch");
        for(const auto& sample:stage.levels[l])Require(sample.size()==ids.size(),"angular level row mismatch");
    }
    struct Stencil {std::vector<size_t> nodes;int span;double noise,allowance,mean;bool passed;int element;bool ratio;};
    std::vector<Stencil> stencils;
    auto check=[&](const std::vector<size_t>& nodes,const std::vector<double>& weights,int span){
        std::vector<size_t> rows;bool positive=true;double scale=1e300;
        for(size_t id:nodes){
            const auto at=std::lower_bound(ids.begin(),ids.end(),id);
            Require(at!=ids.end()&&*at==id,"angular support absent from production");
            const size_t row=at-ids.begin();rows.push_back(row);
            const double m11=stats.mean[row][0];positive=positive&&std::isfinite(m11)&&m11>0;
            scale=std::min(scale,m11);
        }
        // Nonpositive intermediate estimates already fail the ordinary gate.
        // Keep them refinable instead of terminating the controller here.
        if(!positive)scale=1.;
        double worst=-1,worstScore=-1;int element=0;bool ratio=false;
        std::array<double,8> contrast{};
        auto value=[&](const Matrix& sample,int j,bool normalized){
            double total=0;
            for(size_t k=0;k<rows.size();++k){
                const auto& v=sample[rows[k]];const auto& m=stats.mean[rows[k]];
                total+=weights[k]*(normalized?(v[j]-(m[j]/m[0])*v[0])/m[0]:v[j]/scale);
            }
            return total;
        };
        for(int j=0;j<16;++j)for(int normalized=0;normalized<=(positive&&e.ratios&&j?1:0);++normalized){
            std::array<double,8> v{};for(int s=0;s<8;++s)v[s]=value(samples[s],j,normalized);
            double noise=AngularNoise(v);
            // A linearized ratio contrast cannot rescue an unresolved
            // denominator. The ordinary Fieller gate remains authoritative.
            if(normalized)for(size_t row:rows)
                if(!std::isfinite(stats.ratioInterval[row][j])||stats.ratioInterval[row][j]>=1e50)noise=1e100;
            result.maximum=std::max(result.maximum,noise);
            const double score=noise/result.allowance[j];
            if(score>worstScore){worst=noise;worstScore=score;element=j;ratio=normalized;contrast=v;}
        }
        if(!positive){worst=1e100;worstScore=worst/result.allowance[element];}
        if(span)++result.shapes;else ++result.pairs;
        result.maximum=std::max(result.maximum,worst);
        result.maximumRatio=std::max(result.maximumRatio,worstScore);
        if(span)result.maximumShapeRatio=std::max(result.maximumShapeRatio,worstScore);
        const bool passed=worstScore<=1;
        if(!passed){
            result.passed=false;if(span)++result.failedShapes;else ++result.failed;
            // Positive covariance identifies levels that contribute noise to
            // the actual combined contrast, including cutoff corrections.
            // Freeze all proposals from the current counts; no recursive
            // propagation of a remote maximum along the support grid.
            double mean=0;for(double v:contrast)mean+=v/8;
            std::vector<double> covariance(stage.levels.size());
            for(size_t l=0;l<stage.levels.size();++l){
                Require(stage.levels[l].size()==8,"angular sampling requires eight replicas");
                std::array<double,8> v{};double average=0;
                for(int s=0;s<8;++s){v[s]=value(stage.levels[l][s],element,ratio);average+=v[s]/8;}
                for(int s=0;s<8;++s)covariance[l]+=(v[s]-average)*(contrast[s]-mean)/7;
            }
            size_t dominant=0;for(size_t l=1;l<covariance.size();++l)if(covariance[l]>covariance[dominant])dominant=l;
            bool synchronized=false;
            for(size_t l=0;l<covariance.size();++l){
                if(covariance[dominant]>0&&covariance[l]<covariance[dominant]/8)continue;
                int largest=0;for(size_t id:nodes)largest=std::max(largest,state.counts[l][id]);
                for(size_t id:nodes){
                    const int n=state.counts[l][id],cap=n>e.maxCount/4?e.maxCount:4*n;
                    const int next=std::min(largest,cap);
                    result.proposed[l][id]=std::max(result.proposed[l][id],next);synchronized=synchronized||next>n;
                }
            }
            if(!synchronized){
                // Matching prefixes already cancel most common fluctuations.
                // If noise remains, grow the complete stencil together.
                for(size_t id:nodes)result.proposed[dominant][id]=
                    std::max(result.proposed[dominant][id],Doubled(e,state.counts[dominant][id]));
            }
            result.priorities.push_back({worstScore*result.tolerance,nodes});
        }
        double mean=0;for(size_t k=0;k<rows.size();++k){const auto& m=stats.mean[rows[k]];
            mean+=weights[k]*(ratio?m[element]/m[0]:m[element]/scale);}
        stencils.push_back({nodes,span,worst,result.allowance[element],mean,passed,element,ratio});
    };
    Require(std::is_sorted(state.knots.begin(),state.knots.end()),"unordered angular supports");
    for(size_t k=1;k<state.knots.size();++k)
        check({state.knots[k-1],state.knots[k]},{-1.,1.},0);
    // Test the noise of a local bend, not its magnitude. Span 2 also resolves
    // a correlated bump spread over several individually acceptable pairs.
    // Two fixed spans keep CPU work and checkpoint diagnostics linear in knots.
    for(size_t span:{size_t(1),size_t(2)})for(size_t k=span;k+span<state.knots.size();++k){
        const size_t a=state.knots[k-span],b=state.knots[k],c=state.knots[k+span];
        Require(e.theta[a]<e.theta[b]&&e.theta[b]<e.theta[c],"degenerate angular shape stencil");
        const double w=(e.theta[b]-e.theta[a])/(e.theta[c]-e.theta[a]);
        check({a,b,c},{-(1-w),1.,-w},int(span));
    }
    // Report the combined proposal after all overlapping stencils have voted.
    for(const auto& s:stencils){
        auto& out=s.span?shapeTable:table;
        for(size_t id:s.nodes)out<<Number(e.theta[id])<<',';
        if(s.span)out<<s.span<<',';
        out<<Number(s.noise)<<','<<Number(s.allowance)<<','<<s.passed<<",M"<<s.element/4+1<<s.element%4+1<<','<<s.ratio;
        if(s.span)out<<','<<Number(s.mean);
        for(size_t l=0;l<state.counts.size();++l){
            for(size_t id:s.nodes)out<<','<<state.counts[l][id];
            for(size_t id:s.nodes)out<<','<<result.proposed[l][id];
        }
        out<<'\n';
    }
    result.table=table.str();result.shapeTable=shapeTable.str();return result;
}
std::string AngularSamplingFields(const AngularSamplingDecision& result){
    return ",\"angular_sampling_policy\":"+Quoted(AngularSamplingPolicy())
        +",\"angular_sampling_enabled\":"+(result.enabled?"true":"false")
        +",\"angular_sampling_pass\":"+(result.passed?"true":"false")
        +",\"angular_sampling_noise95\":"+Number(result.maximum)
        +",\"angular_sampling_tolerance\":"+Number(result.tolerance)
        +",\"angular_sampling_max_allowance_ratio\":"+Number(result.maximumRatio)
        +",\"angular_sampling_failed_pairs\":"+Number(result.failed)
        +",\"angular_shape_checked_stencils\":"+Number(result.shapes)
        +",\"angular_shape_failed_stencils\":"+Number(result.failedShapes)
        +",\"angular_shape_max_allowance_ratio\":"+Number(result.maximumShapeRatio);
}
void SaveAngularSampling(Engine& e,const AngularSamplingDecision& result,const std::string& folder){
    Atomic(folder+"/angular_sampling.csv",result.table);
    Atomic(folder+"/angular_shape.csv",result.shapeTable);
    std::string scales="[",limits="[";
    for(int j=0;j<16;++j){if(j){scales+=",";limits+=",";}scales+=Number(result.envelope[j]);limits+=Number(result.allowance[j]);}
    scales+="]";limits+="]";
    const auto report="{\"schema\":3"+AngularSamplingFields(result)
        +",\"component_envelopes\":"+scales+",\"component_allowances\":"+limits
        +",\"component_amplitude_floor\":"+Number(e.eps)
        +",\"component_angular_fraction\":"+Number(std::min(result.tolerance,e.eps/4))
        +",\"shape_support_spans\":[1,2],\"physical_curvature_constrained\":false"
        +",\"checked_pairs\":"+Number(result.pairs)
        +",\"replicas\":8,\"physical_slope_constrained\":false,\"matrices_modified\":false"
         ",\"scope\":\"noise of support increments and centered three-point contrasts at two spans; ratios use linearized influence with an unresolved-Fieller guard; not a simultaneous confidence bound\"}";
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
        +",\"component_rule\":\"min(interpolation_tolerance, min(interpolation_tolerance,epsilon/4)*max(epsilon,peak_abs_Mij_over_M11))\""
        +",\"shape_support_spans\":[1,2],\"physical_curvature_constrained\":false"
        +",\"native_blocks_deleted\":0,\"sampling_noise_only\":true,\"all16\":true,\"ratios_follow_fullauto_ratios\":true}");
}
// Saved diagnostics are scheduling hints only. Never restore acceptance from
// them. The current samples, means and accuracy gates are recomputed on resume.
std::vector<RoundOutputPriority> SavedAngularPriorities(const Engine& e,const AdaptiveState& state,
        const std::string& folder){
    if(e.go||e.coherent||!e.local||!Exists(folder+"/angular_sampling.json")
            ||Read(folder+"/angular_sampling.json").find(AngularSamplingPolicy())==std::string::npos
            ||!Exists(folder+"/angular_sampling.csv")||!Exists(folder+"/angular_shape.csv"))return {};
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
        priorities.push_back({noise/tolerance*e.interpolation,{i,j}});
    }
    std::istringstream shapes(Read(folder+"/angular_shape.csv"));
    if(!std::getline(shapes,line)||line.find("theta_left_deg,theta_center_deg,theta_right_deg,span_support_steps,noise95,allowance,passed,element,")!=0)return {};
    while(std::getline(shapes,line)){
        std::replace(line.begin(),line.end(),',',' ');std::istringstream values(line);
        double a,b,c,noise,tolerance;int span,passed;
        if(!(values>>a>>b>>c>>span>>noise>>tolerance>>passed)||!std::isfinite(a)||!std::isfinite(b)||!std::isfinite(c)
                ||!(a<b&&b<c)||(span!=1&&span!=2)||!std::isfinite(noise)||noise<0
                ||!std::isfinite(tolerance)||!(tolerance>0)||(passed!=0&&passed!=1))return {};
        if(passed)continue;
        std::vector<size_t> nodes,positions;bool valid=true;
        for(double angle:{a,b,c}){
            const auto at=std::lower_bound(e.theta.begin(),e.theta.end(),angle);
            if(at==e.theta.end()||*at!=angle){valid=false;break;}
            const size_t id=at-e.theta.begin();const auto knot=std::lower_bound(state.knots.begin(),state.knots.end(),id);
            if(knot==state.knots.end()||*knot!=id){valid=false;break;}
            nodes.push_back(id);positions.push_back(knot-state.knots.begin());
        }
        if(!valid||positions[1]-positions[0]!=size_t(span)||positions[2]-positions[1]!=size_t(span))continue;
        priorities.push_back({noise/tolerance*e.interpolation,nodes});
    }
    return priorities;
}
