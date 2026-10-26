// A frozen, production-trained control variate in theta. Every audit replica
// uses fresh Owen seeds: I(mean_N F(knots)) + mean_n(F(theta)-I(F(knots))).
// I is the existing linear local cubic interpolator. All covariance is kept by
// assembling eight complete replicas before Estimate/Fieller, never by adding
// separately reported confidence intervals or interpolating ratios.
struct AngularAuditPlan {
    bool enabled=false;Budgets anchors,residuals;std::string identity;
    double nativeSeconds=0,controlSeconds=0;size_t intervals=0;
};
std::vector<Engine::WorkRequest> AngularAuditRequests(const Engine& e,const std::vector<int>& depths,
        const AngularAuditPlan& plan,const std::vector<size_t>& ids,const std::vector<size_t>& knots){
    auto result=BudgetRequests(e,depths,plan.anchors,knots);
    for(const auto& group:FullAutoAuditPlan::Stencils(knots,plan.residuals,ids)){
        std::vector<size_t> rows(group.second.rows.begin(),group.second.rows.end());
        bool interior=false;for(size_t id:rows)interior=interior||!std::binary_search(knots.begin(),knots.end(),id);
        if(!interior)continue;
        Budgets ns;for(int n:group.first)ns.push_back(std::vector<int>(e.theta.size(),n));
        const auto requests=BudgetRequests(e,depths,ns,rows);result.insert(result.end(),requests.begin(),requests.end());
    }return result;
}
AngularAuditPlan PlanAngularAudit(Engine& e,const AdaptiveState& state,const std::vector<int>& depths,
        const AdaptiveData& production,const std::vector<size_t>& productionIds,const std::vector<size_t>& ids,
        double fraction,bool local){
    AngularAuditPlan plan;
    if(e.angularAuditPolicy=="off"||!e.angularCache||fraction!=1||state.knots.size()>=ids.size())return plan;
    // The residual cost/variance trainer currently assumes one common phi
    // and one N per level. Mixed rules use the complete direct holdout.
    if(!e.thetaPhi.empty()){e.Log("angular_audit_direct_fallback","{\"reason\":\"mixed theta phi rules require separate residual training\"}");return plan;}
    plan.anchors=FullAutoAuditPlan::Scale(state.counts,ids,fraction,e.initial,e.minimum,local);
    plan.residuals=plan.anchors;for(auto& level:plan.residuals)for(int& n:level)n=std::min(n,256);
    const auto productionMean=Mean(Combine(production.levels));
    const size_t width=std::min(size_t(4),state.knots.size());
    std::ostringstream observations;observations<<"interval,level,source_count,residual_count,estimated95_at_source,bias,total_interval_budget\n";
    // These common-prefix production groups have just been computed by the
    // unchanged interpolation gate. No holdout data, timings or cache coverage
    // can train this allocation. The native forecast ignores saved coverage.
    for(const auto& group:FullAutoAuditPlan::Stencils(state.knots,state.counts,productionIds)){
        const std::vector<size_t> rows(group.second.rows.begin(),group.second.rows.end());
        Budgets ns;for(int n:group.first)ns.push_back(std::vector<int>(e.theta.size(),n));
        const auto data=EvaluateBudget(e,depths,ns,rows,e.seeds.data());
        for(size_t interval:group.second.intervals){
            const size_t left=std::min(interval?interval-1:0,state.knots.size()-width),right=std::min(interval,state.knots.size()-width);
            std::vector<double> knots;std::vector<size_t> supportRows;
            for(size_t k=left;k<right+width;++k){knots.push_back(e.theta[state.knots[k]]);supportRows.push_back(std::lower_bound(rows.begin(),rows.end(),state.knots[k])-rows.begin());}
            std::vector<size_t> queryRows;std::vector<double> query;
            for(size_t r=0;r<rows.size();++r)if(rows[r]>state.knots[interval]&&rows[r]<state.knots[interval+1]){queryRows.push_back(r);query.push_back(e.theta[rows[r]]);}
            std::vector<double> ci(depths.size(),0);double bias=0;
            std::vector<Samples> deltas(depths.size(),Samples(8,Matrix(query.size())));
            for(size_t l=0;l<depths.size();++l)for(int s=0;s<8;++s){Matrix support;for(size_t r:supportRows)support.push_back(data.levels[l][s][r]);
                const auto predicted=LocalCubicInterpolate(knots,support,query);
                for(size_t r=0;r<query.size();++r)for(int j=0;j<16;++j)deltas[l][s][r][j]=data.levels[l][s][queryRows[r]][j]-predicted[r][j];}
            for(size_t r=0;r<query.size();++r){const size_t source=std::lower_bound(productionIds.begin(),productionIds.end(),rows[queryRows[r]])-productionIds.begin();
                const auto& reference=productionMean[source];Require(reference[0]>0,"nonpositive production angular-control scale");
                for(int j=0;j<16;++j){double totalBias=0;const double ratio=j?reference[j]/reference[0]:0;
                    for(size_t l=0;l<depths.size();++l){double mean=0,var=0,ratioMean=0,ratioVar=0;
                        for(int s=0;s<8;++s){const auto& z=deltas[l][s][r];mean+=z[j]/8;ratioMean+=(z[j]-ratio*z[0])/8;}
                        for(int s=0;s<8;++s){const auto& z=deltas[l][s][r];var+=std::pow(z[j]-mean,2)/7;ratioVar+=std::pow(z[j]-ratio*z[0]-ratioMean,2)/7;}
                        ci[l]=std::max(ci[l],t95*std::sqrt(var/8)/reference[0]);
                        if(e.ratios&&j)ci[l]=std::max(ci[l],t95*std::sqrt(ratioVar/8)/reference[0]);
                        totalBias+=mean;
                    }bias=std::max(bias,std::fabs(totalBias)/reference[0]);
                }
            }
            // A missing interior production node is not evidence of zero
            // variance: retain the complete allocation for that interval.
            const double allowance=std::min(e.statistical*.25,std::max(e.interpolation*.05,(e.interpolation-bias)*.5))/depths.size();
            for(size_t l=0;l<depths.size();++l){const int original=group.first[l];
                int n=query.empty()?original:Power2(original*std::pow(ci[l]/std::max(1e-30,allowance),2));
                n=std::min(original,std::max(std::min(original,256),n));
                for(size_t id:ids)if(id>state.knots[interval]&&id<state.knots[interval+1])plan.residuals[l][id]=std::max(plan.residuals[l][id],n);
                observations<<interval<<','<<l<<','<<original<<','<<n<<','<<Number(ci[l])<<','<<Number(bias)<<','<<Number(allowance)<<'\n';
            }++plan.intervals;
        }
    }
    std::ostringstream frozen;frozen<<"angular-residual-v1\n"<<production.sources<<Number(fraction)<<'\n';
    for(size_t id:ids){frozen<<id;for(size_t l=0;l<depths.size();++l)frozen<<' '<<plan.anchors[l][id]<<' '<<plan.residuals[l][id];frozen<<'\n';}
    plan.identity=Hash(frozen.str());
    AdditionalWork work(e,e.seeds.data(),true,true);
    const auto point=work.Forecast(AuditRequests(e,depths,plan.anchors,ids,state.knots));
    const auto controlled=work.Forecast(AngularAuditRequests(e,depths,plan,ids,state.knots));
    plan.nativeSeconds=point.seconds;plan.controlSeconds=controlled.seconds;
    plan.enabled=e.angularAuditPolicy=="on"||(point.profiled&&controlled.profiled&&controlled.seconds<point.seconds*.90);
    // A resumed holdout keeps its scheduling decision as well as its counts.
    // New cached aliases/timing observations must not change the selected
    // estimator after any independent replicas have been observed.
    if(Exists(e.root+"/angular_audit_plan.json")){
        const auto saved=CompactJson(Read(e.root+"/angular_audit_plan.json"));
        if(saved.find("\"allocation_fnv1a64\":"+Quoted(plan.identity))!=std::string::npos)
            plan.enabled=saved.find("\"enabled\":true")!=std::string::npos;
    }
    e.Record(e.root+"/angular_audit_training.csv",observations.str());
    std::ostringstream budget;budget<<"theta_deg,is_support";for(size_t l=0;l<depths.size();++l)budget<<",anchor_N"<<l<<",residual_N"<<l;budget<<'\n';
    for(size_t id:ids){budget<<Number(e.theta[id])<<','<<std::binary_search(state.knots.begin(),state.knots.end(),id);for(size_t l=0;l<depths.size();++l)budget<<','<<plan.anchors[l][id]<<','<<plan.residuals[l][id];budget<<'\n';}
    e.Record(e.root+"/angular_audit_budget.csv",budget.str());
    e.Record(e.root+"/angular_audit_plan.json","{\"schema\":1,\"strategy\":\"independent anchors plus common-orientation angular residual\",\"allocation_fnv1a64\":"+Quoted(plan.identity)
        +",\"enabled\":"+(plan.enabled?"true":"false")+",\"support_theta\":"+Number(state.knots.size())+",\"checked_theta\":"+Number(ids.size())
        +",\"predicted_point_worker_seconds\":"+Number(point.seconds)+",\"predicted_control_worker_seconds\":"+Number(controlled.seconds)
        +",\"profiled\":"+(point.profiled&&controlled.profiled?"true":"false")+",\"training\":\"production replicas only; frozen before holdout\",\"acceptance_gates_unchanged\":true}");
    return plan;
}
AdaptiveData EvaluateAngularAudit(Engine& e,const std::vector<int>& depths,const AngularAuditPlan& plan,
        const std::vector<size_t>& ids,const std::vector<size_t>& knots,const unsigned* seeds){
    Require(plan.enabled,"disabled angular audit evaluated");
    e.Record(e.root+"/independent_audit_progress.json","{\"phase\":\"native prefetch for frozen angular control\",\"allocation_fnv1a64\":"+Quoted(plan.identity)+"}");
    const auto residualGroups=FullAutoAuditPlan::Stencils(knots,plan.residuals,ids);
    std::ostringstream metadata;metadata<<"{\"schema\":1,\"groups\":[";size_t metadataGroup=0;
    for(const auto& group:residualGroups){metadata<<(metadataGroup?",":"")<<"{\"id\":"<<metadataGroup<<",\"counts\":[";
        ++metadataGroup;
        for(size_t l=0;l<group.first.size();++l)metadata<<(l?",":"")<<group.first[l];
        metadata<<"],\"theta_indices\":[";
        bool first=true;
        for(size_t id:group.second.rows){metadata<<(first?"":",")<<id;first=false;}metadata<<"],\"intervals\":[";
        for(size_t i=0;i<group.second.intervals.size();++i)metadata<<(i?",":"")<<group.second.intervals[i];
        metadata<<"]}";
        }
    metadata<<"]}";e.Record(e.root+"/angular_audit_groups.json",metadata.str());WriteGrid(e.root+"/angular_audit_support.csv",Angles(e.theta,knots));
    PrefetchBudget(e,AngularAuditRequests(e,depths,plan,ids,knots),seeds);
    e.Record(e.root+"/independent_audit_progress.json","{\"phase\":\"independent anchor means\",\"support_theta\":"+Number(knots.size())+"}");
    const auto anchors=EvaluateBudget(e,depths,plan.anchors,knots,seeds);
    AdaptiveData result;result.levels.resize(depths.size());result.costs.assign(depths.size(),std::vector<double>(ids.size()));result.analytic.assign(ids.size(),0);
    const auto x=Angles(e.theta,ids),sx=Angles(e.theta,knots);std::ostringstream sources,roles;sources<<anchors.sources;
    {std::istringstream records(anchors.sources);std::string line;while(std::getline(records,line))roles<<"anchor -1 "<<line<<'\n';}
    for(size_t l=0;l<depths.size();++l){result.levels[l].resize(8);for(int s=0;s<8;++s)result.levels[l][s]=LocalCubicInterpolate(sx,anchors.levels[l][s],x);}
    const size_t width=std::min(size_t(4),knots.size());std::vector<bool> assigned(ids.size(),false);
    for(size_t r=0;r<ids.size();++r)assigned[r]=std::binary_search(knots.begin(),knots.end(),ids[r]);
    size_t completedGroups=0;
    for(const auto& group:FullAutoAuditPlan::Stencils(knots,plan.residuals,ids)){
        const size_t groupId=completedGroups++;e.Record(e.root+"/independent_audit_progress.json","{\"phase\":\"common-orientation angular residuals\",\"group\":"+Number(groupId)+"}");
        std::vector<size_t> rows(group.second.rows.begin(),group.second.rows.end());bool interior=false;
        for(size_t id:rows)interior=interior||!std::binary_search(knots.begin(),knots.end(),id);
        if(!interior)continue;
        Budgets ns;for(int n:group.first)ns.push_back(std::vector<int>(e.theta.size(),n));
        const auto data=EvaluateBudget(e,depths,ns,rows,seeds);sources<<data.sources;
        {std::istringstream records(data.sources);std::string line;while(std::getline(records,line))roles<<"residual "<<groupId<<' '<<line<<'\n';}
        for(size_t interval:group.second.intervals){const size_t left=std::min(interval?interval-1:0,knots.size()-width),right=std::min(interval,knots.size()-width);
            std::vector<double> supportX;std::vector<size_t> supports;
            for(size_t k=left;k<right+width;++k){supportX.push_back(e.theta[knots[k]]);supports.push_back(std::lower_bound(rows.begin(),rows.end(),knots[k])-rows.begin());}
            std::vector<double> query;std::vector<size_t> positions,rawRows;
            for(size_t r=0;r<ids.size();++r)if(ids[r]>knots[interval]&&ids[r]<knots[interval+1]){
                Require(!assigned[r],"angular residual added twice");assigned[r]=true;positions.push_back(r);rawRows.push_back(std::lower_bound(rows.begin(),rows.end(),ids[r])-rows.begin());query.push_back(e.theta[ids[r]]);}
            for(size_t l=0;l<depths.size();++l)for(int s=0;s<8;++s){Matrix supportsM;for(size_t r:supports)supportsM.push_back(data.levels[l][s][r]);
                const auto predicted=LocalCubicInterpolate(supportX,supportsM,query);
                for(size_t r=0;r<query.size();++r)for(int j=0;j<16;++j)result.levels[l][s][positions[r]][j]+=data.levels[l][s][rawRows[r]][j]-predicted[r][j];}
        }
    }
    Require(std::find(assigned.begin(),assigned.end(),false)==assigned.end(),"angular audit skipped a checked node");
    e.Record(e.root+"/angular_audit_sources.tsv",roles.str());result.sources=sources.str();return result;
}
