// Pointwise numerical error budget. Physical/execution reserves remain fixed.
// Unchecked output angles retain the full interpolation allowance: no assumed
// envelope between guards, and no reuse of pilot depth errors as final bounds.
const char* BudgetPolicy(){return "pointwise-stat-plus-interpolation-v1";}
double NumericalBudget(const Engine& e){return e.statistical+e.interpolation;}
double FixedBudgetReserve(const Engine& e){return e.eps-NumericalBudget(e);}
struct BudgetResidual {Mueller matrix{},ratio{};bool measured=false;};
BudgetResidual ResidualAt(const Engine& e,const InterpolationResult& interp,const std::vector<size_t>& ids,size_t id){
    BudgetResidual out;out.matrix.fill(e.interpolation);out.ratio.fill(e.interpolation);out.ratio[0]=0;
    const auto at=std::lower_bound(ids.begin(),ids.end(),id);
    if(at==ids.end()||*at!=id)return out;
    const size_t r=at-ids.begin();
    Require(interp.matrixErrors.size()==ids.size()&&interp.ratioErrors.size()==ids.size(),"missing component interpolation estimates");
    out.matrix=interp.matrixErrors[r];out.ratio=interp.ratioErrors[r];out.measured=true;return out;
}
double RemainingStatistics(const Engine& e,const BudgetResidual& residual){
    double used=0;for(int j=0;j<16;++j){used=std::max(used,residual.matrix[j]);if(e.ratios&&j)used=std::max(used,residual.ratio[j]);}
    // Only a scheduling target. Interpolation failures first distinguish
    // noisy residuals from angular structure. Acceptance tests each component below.
    return std::max(e.statistical,NumericalBudget(e)-used);
}
double BudgetUsed(double stat,double residual){
    if(!(stat>=0)||!(residual>=0)||!std::isfinite(stat)||!std::isfinite(residual))return 1e100;
    return stat+residual;
}
double BudgetSum(const Engine& e,double stat,double residual){return BudgetUsed(stat,residual)+FixedBudgetReserve(e);}
ErrorMetric BudgetMetric(const Engine& e,const Statistics& stats,size_t row,const BudgetResidual& residual){
    ErrorMetric metric;double maximum=-1;
    for(int j=0;j<16;++j){
        const double raw=BudgetUsed(stats.interval[row][j],residual.matrix[j]);
        if(raw>maximum){maximum=raw;metric=ErrorMetric(j,false,stats.interval[row][j]);}
        if(e.ratios&&j){const double ratio=BudgetUsed(stats.ratioInterval[row][j],residual.ratio[j]);
            if(ratio>maximum){maximum=ratio;metric=ErrorMetric(j,true,stats.ratioInterval[row][j]);}}
    }return metric;
}
double MetricTotal(const Engine& e,const ErrorMetric& metric,const BudgetResidual& residual){
    return BudgetSum(e,metric.error,metric.ratio?residual.ratio[metric.element]:residual.matrix[metric.element]);
}
bool MetricWithinBudget(const Engine& e,const ErrorMetric& metric,const BudgetResidual& residual){
    return BudgetUsed(metric.error,metric.ratio?residual.ratio[metric.element]:residual.matrix[metric.element])<=NumericalBudget(e);
}
struct PointwiseDecision {double direct=0,output=0;bool directPass=false,outputPass=false;};
PointwiseDecision CheckPointwiseBudget(Engine& e,const Statistics& stats,const Statistics& output,
        const InterpolationResult& interp,const std::vector<size_t>& ids,const std::string& label){
    Require(stats.interval.size()==ids.size()&&output.interval.size()==e.theta.size(),"pointwise budget dimensions differ");
    PointwiseDecision result;result.directPass=result.outputPass=true;std::ostringstream rows;
    rows<<"theta_deg,checked_interpolation,direct_total95,output_total95,remaining_statistical_min";
    for(int j=0;j<16;++j)rows<<",output_total95_M"<<j/4+1<<j%4+1;
    for(int j=1;j<16;++j)rows<<",output_total95_M"<<j/4+1<<j%4+1<<"_over_M11";rows<<'\n';
    for(size_t q=0;q<e.theta.size();++q){const auto residual=ResidualAt(e,interp,ids,q);
        const auto outMetric=BudgetMetric(e,output,q,residual);result.outputPass=result.outputPass&&MetricWithinBudget(e,outMetric,residual);
        const double out=MetricTotal(e,outMetric,residual);result.output=std::max(result.output,out);
        rows<<Number(e.theta[q])<<','<<residual.measured<<',';
        if(residual.measured){const size_t r=std::lower_bound(ids.begin(),ids.end(),q)-ids.begin();
            const auto rawMetric=BudgetMetric(e,stats,r,residual);result.directPass=result.directPass&&MetricWithinBudget(e,rawMetric,residual);
            const double raw=MetricTotal(e,rawMetric,residual);result.direct=std::max(result.direct,raw);rows<<Number(raw);}
        rows<<','<<Number(out)<<','<<Number(RemainingStatistics(e,residual));
        for(int j=0;j<16;++j)rows<<','<<Number(BudgetSum(e,output.interval[q][j],residual.matrix[j]));
        for(int j=1;j<16;++j)rows<<','<<(e.ratios?Number(BudgetSum(e,output.ratioInterval[q][j],residual.ratio[j])):"");rows<<'\n';
    }
    if(label=="production")e.lastCombined=std::max(result.direct,result.output);
    else if(label.find("independent_serial")==0)e.lastIndependentCombined=std::max(result.direct,result.output);
    Atomic(e.root+"/pointwise_budget_"+label+".csv",rows.str());
    e.Log("pointwise_budget","{\"policy\":"+Quoted(BudgetPolicy())+",\"label\":"+Quoted(label)
        +",\"target\":"+Number(e.eps)+",\"fixed_reserve\":"+Number(FixedBudgetReserve(e))+",\"depth_reserve\":"+Number(e.depthBudget)
        +",\"direct_total95\":"+Number(result.direct)+",\"output_total95\":"+Number(result.output)
        +",\"direct_pass\":"+(result.directPass?"true":"false")+",\"output_pass\":"+(result.outputPass?"true":"false")
        +",\"unchecked_output_interpolation_reserve\":"+Number(e.interpolation)+",\"global_simultaneous_confidence\":false}");
    return result;
}

void PreparePointwisePolicy(Engine& e,AdaptiveState& state,bool resumed){
    const std::string path=e.root+"/error_budget_policy.json";
    const bool changed=!Exists(path)||Read(path).find(BudgetPolicy())==std::string::npos;
    e.pendingLegacyAuditSerial=-1; // Unversioned holdouts cannot certify this policy.
    if(changed&&resumed){state.streak=0;e.independentPassed=false;SaveAdaptive(e,state);
        e.Log("error_budget_policy_migration","{\"production_blocks_retained\":true,\"fresh_independent_audit_required\":true,\"pass_streak_reset\":true}");}
    e.Record(path,"{\"policy\":"+Quoted(BudgetPolicy())+",\"interpolation_cap\":"+Number(e.interpolation)+",\"depth_reserve_retained\":true,\"unchecked_output_uses_full_interpolation_reserve\":true,\"global_simultaneous_confidence\":false}");
}
