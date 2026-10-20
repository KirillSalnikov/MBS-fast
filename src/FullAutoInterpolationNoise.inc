// Scheduling only: the existing |mean residual| + t95 SE acceptance tests
// remain authoritative. Means and noise below always refer to the same
// component, stencil, counts and eight replicas.
void InitInterpolationNoise(InterpolationResult& result,size_t rows){
    result.matrixMeans.resize(rows);result.matrixNoise.resize(rows);
    result.ratioMeans.resize(rows);result.ratioNoise.resize(rows);result.growth.resize(rows);
}
void ObserveInterpolationNoise(const Engine& e,InterpolationResult& result,size_t row,int element,
        bool ratio,double upper,double mean,double noise){
    auto& means=ratio?result.ratioMeans:result.matrixMeans;
    auto& noises=ratio?result.ratioNoise:result.matrixNoise;
    // Retain a matched decomposition, never maxima from different stencils.
    if(mean+noise>=means[row][element]+noises[row][element]){
        means[row][element]=mean;noises[row][element]=noise;}
    auto& plan=result.growth[row];
    if(upper>e.interpolation&&mean>e.interpolation&&mean>noise){
        plan.biasDominated=true;plan.level=-1;plan.phi=0;plan.seconds=0;plan.reduction=0;
        plan.counts.clear();plan.observedCounts.clear();}
    if(upper>plan.upper){plan.upper=upper;plan.mean=mean;plan.noise=noise;
        plan.element=element;plan.ratio=ratio;plan.level=-1;plan.phi=0;plan.seconds=0;plan.reduction=0;
        plan.counts.clear();plan.observedCounts.clear();}
}
void PlanInterpolationNoise(const Engine& e,InterpolationGrowth& plan,const AdaptiveData& stage,
        const Budgets& observed,const std::vector<size_t>& ids,
        const std::vector<std::pair<size_t,double>>& weights,size_t guard,
        const Mueller& prediction,const Mueller& reference){
    plan.level=-1;plan.phi=0;plan.seconds=0;plan.reduction=0;plan.counts.clear();plan.observedCounts.clear();
    if(plan.biasDominated||plan.upper<=e.interpolation||!(reference[0]>0)
            ||!(prediction[0]>0)||!std::isfinite(plan.upper))return;
    // Exact first-order influence of the ratio of two means. Prediction and
    // direct value have separate denominators and polarization ratios.
    std::map<size_t,std::pair<double,double>> coefficients;
    for(const auto& w:weights){auto& c=coefficients[w.first];
        c.first+=w.second/(plan.ratio?prediction[0]:reference[0]);
        if(plan.ratio)c.second-=w.second*prediction[plan.element]/(prediction[0]*prediction[0]);}
    coefficients[guard].first-=1/reference[0];
    if(plan.ratio)coefficients[guard].second+=reference[plan.element]/(reference[0]*reference[0]);
    for(auto at=coefficients.begin();at!=coefficients.end();){
        if(at->second.first==0&&at->second.second==0)at=coefficients.erase(at);else ++at;}
    if(coefficients.empty())return;
    std::vector<std::map<int,std::array<double,8>>> centered(stage.levels.size());
    std::array<double,8> total{};
    for(size_t l=0;l<stage.levels.size();++l){
        for(const auto& c:coefficients){const auto at=std::lower_bound(ids.begin(),ids.end(),c.first);
            Require(at!=ids.end()&&*at==c.first,"missing interpolation influence row");const size_t row=at-ids.begin();
            auto& values=centered[l][e.PhiForTheta(c.first)];
            for(int s=0;s<8;++s)values[s]+=c.second.first*stage.levels[l][s][row][plan.element]
                +c.second.second*stage.levels[l][s][row][0];}
        for(auto& group:centered[l]){double mean=0;for(double value:group.second)mean+=value/8;
            for(int s=0;s<8;++s){group.second[s]-=mean;total[s]+=group.second[s];}}}
    struct Candidate {bool eligible=false;int level=0,phi=0;double variance=0,reduction=0,cost=0;
        std::vector<std::pair<size_t,int>> counts;std::vector<int> observed;};
    std::vector<Candidate> choices;bool positive=false;
    for(size_t l=0;l<stage.levels.size();++l){
        // Compare one phi cohort with the joint cohort. Quiet old-phi rows
        // need not grow just because a newly sampled phi is noisy. The joint
        // alternative retains useful cancellation between the phi groups.
        std::vector<int> phis{0};if(e.local&&centered[l].size()>1)for(const auto& group:centered[l])phis.push_back(group.first);
        for(int phi:phis){Candidate choice;choice.level=int(l);choice.phi=phi;bool capped=false;
            std::map<int,std::pair<double,int>> costs;
            std::vector<size_t> nodes;if(e.local){for(const auto& c:coefficients)if(!phi||e.PhiForTheta(c.first)==phi)nodes.push_back(c.first);}else nodes=ids;
            int peak=0;if(!e.local)for(size_t id:nodes)peak=std::max(peak,observed[l][id]);
            for(size_t id:nodes){const int n=observed[l][id],base=e.local?n:peak,next=base<=e.maxCount/2?2*base:n;
                if(next<=n){capped=true;break;}
                // Uniform-mode cost includes all evaluated rows, but each
                // guard retains only its small signed stencil. Applying the
                // proposal expands the selected level to the global peak.
                // Storing all rows per guard would be quadratic in grid size.
                if(e.local||coefficients.count(id)){
                    choice.counts.push_back({id,next});choice.observed.push_back(n);}
                const size_t row=std::lower_bound(ids.begin(),ids.end(),id)-ids.begin();
                auto& cost=costs[e.PhiForTheta(id)];cost.first=std::max(cost.first,stage.costs[l][row]);
                cost.second=std::max(cost.second,next-n);}
            if(capped)continue;
            std::array<double,8> values{};
            for(const auto& group:centered[l])if(!phi||group.first==phi)for(int s=0;s<8;++s)values[s]+=group.second[s];
            // Local cohorts double every participating observed prefix. For
            // uniform legacy mixed counts the forecast uses the least shrink;
            // actual acceptance is always recomputed, never forecast.
            const double change=1-std::sqrt(.5);double covariance=0;
            for(int s=0;s<8;++s){choice.variance+=values[s]*values[s]/7;covariance+=values[s]*total[s]/7;}
            choice.reduction=2*change*covariance-change*change*choice.variance;
            for(const auto& cost:costs)choice.cost+=8*cost.second.first*cost.second.second;
            choice.eligible=choice.variance>0&&std::isfinite(choice.variance)&&std::isfinite(choice.cost);
            positive=positive||(choice.eligible&&choice.reduction>0);choices.push_back(choice);}}
    // If every available cohort would destroy useful cross-level
    // cancellation, retain angular refinement rather than deliberately
    // increasing the predicted noise. Joint depth growth is not forecast here.
    if(!positive)return;
    double best=-1;
    for(const auto& c:choices){if(!c.eligible)continue;
        const double gain=std::max(0.,c.reduction);if(!(gain>0))continue;
        const double score=gain/std::max(1e-12,c.cost);if(score<=best)continue;
        best=score;plan.level=c.level;plan.phi=c.phi;plan.seconds=c.cost;plan.reduction=c.reduction;plan.counts=c.counts;plan.observedCounts=c.observed;}

}
void CopyInterpolationCheck(InterpolationResult& into,const InterpolationResult& from){
    // Keep the original output interpolant. Common-prefix checks replace only
    // its validation values and the scheduling diagnostics derived from them.
    into.errors=from.errors;into.matrixErrors=from.matrixErrors;into.ratioErrors=from.ratioErrors;
    into.maximum=from.maximum;into.bad=from.bad;into.matrixMeans=from.matrixMeans;
    into.matrixNoise=from.matrixNoise;into.ratioMeans=from.ratioMeans;into.ratioNoise=from.ratioNoise;into.growth=from.growth;
}
void WriteInterpolationNoise(Engine& e,const InterpolationResult& result,const std::vector<size_t>& ids,const std::string& folder){
    std::ostringstream out;out<<"theta_deg,metric,absolute_mean,noise95,upper95,decision,selected_level,selected_phi_0_means_joint,worker_cost_proxy\n";
    for(size_t r=0;r<ids.size();++r)for(int j=0;j<16;++j)for(int ratio=0;ratio<=(e.ratios&&j?1:0);++ratio){
        const auto& p=result.growth[r];const double upper=(ratio?result.ratioErrors:result.matrixErrors)[r][j];
        const char* decision=result.errors[r]<=e.interpolation?"passed":p.biasDominated?"refine_theta":p.level>=0?"sample_residual":"refine_theta_no_growth";
        out<<Number(e.theta[ids[r]])<<",M"<<j/4+1<<j%4+1<<(ratio?"/M11":"/scale(M11)")<<','
            <<Number((ratio?result.ratioMeans:result.matrixMeans)[r][j])<<','
            <<Number((ratio?result.ratioNoise:result.matrixNoise)[r][j])<<','<<Number(upper)<<','
            <<decision<<','<<p.level<<','<<p.phi<<','<<Number(p.seconds)<<'\n';}
    Atomic(folder+"/interpolation_noise.csv",out.str());
}
void RefineInterpolationFailures(Engine& e,AdaptiveState& state,const std::vector<size_t>& ids,
        const InterpolationResult& result,bool positive){
    const Budgets before=state.counts;Budgets proposed=before;
    std::set<size_t> knots(state.knots.begin(),state.knots.end());size_t sampled=0,refined=0;
    std::ostringstream rows;rows<<"theta_deg,decision,metric,absolute_mean,noise95,upper95,level,selected_phi_0_means_joint,node_theta,observed_count,next_count\n";
    for(size_t id:result.bad){const size_t r=std::lower_bound(ids.begin(),ids.end(),id)-ids.begin();
        const auto& p=result.growth[r];bool canGrow=positive&&!p.biasDominated&&p.level>=0&&!p.counts.empty();
        // A proposal for an overlapping stencil may already cover this one.
        // It is still progress relative to the observed completed prefix.
        if(canGrow){bool grows=false;for(const auto& n:p.counts)grows=grows||n.second>before[p.level][n.first];canGrow=grows;}
        const char* action=canGrow?"sample_residual":p.biasDominated?"refine_theta":"refine_theta_no_growth";
        auto log=[&](int level,size_t node,int observed,int next){
            rows<<Number(e.theta[id])<<','<<action<<",M"<<p.element/4+1<<p.element%4+1<<(p.ratio?"/M11":"/scale(M11)")
                <<','<<Number(p.mean)<<','<<Number(p.noise)<<','<<Number(p.upper)<<','<<level<<','<<p.phi<<','
                <<Number(e.theta[node])<<','<<observed<<','<<next<<'\n';};
        if(canGrow){++sampled;for(size_t k=0;k<p.counts.size();++k){const auto& n=p.counts[k];
                proposed[p.level][n.first]=std::max(proposed[p.level][n.first],n.second);
                log(p.level,n.first,p.observedCounts.at(k),n.second);}}
        else{refined+=knots.insert(id).second;log(-1,id,0,0);}
    }
    if(!positive){knots.insert(ids.begin(),ids.end());}
    state.knots.assign(knots.begin(),knots.end());
    ExtendPhiBudgets(e,proposed,ids,AdaptiveGuards(state.knots,e.theta.size()));
    if(!e.local){const auto nextIds=AdaptiveGuards(state.knots,e.theta.size());
        for(auto& level:proposed){int maximum=0;for(size_t id:nextIds)maximum=std::max(maximum,level[id]);for(size_t id:nextIds)level[id]=maximum;}}
    state.counts=proposed;state.streak=0;
    const auto filename="interpolation_refinement_round"+Number(state.round)+"_"+Hash(rows.str())+".csv";
    Atomic(e.root+"/"+filename,rows.str());
    e.Log("interpolation_refinement","{\"policy\":\"noise before angular refinement\",\"sampled_guards\":"+Number(sampled)
        +",\"added_theta_supports\":"+Number(refined)+",\"details\":"+Quoted(filename)
        +",\"growth_per_selected_stencil_level\":2,\"allocation_uses_observed_stencil_prefix\":true"
        ",\"accuracy_gates_unchanged\":true,\"variance_cost_is_scheduling_forecast\":true}");
    // No per-row cost limiter may cut this cohort apart: its noise forecast
    // relies on the chosen signed phi cohort receiving a larger common prefix.
    SaveAdaptive(e,state);
}
