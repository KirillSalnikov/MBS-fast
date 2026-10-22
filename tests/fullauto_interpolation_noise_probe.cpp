// Decision-level regressions. Synthetic replicas isolate allocation failure
// modes; these fixtures are not physical accuracy or campaign speed claims.
#include "../src/FullAuto.cpp"
using namespace FullAuto;
namespace {
const double signs[8]={-1,-1,-1,-1,1,1,1,1};
const double orthogonal[8]={-1,-1,1,1,-1,-1,1,1};
void Close(double actual,double expected,const char* message){
    Require(std::isfinite(actual)&&std::isfinite(expected)&&
        std::fabs(actual-expected)<=2e-13*std::max(1.,std::fabs(expected)),message);
}
void Init(Engine& e,AdaptiveState& state,const std::string& root,size_t levels=2){
    e.root=root;e.signature="interpolation_noise_probe";e.controls=false;
    e.local=true;e.localPolicy="off";e.depthSharing=false;e.ratios=true;
    e.eps=.02;e.statistical=.015;e.interpolation=.005;
    e.initial=512;e.minimum=64;e.maxCount=16777216;e.phi=64;
    e.minPhi=1;e.maxPhi=64;e.fixedPhi=0;e.seeds.assign(production,production+8);
    e.theta={170,172.5,175,177.5,180};
    state.knots={0,2,4};state.counts=Budgets(levels,std::vector<int>(5,256));
    state.counts[0].assign(5,524288);state.previous.resize(5);
    state.valid.assign(5,false);state.coarse.assign(5,false);
    Directory(root);Directory(root+"/jobs");
}
AdaptiveData Quiet(size_t levels=2){
    AdaptiveData data;data.levels.assign(levels,Samples(8,Matrix(5)));
    data.costs.assign(levels,std::vector<double>(5,.01));
    data.costs[0].assign(5,10.);
    for(int s=0;s<8;++s)for(size_t r=0;r<5;++r){
        data.levels[0][s][r][0]=100;
        for(int j=1;j<16;++j)data.levels[0][s][r][j]=.1*j+r*.02;
    }
    return data;
}
void CheckUnchangedGates(const Engine& e,const AdaptiveState& state,const AdaptiveData& data,
                         const InterpolationResult& result){
    const auto samples=Combine(data.levels);
    const std::vector<size_t> ids={0,1,2,3,4};
    const auto supports=Support(samples,ids,state.knots);
    Samples predicted(8);for(int s=0;s<8;++s)
        predicted[s]=LocalCubicInterpolate(Angles(e.theta,state.knots),supports[s],e.theta);
    double maximum=0;std::vector<size_t> bad;
    for(size_t row=0;row<5;++row){
        Mueller mean{},reference{},prediction{};
        for(int s=0;s<8;++s)for(int j=0;j<16;++j){
            mean[j]+=(predicted[s][row][j]-samples[s][row][j])/8;
            reference[j]+=samples[s][row][j]/8;prediction[j]+=predicted[s][row][j]/8;
        }
        double worst=0;
        for(int j=0;j<16;++j){
            double variance=0;for(int s=0;s<8;++s)
                variance+=std::pow(predicted[s][row][j]-samples[s][row][j]-mean[j],2)/7;
            const double noise=t95*std::sqrt(variance/8)/reference[0];
            const double absMean=std::fabs(mean[j])/reference[0],upper=absMean+noise;
            Close(result.matrixErrors[row][j],upper,"raw acceptance gate changed");
            Close(result.matrixMeans[row][j],absMean,"raw mean diagnostic wrong");
            Close(result.matrixNoise[row][j],noise,"raw noise diagnostic wrong");
            worst=std::max(worst,upper);
            if(e.ratios&&j){
                const double a=prediction[j]/prediction[0],b=reference[j]/reference[0];
                double influenceVariance=0;
                for(int s=0;s<8;++s){
                    const double influence=(predicted[s][row][j]-a*predicted[s][row][0])/prediction[0]
                        -(samples[s][row][j]-b*samples[s][row][0])/reference[0];
                    influenceVariance+=influence*influence/7;
                }
                const double ratioNoise=t95*std::sqrt(influenceVariance/8);
                const double ratioUpper=std::fabs(a-b)+ratioNoise;
                Close(result.ratioErrors[row][j],ratioUpper,"ratio acceptance gate changed");
                Close(result.ratioMeans[row][j],std::fabs(a-b),"ratio mean diagnostic wrong");
                Close(result.ratioNoise[row][j],ratioNoise,"ratio noise diagnostic wrong");
                worst=std::max(worst,ratioUpper);
            }
        }
        Close(result.errors[row],worst,"pointwise maximum changed");
        maximum=std::max(maximum,worst);if(worst>e.interpolation)bad.push_back(row);
    }
    Close(result.maximum,maximum,"maximum interpolation gate changed");
    Require(result.bad==bad,"bad-node acceptance set changed");
}
void ReplayMeasuredGroup(const std::string& fixture,const std::string& root){
    std::ifstream in(fixture);size_t fullRows,rows,levels,replicas,supportCount,query,weightCount;
    Require(bool(in>>fullRows>>rows>>levels>>replicas>>supportCount>>query>>weightCount)&&
        fullRows==4296&&rows==20&&levels==3&&replicas==8,"invalid measured replay header");
    ArgPP args;RunConfig config;Engine e(args,config);AdaptiveState state;
    Init(e,state,root+"/real794",levels);e.theta.resize(fullRows);
    for(double& theta:e.theta)Require(bool(in>>theta),"missing measured theta");
    std::vector<size_t> ids(rows);Budgets observed(levels,std::vector<int>(fullRows));
    for(size_t r=0;r<rows;++r){
        int phi;Require(bool(in>>ids[r]>>phi),"missing measured row/phi");
        e.thetaPhi[ids[r]]=phi;
        for(auto& level:observed)Require(bool(in>>level[ids[r]])&&level[ids[r]]>0,"missing measured prefix");
    }
    state.knots.resize(supportCount);for(size_t& id:state.knots)Require(bool(in>>id),"missing measured support");
    AdaptiveData data;data.levels.assign(levels,Samples(8,Matrix(rows)));
    data.costs.assign(levels,std::vector<double>(rows));
    for(size_t l=0;l<levels;++l)for(size_t r=0;r<rows;++r){
        Require(bool(in>>data.costs[l][r]),"missing measured cost");
        for(int s=0;s<8;++s)for(double& value:data.levels[l][s][r])
            Require(bool(in>>value)&&std::isfinite(value),"missing measured matrix");
    }
    std::vector<std::pair<size_t,double>> weights(weightCount);
    for(auto& w:weights)Require(bool(in>>w.first>>w.second),"missing measured interpolation weight");
    std::string extra;Require(!(in>>extra),"extra measured fixture data");
    const auto samples=Combine(data.levels);
    const size_t qrow=std::lower_bound(ids.begin(),ids.end(),query)-ids.begin();
    Require(qrow<rows&&ids[qrow]==query,"missing measured query");
    Samples predictionSamples(8,Matrix(1));Mueller prediction{},reference{},delta{};
    for(const auto& w:weights){
        const size_t row=std::lower_bound(ids.begin(),ids.end(),w.first)-ids.begin();
        Require(row<rows&&ids[row]==w.first,"missing measured weight row");
        for(int s=0;s<8;++s)for(int j=0;j<16;++j)
            predictionSamples[s][0][j]+=w.second*samples[s][row][j];
    }
    for(int s=0;s<8;++s)for(int j=0;j<16;++j){
        prediction[j]+=predictionSamples[s][0][j]/8;
        reference[j]+=samples[s][qrow][j]/8;
        delta[j]+=(predictionSamples[s][0][j]-samples[s][qrow][j])/8;
    }
    InterpolationResult result;InitInterpolationNoise(result,1);
    for(int j=0;j<16;++j){
        double variance=0;for(int s=0;s<8;++s)
            variance+=std::pow(predictionSamples[s][0][j]-samples[s][qrow][j]-delta[j],2)/7;
        const double mean=std::fabs(delta[j])/reference[0],noise=t95*std::sqrt(variance/8)/reference[0];
        ObserveInterpolationNoise(e,result,0,j,false,mean+noise,mean,noise);
        if(j){
            const double a=prediction[j]/prediction[0],b=reference[j]/reference[0];
            double v=0;for(int s=0;s<8;++s){
                const double influence=(predictionSamples[s][0][j]-a*predictionSamples[s][0][0])/prediction[0]
                    -(samples[s][qrow][j]-b*samples[s][qrow][0])/reference[0];
                v+=influence*influence/7;
            }
            const double ratioMean=std::fabs(a-b),ratioNoise=t95*std::sqrt(v/8);
            ObserveInterpolationNoise(e,result,0,j,true,ratioMean+ratioNoise,ratioMean,ratioNoise);
        }
    }
    auto& plan=result.growth[0];
    Close(plan.upper,.5161183898529984,"measured residual bound changed");
    Close(plan.mean,.1501735629637896,"measured residual mean changed");
    Close(plan.noise,.3659448268892089,"measured residual noise changed");
    PlanInterpolationNoise(e,plan,data,observed,ids,weights,query,prediction,reference);
    Require(!plan.biasDominated&&plan.level>=0&&plan.phi==64&&!plan.counts.empty(),
        "measured noisy new-phi residual did not select phi64 growth");
    // The direct estimate must be scheduled even while this guard fails.
    // Its measured matrices and original sampling counts remain unchanged.
    state.counts=observed;
    InterpolationResult directResidual;directResidual.matrixErrors.resize(rows);
    directResidual.ratioErrors.resize(rows);
    const auto directStats=Estimate(samples);
    const auto directPlan=PlanPointwiseStatistics(e,state,data,ids,directStats,directResidual);
    Require(directPlan[1][query]>observed[1][query]||directPlan[2][query]>observed[2][query],
        "measured 794 noisy correction was starved by interpolation failure");
    Require(directPlan[0][query]==observed[0][query],
        "measured 794 correction unnecessarily grew the large base");
    std::ostringstream report;report<<"{\"fixture\":"<<Quoted(fixture)
        <<",\"theta_deg\":"<<Number(e.theta[query])<<",\"absolute_mean\":"<<Number(plan.mean)
        <<",\"noise95\":"<<Number(plan.noise)<<",\"upper95\":"<<Number(plan.upper)
        <<",\"metric\":\"M"<<plan.element/4+1<<plan.element%4+1<<(plan.ratio?"/M11":"/scale(M11)")
        <<"\",\"decision\":\"sample_residual\",\"selected_level\":"<<plan.level<<",\"selected_phi\":"<<plan.phi
        <<",\"worker_cost_proxy_seconds\":"<<Number(plan.seconds)
        <<",\"predicted_variance_reduction\":"<<Number(plan.reduction)
        <<",\"cost_is_measured_speedup\":false,\"growth\":[";
    for(size_t k=0;k<plan.counts.size();++k){
        const auto& node=plan.counts[k];
        Require(e.PhiForTheta(node.first)==64&&node.second==2*observed[plan.level][node.first],
            "measured proposal changed old phi or missed observed prefix");
        report<<(k?",":"")<<"{\"theta_deg\":"<<Number(e.theta[node.first])
            <<",\"observed_count\":"<<observed[plan.level][node.first]<<",\"next_count\":"<<node.second<<"}";
    }
    report<<"]}";Atomic(root+"/real794_residual_plan.json",report.str());
    std::cout<<"measured 794 replay: upper="<<Number(plan.upper)<<", mean="<<Number(plan.mean)
        <<", noise="<<Number(plan.noise)<<", selected level="<<plan.level<<", phi="<<plan.phi
        <<", growing nodes="<<plan.counts.size()<<"\n";
}
}
int main(int argc,char** argv){try{
    Require(argc==2||argc==3,"expected private output directory and optional measured fixture");
    const std::string root=argv[1];Directory(root);
    ArgPP args;RunConfig config;const std::vector<size_t> ids={0,1,2,3,4};
    size_t passed=0;
    // Every Mueller component must retain its acceptance bound while noisy
    // residuals grow the cheap correction rather than the quiet large base.
    for(int j=0;j<16;++j){
        Engine e(args,config);AdaptiveState state;Init(e,state,root+"/raw"+Number(j));
        e.ratios=false;auto data=Quiet();
        for(int s=0;s<8;++s)data.levels[1][s][1][j]=2*signs[s];
        const auto before=state.counts;const auto knots=state.knots;
        const auto result=InterpolateBudget(e,state.knots,ids,Combine(data.levels),&data,&state.counts);
        CheckUnchangedGates(e,state,data,result);
        Require(result.bad==std::vector<size_t>{1},"noise fixture affected another guard");
        Require(result.growth[1].level==1&&!result.growth[1].biasDominated&&
            result.growth[1].element==j,"noisy cheap correction was not selected");
        RefineInterpolationFailures(e,state,ids,result,true);
        Require(state.knots==knots&&state.counts[0]==before[0],"noise added theta or base work");
        for(size_t id:{size_t(0),size_t(1),size_t(2),size_t(4)})
            Require(state.counts[1][id]==512,"signed stencil did not grow together");
        Require(state.counts[1][3]==256,"unrelated guard multiplied");
        ++passed;
    }
    // Deterministic curvature calls for a support node, not more samples.
    {
        Engine e(args,config);AdaptiveState state;Init(e,state,root+"/bias");
        auto data=Quiet();for(int s=0;s<8;++s)data.levels[1][s][1][15]=2.;
        const auto before=state.counts;
        const auto result=InterpolateBudget(e,state.knots,ids,Combine(data.levels),&data,&state.counts);
        CheckUnchangedGates(e,state,data,result);
        Require(result.growth[1].biasDominated&&result.growth[1].level<0,"resolved bias sampled");
        RefineInterpolationFailures(e,state,ids,result,true);
        Require(std::binary_search(state.knots.begin(),state.knots.end(),size_t(1))&&
            state.counts==before,"deterministic bias did not add theta alone");
        ++passed;
    }
    // A smaller resolved bias cannot be hidden by another element's larger
    // statistical bound. Classification must inspect all failing metrics.
    {
        Engine e(args,config);AdaptiveState state;Init(e,state,root+"/mixed_metrics");
        auto data=Quiet();for(int s=0;s<8;++s){
            data.levels[1][s][1][1]=10*signs[s];
            data.levels[1][s][1][15]=1;
        }
        const auto result=InterpolateBudget(e,state.knots,ids,Combine(data.levels),&data,&state.counts);
        CheckUnchangedGates(e,state,data,result);
        Require(result.growth[1].element==1&&result.growth[1].biasDominated,
            "noise maximum hid resolved bias in another component");
        ++passed;
    }
    // Ratio noise can require a different level than raw intensity noise.
    // Level zero has large cost advantage, but its M12/M11 influence is zero.
    {
        Engine e(args,config);AdaptiveState state;Init(e,state,root+"/ratio");
        auto data=Quiet();data.costs[0].assign(5,.01);data.costs[1].assign(5,100.);
        for(int s=0;s<8;++s)for(size_t r=0;r<5;++r){
            const double base=100+(r==1?3*signs[s]:0);
            data.levels[0][s][r][0]=base;data.levels[0][s][r][1]=.8*base;
            data.levels[1][s][r][0]=r==1?4*orthogonal[s]:0;
            data.levels[1][s][r][1]=r==1?-4*orthogonal[s]:0;
        }
        const auto result=InterpolateBudget(e,state.knots,ids,Combine(data.levels),&data,&state.counts);
        CheckUnchangedGates(e,state,data,result);
        Require(result.growth[1].ratio&&result.growth[1].element==1&&result.growth[1].level==1,
            "ratio influence allocated the quiet polarization level");
        ++passed;
    }
    // Ratio influence needs separate prediction and guard denominators and
    // ratios. This tests the allocation primitive independently of bias choice.
    {
        Engine e(args,config);AdaptiveState state;Init(e,state,root+"/ratio_denominators");
        auto data=Quiet();data.costs[0].assign(5,.001);data.costs[1].assign(5,100.);
        for(int s=0;s<8;++s)for(size_t r=0;r<5;++r){
            const bool guard=r==1;
            const double base=(guard?100:120)+20*signs[s];
            data.levels[0][s][r][0]=base;data.levels[0][s][r][1]=(guard?.8:.5)*base;
            data.levels[1][s][r][0]=guard?4*orthogonal[s]:0;
            data.levels[1][s][r][1]=guard?-4*orthogonal[s]:0;
        }
        InterpolationGrowth plan;plan.ratio=true;plan.element=1;plan.upper=1;
        Mueller prediction{},reference{};prediction[0]=120;prediction[1]=60;
        reference[0]=100;reference[1]=80;
        PlanInterpolationNoise(e,plan,data,state.counts,ids,
            QueryWeights(e,state.knots,e.theta[1]),1,prediction,reference);
        Require(plan.level==1,"ratio allocator merged different denominators/ratios");
        ++passed;
    }
    // A large observed mean with still larger uncertainty is unresolved,
    // while another noisy mean below tolerance must also get more samples.
    for(double mean:{.3,3.}){
        Engine e(args,config);AdaptiveState state;Init(e,state,root+"/uncertain"+Number(mean));
        auto data=Quiet();for(int s=0;s<8;++s)data.levels[1][s][1][15]=mean+10*signs[s];
        const auto knots=state.knots;
        const auto result=InterpolateBudget(e,state.knots,ids,Combine(data.levels),&data,&state.counts);
        CheckUnchangedGates(e,state,data,result);
        Require(!result.growth[1].biasDominated&&result.growth[1].level==1,
            "unresolved residual mean was classified as interpolation bias");
        RefineInterpolationFailures(e,state,ids,result,true);
        Require(state.knots==knots,"unresolved residual inserted theta");
        ++passed;
    }
    // Common-prefix data have N=4096 even when production contains N=256.
    // Mixed phi rules keep their own prefixes; neither can be silently reset.
    {
        Engine e(args,config);AdaptiveState state;Init(e,state,root+"/actual_prefix");
        e.thetaPhi={{0,4},{2,4},{4,4}};auto data=Quiet();
        for(int s=0;s<8;++s)data.levels[1][s][1][15]=2*signs[s];
        state.counts[1]={8192,256,8192,4096,8192};
        Budgets observed=state.counts;observed[1]={8192,4096,8192,4096,8192};
        const auto result=InterpolateBudget(e,state.knots,ids,Combine(data.levels),&data,&observed);
        Require(result.growth[1].level==1,"actual-prefix fixture lost growth");
        const auto phi=e.thetaPhi;RefineInterpolationFailures(e,state,ids,result,true);
        Require(state.counts[1]==std::vector<int>({8192,8192,8192,4096,8192})&&e.thetaPhi==phi,
            "growth used production count or multiplied the quiet old phi prefix");
        AdaptiveState restored;restored.counts=Budgets(2,std::vector<int>(5,256));
        restored.previous.resize(5);restored.valid.assign(5,false);restored.coarse.assign(5,false);
        Require(RestoreAdaptive(e,restored)&&restored.counts==state.counts&&
            restored.knots==state.knots&&e.thetaPhi==phi,"resume lost residual allocation or phi");
        ++passed;
    }
    // A cancelling pair of phi rules may need joint growth: each isolated
    // shrink would increase the residual variance, whereas both together help.
    {
        Engine e(args,config);AdaptiveState state;Init(e,state,root+"/phi_cancellation");
        e.thetaPhi={{0,4},{2,4},{4,4}};auto data=Quiet();
        for(int s=0;s<8;++s)for(size_t r=0;r<5;++r)
            data.levels[1][s][r][15]=(r%2?9:10)*signs[s];
        const auto result=InterpolateBudget(e,state.knots,ids,Combine(data.levels),&data,&state.counts);
        CheckUnchangedGates(e,state,data,result);
        const auto& plan=result.growth[1];std::set<int> rules;
        for(const auto& node:plan.counts)rules.insert(e.PhiForTheta(node.first));
        Require(plan.level==1&&plan.reduction>0&&rules==std::set<int>({4,64}),
            "cross-phi cancellation was split despite negative isolated gains");
        ++passed;
    }
    // Covariance ranking must not choose the level that currently cancels
    // another level: growing that level alone would increase residual noise.
    {
        Engine e(args,config);AdaptiveState state;Init(e,state,root+"/cancellation");
        auto data=Quiet();data.costs[0].assign(5,1.);data.costs[1].assign(5,1e-10);
        for(int s=0;s<8;++s){
            data.levels[0][s][1][15]+=10*signs[s];
            data.levels[1][s][1][15]=-8*signs[s];
        }
        const auto result=InterpolateBudget(e,state.knots,ids,Combine(data.levels),&data,&state.counts);
        CheckUnchangedGates(e,state,data,result);
        Require(result.growth[1].level==0&&result.growth[1].reduction>0,
            "anti-correlated cheap level grew despite negative total gain");
        ++passed;
    }
    // If every single-level doubling would break cross-level cancellation,
    // the limited allocator must fall back rather than forecast improvement.
    {
        Engine e(args,config);AdaptiveState state;Init(e,state,root+"/no_positive_gain");
        auto data=Quiet();
        for(int s=0;s<8;++s){
            data.levels[0][s][1][15]+=10*signs[s];
            data.levels[1][s][1][15]=-9*signs[s];
        }
        const auto before=state.counts;
        const auto result=InterpolateBudget(e,state.knots,ids,Combine(data.levels),&data,&state.counts);
        Require(!result.growth[1].biasDominated&&result.growth[1].level<0,
            "negative covariance gain was scheduled as noise reduction");
        RefineInterpolationFailures(e,state,ids,result,true);
        Require(state.counts==before&&std::binary_search(state.knots.begin(),state.knots.end(),size_t(1)),
            "no-positive-gain fallback grew a harmful isolated level");
        ++passed;
    }
    // No eligible complete prefix remains at maxCount: insert the guard.
    {
        Engine e(args,config);AdaptiveState state;Init(e,state,root+"/cap");
        for(auto& level:state.counts)level.assign(5,e.maxCount);
        auto data=Quiet();for(int s=0;s<8;++s)data.levels[1][s][1][15]=2*signs[s];
        const auto result=InterpolateBudget(e,state.knots,ids,Combine(data.levels),&data,&state.counts);
        Require(result.growth[1].level<0,"capped prefix produced a growth plan");
        RefineInterpolationFailures(e,state,ids,result,true);
        Require(std::binary_search(state.knots.begin(),state.knots.end(),size_t(1)),
            "maxCount left an endlessly deferred guard");
        ++passed;
    }
    // A user maximum need not be a power of two. Never create a checkpoint
    // count at that fractional cap; fallback if a full doubling cannot fit.
    {
        Engine e(args,config);AdaptiveState state;Init(e,state,root+"/fractional_cap");
        e.maxCount=1000;state.counts[0].assign(5,512);state.counts[1].assign(5,512);
        auto data=Quiet();for(int s=0;s<8;++s)data.levels[1][s][1][15]=2*signs[s];
        const auto result=InterpolateBudget(e,state.knots,ids,Combine(data.levels),&data,&state.counts);
        Require(result.growth[1].level<0,"fractional cap admitted a non-power-of-two prefix");
        RefineInterpolationFailures(e,state,ids,result,true);
        Require(state.counts[0]==std::vector<int>(5,512)&&state.counts[1]==std::vector<int>(5,512),
            "fallback changed prefix at fractional cap");
        ++passed;
    }
    // Explicit uniform mode remains uniform, including unrelated guard rows.
    {
        Engine e(args,config);AdaptiveState state;Init(e,state,root+"/uniform");
        e.local=false;auto data=Quiet();
        for(int s=0;s<8;++s)data.levels[1][s][1][15]=2*signs[s];
        const auto knots=state.knots;
        const auto result=InterpolateBudget(e,state.knots,ids,Combine(data.levels),&data,&state.counts);
        RefineInterpolationFailures(e,state,ids,result,true);
        Require(state.knots==knots&&state.counts[1]==std::vector<int>(5,512)&&
            state.counts[0]==std::vector<int>(5,524288),"uniform allocation was silently localized");
        ++passed;
    }
    // Uniform work applies to the whole grid, but retained per-guard plans
    // must contain only their bounded stencil, not a quadratic grid copy.
    {
        Engine e(args,config);AdaptiveState state;Init(e,state,root+"/uniform_storage");
        e.local=false;const size_t n=257;e.theta.resize(n);std::vector<size_t> all(n);
        state.knots.clear();for(size_t r=0;r<n;++r){
            all[r]=r;e.theta[r]=170+10.*r/(n-1);if(r%2==0)state.knots.push_back(r);
        }
        state.counts=Budgets(2,std::vector<int>(n,256));state.counts[0].assign(n,524288);
        state.previous.resize(n);state.valid.assign(n,false);state.coarse.assign(n,false);
        AdaptiveData data;data.levels.assign(2,Samples(8,Matrix(n)));
        data.costs.assign(2,std::vector<double>(n,.01));data.costs[0].assign(n,10);
        for(int s=0;s<8;++s)for(size_t r=0;r<n;++r){
            data.levels[0][s][r][0]=100;
            if(r%2)data.levels[1][s][r][15]=2*signs[s];
        }
        const auto knots=state.knots;
        const auto result=InterpolateBudget(e,state.knots,all,Combine(data.levels),&data,&state.counts);
        Require(result.bad.size()==128,"large uniform storage fixture not exercised");
        size_t saved=0;for(size_t id:result.bad){
            const auto& plan=result.growth[id];
            Require(plan.level==1&&plan.counts.size()<=6&&
                plan.observedCounts.size()==plan.counts.size(),"uniform plan retained a grid-wide prefix copy");
            saved+=plan.counts.size();
        }
        Require(saved<=6*result.bad.size(),"uniform plan memory became quadratic");
        RefineInterpolationFailures(e,state,all,result,true);
        Require(state.knots==knots&&state.counts[1]==std::vector<int>(n,512)&&
            state.counts[0]==std::vector<int>(n,524288),"compact uniform plan changed actual allocation");
        ++passed;
    }
    // Complete cross-angle cancellation must not produce a noisy guard.
    {
        Engine e(args,config);AdaptiveState state;Init(e,state,root+"/cancel_angles");
        auto data=Quiet();for(int s=0;s<8;++s)for(size_t r=0;r<5;++r)
            data.levels[1][s][r][15]=20*signs[s]*(1+.1*r);
        const auto result=InterpolateBudget(e,state.knots,ids,Combine(data.levels),&data,&state.counts);
        CheckUnchangedGates(e,state,data,result);
        Require(result.bad.empty()&&result.maximum<1e-13,"paired angular cancellation lost");
        ++passed;
    }
    // Noisy support values can cancel completely from the interpolation
    // residual. A separate angular failure must not freeze their statistics.
    for(int j=0;j<16;++j){
        Engine e(args,config);AdaptiveState state;Init(e,state,root+"/starved"+Number(j));
        e.ratios=false;auto data=Quiet();
        for(int s=0;s<8;++s)for(size_t r=0;r<5;++r)data.levels[1][s][r][j]=5*signs[s];
        for(int s=0;s<8;++s)data.levels[0][s][1][15]+=3;
        const auto samples=Combine(data.levels);const auto stats=Estimate(samples);
        const auto interp=InterpolateBudget(e,state.knots,ids,samples,&data,&state.counts);
        Require(!interp.bad.empty()&&interp.errors[0]<1e-12,"starvation fixture needs a quiet support and failing guard");
        const auto before=state.counts;
        const auto proposal=PlanPointwiseStatistics(e,state,data,ids,stats,interp);
        Require(proposal[1][0]>before[1][0]&&proposal[0]==before[0],"noisy support did not grow its cheap correction");
        RefineInterpolationFailures(e,state,ids,interp,true);const auto interpolationCounts=state.counts;
        MergeGrowingBudgets(state.counts,proposal);
        Require(state.counts[1][0]>before[1][0],"interpolation refinement erased the statistical proposal");
        for(size_t l=0;l<before.size();++l)for(size_t r=0;r<ids.size();++r)
            Require(state.counts[l][r]>=interpolationCounts[l][r],"statistical merge reduced angular work");
        ++passed;
    }
    if(argc==3){ReplayMeasuredGroup(argv[2],root);++passed;}
    std::cout<<"passed "<<passed<<" interpolation-noise decision fixtures; unchanged all-element gates\n";
    return 0;
}catch(const std::exception& error){std::cerr<<error.what()<<'\n';return 1;}}
