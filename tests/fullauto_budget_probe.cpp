// Real saved confidence arrays can be replayed without running optics again.
// Adversarial cases exercise missing guards, all components, ratios and NaNs.
#include "../src/FullAuto.cpp"
using namespace FullAuto;
void NearBudget(double a,double b){Require(std::fabs(a-b)<1e-11,"pointwise replay differs from independent NumPy result");}
int main(int argc,char** argv){try{
    Require(argc==3,"expected replay fixture and output directory");ArgPP args;RunConfig config;Engine e(args,config);
    e.root=argv[2];Directory(e.root);e.eps=.02;e.statistical=.00951171875;e.interpolation=.005;e.depthBudget=.005;e.ratios=true;
    BudgetResidual residual;residual.matrix[0]=.001;
    Statistics sample;sample.interval.resize(1);sample.ratioInterval.resize(1);sample.mean.resize(1);sample.mean[0][0]=1;
    sample.interval[0][0]=.011;
    Require(MetricTotal(e,BudgetMetric(e,sample,0,residual),residual)<e.eps,"feasible combined error rejected by old fixed split");
    residual.matrix[0]=.004;Require(MetricTotal(e,BudgetMetric(e,sample,0,residual),residual)>e.eps,"sum above target accepted");
    residual.matrix[0]=0;sample.ratioInterval[0][5]=.02;
    Require(BudgetMetric(e,sample,0,residual).ratio,"ratio failure hidden by raw matrix");
    e.ratios=false;Require(MetricTotal(e,BudgetMetric(e,sample,0,residual),residual)<e.eps,"explicit ratio-off ignored");e.ratios=true;
    sample.ratioInterval[0][5]=0;
    for(int j=0;j<16;++j){sample.interval[0][j]=.02;Require(MetricTotal(e,BudgetMetric(e,sample,0,residual),residual)>e.eps,"matrix component skipped");sample.interval[0][j]=0;}
    sample.interval[0][7]=std::numeric_limits<double>::quiet_NaN();Require(MetricTotal(e,BudgetMetric(e,sample,0,residual),residual)>e.eps,"NaN accepted");
    NearBudget(FixedBudgetReserve(e),.00548828125);
    // Tiny budgets must not disappear when added to a much larger reserve.
    {const double saved=e.statistical;e.statistical=1e-30;const double interpolation=e.interpolation;e.interpolation=1e-30;
        sample.interval[0].fill(0);sample.interval[0][0]=1e-20;
        Require(!MetricWithinBudget(e,BudgetMetric(e,sample,0,residual),residual),"roundoff erased a failing tiny numerical budget");
        e.statistical=saved;e.interpolation=interpolation;}
    InterpolationResult empty;const auto missing=ResidualAt(e,empty,{},0);NearBudget(missing.matrix[3],.005);Require(!missing.measured,"unchecked theta treated as measured");
    {e.signature="budget_migration_test";e.theta={0,1};Directory(e.root+"/jobs");AdaptiveState state;
        state.knots={0,1};state.counts=Budgets(1,std::vector<int>(2,512));state.previous.resize(2);state.valid.assign(2,true);state.coarse.assign(2,true);
        state.streak=2;state.auditSerial=7;e.pendingLegacyAuditSerial=6;
        Atomic(e.root+"/error_budget_policy.json","{\"policy\":\"fixed-split-legacy\"}");
        PreparePointwisePolicy(e,state,true);
        Require(state.streak==0&&state.auditSerial==7&&e.pendingLegacyAuditSerial<0,"old holdout or pass streak survived changed acceptance policy");
        Require(state.counts[0][0]==512&&state.valid[0]&&state.coarse[0],"policy migration discarded working samples");
        const auto checkpoint=Read(e.root+"/adaptive_checkpoint.tsv");Require(checkpoint.find("budget_migration_test active 0 0 ")!=std::string::npos,"policy reset was not durable");
        state.streak=1;PreparePointwisePolicy(e,state,true);Require(state.streak==1&&state.auditSerial==7,"same-policy resume reset progress");}
    std::ifstream in(argv[1]);size_t count,checked;Require(bool(in>>count>>checked),"invalid replay header");
    e.theta.resize(count);for(auto& x:e.theta)in>>x;std::vector<size_t> ids(checked);for(auto& id:ids)in>>id;
    Statistics direct,output;direct.interval.resize(checked);direct.ratioInterval.resize(checked);output.interval.resize(count);output.ratioInterval.resize(count);
    InterpolationResult interp;interp.matrixErrors.resize(checked);interp.ratioErrors.resize(checked);
    for(auto* matrix:{&direct.interval,&direct.ratioInterval,&output.interval,&output.ratioInterval,&interp.matrixErrors,&interp.ratioErrors})for(auto& row:*matrix)for(auto& v:row)in>>v;
    double expectedDirect,expectedOutput;size_t expectedFailures;Require(bool(in>>expectedDirect>>expectedOutput>>expectedFailures),"incomplete replay fixture");
    const auto decision=CheckPointwiseBudget(e,direct,output,interp,ids,"replay");NearBudget(decision.direct,expectedDirect);NearBudget(decision.output,expectedOutput);
    size_t failures=0;for(size_t r=0;r<ids.size();++r){const auto residual=ResidualAt(e,interp,ids,ids[r]);
        const auto raw=MetricTotal(e,BudgetMetric(e,direct,r,residual),residual),out=MetricTotal(e,BudgetMetric(e,output,ids[r],residual),residual);
        failures+=std::max(raw,out)>e.eps;}
    Require(failures==expectedFailures,"wrong failed-angle count");
    std::cout<<"PASS: all components, ratios, unchecked guards, fixed depth reserve, NaN; replay "<<checked<<" checked / "<<count<<" output angles, "<<failures<<" failed; direct="<<Number(decision.direct)<<" output="<<Number(decision.output)<<'\n';return 0;
}catch(const std::exception& error){std::cerr<<error.what()<<'\n';return 1;}}
