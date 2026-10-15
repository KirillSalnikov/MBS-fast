// Frozen allocation; the worker computes real optical matrices. This harness
// isolates final validation from timing-dependent production adaptation.
#include "../src/FullAuto.cpp"
int main(int argc,const char*argv[]){try{
    using namespace FullAuto;
    ArgPP args;RegisterCliOptions(args);args.Parse(argc,argv);
    RunConfig config=RunConfig::FromCommandLine(args,false,false);
    Engine e(args,config);e.Init(argc,argv);const char*worker=getenv("TASK_MBS_WORKER");Require(worker,"missing worker");
    e.exe=worker;e.binaryHash=FileHash(e.exe);e.controls=false;e.depthSharing=false;e.workerPolicy="process";
    e.phi=args.GetIntValue("nphi");e.kernel.pipeline=false;e.kernel.block=64;e.kernel.warp="auto";
    e.seeds.assign(production,production+8);e.weightTheta={e.theta.front(),e.theta.back()};
    Directory(e.root+"/jobs");WriteGrid(e.root+"/requested_theta.csv",e.theta);
    // Real native short-run fixture with a fixed heterogeneous allocation.
    const std::vector<int> depths={1,2,config.maxReflections};AdaptiveState state;
    state.counts=Budgets(3,std::vector<int>(e.theta.size(),128));state.counts[0].assign(e.theta.size(),512);state.counts[1].assign(e.theta.size(),256);
    for(size_t i=0;i<e.theta.size();++i){if(i<e.theta.size()/3)state.counts[0][i]=1024;if(i==e.theta.size()/2)state.counts[0][i]=4096;}
    for(size_t i=0;i<e.theta.size();i+=8)state.knots.push_back(i);if(state.knots.back()!=e.theta.size()-1)state.knots.push_back(e.theta.size()-1);
    state.valid.assign(e.theta.size(),true);state.coarse.assign(e.theta.size(),true);state.previous.resize(e.theta.size());
    if(getenv("TASK_MBS_WARM_CACHE"))RestoreAdaptive(e,state);
    const auto ids=AdaptiveGuards(state.knots,e.theta.size());
    auto productionData=EvaluateBudget(e,depths,state.counts,ids,e.seeds.data());
    WriteMatrix(e.root+"/frozen_production.dat",Angles(e.theta,ids),Mean(Combine(productionData.levels)));
    // Establish explicit, fixed tolerances before reading any holdout.
    e.eps=.10;e.statistical=.03;e.interpolation=.025;e.finalAudit=true;e.ratios=true;e.local=true;e.localPolicy="auto";
    const char*policy=getenv("TASK_MBS_AUDIT_POLICY");if(policy&&std::string(policy)=="uniform"){e.local=false;e.localPolicy="off";}
    else if(policy&&std::string(policy)=="local")e.localPolicy="on";
    const double started=Seconds();const bool pass=IndependentAudit(e,state,depths,productionData,ids);
    Require(pass,"short physical independent audit failed");
    WriteGrid(e.root+"/theta_support.csv",Angles(e.theta,state.knots));
    e.State("audit_test_passed",0,PeakBudgets(state.counts),ids.size(),0,0,0);
    std::ostringstream polar;polar<<"theta_deg";for(int i=0;i<32;++i)polar<<" c"<<i;polar<<'\n';
    for(double theta:e.theta){polar<<Number(theta);for(int i=0;i<32;++i)polar<<" 0";polar<<'\n';}
    Atomic(e.root+"/polarization_weights.tsv",polar.str());
    Atomic(e.root+"/audit_seconds.json","{\"seconds\":"+Number(Seconds()-started)+",\"passed\":true}");
    return 0;
}catch(const std::exception&error){std::cerr<<error.what()<<'\n';return 1;}}
