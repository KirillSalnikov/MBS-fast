// Real optical integration on a frozen grid; all 16 Mueller components and
// ratios, the equal-budget identity, estimator covariance and planning isolation.
#include "../src/FullAuto.cpp"
int main(int argc,const char* argv[]){try{
    using namespace FullAuto;ArgPP args;RegisterCliOptions(args);args.Parse(argc,argv);
#ifdef USE_CUDA
    const bool cuda=true;
#else
    const bool cuda=false;
#endif
    RunConfig config=RunConfig::FromCommandLine(args,false,cuda);Engine e(args,config);e.Init(argc,argv);
    const char* worker=getenv("TASK_MBS_WORKER");Require(worker,"TASK_MBS_WORKER required");e.exe=worker;e.binaryHash=FileHash(worker);
    WriteGrid(e.root+"/requested_theta.csv",e.theta);Directory(e.root+"/jobs");signal(SIGTERM,StopWorkers);signal(SIGINT,StopWorkers);
    e.phi=args.GetIntValue("nphi");e.kernel.pipeline=args.IsCatched("orientation_pipeline");e.depthSharing=false;e.SetMirror(false);
    if(getenv("TASK_MBS_WARM_CACHE"))e.WarmCache();
    AdaptiveState state;for(size_t i=0;i<e.theta.size();i+=4)state.knots.push_back(i);if(state.knots.back()!=e.theta.size()-1)state.knots.push_back(e.theta.size()-1);
    std::vector<size_t> ids;for(size_t i=0;i<e.theta.size();++i)ids.push_back(i);
    const std::vector<int> depths={std::max(1,config.maxReflections-2),config.maxReflections};
    state.counts=Budgets(depths.size(),std::vector<int>(e.theta.size(),e.initial));state.counts[1].assign(e.theta.size(),std::max(16,e.initial/2));
    const auto productionData=EvaluateBudget(e,depths,state.counts,ids,e.seeds.data());
    const auto trained=PlanAngularAudit(e,state,depths,productionData,ids,ids,1,true);
    Require(trained.enabled,"forced angular audit did not train");
    const unsigned fresh[8]={1000000000,1000000017,1000000034,1000000051,1000000068,1000000085,1000000102,1000000119};
    AngularAuditPlan exact=trained;exact.residuals=state.counts;exact.anchors=state.counts;
    const auto direct=EvaluateBudget(e,depths,state.counts,ids,fresh);
    const auto assembled=EvaluateAngularAudit(e,depths,exact,ids,state.knots,fresh);
    const auto a=Combine(direct.levels),b=Combine(assembled.levels);double all16=0,ratio=0;
    for(int s=0;s<8;++s)for(size_t r=0;r<ids.size();++r){const double scale=std::max(1e-100,std::fabs(a[s][r][0]));for(int j=0;j<16;++j)all16=std::max(all16,std::fabs(a[s][r][j]-b[s][r][j])/scale);
        if(a[s][r][0]>0&&b[s][r][0]>0)for(int j=1;j<16;++j)ratio=std::max(ratio,std::fabs(a[s][r][j]/a[s][r][0]-b[s][r][j]/b[s][r][0]));}
    Require(std::max(all16,ratio)<1e-11,"equal-prefix angular estimator differs from direct optics");
    const auto repeated=PlanAngularAudit(e,state,depths,productionData,ids,ids,1,true);
    Require(trained.identity==repeated.identity&&trained.anchors==repeated.anchors&&trained.residuals==repeated.residuals,"holdout changed frozen production allocation");
    const auto reduced=EvaluateAngularAudit(e,depths,trained,ids,state.knots,fresh);const auto z=Combine(reduced.levels);const auto stats=Estimate(z);
    WriteMatrix(e.root+"/angular_matrix.dat",e.theta,Mean(z));WriteMatrix(e.root+"/direct_matrix.dat",e.theta,Mean(a));
    SaveReplicaMatrices(e,"angular",z,e.theta,fresh);SaveReplicaMatrices(e,"direct",a,e.theta,fresh);
    Atomic(e.root+"/angular_probe.json","{\"equal_prefix_all16_scaled_M11\":"+Number(all16)+",\"equal_prefix_ratio_difference\":"+Number(ratio)
        +",\"frozen_plan_unchanged_after_holdout\":true,\"support_theta\":"+Number(state.knots.size())+",\"checked_theta\":"+Number(ids.size())+",\"reduced_statistical95\":"+Number(StatisticalError(e,stats))+"}\n");
    if(getenv("TASK_FORCE_AUDIT_REJECTION")){
        state.previous.assign(e.theta.size(),Mueller{});state.valid.assign(e.theta.size(),false);state.coarse.assign(e.theta.size(),false);state.auditSerial=32;
        e.statistical=1e-30;e.interpolation=1e-30;
        Require(!IndependentAudit(e,state,depths,productionData,ids),"impossible tolerance unexpectedly passed");
        Require(state.auditSerial==34,"failed angular audit did not reserve a fresh direct seed set");
        const auto first=Read(e.root+"/independent_final_validation_serial32.json"),second=Read(e.root+"/independent_final_validation_serial33.json");
        Require(first.find("angular_residual_v1")!=std::string::npos&&second.find("direct")!=std::string::npos,"angular rejection did not use direct fallback");
        Atomic(e.root+"/fallback_probe.json","{\"passed\":true,\"first_serial\":32,\"fresh_direct_serial\":33}\n");
    }
    std::cout<<Read(e.root+"/angular_probe.json");return 0;
}catch(const std::exception& error){std::cerr<<error.what()<<'\n';return 1;}}
