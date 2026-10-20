// Inspect copied production matrices only. No native optical work is allowed.
#include "../src/FullAuto.cpp"
int main(int argc,const char* argv[]){try{
    using namespace FullAuto;const char* binary=getenv("TASK_NEW_BINARY");Require(binary,"TASK_NEW_BINARY required");
    if(argc==2&&std::string(argv[1])=="--version"){execl(binary,binary,"--version",static_cast<char*>(nullptr));return 127;}
    ArgPP args;RegisterCliOptions(args);args.Parse(argc,argv);RunConfig config=RunConfig::FromCommandLine(args,false,true);Engine e(args,config);e.Init(argc,argv);
    Require(Read(e.root+"/configuration.txt")==e.signature+"\n","production command fingerprint differs");
    e.root=getenv("TASK_CHECKPOINT_COPY");e.exe=binary;e.binaryHash=FileHash(binary);e.forbidNativeWork=true;
    e.warmedCacheControls.insert(false);e.warmedCacheControls.insert(true);Require(e.RestoreTuning(),"missing copied tuning");SelectAutomaticDepth(e);
    std::vector<int> depths;for(int d:{std::max(1,e.target-4),std::max(1,e.target-2),e.target})if(depths.empty()||depths.back()!=d)depths.push_back(d);
    AdaptiveState state;state.counts=Budgets(depths.size(),std::vector<int>(e.theta.size(),e.minimum));state.counts[0].assign(e.theta.size(),e.initial);state.previous=Matrix(e.theta.size());state.valid.assign(e.theta.size(),false);state.coarse.assign(e.theta.size(),false);Require(RestoreAdaptive(e,state),"copied state not restored");
    std::istringstream list(Read(e.root+"/selected_jobs.txt"));std::string folder;size_t count=0;std::set<std::string> visited;
    while(std::getline(list,folder)){Job j;const auto path=SamplingIndexPath(folder);std::istringstream index(DecodeSamplingIndex(path,Read(path)));std::string signature;
        Require(bool(index>>signature>>j.identity>>j.count>>j.offset>>j.seed>>j.depth>>j.phi>>j.block>>j.pipeline>>j.mirrored>>j.paired>>j.assembled>>j.hasControls>>j.pairedDepth>>j.warp>>j.weightHash>>j.executionHash)&&signature==e.signature,"copied matrix sampling differs");
        if(j.hasControls!=e.controls)continue;double x;while(index>>x)j.inputAngles.push_back(x);j.folder=folder;j.data=folder+"/result/result.dat";j.summary=folder+"/result/result_analytic_facets.tsv";
        ValidateStoredJob(folder,visited);std::istringstream meta(Read(folder+"/complete.txt"));std::string identity;meta>>identity>>j.seconds;Require(identity==j.identity,"copied identity differs");e.CacheJob(e.Load(j,false));++count;
    }
    const auto productionIds=AdaptiveGuards(state.knots,e.theta.size());auto ids=AdaptiveGuards(state.knots,e.theta.size(),true);if(e.theta.size()<=4097){ids.clear();for(size_t i=0;i<e.theta.size();++i)ids.push_back(i);}
    FullAutoAuditPlan::Extend(state.counts,productionIds,ids);const std::string queue=e.workQueuePolicy;e.workQueuePolicy="off";e.angularAuditPolicy="on";const bool sharing=e.depthSharing;e.depthSharing=false;
    const auto data=EvaluateBudget(e,depths,state.counts,productionIds,e.seeds.data());e.workQueuePolicy=queue;e.depthSharing=sharing;
    const bool local=UseLocalAudit(e,state,depths,ids,1);e.workQueuePolicy="off";e.depthSharing=false;
    const auto plan=PlanAngularAudit(e,state,depths,data,productionIds,ids,1,local);e.workQueuePolicy=queue;e.depthSharing=sharing;
    AdditionalWork forecast(e,e.seeds.data(),true,true);const auto full=forecast.Forecast(AuditRequests(e,depths,plan.anchors,ids,state.knots)),reduced=forecast.Forecast(AngularAuditRequests(e,depths,plan,ids,state.knots));
    Atomic(e.root+"/measured_plan.json","{\"production_matrices_verified\":"+Number(count)+",\"new_optical_workers\":0,\"support_theta\":"+Number(state.knots.size())+",\"checked_theta\":"+Number(ids.size())+",\"profiled\":"+(full.profiled&&reduced.profiled?"true":"false")+",\"direct_worker_seconds_forecast\":"+Number(full.seconds)+",\"residual_worker_seconds_forecast\":"+Number(reduced.seconds)+",\"forecast_ratio\":"+Number(reduced.seconds/std::max(1e-30,full.seconds))+"}\n");
    std::cout<<Read(e.root+"/measured_plan.json");return 0;
}catch(const std::exception& error){std::cerr<<error.what()<<'\n';return 1;}}
