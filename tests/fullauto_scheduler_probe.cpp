// Decision regressions: measured saved matrices plus adversarial numerical
// and cache cases. These tests do not label synthetic/replayed data an audit.
#include "../src/FullAuto.cpp"
using namespace FullAuto;
void Near(double a,double b,const char* message){Require(std::fabs(a-b)<=1e-10*std::max(1.,std::fabs(b)),message);}
void InitTest(Engine& e,const std::string& root,size_t rows,size_t levels,AdaptiveState& state){
    e.root=root;e.signature="scheduler_test";e.controls=false;e.local=true;e.localPolicy="off";e.depthSharing=false;
    e.eps=.02;e.statistical=.0074;e.interpolation=.005;e.initial=512;e.minimum=64;
    e.seeds.assign(production,production+8);e.theta.resize(rows);for(size_t i=0;i<rows;++i)e.theta[i]=170+double(i)*10/(rows-1);
    state.counts=Budgets(levels,std::vector<int>(rows,64));state.counts[0].assign(rows,512);
    state.valid.assign(rows,false);state.coarse.assign(rows,false);state.previous.resize(rows);
    for(size_t i=0;i<rows;++i)state.knots.push_back(i);Directory(root);Directory(root+"/jobs");
}
void ProfileJob(Engine& e,int count,int offset,int depth,int paired,unsigned seed,const std::vector<double>& angles,const std::string& id){
    Job j;e.SetSampling(j,count,seed,depth,angles,false);j.pairedDepth=paired;j.offset=offset;j.identity=id;j.folder=id;
    j.traceSeconds=count;j.diffractionSeconds=count*angles.size();j.seconds=j.traceSeconds+j.diffractionSeconds;e.CacheJob(j);
}
int main(int argc,char** argv){try{
    Require(argc==3,"expected output directory and saved level fixture");Directory(argv[1]);
    // Actual epyc1 replay: adding overlapping phases (2.2616 s) to a
    // 2.1188 s worker incorrectly reported zero launch/output overhead.
    Job timed;timed.pipeline=true;timed.seconds=2.1187703609466553;
    ReadPhaseTimings(timed,"phase\tseconds\tclock\ntrace_prepare\t1.733824431\tstage_wall_sum\n"
        "diffraction_and_controls\t0.527801716\tstage_wall_sum\nprocessing\t1.770477349\twall\n");
    Near(timed.NativeWall(),1.770477349,"overlapping phases substituted for processing wall");
    Near(timed.OutsideProcessing(),.3482930119466552,"lost launch/output overhead under overlap");
    Near(timed.traceSeconds+timed.diffractionSeconds-timed.OverlapFraction()*std::min(timed.traceSeconds,timed.diffractionSeconds),
        timed.NativeWall(),"overlap model does not reproduce measured work");
    Job legacy;ReadPhaseTimings(legacy,"phase seconds\ntrace_prepare 2\ndiffraction_and_controls 3\n");
    Near(legacy.NativeWall(),5,"legacy phase logs no longer load");
    ReadPhaseTimings(legacy,"phase seconds clock\nprocessing 9 summed_cpu_wall\nprocessing -1 wall\n");
    Require(legacy.processingSeconds<0,"invalid or CPU-summed processing clock accepted");
    for(size_t n:{size_t(2),size_t(3),size_t(4),size_t(7)}){std::vector<double>x(n);Matrix y(n);
        for(size_t i=0;i<n;++i){x[i]=i*i+.7*i;for(int j=0;j<16;++j)y[i][j]=std::sin(3*i+.4*j);}
        for(int r=0;r<=100;++r){const double q=x.back()*r/100;const auto ws=LocalCubicWeights(x,q);Mueller value{};double sum=0;
            for(const auto& w:ws){sum+=w.second;for(int j=0;j<16;++j)value[j]+=w.second*y[w.first][j];}
            Near(sum,1,"spline weights do not preserve constants");const auto reference=LocalCubicInterpolate(x,y,{q});for(int j=0;j<16;++j)Near(value[j],reference[0][j],"spline covariance uses a different interpolant");}
        if(n==7)Require(LocalCubicWeights(x,(x[2]+x[3])*.5).size()==5,"lost fifth support node");}
    ArgPP args;RunConfig config;Engine e(args,config);AdaptiveState state;InitTest(e,std::string(argv[1])+"/work",3,3,state);
    {Engine clockEngine(args,config);AdaptiveState clockState;InitTest(clockEngine,std::string(argv[1])+"/clock",3,1,clockState);
        clockEngine.kernel.pipeline=true;clockEngine.workQueuePolicy="on";clockEngine.orientationBlock=4096;
        clockEngine.SetSampling(timed,4096,999,14,{170,175},false);timed.identity="measured_pipeline";timed.folder="measured_pipeline";clockEngine.CacheJob(timed);
        AdditionalWork cost(clockEngine,clockEngine.seeds.data());const auto work=cost.Forecast({{4096,0,14,0,{170,175}}});
        Require(work.profiled&&work.jobs==8,"pipeline profile not used");Near(work.seconds,8*timed.seconds,"forecast double-counted overlap or lost outside-processing time");}
    e.workQueuePolicy="on";e.orientationBlock=64;
    ProfileJob(e,64,0,6,0,999,{170,175},"profile");
    std::vector<Engine::WorkRequest> requests={{128,0,6,0,{170,175}},{64,0,6,0,{175,180}}};
    {AdditionalWork cost(e,e.seeds.data());const auto cold=cost.Forecast(requests);Require(cold.profiled&&cold.jobs==16,"wrong cold job count");Near(cold.seconds,3584,"common trace counted twice");}
    for(int s=0;s<8;++s)ProfileJob(e,64,0,6,0,e.seeds[s],{170,175,180},"saved"+Number(s));
    {AdditionalWork cost(e,e.seeds.data());const auto warm=cost.Forecast(requests);Require(warm.jobs==8,"saved blocks counted as new jobs");Near(warm.seconds,1536,"saved orientations counted as new work");}
    for(int s=0;s<7;++s)ProfileJob(e,64,64,6,0,e.seeds[s],{170,175},"partial"+Number(s));
    {AdditionalWork cost(e,e.seeds.data());const auto partial=cost.Forecast(requests);Require(partial.jobs==1,"partial seed completion ignored");Near(partial.seconds,192,"partial restart cost wrong");}
    e.base={"changed_optical_recipe"};{AdditionalWork cost(e,e.seeds.data());const auto incompatible=cost.Forecast(requests);Require(!incompatible.profiled&&incompatible.jobs==16,"incompatible optical cache/profile accepted");}e.base.clear();
    ProfileJob(e,64,0,8,6,999,{170,175},"paired_profile");
    {AdditionalWork cost(e,e.seeds.data());auto paired=cost.Forecast({{64,0,8,6,{170,175}}});Require(paired.jobs==8,"paired depths counted twice");Near(paired.seconds,1536,"paired trace cost wrong");}
    // Dominant raw variance is not necessarily dominant ratio uncertainty.
    Engine ratioEngine(args,config);AdaptiveState ratioState;InitTest(ratioEngine,std::string(argv[1])+"/ratio",3,2,ratioState);
    AdaptiveData ratioData;ratioData.levels.assign(2,Samples(8,Matrix(3)));ratioData.costs.assign(2,std::vector<double>(3,1));
    const double a[8]={-1,-1,-1,-1,1,1,1,1},b[8]={-1,-1,1,1,-1,-1,1,1};
    for(int s=0;s<8;++s)for(int r=0;r<3;++r){ratioData.levels[0][s][r][0]=10+a[s];ratioData.levels[0][s][r][1]=.8*(10+a[s]);
        ratioData.levels[1][s][r][0]=.8*b[s];ratioData.levels[1][s][r][1]=-.8*b[s];}
    const auto rs=Estimate(Combine(ratioData.levels));const auto metric=WorstMetric(rs,1,true);Require(metric.ratio,"fixture does not isolate ratio error");
    Budgets proposed=ratioState.counts;AdditionalWork ratioCost(ratioEngine,ratioEngine.seeds.data());std::ostringstream log;log<<RefinementHeader();
    RefineError(ratioEngine,ratioCost,{6,8},ratioState.counts,proposed,{0,1,2},ratioData,{{1,1}},metric,rs.mean[1],175,"ratio_regression",log);
    Require(proposed[0]==ratioState.counts[0]&&proposed[1][1]==128,"selected raw variance instead of ratio influence");
    // Failed audit at one angle must not double other angles or quiet levels.
    Engine auditEngine(args,config);AdaptiveState auditState;InitTest(auditEngine,std::string(argv[1])+"/audit",3,3,auditState);
    AdaptiveData auditData;auditData.levels.assign(3,Samples(8,Matrix(3)));auditData.costs.assign(3,std::vector<double>(3,1));
    for(int s=0;s<8;++s)for(int r=0;r<3;++r){auditData.levels[0][s][r][0]=100;auditData.levels[1][s][r][0]=r==1?10*a[s]:0;}
    const auto stats=Estimate(Combine(auditData.levels));const auto interp=InterpolateBudget(auditEngine,auditState.knots,{0,1,2},Combine(auditData.levels));
    RefineFailedAudit(auditEngine,auditState,{6,8,10},auditData,{0,1,2},auditData,{0,1,2},Budgets({{512,512,512},{64,64,64},{64,64,64}}),stats,stats,interp);
    Require(auditState.counts==Budgets({{512,512,512},{64,128,64},{64,64,64}}),"audit failure multiplied unrelated work");
    const auto first=ReserveAuditSeeds(auditEngine,auditState);AdaptiveState restored;restored.counts=Budgets(3,std::vector<int>(3,64));restored.previous.resize(3);restored.valid.assign(3,false);restored.coarse.assign(3,false);
    Require(RestoreAdaptive(auditEngine,restored)&&restored.counts==auditState.counts,"refined budgets lost on resume");
    const auto second=ReserveAuditSeeds(auditEngine,restored);for(unsigned seed:first)Require(std::find(second.begin(),second.end(),seed)==second.end(),"resumed audit reused an exposed holdout");
    // An unexposed pending holdout is resumed with its frozen allocation;
    // a completed/failed report, changed production or corruption cannot reuse it.
    Engine pending(args,config);AdaptiveState pendingState;InitTest(pending,std::string(argv[1])+"/pending",3,3,pendingState);
    pending.orientationBlock=64;const auto reserved=ReserveAuditSeeds(pending,pendingState);
    pending.pendingLegacyAuditSerial=0;pending.pendingLegacyAuditRound=0;pending.pendingLegacyAuditSources=Hash(auditData.sources);
    int attempt=0;double fraction=.5;
    const auto imported=ReserveOrResumeAudit(pending,pendingState,auditData,{0,1,2},attempt,fraction);
    Require(imported.resumed&&imported.seeds==reserved&&imported.legacyFloor==64,"unexposed legacy audit was discarded");
    const auto again=ReserveOrResumeAudit(pending,pendingState,auditData,{0,1,2},attempt,fraction);
    Require(again.resumed&&again.seeds==reserved&&pendingState.auditSerial==1,"pending audit reseeded on restart");
    Atomic(pending.root+"/independent_final_validation_serial0.json","{\"passed\":false}");attempt=1;fraction=1;
    const auto fallback=ReserveOrResumeAudit(pending,pendingState,auditData,{0,1,2},attempt,fraction);
    Require(!fallback.resumed&&fallback.seeds!=reserved&&pendingState.auditSerial==2,"exposed holdout reused");
    attempt=0;fraction=.5;const auto fallbackResume=ReserveOrResumeAudit(pending,pendingState,auditData,{0,1,2},attempt,fraction);
    Require(fallbackResume.resumed&&attempt==1&&fraction==1&&fallbackResume.seeds==fallback.seeds,"fallback plan changed on resume");
    const auto reservationPath=pending.root+"/independent_audit_reservation.tsv";const auto payload=Read(reservationPath);Atomic(reservationPath,payload+"changed");
    bool corrupt=false;try{ReserveOrResumeAudit(pending,pendingState,auditData,{0,1,2},attempt,fraction);}catch(const std::exception&){corrupt=true;}
    Require(corrupt,"corrupt frozen audit allocation accepted");Atomic(reservationPath,payload);
    Engine warmed(args,config);AdaptiveState warmState;InitTest(warmed,std::string(argv[1])+"/warm_once",3,1,warmState);
    warmed.WarmCache();warmed.WarmCache();Require(warmed.warmedCacheControls.size()==1,"same cache domain scanned twice");
    warmed.controls=true;warmed.WarmCache();warmed.controls=false;warmed.WarmCache();
    Require(warmed.warmedCacheControls.size()==2,"distinct controls domain not verified");
    Engine decisions(args,config);AdaptiveState decisionState;InitTest(decisions,std::string(argv[1])+"/idempotent",3,1,decisionState);
    decisions.stablePasses=2;AdaptiveData decided;decided.sources="prefix_A";
    UpdatePassStreak(decisions,decisionState,true,decided,"snapshot_A");Require(decisionState.streak==1,"first decision not counted");
    UpdatePassStreak(decisions,decisionState,true,decided,"snapshot_A");Require(decisionState.streak==1,"same estimator counted twice");
    SaveAdaptive(decisions,decisionState);AdaptiveState reload;InitTest(decisions,decisions.root,3,1,reload);RestoreAdaptive(decisions,reload);decisions.statistical=.0074;
    UpdatePassStreak(decisions,reload,true,decided,"snapshot_A");Require(reload.streak==1,"resume inflated validation count");
    decided.sources="prefix_B";UpdatePassStreak(decisions,reload,true,decided,"snapshot_B");Require(reload.streak==2,"fresh prefix did not count");
    UpdatePassStreak(decisions,reload,true,decided,"snapshot_B");Require(reload.streak==2,"stable count exceeded requirement");
    UpdatePassStreak(decisions,reload,false,decided,"snapshot_B");Require(reload.streak==0,"failed decision did not reset count");
    const auto decisionPath=decisions.root+"/production_validation_receipt.tsv",decisionText=Read(decisionPath);Atomic(decisionPath,decisionText+"changed");
    bool decisionCorrupt=false;try{UpdatePassStreak(decisions,reload,true,decided,"snapshot_B");}catch(const std::exception&){decisionCorrupt=true;}
    Require(decisionCorrupt,"corrupt production decision receipt accepted");Atomic(decisionPath,decisionText);
    Engine legacyDecision(args,config);AdaptiveState old;InitTest(legacyDecision,std::string(argv[1])+"/legacy_decision",3,1,old);legacyDecision.stablePasses=2;old.streak=4;
    Atomic(legacyDecision.root+"/intermediate_current.json","{\"snapshot_folder\":\"legacy_snapshot\",\"numerical_round_passed\":true}");
    UpdatePassStreak(legacyDecision,old,true,decided,"legacy_snapshot");Require(old.streak==2,"legacy repeated passes not capped");
    RunConfig gpu=config;gpu.useGpu=true;Engine preparation(args,gpu);AdaptiveState unused;InitTest(preparation,std::string(argv[1])+"/preparation",3,1,unused);
    preparation.prepareThreadCap=8;preparation.workerThreads=1;Job trace;preparation.SetSampling(trace,1024,1,6,{0.,1.,2.,3.,4.},false);trace.pairedDepth=4;trace.identity="trace_profile";trace.folder=trace.identity;trace.traceSeconds=10;trace.diffractionSeconds=1;preparation.CacheJob(trace);
    Require(preparation.AutoPrepareThreads(16384,6,4,20)==8,"CPU-dominant job left serial");
    Require(preparation.AutoPrepareThreads(16384,6,4,1000)==0&&preparation.AutoPrepareThreads(16,6,4,20)==0,"GPU-bound or tiny job over-threaded");
    preparation.resourcePolicy="off";Require(preparation.AutoPrepareThreads(16384,6,4,20)==0,"resources off overridden");
    ArgPP fixed;RegisterCliOptions(fixed);const char* cli[]={"test","--threads","1"};fixed.Parse(3,cli);Engine pinned(fixed,gpu);pinned.prepareThreadCap=8;Require(pinned.AutoPrepareThreads(16384,6,4,20)==0,"explicit threads overridden");
    // Use real large-particle backscatter samples as failure-scheduling input.
    std::ifstream input(argv[2]);size_t rows;double tolerance;Require(bool(input>>rows>>tolerance),"invalid measured fixture");std::ostringstream replay;replay<<"theta,level,before,after\n";
    for(size_t r=0;r<rows;++r){Engine real(args,config);AdaptiveState saved;InitTest(real,std::string(argv[1])+"/real"+Number(r),2,3,saved);double theta;input>>theta;real.statistical=tolerance;
        AdaptiveData data;data.levels.assign(3,Samples(8,Matrix(2)));data.costs.assign(3,std::vector<double>(2));std::vector<int> ns;
        for(int l=0;l<3;++l){int n;double seconds;input>>n>>seconds;ns.push_back(n);saved.counts[l].assign(2,n);data.costs[l].assign(2,seconds);
            for(int s=0;s<8;++s)for(int j=0;j<16;++j){Require(bool(input>>data.levels[l][s][0][j]),"incomplete real fixture");data.levels[l][s][1][j]=data.levels[l][s][0][j];}}
        const auto st=Estimate(Combine(data.levels));const auto in=InterpolateBudget(real,saved.knots,{0,1},Combine(data.levels));
        RefineFailedAudit(real,saved,{14,16,18},data,{0,1},data,{0,1},Budgets({{ns[0],ns[0]},{ns[1],ns[1]},{ns[2],ns[2]}}),st,st,in);
        Require(saved.counts[0][0]==ns[0]&&saved.counts[1][0]==2*ns[1]&&saved.counts[2][0]==ns[2],"measured audit replay did not prioritize first correction");
        for(size_t l=0;l<3;++l)replay<<Number(theta)<<','<<l<<','<<ns[l]<<','<<saved.counts[l][0]<<'\n';}
    Atomic(std::string(argv[1])+"/real_failure_replay.csv",replay.str());
    std::cout<<"PASS: true spline stencil; ratio covariance; cached/paired/shared additional work; targeted audit; fresh seeds after resume; real backscatter scheduling\n";return 0;
}catch(const std::exception& error){std::cerr<<error.what()<<'\n';return 1;}}
