#define FULLAUTO_STALL_HELPERS_ONLY
#include "fullauto_stall_probe.cpp"
void SchedulingMigration(const std::string& root){
    ArgPP args;RunConfig config;Engine e(args,config);SetupMemory(e,root,"/bin/true",false);
    e.theta={170,171,172,173,174};e.phi=2;e.initial=128;e.minimum=64;e.maxCount=65536;e.local=true;e.interpolation=.005;
    AdaptiveState s;s.knots={0,4};s.counts=Budgets(1,std::vector<int>(5,128));s.counts[0][0]=8192;s.counts[0][2]=8192;s.counts[0][4]=8192;
    s.previous=Matrix(5);s.valid.assign(5,false);s.coarse.assign(5,false);AdaptiveData data;data.sources="private-scheduling-fixture";
    Statistics stats;stats.mean=Matrix(3);stats.interval.resize(3);stats.ratioInterval.resize(3);for(auto& row:stats.interval)row[0]=.01;
    auto folder=SaveRoundDiagnostics(e,s,data,stats,{0,2,4},{0,.01,0},{170,172,174});
    Atomic(root+"/intermediate_current.json","{\"snapshot_folder\":"+Quoted(folder)+",\"numerical_round_passed\":false,\"interpolation_residual95\":0.01}");
    s.knots={0,2,4};for(auto& n:s.counts[0])n=16384;PrepareResidualBudget(e,s);s.interpolationCounts[0][2]*=2;
    const auto pending=s.interpolationCounts;Budgets baseline;Statistics priority;std::string path;bool passed=false,residual=false;
    Require(!PublishedRoundBudget(e,s,baseline,priority,path,passed),"strict acceptance hint allowed changed support");
    Require(PublishedSchedulingBudget(e,s,baseline,priority,path,passed,residual)&&residual&&!passed,"changed-support scheduling hint rejected");
    Require(baseline[0]==std::vector<int>({8192,128,8192,128,8192}),"new guards inherited expensive neighbours");
    // Completed eight-seed work above the new-row startup remains mandatory.
    for(unsigned seed:e.seeds){Job j;e.SetSampling(j,512,seed,6,{171,173},false);j.folder=root+"/metadata_"+Number(seed);j.identity="cached_"+Number(seed);e.CacheJob(j);}
    e.forbidNativeWork=true;BoundResumedRound(e,s,{6});
    Require(s.counts[0]==std::vector<int>({8192,512,8192,512,8192})&&s.interpolationCounts==pending,"resume discarded completed work or pending residual growth");
    Require(!Exists(root+"/final_outputs.fnv"),"scheduling diagnostics accepted a matrix");
    e.thetaPhi[2]=4;Require(PublishedSchedulingBudget(e,s,baseline,priority,path,passed,residual)&&baseline[0][2]==128,"old-phi count treated as current-phi completed work");
}
void SamplingOrder(const std::string& root){
    ArgPP args;RunConfig config;Engine e(args,config);SetupMemory(e,root,"/bin/true",false);
    e.theta={170,172,174};e.phi=2;e.initial=512;e.minimum=64;e.interpolation=.005;e.statistical=.01;e.local=true;
    AdaptiveState state;state.knots={0,2};state.counts=Budgets(1,std::vector<int>(3,8192));PrepareResidualBudget(e,state);
    Statistics stats;stats.mean=Matrix(3);stats.interval.resize(3);stats.ratioInterval.resize(3);for(auto& row:stats.interval)row[0]=.012;
    InterpolationResult interp;interp.maximum=.02;interp.matrixErrors.resize(3);interp.ratioErrors.resize(3);for(auto& row:interp.matrixErrors)row.fill(.02);
    AdaptiveData data;const auto proposal=PlanPointwiseStatistics(e,state,data,{0,1,2},stats,interp);
    Require(proposal==state.counts,"unresolved residual caused unnecessary production growth");
    InitInterpolationNoise(interp,3);ObserveInterpolationNoise(e,interp,1,0,false,.00852,.00561,.00291,true);
    Require(!interp.growth[1].biasDominated,"uncertain near-threshold mean triggered expensive angular refinement");
    ObserveInterpolationNoise(e,interp,2,0,false,.012,.01,.002,true);Require(interp.growth[2].biasDominated,"resolved angular bias was not refined");
    Budgets counts(2,std::vector<int>(5,32768));e.theta={170,171,172,173,174};
    ExtendProductionGuards(e,counts,{0,2,4},{0,1,2,3,4});
    Require(counts[0][1]==512&&counts[1][1]==64&&counts[0][0]==32768,"new guard startup corrupted old budgets");
    e.local=false;counts=Budgets(2,std::vector<int>(5,32768));ExtendProductionGuards(e,counts,{0,2,4},{0,1,2,3,4});
    Require(counts[0][1]==32768,"explicit uniform sampling was silently changed");
}
void CorrectionResume(const std::string& root){
    ArgPP args;RunConfig config;Engine e(args,config);SetupMemory(e,root,"/bin/true",false);
    e.theta={170,172,174};e.phi=2;e.initial=128;e.minimum=64;e.maxCount=65536;e.local=true;
    AdaptiveState s;s.knots={0,2};s.counts=Budgets(2,std::vector<int>(3,128));s.previous=Matrix(3);s.valid.assign(3,false);s.coarse.assign(3,false);
    Statistics stats;stats.mean=Matrix(3);stats.interval.resize(3);stats.ratioInterval.resize(3);
    AdaptiveData data;data.sources="unchanged-production";const auto folder=SaveRoundDiagnostics(e,s,data,stats,{0,1,2},{0,.02,0},e.theta);
    Atomic(root+"/intermediate_current.json","{\"snapshot_folder\":"+Quoted(folder)+",\"numerical_round_passed\":false,\"interpolation_residual95\":0.02}");
    PrepareResidualBudget(e,s);s.counts[0][0]=512;
    CorrectionProbeReceipt probe;probe.attempt=1;probe.sources=Hash(data.sources);probe.rows.push_back({1,1,128,256});SaveCorrectionProbe(e,s,probe);
    const auto pending=s.counts;e.forbidNativeWork=true;BoundResumedRound(e,s,{4,6});Require(s.counts==pending,"pending correction receipt production changed");
    probe.complete=true;SaveCorrectionProbe(e,s,probe);RestoreCorrectionProbe(e,s);BoundResumedRound(e,s,{4,6});
    Require(s.counts[1][1]>=256&&!RestoreCorrectionProbe(e,s),"completed correction lost during resume bound");
    // The widened subround can contain two correction depths and 4x growth.
    // Recovery commits each allocation once, without extending the base.
    e.theta.clear();for(size_t i=0;i<128;++i)e.theta.push_back(170+i/16.);
    s.round=17;s.counts=Budgets(3,std::vector<int>(128,512));s.previous=Matrix(128);s.valid.assign(128,false);s.coarse.assign(128,false);
    s.knots={0,127};s.interpolationCounts.clear();probe.rows.clear();probe.complete=true;
    for(size_t i=0;i<128;++i)probe.rows.push_back({i,1+i%2,512,2048});
    SaveCorrectionProbe(e,s,probe);const auto restored=ReadCorrectionProbe(e,s);
    Require(restored.rows.size()==128&&CorrectionProbeRequests(e,{8,10,12},restored).size()==2,"wide correction lost depth grouping");
    Require(RestoreCorrectionProbe(e,s)&&!RestoreCorrectionProbe(e,s),"wide correction replayed or failed recovery");
    for(size_t i=0;i<128;++i)Require(s.counts[0][i]==512&&s.counts[1+i%2][i]==2048&&s.counts[2-i%2][i]==512,"wide correction modified another level");
}
void CommonPhiScopeTest(const std::string& root){
    ArgPP args;RunConfig config;Engine e(args,config);SetupMemory(e,root,"/bin/true",false);
    e.theta={170,175,180};e.phi=32;e.thetaPhi[0]=4;e.thetaPhi[2]=64;const auto saved=e.thetaPhi;
    Require(ResidualPhi(e)==64&&MixedResidualPhi(e),"residual phi selection wrong");
    try{ResidualPhiScope scope(e,true);Require(e.thetaPhi.empty()&&e.phi==64,"mixed rules survived inside residual");throw std::runtime_error("test unwind");}catch(const std::runtime_error&){}
    Require(e.phi==32&&e.thetaPhi==saved,"residual phi scope leaked into production");
    {ResidualPhiScope scope(e,false);Require(e.phi==32&&e.thetaPhi==saved,"legacy audit changed phi");}
    e.local=true;e.minimum=256;e.initial=524288;
    AdaptiveState state;state.knots={0,2};state.counts=Budgets(2,std::vector<int>(3,1048576));
    state.interpolationCounts=Budgets(2,std::vector<int>(3,65536));
    state.previous=Matrix(3);state.valid.assign(3,false);state.coarse.assign(3,false);const auto production=state.counts;
    PrepareResidualBudget(e,state);
    Require(state.interpolationPhi==64&&PeakBudgets(state.interpolationCounts)==std::vector<int>({512,512})&&state.counts==production,
        "mixed-phi legacy noise budget infected new common residual");
    state.interpolationCounts[0][1]=1024;PrepareResidualBudget(e,state);Require(state.interpolationCounts[0][1]==1024,"new residual growth repeatedly reset");
    SaveAdaptive(e,state);AdaptiveState copy;copy.counts=production;copy.previous=Matrix(3);copy.valid.assign(3,false);copy.coarse.assign(3,false);
    Require(RestoreAdaptive(e,copy,false)&&copy.interpolationPhi==64,"schema6 lost residual phi");
    PrepareResidualBudget(e,copy);Require(copy.interpolationCounts==state.interpolationCounts,"resume reset completed common-phi growth");
    state.interpolationPhi=0;SaveAdaptive(e,state);Require(Read(root+"/adaptive_checkpoint.tsv").substr(0,2)=="5 ","legacy v5 fixture missing");
    Require(RestoreAdaptive(e,copy,false)&&copy.interpolationPhi==0,"legacy v5 restore failed");PrepareResidualBudget(e,copy);
    Require(copy.interpolationPhi==64&&copy.interpolationCounts[0][1]==512,"legacy v5 common-phi migration failed");
}
void LongMixedRound(const std::string& root){
    ArgPP args;RunConfig config;Engine e(args,config);SetupMemory(e,root,"/bin/true",false);e.orientationBlock=256;e.maxCount=65536;e.local=true;e.depthSharing=false;
    std::vector<size_t> ids;e.theta.clear();for(size_t i=0;i<32;++i){ids.push_back(i);e.theta.push_back(170+i*.25);}const std::vector<double> half(e.theta.begin(),e.theta.begin()+16);
    for(int offset=0;offset<16384;offset+=256)for(unsigned seed:e.seeds){const auto& grid=offset<8192?e.theta:half;Job j;e.SetSampling(j,256,seed,6,grid,false);
        j.offset=offset;j.identity=Number(offset)+"_"+Number(seed);j.folder=root+"/fake_"+j.identity;j.traceSeconds=40;j.diffractionSeconds=60*grid.size()/32.;j.seconds=j.traceSeconds+j.diffractionSeconds;e.CacheJob(j);}
    Budgets current(1,std::vector<int>(32,8192)),proposed(1,std::vector<int>(32,16384));Statistics stats;stats.interval.resize(32);stats.ratioInterval.resize(32);for(size_t r=0;r<32;++r)stats.interval[r][0]=.1+r*.01;
    AdditionalWork work(e,e.seeds.data());const auto original=work.Forecast(BudgetRequests(e,{6},proposed,ids));const auto full=proposed;
    BoundRoundWork(e,work,current,proposed,{6},ids,stats);const auto next=work.Forecast(BudgetRequests(e,{6},proposed,ids));
    Require(proposed!=full&&proposed!=current&&next.seconds<.7*original.seconds,"trace-sharing gate retained a multi-hour round");
    std::cout<<"PASS mixed tracing/diffraction: forecast "<<original.seconds<<" -> "<<next.seconds<<" worker seconds; complete power-of-two prefixes retained\n";
}
void InitialProduction(const std::string& root){
    ArgPP args;RegisterCliOptions(args);RunConfig config;Engine e(args,config);SetupMemory(e,root,"/bin/true",false);
    e.theta={0,5,10,15,20,25};e.initial=524288;e.minimum=256;e.orientationBlock=4096;e.maxCount=1048576;e.local=true;e.depthSharing=false;
    AdaptiveState state;state.knots={0,5};state.counts=Budgets(1,std::vector<int>(6,e.initial));state.valid.assign(6,false);
    const auto ids=AdaptiveGuards(state.knots,e.theta.size());const auto angles=Angles(e.theta,ids);
    // One expensive wide-grid calibration: the 524288-start must not be run
    // before the first estimate, and unprofiled levels must not look free.
    Job timing;e.SetSampling(timing,4096,9999,6,angles,false);timing.identity="timing";timing.folder=root+"/timing";
    timing.seconds=100;timing.traceSeconds=10;timing.diffractionSeconds=90;e.CacheJob(timing);
    e.forbidNativeWork=true;BoundInitialProduction(e,state,{6});
    Require(PeakBudgets(state.counts)[0]<e.initial&&PeakBudgets(state.counts)[0]>=512,"first-round forecast did not bound expensive start");
    state.counts[0].assign(6,e.initial);config.adaptive.loadedFromFile=true;
    BoundInitialProduction(e,state,{6});Require(PeakBudgets(state.counts)[0]==e.initial,"adaptive-config start replaced");config.adaptive.loadedFromFile=false;
    e.local=false;BoundInitialProduction(e,state,{6});Require(PeakBudgets(state.counts)[0]==e.initial,"explicit uniform start replaced");e.local=true;
    const char* av[]={"probe","--fullauto-initial","524288"};args.Parse(3,av);
    BoundInitialProduction(e,state,{6});Require(PeakBudgets(state.counts)[0]==e.initial,"explicit start replaced");
    ArgPP automatic;RunConfig other;Engine resumed(automatic,other);SetupMemory(resumed,root+"/resume","/bin/true",false);
    resumed.theta=e.theta;resumed.initial=e.initial;resumed.minimum=256;resumed.orientationBlock=4096;resumed.maxCount=e.maxCount;resumed.local=true;resumed.depthSharing=false;
    // Already completed native leaves have not yet acquired a parent prefix.
    for(int offset=0;offset<16384;offset+=4096)for(unsigned seed:resumed.seeds){Job j;
        resumed.SetSampling(j,4096,seed,6,angles,false);j.offset=offset;j.identity=Number(offset)+"_"+Number(seed);j.folder=root+"/leaf_"+j.identity;
        j.seconds=100;j.traceSeconds=10;j.diffractionSeconds=90;resumed.CacheJob(j);}
    state.counts[0].assign(6,resumed.initial);BoundInitialProduction(resumed,state,{6});
    for(size_t id:ids)Require(state.counts[0][id]>=16384,"complete eight-seed native leaves were discarded");
    Require(PeakBudgets(state.counts)[0]<resumed.initial,"unfinished initial round remained unbounded");
    state.round=1;state.counts[0].assign(6,resumed.initial);BoundInitialProduction(resumed,state,{6});
    Require(PeakBudgets(state.counts)[0]==resumed.initial,"ordinary production allocation replaced by startup pilot");
    Require(!Exists(root+"/final_outputs.fnv"),"startup scheduling accepted accuracy");
    std::cout<<"PASS first production cost bound, explicit/config starts, native-leaf eight-seed resume floor and unchanged later rounds\n";
}
int main(int argc,char** argv){try{Require(argc==2,"OUTPUT");const std::string r=argv[1];SchedulingMigration(r+"/migration");SamplingOrder(r+"/order");CorrectionResume(r+"/correction");CommonPhiScopeTest(r+"/phi");LongMixedRound(r+"/mixed");InitialProduction(r+"/initial");
    PartitionWork(r+"/pure_trace",true);PartitionWork(r+"/diffraction",false);
    std::cout<<"PASS residual-first scheduling, changed-support resume, eight-seed complete-prefix floor, explicit uniform mode and confirmed-bias refinement\n";return 0;
}catch(const std::exception& error){std::cerr<<error.what()<<'\n';return 1;}}
