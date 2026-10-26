#define FULLAUTO_STALL_HELPERS_ONLY
#include "fullauto_stall_probe.cpp"
void StopPolicy(){
    for(bool final:{false,true})for(bool loaded:{false,true})for(bool explicitFlag:{false,true}){
        ArgPP args;RegisterCliOptions(args);const char* av[]={"probe","--convergence-passes","3"};if(explicitFlag)args.Parse(3,av);
        RunConfig c;c.adaptive.stablePasses=explicitFlag?3:2;c.adaptive.loadedFromFile=loaded;
        Engine e(args,c);e.finalChecks=final;
        Require(e.RequiredProductionPasses()==(!final&&!loaded&&!explicitFlag?1:c.adaptive.stablePasses),"explicit/config pass policy lost");
    }
    ArgPP args;RunConfig c;c.adaptive.loadedFromFile=true;c.adaptive.stablePasses=0;Engine e(args,c);bool rejected=false;
    try{e.RequiredProductionPasses();}catch(const std::exception&){rejected=true;}Require(rejected,"invalid pass count accepted");
    e.gpus={"0","1"};e.dispatchGpus={"1"};Require(RoundPhysicalGpus(e)==1,"original GPU selection overrode dispatch");
    e.dispatchGpus={"1","1"};Require(RoundPhysicalGpus(e)==1,"duplicate workers counted as two physical GPUs");
    e.dispatchGpus={"0","1","1"};Require(RoundPhysicalGpus(e)==2,"distinct physical GPUs lost");
}
void SavedPass(const std::string& root){
    ArgPP args;RunConfig config;Engine e(args,config);SetupMemory(e,root,"/bin/true",false);
    e.theta={170,175,180};e.initial=64;e.minimum=32;e.maxCount=4096;e.finalChecks=false;e.stablePasses=1;e.phi=2;
    AdaptiveState state;state.knots={0,2};state.counts=Budgets(2,std::vector<int>(3,128));
    state.valid.assign(3,true);state.coarse.assign(3,true);state.previous=Matrix(3);state.streak=1;
    AdaptiveData stage;stage.sources="saved-source-graph\n";Statistics stats;stats.mean=Matrix(3);stats.interval.resize(3);stats.ratioInterval.resize(3);
    const auto folder=SaveRoundDiagnostics(e,state,stage,stats,{0,1,2},{0,0,0},e.theta);
    Atomic(root+"/intermediate_current.json","{\"snapshot_folder\":"+Quoted(folder)+",\"numerical_round_passed\":true}");
    const auto saved=state.counts;for(auto& l:state.counts)for(auto& n:l)n*=2;
    Budgets loaded;Statistics priorities;std::string path;bool passed=false;
    Require(PublishedRoundBudget(e,state,loaded,priorities,path,passed)&&passed&&loaded==saved,"last published allocation not read");
    const auto previous=state.previous;RevisitPublishedPass(e,state);
    Require(state.counts==saved&&state.previous==previous&&state.streak==1,"saved pass was forged or allocation not revisited");
    Require(!Exists(root+"/final_outputs.fnv"),"diagnostics accepted a matrix");
    for(auto& l:state.counts)for(auto& n:l)n*=2;const auto original=state.counts;
    e.stablePasses=2;RevisitPublishedPass(e,state);Require(state.counts==original,"explicit repeated production pass bypassed");e.stablePasses=1;
    e.thetaPhi[1]=4;Require(!PublishedRoundBudget(e,state,loaded,priorities,path,passed),"changed phi reused snapshot");e.thetaPhi.clear();
    state.knots={0,1,2};Require(!PublishedRoundBudget(e,state,loaded,priorities,path,passed),"changed support reused snapshot");state.knots={0,2};
    Atomic(root+"/intermediate_current.json","{\"snapshot_folder\":"+Quoted(folder)+",\"numerical_round_passed\":false}");
    RevisitPublishedPass(e,state);Require(state.counts==original,"failed published round selected for early acceptance");
}
void InterpolationEnvelope(const std::string& root){
    ArgPP args;RunConfig config;Engine e(args,config);SetupMemory(e,root,"/bin/true",false);
    e.theta={170,172.5,175,177.5,180};e.maxCount=65536;e.orientationBlock=256;e.depthSharing=false;e.local=true;
    const std::vector<size_t> ids={0,1,2,3,4},knots={0,2,4};Budgets counts(1,std::vector<int>(5,256));counts[0][0]=8192;
    const auto requests=AuditRequests(e,{6},counts,ids,knots);bool extended=false;
    for(const auto& r:requests)if(r.count==8192&&std::find(r.angles.begin(),r.angles.end(),172.5)!=r.angles.end())extended=true;
    Require(extended,"interpolation common-prefix work omitted");
}
void CompleteResidualCohorts(const std::string& root){
    ArgPP args;RunConfig config;Engine e(args,config);SetupMemory(e,root,"/bin/true",false);
    e.theta.clear();std::vector<size_t> ids;for(size_t i=0;i<32;++i){e.theta.push_back(170+i*.1);ids.push_back(i);}
    e.maxCount=65536;e.orientationBlock=256;e.depthSharing=false;e.local=true;
    Job j;e.SetSampling(j,256,9999,6,e.theta,false);j.identity="timing";j.folder=root+"/timing";j.seconds=100;j.traceSeconds=.01;j.diffractionSeconds=99.99;e.CacheJob(j);
    Budgets current(1,std::vector<int>(32,256)),proposed(1,std::vector<int>(32,8192));const auto wanted=proposed;
    Statistics stats;stats.interval.resize(32);stats.ratioInterval.resize(32);
    for(size_t r=0;r<32;++r)stats.interval[r][0]=.1+r*.01;
    std::vector<std::vector<size_t>> cohorts;for(size_t id=0;id<32;id+=4)cohorts.push_back({id,id+1,id+2,id+3});
    AdditionalWork work(e,e.seeds.data());BoundRoundWork(e,work,current,proposed,{6},ids,stats,ids,cohorts);
    Require(proposed!=wanted&&proposed!=current,"residual cohort scheduling made no bounded progress");
    for(const auto& group:cohorts)for(size_t id:group){
        Require(proposed[0][id]==proposed[0][group.front()],"signed residual cohort split");
        Require(proposed[0][id]==256||proposed[0][id]==8192,"cohort common-prefix prescription truncated");}
}
int main(int argc,char** argv){try{
    Require(argc==2,"OUTPUT");StopPolicy();SavedPass(std::string(argv[1])+"/saved");InterpolationEnvelope(std::string(argv[1])+"/envelope");
    CompleteResidualCohorts(std::string(argv[1])+"/cohorts");
    PartitionWork(std::string(argv[1])+"/diffraction",false);PartitionWork(std::string(argv[1])+"/trace",true);
    std::cout<<"PASS: stop policy priority, physical GPU capacity, saved allocation without acceptance, interpolation envelope, bounded work and trace sharing\n";return 0;
}catch(const std::exception& error){std::cerr<<error.what()<<'\n';return 1;}}
