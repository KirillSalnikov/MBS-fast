// Regressions for the October 28 phi promotion and corrected-cutoff round
// scheduling failures. Metadata fixtures exercise allocation, not accuracy.
#define FULLAUTO_STALL_HELPERS_ONLY
#include "fullauto_stall_probe.cpp"

void Timing(Engine& e,int n,int depth,int pair,int role,double trace,double diffraction,double setup,const std::string& id){
    Job j;e.SetSampling(j,n,9999,depth,e.theta,false);j.pairedDepth=pair;
    if(role)j.executionHash=Hash(e.ExecutionHash()+(role==1?" balanced-child-v1":" off-child-v1"));
    j.identity=id;j.folder=e.root+"/"+id;j.traceSeconds=trace;j.diffractionSeconds=diffraction;
    j.seconds=trace+diffraction+setup;e.CacheJob(j);
}
void CorrectedProfiles(const std::string& root){
    ArgPP args;RunConfig config;Engine e(args,config);SetupMemory(e,root,"/bin/true",false);
    e.theta={170,175,180};e.orientationBlock=4096;e.balancedCorrection=true;e.depthSharing=true;
    Timing(e,4096,14,0,1,4,32,1,"balanced");Timing(e,512,14,0,2,2,12,2,"off");
    Job parent;e.SetSampling(parent,4096,9999,14,e.theta,false);parent.assembled=true;parent.seconds=1e9;parent.identity="parent";parent.folder=root+"/parent";e.CacheJob(parent);
    AdditionalWork cost(e,e.seeds.data());const auto forecast=cost.Forecast({{4096,0,14,0,e.theta}});
    // B4096=37, B512=5.5, O512=16; no fake direct n14 profile exists.
    Require(forecast.profiled&&forecast.jobs==24&&std::fabs(forecast.seconds-468)<1e-10,"corrected base timing omitted a child or counted assembled parent");
    const auto small=cost.Forecast({{256,0,14,0,e.theta}});
    Require(small.profiled&&small.jobs==8&&small.seconds==72,"small off block did not use compatible off-child timing");
    Timing(e,4096,16,14,1,5,48,1,"balanced_pair");Timing(e,512,16,14,2,3,18,2,"off_pair");
    AdditionalWork paired(e,e.seeds.data());const auto pf=paired.Forecast({{4096,0,16,14,e.theta}});
    Require(pf.profiled&&pf.jobs==24&&pf.seconds==677,"paired correction did not account for three native legs");
    Timing(e,4096,16,0,1,5,48,1,"legacy_balanced_lower");Timing(e,4096,18,0,1,6,60,1,"legacy_balanced_upper");
    Timing(e,512,18,16,2,4,24,2,"legacy_off_pair");AdditionalWork legacy(e,e.seeds.data());
    const auto five=legacy.Forecast({{4096,0,18,16,e.theta}});Require(five.profiled&&five.jobs==40,"legacy five-leg correction missing a balanced depth");
    e.phi=2;AdditionalWork wrongPhi(e,e.seeds.data());Require(!wrongPhi.Forecast({{4096,0,14,0,e.theta}}).profiled,"other phi timing treated as measured");e.phi=1;
    e.base.push_back("different_optics");AdditionalWork wrongOptics(e,e.seeds.data());Require(!wrongOptics.Forecast({{4096,0,14,0,e.theta}}).profiled,"other recipe timing accepted");e.base.pop_back();
    for(unsigned seed:e.seeds){Job j=parent;j.seed=seed;j.identity="complete_"+Number(seed);j.folder=root+"/"+j.identity;e.CacheJob(j);}
    AdditionalWork cached(e,e.seeds.data());auto saved=cached.Forecast({{4096,0,14,0,e.theta}});
    Require(saved.jobs==0&&saved.seconds==0&&saved.profiled,"complete corrected parents scheduled again");
    // A child is a timing observation, never a completed physical target.
    AdditionalWork children(e,e.seeds.data(),true);Require(children.Forecast({{4096,0,14,0,e.theta}}).jobs==24,"child hashes contaminated physical coverage");
    std::cout<<"PASS corrected profiles: base, paired, small off, recipe isolation, assembled/cached parents\n";
}
void UnknownCost(const std::string& root){
    ArgPP args;RunConfig config;Engine e(args,config);SetupMemory(e,root,"/bin/true",false);e.local=true;e.maxCount=1<<22;e.orientationBlock=4096;
    std::vector<size_t> ids;e.theta.clear();for(size_t i=0;i<32;++i){e.theta.push_back(170+i*.25);ids.push_back(i);}
    Budgets before(1,std::vector<int>(32,524288)),wanted(1,std::vector<int>(32,2097152));
    Statistics stats;stats.interval.resize(32);stats.ratioInterval.resize(32);for(size_t i=0;i<32;++i)stats.interval[i][0]=.1+i*.01;
    std::vector<std::vector<size_t>> cohorts;for(size_t i=0;i<32;i+=4)cohorts.push_back({i,i+1,i+2,i+3});
    AdditionalWork work(e,e.seeds.data());auto selected=wanted;BoundRoundWork(e,work,before,selected,{14},ids,stats,{},cohorts);
    for(size_t i=0;i<32;++i)Require(selected[0][i]==(i>=28?2097152:524288),"unknown timing did not select one intact highest-priority cohort");
    const auto log=Read(root+"/round_work_limit.json");Require(log.find("\"unknown_cost_fallback\":true")!=std::string::npos&&log.find("\"selected_worker_seconds\":null")!=std::string::npos,
        "unknown work reported as zero seconds or silently bypassed limiter");
    std::cout<<"PASS unknown cost: one complete cohort, eight seeds unchanged, missing time is null\n";
}
void PhiPromotion(const std::string& root){
    ArgPP args;RunConfig config;Engine e(args,config);SetupMemory(e,root,"/bin/true",false);
    e.local=true;e.initial=524288;e.minimum=256;e.maxCount=1<<22;e.orientationBlock=4096;e.target=18;e.phi=64;e.balancedCorrection=true;e.depthSharing=true;
    AdaptiveState state;state.round=27;e.theta.clear();for(size_t i=0;i<1507;++i){e.theta.push_back(170+10.*i/1506);state.knots.push_back(i);}
    for(int l=0;l<3;++l){Timing(e,4096,14+2*l,l?12+2*l:0,1,.1,120,1,"balanced"+Number(l));Timing(e,512,14+2*l,l?12+2*l:0,2,.1,15,1,"off"+Number(l));}
    state.counts=Budgets(3,std::vector<int>(1507,8192));state.counts[0].assign(1507,524288);
    state.valid.assign(1507,true);state.coarse=state.valid;state.previous=Matrix(1507);state.auditSerial=7;
    WorkingPhiPlan plan;plan.phase=2;plan.oldPhi=16;plan.newPhi=64;plan.count=4096;plan.counts=state.counts;plan.sourceHash="frozen_pilot";
    e.phi=16;SaveWorkingPhi(e,plan);e.forbidNativeWork=true;ApplyWorkingPhi(e,state,plan);
    Require(state.round==27&&PeakBudgets(state.counts)[0]<524288&&PeakBudgets(state.counts)[0]>=512,"late phi promotion bypassed first-production bound");
    AdditionalWork cost(e,e.seeds.data());const auto estimate=cost.Forecast(BudgetRequests(e,{14,16,18},state.counts,state.knots));
    Require(estimate.profiled&&estimate.seconds<=1800,"measured phi restart exceeds scheduling budget");
    Require(state.auditSerial==7&&!e.independentPassed&&std::find(state.valid.begin(),state.valid.end(),true)==state.valid.end(),"phi scheduling retained invalid accuracy");
    const auto selected=state.counts;plan.phase=2;SaveWorkingPhi(e,plan);RestoreWorkingPhi(e,state);Require(state.counts==selected,"interrupted promotion changed frozen bounded allocation");
    // An older phase3 checkpoint has the huge allocation but no new-phi result.
    state.counts=plan.counts;RestoreWorkingPhi(e,state);Require(state.counts==selected,"legacy committed phi promotion left unbounded on resume");
    state.valid.front()=true;state.counts=plan.counts;RestoreWorkingPhi(e,state);Require(state.counts==plan.counts,"published new-phi estimator reset as a cold start");
    std::cout<<"PASS late phi 16->64, 1507 rows: 524288 -> "<<PeakBudgets(selected)[0]<<" base/seed; "<<estimate.seconds<<" forecast worker seconds; interrupted/legacy resume\n";
}
void PrefixRecovery(const std::string& root){
    ArgPP args;RunConfig config;Engine e(args,config);SetupMemory(e,root,"/bin/true",true);
    e.theta={170,175,180};e.phi=64;e.orientationBlock=512;e.depthSharing=true;e.maxCount=65536;e.forbidNativeWork=true;
    for(int offset=0;offset<8192;offset+=512)for(size_t s=0;s<8;++s){
        std::vector<double> x={170,175};if(offset<4096&&!(offset==3584&&s==7))x.push_back(180);
        Job j;e.SetSampling(j,512,e.seeds[s],14,x,false);j.offset=offset;j.assembled=true;j.identity=Number(offset)+"_"+Number(s);j.folder=root+"/"+j.identity;e.CacheJob(j);
    }
    Budgets wanted(1,std::vector<int>(3,65536)),floor(1,std::vector<int>(3,512));
    AdditionalWork work(e,e.seeds.data(),false,false,false,0);RetainScheduledPrefixes(e,work,{14},wanted,{0,1,2},floor);
    Require(floor[0]==std::vector<int>({8192,8192,2048}),"resume lost heterogeneous complete leaf prefixes or accepted incomplete eighth seed");
    Require(work.CoverageIndexBytes()==0&&!Exists(root+"/final_outputs.fnv"),"metadata recovery read matrices/accepted accuracy");
    for(int offset=0;offset<1024;offset+=512)for(unsigned seed:e.seeds){Job j;e.SetSampling(j,512,seed,16,e.theta,false);
        j.pairedDepth=14;j.offset=offset;j.identity="paired_"+Number(offset)+"_"+Number(seed);j.folder=root+"/"+j.identity;e.CacheJob(j);}
    wanted.push_back(std::vector<int>(3,65536));floor.push_back(std::vector<int>(3,512));AdditionalWork paired(e,e.seeds.data());
    RetainScheduledPrefixes(e,paired,{14,16},wanted,{0,1,2},floor);Require(floor[1]==std::vector<int>(3,1024),"complete paired-depth correction leaves lost on recovery");
    std::cout<<"PASS recovery: unassembled eight-seed prefixes, heterogeneous angles, partial tail, zero-budget index fallback\n";
}
void MeasuredReplay(const std::string& root,const std::string& fixture){
    struct Row {int role,count,depth,pair,rows,phi;double seconds,trace,diffraction,processing;};
    std::map<std::string,std::vector<Row>> data;std::istringstream input(Read(fixture));std::string line;
    while(std::getline(input,line)){if(line.empty()||line[0]=='#')continue;std::istringstream stream(line);Row row;std::string host;
        Require(bool(stream>>host>>row.role>>row.count>>row.depth>>row.pair>>row.rows>>row.phi>>row.seconds>>row.trace>>row.diffraction>>row.processing),"invalid measured timing fixture");data[host].push_back(row);}
    for(const auto& host:data){ArgPP args;RunConfig config;Engine e(args,config);SetupMemory(e,root+"/"+host.first,"/bin/true",host.first=="gpu1");
        e.balancedCorrection=true;e.phi=host.second.front().phi;e.kernel.pipeline=true;e.orientationBlock=4096;e.local=true;e.maxCount=1<<24;e.target=18;
        e.theta.clear();std::vector<size_t> ids;for(int i=0;i<host.second.front().rows;++i){e.theta.push_back(170+10.*i/(host.second.front().rows-1));ids.push_back(i);}
        double measured=0;int serial=0;
        for(const auto& row:host.second){Job j;e.SetSampling(j,row.count,9999,row.depth,e.theta,false);j.pairedDepth=row.pair;
            j.executionHash=Hash(e.ExecutionHash()+(row.role==1?" balanced-child-v1":" off-child-v1"));j.identity="measured_"+Number(serial++);j.folder=e.root+"/"+j.identity;
            j.seconds=row.seconds;j.traceSeconds=row.trace;j.diffractionSeconds=row.diffraction;j.processingSeconds=row.processing;measured+=j.seconds;e.CacheJob(j);}
        AdditionalWork work(e,e.seeds.data());const auto leaf=work.Forecast({{4096,0,14,0,e.theta}});
        Require(leaf.profiled&&std::fabs(leaf.seconds/8-measured)<.05*measured,"three-child model failed measured block replay");
        const int requested=host.first=="gpu1"?524288:2097152;
        const auto original=work.Forecast({{requested,0,14,0,e.theta}});
        if(host.first=="gpu1"){
            e.initial=524288;e.minimum=256;e.depthSharing=true;
            AdaptiveState state;state.round=27;state.knots=ids;state.counts=Budgets(3,std::vector<int>(e.theta.size(),8192));state.counts[0].assign(e.theta.size(),requested);
            // Missing correction-depth profiles force a small metadata-only
            // first wave. Only the base cost is measured in this fixture.
            BoundPhiProduction(e,state);Require(PeakBudgets(state.counts)[0]==512,"unknown correction timings bypassed phi restart pilot");
            auto base=work.Forecast(BudgetRequests(e,{14},Budgets{state.counts[0]},ids));
            std::cout<<host.first<<" measured base replay: "<<original.seconds<<" -> "<<base.seconds<<" worker seconds; corrections unprofiled, conservative pilot selected\n";
        }else{
            // Replay the n14 growth, with an already completed common prefix.
            // These are coverage metadata; no synthetic matrices enter tests.
            const int previous=requested/2;
            for(int offset=0;offset<previous;offset+=4096)for(unsigned seed:e.seeds){Job j;e.SetSampling(j,4096,seed,14,e.theta,false);
                j.offset=offset;j.assembled=true;j.identity=Number(offset)+"_"+Number(seed);j.folder=e.root+"/"+j.identity;e.CacheJob(j);}
            AdditionalWork warmed(e,e.seeds.data());Budgets current(1,std::vector<int>(e.theta.size(),previous)),next(1,std::vector<int>(e.theta.size(),requested));
            Statistics stats;stats.interval.resize(ids.size());stats.ratioInterval.resize(ids.size());for(size_t i:ids)stats.interval[i][0]=.22;
            const auto before=warmed.Forecast(BudgetRequests(e,{14},next,ids));BoundRoundWork(e,warmed,current,next,{14},ids,stats);
            const auto after=warmed.Forecast(BudgetRequests(e,{14},next,ids));size_t rows=0;for(size_t i:ids)rows+=next[0][i]>previous;
            Require(before.profiled&&after.profiled&&rows>0&&rows<ids.size()&&after.seconds<before.seconds/10,"measured corrected base still admitted whole multi-day round");
            std::cout<<host.first<<" measured n14 growth replay: "<<before.seconds<<" -> "<<after.seconds<<" worker seconds, "<<rows<<"/"<<ids.size()<<" growing rows; eight-seed prefixes preserved\n";
        }
    }
}
int main(int argc,char** argv){try{
    Require(argc==3,"OUTPUT MEASURED_TIMINGS");const std::string root=argv[1];CorrectedProfiles(root+"/profiles");UnknownCost(root+"/unknown");PhiPromotion(root+"/phi");PrefixRecovery(root+"/resume");MeasuredReplay(root+"/measured",argv[2]);
    return 0;
}catch(const std::exception& error){std::cerr<<error.what()<<'\n';return 1;}}
