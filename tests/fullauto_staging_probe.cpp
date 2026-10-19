// Numerical scheduling tests plus real native production subround recovery.
#include "../src/FullAuto.cpp"
using namespace FullAuto;
void Close(double a,double b,const char* message){Require(std::fabs(a-b)<1e-10*std::max(1.,std::fabs(b)),message);}
void PhiDecisionBranches(const std::string& root,bool rejectConfirmation,bool resolveOnRetry=false,bool resumeInconclusive=false,bool retainFirst=false){
    // Fault-injected statistical decisions, not optical measurements. Real
    // native phi matrices and production accuracy are checked separately.
    ArgPP args;RunConfig config;Engine e(args,config);e.root=root;e.signature="phi_decision_fixture";
    Directory(root);e.exe="/bin/true";e.binaryHash=FileHash(e.exe);e.controls=false;e.depthSharing=false;
    e.phi=1;e.minPhi=1;e.maxPhi=4;e.target=4;e.initial=512;e.minimum=64;e.maxCount=1<<20;e.statistical=.001;
    e.gpus={""};e.seeds.assign(production,production+8);for(int r=0;r<9;++r)e.theta.push_back(170+1.25*r);
    AdaptiveState state;state.round=2;state.auditSerial=42;state.knots={0,1,2,3,4,5,6,7,8};
    state.counts=Budgets({std::vector<int>(9,4096),std::vector<int>(9,64),std::vector<int>(9,64)});
    state.valid.assign(9,true);state.coarse.assign(9,true);
    const int depths[3]={1,2,4};const double a[8]={-1,-1,-1,-1,1,1,1,1},b[8]={-1,-1,1,1,-1,-1,1,1},c[8]={-1,1,-1,1,-1,1,-1,1};
    AdaptiveData data;data.sources="fault_injected_decision_data";data.levels.assign(3,Samples(8,Matrix(9)));
    for(int seed=0;seed<8;++seed)for(int row=0;row<9;++row){data.levels[0][seed][row][0]=100+20*a[seed];data.levels[1][seed][row][0]=5*b[seed];data.levels[2][seed][row][0]=2*c[seed];}
    if(retainFirst)for(int seed=0;seed<8;++seed){data.levels[0][seed][0][0]=100;data.levels[1][seed][0][0]=data.levels[2][seed][0][0]=0;}
    const auto stats=Estimate(Combine(data.levels));state.previous=stats.mean;SaveAdaptive(e,state);
    for(int attempt=0;attempt<3;++attempt)for(bool confirm:{false,true})for(int phi:{1,2,4})for(int seed=0;seed<8;++seed)for(int depth:depths){
        const int count=attempt==0?2048:(attempt==1?8192:16384);
        e.phi=phi;const unsigned scramble=(confirm?1794000000u:1694000000u)+100000*attempt+17*seed;
        const bool reject=rejectConfirmation&&(!resolveOnRetry||attempt==0);
        const double scale=phi==1?1.:(phi==2?.5:(confirm&&reject?2.:.01));
        Job job;e.SetSampling(job,count,scramble,depth,e.theta,false);job.weightHash=Hash("");job.angles=e.theta;job.matrix=Matrix(9);
        for(int row=0;row<9;++row)job.matrix[row][0]=100+scale*(20*a[seed]+(depth>=2?5*b[seed]:0)+(depth>=4?2*c[seed]:0));
        job.identity=Number(phi)+"_"+Number(scramble)+"_"+Number(depth);job.folder=root+"/fake_"+job.identity;job.seconds=count*1e-5*(1+.1*phi);e.CacheJob(job);
    }
    e.phi=1;const auto counts=state.counts;
    if(resumeInconclusive){WorkingPhiPlan pending;pending.phase=4;pending.oldPhi=1;pending.newPhi=1;pending.count=2048;
        pending.peak=4096;pending.rows=state.knots;pending.sourceHash=Hash(data.sources);SaveWorkingPhi(e,pending);
        Require(ReadWorkingPhi(e).phase==4,"inconclusive receipt not recoverable");
        const auto path=WorkingPhiPath(e),saved=Read(path);Atomic(path,saved+" ");bool rejected=false;
        try{ReadWorkingPhi(e);}catch(const std::exception&){rejected=true;}Require(rejected,"corrupted retry receipt accepted");Atomic(path,saved);}
    InterpolationResult interpolation;interpolation.errors.resize(9);interpolation.matrixErrors.resize(9);interpolation.ratioErrors.resize(9);
    const bool promoted=TryWorkingPhi(e,state,{1,2,4},data,state.knots,stats,0,retainFirst?&interpolation:nullptr);
    Require(promoted==(!rejectConfirmation||resolveOnRetry),"independent phi confirmation branch selected incorrectly");
    const auto report=Read(root+"/working_phi_n4.json");Require(report.find("\"independent_confirmation\":true")!=std::string::npos,"phi decision skipped independent confirmation");
    Require(state.auditSerial==42,"phi decision rewound audit seed serial");
    if(rejectConfirmation&&!resolveOnRetry)Require(e.phi==1&&state.counts==counts&&state.valid[0]&&ReadWorkingPhi(e).attempt==2,"rejected confirmation changed production or retries unbounded");
    else Require(e.phi==4&&state.valid[0]==retainFirst&&!state.valid[1]&&ReadWorkingPhi(e).phase==3,"accepted confirmation did not commit the selected rows");
    if(retainFirst){Require(ReadWorkingPhi(e).retained==std::vector<size_t>({0})&&e.PhiForTheta(0)==1,"converged theta not retained by automatic rescue");
        for(size_t l=0;l<counts.size();++l)Require(state.counts[l][0]==counts[l][0],"retained theta lost its completed counts");}
    if(resolveOnRetry)Require(ReadWorkingPhi(e).attempt==1&&state.round==2,"contradictory pilot waited for production growth");
}
int Decisions(const char* fixture,bool strict=true){
    const double a[8]={-1,-1,-1,-1,1,1,1,1},b[8]={-1,-1,1,1,-1,-1,1,1};
    std::vector<Samples> levels(3,Samples(8,Matrix(1)));
    for(int s=0;s<8;++s){levels[0][s][0][0]=100+20*a[s];levels[1][s][0][0]=-20*a[s];levels[2][s][0][0]=5*b[s];}
    const auto mean=Mean(Combine(levels));AllocationForecast forecast(levels,{512,512,64},mean,true);
    const auto actual=Estimate(Combine(levels)),predicted=forecast.Predict({512,512,64});
    for(int j=0;j<16;++j){Close(actual.interval[0][j],predicted.interval[0][j],"covariance forecast differs at current N");Close(actual.ratioInterval[0][j],predicted.ratioInterval[0][j],"ratio forecast differs at current N");}
    const auto ns=AllocateCovariant(levels,{1.,1.,1.},{512,512,64},mean,.02,true,1<<20);
    Require(ns[0]==512&&ns[1]==512&&ns[2]>64,"allocator destroyed a cancelling pair");
    levels[2][0][0][1]=20;levels[2][1][0][1]=-20;
    AllocationForecast ratios(levels,{512,512,64},Mean(Combine(levels)),true);
    const auto rp=ratios.Predict({1024,1024,128});
    Samples reference(8,Mean(Combine(levels)));
    for(size_t l=0;l<levels.size();++l){const auto m=Mean(levels[l]);for(int s=0;s<8;++s)for(int j=0;j<16;++j)reference[s][0][j]+=(levels[l][s][0][j]-m[0][j])/std::sqrt(2.);}
    const auto expected=Estimate(reference);for(int j=0;j<16;++j){Close(rp.interval[0][j],expected.interval[0][j],"joint covariance lost");Close(rp.ratioInterval[0][j],expected.ratioInterval[0][j],"Fieller denominator covariance lost");}
    WorkingPhiPilot oldPhi,newPhi;oldPhi.error=.2;newPhi.error=.05;
    Require(WorkingPhiImproves(oldPhi,newPhi,1000,100,50),"useful independent phi confirmation rejected");
    Require(!WorkingPhiImproves(oldPhi,newPhi,1000,600,50),"phi ignored cost of replacing existing samples");
    oldPhi.error=1e100;Require(!WorkingPhiImproves(oldPhi,newPhi,1e200,100,50),"unresolved baseline denominator accepted");
    std::ifstream input(fixture);size_t rows;double tolerance;Require(bool(input>>rows>>tolerance),"invalid saved fixture");
    std::cout<<std::setprecision(17)<<"theta,current0,current1,current2,independent0,independent1,independent2,covariant0,covariant1,covariant2\n";
    for(size_t row=0;row<rows;++row){double theta;input>>theta;std::vector<int> counts(3);std::vector<double> costs(3);std::vector<Samples> real(3,Samples(8,Matrix(1)));
        for(int l=0;l<3;++l){input>>counts[l]>>costs[l];for(int s=0;s<8;++s)for(int j=0;j<16;++j)Require(bool(input>>real[l][s][0][j]),"incomplete saved level data");}
        const auto saved=real;const auto stats=Estimate(Combine(real));auto legacy=Allocate(real,costs,counts,stats.mean,tolerance*.85);
        auto proposed=std::max(stats.maxInterval,stats.maxRatioInterval)>tolerance?AllocateCovariant(real,costs,counts,stats.mean,tolerance*.85,true,1<<30,tolerance):counts;
        auto staged=PrioritizeAllocation(real,costs,counts,proposed,tolerance,true);Require(real==saved,"allocation changed observed matrices");
        if(strict)Require(staged[1]>counts[1],"real dominant correction not refined");
        std::cout<<theta;for(const auto& values:{counts,legacy,staged})for(int n:values)std::cout<<','<<n;std::cout<<'\n';
    }
    return 0;
}
int main(int argc,const char* argv[]){try{
    if(argc==3&&std::string(argv[1])=="--decisions")return Decisions(argv[2]);
    if(argc==3&&std::string(argv[1])=="--bulk-decisions")return Decisions(argv[2],false);
    if(argc==3&&std::string(argv[1])=="--phi-decisions"){
        {ArgPP args;RunConfig config;Engine e(args,config);e.initial=512;e.minimum=64;e.maxCount=8192;e.statistical=.01;
            WorkingPhiPilot high,low;high.error=.2;low.error=.05;high.cost=low.cost={1.};
            const Budgets counts={{8192}};const std::vector<size_t> ids={0};
            const double oldCost=WorkingPhiCost(e,high,counts,ids,true),newCost=WorkingPhiCost(e,low,counts,ids,false);
            Require(oldCost>10*newCost,"launch cap made an unconverged estimator look cheap");
            e.maxCount=1<<28;Close(WorkingPhiCost(e,high,counts,ids,true),oldCost,"forecast depends on launch cap instead of target error");}
        Close(WorkingPhiBaselineError(1e100,.2),.2,"unbounded short pilot hides measured production baseline");
        Close(WorkingPhiBaselineError(.4,.2),.4,"larger measured pilot error lost");
        Close(WorkingPhiBaselineError(.1,.2),.2,"production floor lost");
        Require(!FullAutoScheduling::Resolved(WorkingPhiBaselineError(1e100,1e100)),"unresolved baseline became finite");
        Require(WorkingPhiCandidates(1,1,64)==std::vector<int>({1,4,8,16,32}),"phi=32 absent from rescue");
        Require(WorkingPhiCandidates(16,1,32)==std::vector<int>({16,32}),"explicit phi upper bound ignored");
        Require(WorkingPhiDisagrees(3.4507133628298288,.091573670029267404),"real 398.1 contradictory pilot not detected");
        Require(!WorkingPhiDisagrees(.12,.09),"ordinary pilot fluctuation triggers retry");
        PhiDecisionBranches(std::string(argv[2])+"/accepted",false);PhiDecisionBranches(std::string(argv[2])+"/rejected",true);
        PhiDecisionBranches(std::string(argv[2])+"/retry",true,true);
        PhiDecisionBranches(std::string(argv[2])+"/resume_retry",true,true,true);
        PhiDecisionBranches(std::string(argv[2])+"/local_retention",false,false,false,true);
        std::cout<<"PASS: independently confirmed, bounded contradictory retries, and rejected phi decisions\n";return 0;}
    ArgPP args;RegisterCliOptions(args);args.Parse(argc,argv);RunConfig config=RunConfig::FromCommandLine(args,false,false);Engine e(args,config);e.Init(argc,argv);
    const char* worker=getenv("TASK_MBS_WORKER");Require(worker,"TASK_MBS_WORKER required");e.exe=worker;e.binaryHash=FileHash(e.exe);e.phi=args.GetIntValue("nphi");e.controls=false;e.SetMirror(false);
    e.workerPolicy="persistent";e.local=true;e.depthSharing=true;e.depthVerified=true;Directory(e.root+"/jobs");WriteGrid(e.root+"/requested_theta.csv",e.theta);
    if(getenv("TASK_TEST_LOCAL_FULLAUTO")){
        e.RecordFinalChecksPolicy();e.PrepareTuning(e.theta);
        std::ostringstream command;command<<"{\"test_driver\":\"fullauto_staging_probe\",\"argv\":[";for(int i=0;i<argc;++i)command<<(i?",":"")<<Quoted(argv[i]);command<<"]}";
        Atomic(e.root+"/controller_command.json",command.str());
    }
    AdaptiveState state;state.counts=Budgets({std::vector<int>(e.theta.size(),256),std::vector<int>(e.theta.size(),64),std::vector<int>(e.theta.size(),64)});
    state.valid.assign(e.theta.size(),true);state.coarse.assign(e.theta.size(),true);state.previous.resize(e.theta.size());state.auditSerial=42;
    std::vector<size_t> ids;for(size_t i=0;i<e.theta.size();++i)ids.push_back(i);state.knots=ids;const std::vector<int> depths={1,2,4};
    auto stage=EvaluateBudget(e,depths,state.counts,ids,e.seeds.data());auto stats=Estimate(Combine(stage.levels));state.previous=stats.mean;
    const auto originalBase=state.counts[0];std::map<std::string,std::string> markers;for(const auto& item:e.cache)markers[item.second.folder]=Read(item.second.folder+"/complete.txt");
    CorrectionProbeReceipt pending;pending.attempt=1;pending.sources=Hash(stage.sources);pending.rows={{ids[ids.size()/2],1,64,128}};SaveCorrectionProbe(e,state,pending);SaveAdaptive(e,state);
    Require(TryCorrectionProbe(e,state,depths,stage,ids,stats,state.counts,0),"pending native correction probe not resumed");
    Require(state.counts[0]==originalBase&&state.counts[1][ids[ids.size()/2]]==128,"correction probe extended the base or lost its sample");
    for(const auto& item:markers)Require(Read(item.first+"/complete.txt")==item.second,"correction probe replaced completed work");
    state.counts[1][ids[ids.size()/2]]=64;Require(RestoreCorrectionProbe(e,state),"committed correction not recovered after lost checkpoint");
    Require(!RestoreCorrectionProbe(e,state),"correction applied twice");
    const auto receiptPath=CorrectionProbePath(e,state),receipt=Read(receiptPath);Atomic(receiptPath,receipt+" ");bool corrupt=false;try{ReadCorrectionProbe(e,state);}catch(const std::exception&){corrupt=true;}Require(corrupt,"corrupt correction receipt accepted");Atomic(receiptPath,receipt);
    std::array<unsigned,8> screen;for(size_t i=0;i<8;++i)screen[i]=1690000000u+17*i;
    const int fixed=e.fixedPhi;Require(fixed>0&&!TryWorkingPhi(e,state,depths,stage,ids,stats,0),"manual phi was overridden");
    const auto reference=stats.mean;std::ostringstream pilot;pilot<<"phi,error,cost\n";
    for(int phi:{1,4}){e.phi=phi;const auto result=EvaluateWorkingPhi(e,depths,state.counts,ids,reference,64,screen.data());double cost=0;for(double v:result.cost)cost+=v;Require(std::isfinite(result.error)&&cost>0,"native phi probe invalid");pilot<<phi<<','<<Number(result.error)<<','<<Number(cost)<<'\n';}
    Atomic(e.root+"/native_phi_probe.csv",pilot.str());
    // Resume a small frozen screening plan through the real implementation.
    // The default large-particle trigger/count is tested separately by policy.
    e.phi=1;e.fixedPhi=0;e.minPhi=1;e.maxPhi=64;
    auto current=EvaluateBudget(e,depths,state.counts,ids,e.seeds.data());const auto currentStats=Estimate(Combine(current.levels));
    WorkingPhiPlan screening;screening.oldPhi=1;screening.count=64;screening.rows=ids;screening.sourceHash=Hash(current.sources);SaveWorkingPhi(e,screening);
    TryWorkingPhi(e,state,depths,current,ids,currentStats,0);Require(ReadWorkingPhi(e).phase!=0,"optional screening plan not closed");
    // Exercise the promotion transaction separately from the numerical choice.
    // This test does not claim that a small-prism probe chose phi=4.
    e.phi=1;e.fixedPhi=0;e.minPhi=1;e.maxPhi=64;WorkingPhiPlan promotion;promotion.phase=2;promotion.oldPhi=1;promotion.newPhi=4;promotion.count=64;promotion.sourceHash=Hash(stage.sources);promotion.rows=ids;promotion.counts=state.counts;
    SaveWorkingPhi(e,promotion);ApplyWorkingPhi(e,state,promotion);
    Require(e.phi==4&&state.auditSerial==42&&state.streak==0&&std::find(state.valid.begin(),state.valid.end(),true)==state.valid.end(),"phi promotion reused old convergence or audit serial");
    promotion.phase=2;SaveWorkingPhi(e,promotion);RestoreWorkingPhi(e,state);Require(e.phi==4&&ReadWorkingPhi(e).phase==3,"interrupted phi promotion not recovered");
    // Real optical rows at two rules, all depths and all eight seeds. Mixed
    // assembly must equal a literal per-row splice, without launching work
    // already cached by the two homogeneous estimators.
    e.phi=1;auto oldRows=EvaluateUniformBudget(e,depths,state.counts,ids,e.seeds.data());
    e.phi=4;auto newRows=EvaluateUniformBudget(e,depths,state.counts,ids,e.seeds.data());
    state.previous=Mean(Combine(oldRows.levels));state.valid.assign(ids.size(),true);state.coarse.assign(ids.size(),true);
    e.phi=1;promotion.phase=2;promotion.retained={ids.front()};promotion.counts=state.counts;
    const auto oldPrevious=state.previous;SaveWorkingPhi(e,promotion);ApplyWorkingPhi(e,state,promotion);
    Require(e.PhiForTheta(ids.front())==1&&e.PhiForTheta(ids.back())==4&&state.valid.front()&&!state.valid.back()&&state.previous.front()==oldPrevious.front(),"local promotion discarded retained rows");
    const size_t computed=e.newRows;auto mixed=EvaluateBudget(e,depths,state.counts,ids,e.seeds.data());
    Require(e.newRows==computed,"mixed assembly recalculated complete cached native rows");
    for(size_t l=0;l<depths.size();++l)for(size_t s=0;s<8;++s)for(size_t r=0;r<ids.size();++r)for(int j=0;j<16;++j)
        Close(mixed.levels[l][s][r][j],(r?newRows:oldRows).levels[l][s][r][j],"mixed phi changed an optical matrix element");
    const auto mixedStats=Estimate(Combine(mixed.levels));Require(std::isfinite(mixedStats.maxInterval),"mixed estimator lost confidence intervals");
    const auto mixedInterp=InterpolateBudget(e,state.knots,ids,Combine(mixed.levels));
    WriteFullAutoIntegrals(e,e.root+"/mixed_integrals",mixedInterp.output,mixed.sources,false);
    SaveAdaptive(e,state);Require(Read(e.root+"/adaptive_checkpoint.tsv").substr(0,2)=="4 ","mixed rule checkpoint not versioned");
    auto restored=state;e.thetaPhi.clear();Require(RestoreAdaptive(e,restored)&&e.PhiForTheta(ids.front())==1&&e.PhiForTheta(ids.back())==4,"mixed rule checkpoint not restored");
    Require(restored.counts==state.counts&&restored.previous==state.previous&&restored.valid==state.valid,"mixed resume changed counts or old means");
    auto again=EvaluateBudget(e,depths,state.counts,ids,e.seeds.data());Require(again.levels==mixed.levels,"mixed resumed matrices differ");
    {PhiScope scope(e,1);Require(e.phi==1,"phi scope failed");}Require(e.phi==4,"phi scope leaked to next worker");
    Require(ReadWorkingPhi(e).retained==promotion.retained,"local promotion receipt lost retained rows");
    {AdditionalWork work(e,e.seeds.data());const auto cost=work.Forecast(BudgetRequests(e,depths,state.counts,ids));Require(cost.jobs==0,"mixed forecast ignored saved phi-specific coverage");}
    if(ids.size()>=5){Budgets unequal=state.counts;unequal[0].assign(e.theta.size(),64);unequal[0][ids.front()]=256;
        const auto groups=FullAutoAuditPlan::Stencils({ids.front(),ids[2],ids.back()},unequal,ids);
        for(const auto& group:groups){auto ns=StencilBudgets(e,unequal,group);for(size_t id:group.second.rows)Require(ns[0][id]==(id==ids.front()?256:64),"old phi prefix inflated neighbouring new phi work");}
        Budgets extending=unequal;ExtendPhiBudgets(e,extending,{ids.front(),ids.back()},{ids.front(),ids[1],ids.back()});
        Require(extending[0][ids[1]]==std::max(e.initial,64),"new theta inherited another phi's count");}
    if(getenv("TASK_TEST_LOCAL_FULLAUTO")){
        Require(RunAdaptive(e)==0,"mixed phi full adaptive cycle did not converge");
        Require(e.thetaPhi.count(ids.front())&&e.PhiForTheta(ids.front())==1,"full adaptive cycle lost retained phi");
        e.completedCheckpoint=false;Require(RunAdaptive(e)==0&&e.completedCheckpoint,"completed mixed fullauto did not resume");
    }
    e.thetaPhi.clear();
    e.fixedPhi=fixed;
    std::cout<<"PASS: correction subround; local phi cache reuse all 16 elements; integral output; mixed checkpoint recovery; retained old rows; explicit phi\n";return 0;
}catch(const std::exception& error){std::cerr<<error.what()<<std::endl;return 1;}}
