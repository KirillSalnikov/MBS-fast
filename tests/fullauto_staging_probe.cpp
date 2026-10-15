// Numerical scheduling tests plus real native production subround recovery.
#include "../src/FullAuto.cpp"
using namespace FullAuto;
void Close(double a,double b,const char* message){Require(std::fabs(a-b)<1e-10*std::max(1.,std::fabs(b)),message);}
void PhiDecisionBranches(const std::string& root,bool rejectConfirmation){
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
    const auto stats=Estimate(Combine(data.levels));state.previous=stats.mean;SaveAdaptive(e,state);
    for(bool confirm:{false,true})for(int phi:{1,2,4})for(int seed=0;seed<8;++seed)for(int depth:depths){
        e.phi=phi;const unsigned scramble=(confirm?1794000000u:1694000000u)+17*seed;
        const double scale=phi==1?1.:(phi==2?.5:(confirm&&rejectConfirmation?2.:.01));
        Job job;e.SetSampling(job,2048,scramble,depth,e.theta,false);job.weightHash=Hash("");job.angles=e.theta;job.matrix=Matrix(9);
        for(int row=0;row<9;++row)job.matrix[row][0]=100+scale*(20*a[seed]+(depth>=2?5*b[seed]:0)+(depth>=4?2*c[seed]:0));
        job.identity=Number(phi)+"_"+Number(scramble)+"_"+Number(depth);job.folder=root+"/fake_"+job.identity;job.seconds=2048e-5*(1+.1*phi);e.CacheJob(job);
    }
    e.phi=1;const auto counts=state.counts;
    const bool promoted=TryWorkingPhi(e,state,{1,2,4},data,state.knots,stats,0);
    Require(promoted!=rejectConfirmation,"independent phi confirmation branch selected incorrectly");
    const auto report=Read(root+"/working_phi_n4.json");Require(report.find("\"independent_confirmation\":true")!=std::string::npos,"phi decision skipped independent confirmation");
    Require(state.auditSerial==42,"phi decision rewound audit seed serial");
    if(rejectConfirmation)Require(e.phi==1&&state.counts==counts&&state.valid[0],"rejected confirmation changed production");
    else Require(e.phi==4&&!state.valid[0]&&ReadWorkingPhi(e).phase==3,"accepted confirmation did not commit a separate estimator");
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
    if(argc==3&&std::string(argv[1])=="--phi-decisions"){PhiDecisionBranches(std::string(argv[2])+"/accepted",false);PhiDecisionBranches(std::string(argv[2])+"/rejected",true);std::cout<<"PASS: independently confirmed and rejected phi decisions\n";return 0;}
    ArgPP args;RegisterCliOptions(args);args.Parse(argc,argv);RunConfig config=RunConfig::FromCommandLine(args,false,false);Engine e(args,config);e.Init(argc,argv);
    const char* worker=getenv("TASK_MBS_WORKER");Require(worker,"TASK_MBS_WORKER required");e.exe=worker;e.binaryHash=FileHash(e.exe);e.phi=args.GetIntValue("nphi");e.controls=false;e.SetMirror(false);
    e.workerPolicy="persistent";e.local=true;e.depthSharing=true;e.depthVerified=true;Directory(e.root+"/jobs");WriteGrid(e.root+"/requested_theta.csv",e.theta);
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
    e.fixedPhi=fixed;
    std::cout<<"PASS: real correction-only subround; unchanged base; durable receipts; fresh native phi probes; explicit phi; recoverable promotion\n";return 0;
}catch(const std::exception& error){std::cerr<<error.what()<<std::endl;return 1;}}
