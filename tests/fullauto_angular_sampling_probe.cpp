// Real eight-replica replay plus decision regressions. Synthetic fixtures below
// isolate controller behavior; they are not physical accuracy/speed claims.
#include "../src/FullAuto.cpp"
using namespace FullAuto;
namespace {
const double signs[8]={-1,-1,-1,-1,1,1,1,1};
void Close(double a,double b,const char* message){
    Require(std::isfinite(a)&&std::fabs(a-b)<1e-10*std::max(1.,std::fabs(b)),message);
}
void Init(Engine& e,AdaptiveState& state,const std::string& root,size_t levels=2){
    e.root=root;e.signature="angular_sampling_probe";e.controls=false;e.coherent=false;
    e.local=true;e.localPolicy="on";e.depthSharing=false;e.ratios=true;
    e.eps=.04;e.statistical=.03;e.interpolation=.01;
    e.initial=512;e.minimum=64;e.maxCount=1<<22;e.phi=32;
    e.minPhi=1;e.maxPhi=64;e.fixedPhi=32;e.seeds.assign(production,production+8);
    e.theta={168.45,168.466,168.479};state.knots={0,1,2};
    state.counts=Budgets(levels,std::vector<int>(3,512));
    state.previous.resize(3);state.valid.assign(3,true);state.coarse.assign(3,true);
    Directory(root);Directory(root+"/jobs");
}
AdaptiveData Quiet(size_t levels=2){
    AdaptiveData data;data.levels.assign(levels,Samples(8,Matrix(3)));
    data.costs.assign(levels,std::vector<double>(3,.01));
    for(int s=0;s<8;++s)for(size_t r=0;r<3;++r){
        data.levels[0][s][r][0]=100;
        for(int j=1;j<16;++j)data.levels[0][s][r][j]=j*.1;
    }
    return data;
}
AngularSamplingDecision Check(const Engine& e,const AdaptiveState& state,const AdaptiveData& data){
    const auto samples=Combine(data.levels);
    return CheckAngularSampling(e,state,data,{0,1,2},samples,Estimate(samples));
}
void RealReplay(const std::string& fixture,const std::string& root){
    std::ifstream in(fixture);size_t rows,levels,replicas;
    Require(bool(in>>rows>>levels>>replicas)&&rows==3&&levels==3&&replicas==8,"bad real fixture header");
    ArgPP args;RunConfig config;Engine e(args,config);AdaptiveState state;Init(e,state,root,levels);
    for(double& angle:e.theta)Require(bool(in>>angle),"missing real angle");
    for(auto& level:state.counts)for(int& count:level)Require(bool(in>>count),"missing real prefix");
    auto data=Quiet(levels);
    for(auto& level:data.levels)for(auto& sample:level)for(auto& row:sample)for(double& value:row)
        Require(bool(in>>value)&&std::isfinite(value),"missing real Mueller element");
    std::string extra;Require(!(in>>extra),"trailing fixture data");
    const auto original=data.levels;const auto counts=state.counts;
    const auto samples=Combine(data.levels);const auto stats=Estimate(samples);
    Close(stats.mean[0][0],1385.5082116425174,"left production M11 changed");
    Close(stats.mean[1][0],1243.038413748135,"middle production M11 changed");
    Close(stats.mean[2][0],1448.5684709154918,"right production M11 changed");
    const auto result=Check(e,state,data);
    Require(result.enabled&&!result.passed&&result.failed==2&&result.priorities.size()==result.failed+result.failedShapes,"real jagged supports were accepted");
    Close(result.maximum,.07982800520590279,"raw paired noise changed");
    // v3 reserves a quarter of the component tolerance for angular noise.
    Close(result.priorities[0].error,4*.2551012778371969,"left component priority changed");
    Close(result.priorities[1].error,4*.20594001091880645,"right component priority changed");
    Require(data.levels==original&&state.counts==counts,"checking altered data or current allocation");
    bool grew=false;
    for(size_t l=0;l<levels;++l)for(size_t r=0;r<rows;++r){
        Require(result.proposed[l][r]>=counts[l][r]&&result.proposed[l][r]<=4*counts[l][r],"unbounded or shrinking real proposal");
        grew=grew||result.proposed[l][r]>counts[l][r];
    }
    Require(grew,"real pair did not schedule refinement");
    SaveAngularSampling(e,result,root);
    const auto saved=SavedAngularPriorities(e,state,root);
    Require(saved.size()==result.priorities.size()&&saved[0].nodes==std::vector<size_t>({0,1})&&saved[1].nodes==std::vector<size_t>({1,2}),"resume split/lost measured stencils");
    // Unknown-cost round limiting must retain one entire selected pair.
    e.orientationBlock=4096;AdditionalWork work(e,e.seeds.data());
    auto proposal=result.proposed;
    BoundRoundWork(e,work,counts,proposal,{8,10,12},{0,1,2},stats,{}, {},nullptr,saved);
    bool intact=false;
    for(const auto& group:saved){bool whole=true,changed=false;
        for(size_t r:group.nodes)for(size_t l=0;l<levels;++l){whole=whole&&proposal[l][r]==result.proposed[l][r];changed=changed||proposal[l][r]>counts[l][r];}
        intact=intact||(whole&&changed);
    }
    Require(intact,"time limiter destroyed all paired growth");
    std::cout<<"PASS measured 562.3: two noisy pairs rejected; saved matrices unchanged; local bounded growth\n";
}
void CompletedReplays(const std::string& fixture,const std::string& root){
    std::ifstream in(fixture);size_t cases;Require(bool(in>>cases)&&cases==7,"bad completed fixture header");
    for(size_t k=0;k<cases;++k){
        std::string name;Require(bool(in>>name),"missing completed job name");
        ArgPP args;RunConfig config;Engine e(args,config);AdaptiveState state;Init(e,state,root+"/"+Number(k),3);
        for(double& angle:e.theta)Require(bool(in>>angle),"missing completed angle");
        for(auto& level:state.counts)for(int& n:level)Require(bool(in>>n),"missing completed prefix");
        auto data=Quiet(3);
        for(auto& level:data.levels)for(auto& sample:level)for(auto& row:sample)for(double& v:row)
            Require(bool(in>>v)&&std::isfinite(v),"missing completed Mueller element");
        const auto original=data.levels;const auto counts=state.counts;
        const auto result=Check(e,state,data);
        Require(result.maximum<e.interpolation,"fixture did not pass the old M11-scale angular rule");
        Require(!result.passed&&result.maximumRatio>1&&result.failed>0,"noisy weak completed element accepted");
        Require(result.allowance[5]<e.interpolation,"weak M22 kept the full M11 allowance");
        Require(data.levels==original&&state.counts==counts,"quality check changed measured matrices");
        bool grew=false;for(size_t l=0;l<counts.size();++l)for(size_t r=0;r<3;++r){
            Require(result.proposed[l][r]>=counts[l][r]&&result.proposed[l][r]<=4*counts[l][r],"unbounded weak-element refinement");
            grew=grew||result.proposed[l][r]>counts[l][r];
        }
        Require(grew,"failed completed stencil did not schedule refinement");
        SaveAngularSampling(e,result,e.root);
        const auto saved=SavedAngularPriorities(e,state,e.root);
        Require(saved.size()==result.priorities.size(),"weak-element priorities lost on resume");
        for(size_t i=0;i<saved.size();++i)Close(saved[i].error,result.priorities[i].error,"weak-element priority scale changed on resume");
        std::cout<<"PASS completed replay "<<name<<" old_max="<<result.maximum<<" allowance_ratio="<<result.maximumRatio<<'\n';
    }
    std::string extra;Require(!(in>>extra),"trailing completed fixture data");
}
void ResidualReplay(const std::string& fixture,const std::string& root){
    // These are the ACTUAL combined estimator replicas after v2 accepted
    // 112.2 um. Representing the combined estimator as one level tests only
    // acceptance; it does not claim to reproduce its multilevel allocation.
    std::ifstream in(fixture);size_t rows,knots,replicas;
    Require(bool(in>>rows>>knots>>replicas)&&rows==165&&knots==83&&replicas==8,"bad residual fixture dimensions");
    ArgPP args;RunConfig config;Engine e(args,config);AdaptiveState state;Init(e,state,root,1);
    e.theta.resize(rows);for(double& angle:e.theta)Require(bool(in>>angle),"missing residual angle");
    state.knots.resize(knots);for(size_t& k:state.knots)Require(bool(in>>k)&&k<rows,"bad residual knot");
    state.counts=Budgets(1,std::vector<int>(rows,512));
    AdaptiveData data;data.levels.assign(1,Samples(8,Matrix(rows)));
    for(auto& sample:data.levels[0])for(auto& row:sample)for(double& v:row)
        Require(bool(in>>v)&&std::isfinite(v),"missing residual replica");
    std::string extra;Require(!(in>>extra),"trailing residual fixture");
    const auto original=data.levels;const auto counts=state.counts;
    std::vector<size_t> ids;for(size_t i=0;i<rows;++i)ids.push_back(i);
    const auto samples=Combine(data.levels);const auto stats=Estimate(samples);
    const auto result=CheckAngularSampling(e,state,data,ids,samples,stats);
    Require(!result.passed&&result.pairs==82&&result.failed==72,"v2 residual pairs incorrectly accepted");
    Require(result.shapes==160&&result.failedShapes>0,"residual shape noise not checked");
    Close(result.allowance[5],.0005793091420364889,"M22 residual allowance changed");
    Require(data.levels==original&&state.counts==counts,"residual replay modified physical values");
    SaveAngularSampling(e,result,root);
    const auto saved=SavedAngularPriorities(e,state,root);
    Require(saved.size()==result.priorities.size(),"residual scheduling priorities lost on resume");
    for(size_t i=0;i<saved.size();++i){
        Require(saved[i].nodes==result.priorities[i].nodes,"residual stencil split on resume");
        Close(saved[i].error,result.priorities[i].error,"residual priority changed on resume");
    }
    std::cout<<"PASS measured v2 residual: "<<result.failed<<"/"<<result.pairs<<" pairs and "
        <<result.failedShapes<<"/"<<result.shapes<<" shapes rejected; saved values unchanged\n";
}
}
int main(int argc,char** argv){try{
    Require(argc==5,"expected old fixture, private output directory, completed and residual fixtures");
    const std::string root=argv[2];Directory(root);RealReplay(argv[1],root+"/real562");
    CompletedReplays(argv[3],root+"/completed");
    ResidualReplay(argv[4],root+"/residual112");
    ArgPP args;RunConfig config;size_t tests=0;
    for(int j=0;j<16;++j){
        Engine e(args,config);AdaptiveState state;Init(e,state,root+"/raw"+Number(j));e.ratios=false;
        auto data=Quiet();for(int s=0;s<8;++s)data.levels[1][s][1][j]=5*signs[s];
        const auto result=Check(e,state,data);
        Require(!result.passed&&result.failed==2&&result.proposed[0]==state.counts[0],"raw noise missed or quiet base grew");
        Require(result.proposed[1]==std::vector<int>({1024,1024,1024}),"equal-prefix pair did not grow together");++tests;
    }
    {
        Engine e(args,config);AdaptiveState state;Init(e,state,root+"/weak-zero");auto data=Quiet();
        for(int s=0;s<8;++s){for(int r=0;r<3;++r)data.levels[0][s][r][5]=2*(r-1);
            data.levels[1][s][1][5]=.2*signs[s];}
        const auto noisy=Check(e,state,data);
        Require(noisy.maximum<e.interpolation&&!noisy.passed,"weak zero-crossing noise was accepted");
        Close(noisy.allowance[5],.01*.04,"zero crossing lost the nonzero amplitude floor");
        for(int s=0;s<8;++s)data.levels[1][s][1][5]=0;
        Require(Check(e,state,data).passed,"resolved physical zero crossing failed");++tests;
    }
    {
        Engine e(args,config);AdaptiveState state;Init(e,state,root+"/unequal");auto data=Quiet();
        state.counts[1]={8192,512,8192};
        for(int s=0;s<8;++s)data.levels[1][s][1][5]=5*signs[s];
        const auto result=Check(e,state,data);
        Require(result.proposed[1]==std::vector<int>({8192,2048,8192})&&result.proposed[0]==state.counts[0],"matching propagated maximum or overshot fourfold cap");
        e.maxCount=512;state.counts[1].assign(3,512);
        const auto capped=Check(e,state,data);Require(!capped.passed&&capped.proposed==state.counts,"maxCount manufactured acceptance");++tests;
    }
    {
        Engine e(args,config);AdaptiveState state;Init(e,state,root+"/sharp");auto data=Quiet();
        for(int s=0;s<8;++s)for(int j=0;j<16;++j)data.levels[0][s][1][j]*=100;
        const auto result=Check(e,state,data);Require(result.passed&&result.maximum<1e-12&&result.proposed==state.counts,"deterministic sharp structure constrained");++tests;
    }
    {
        Engine e(args,config);AdaptiveState state;Init(e,state,root+"/common");auto data=Quiet();
        for(int s=0;s<8;++s)for(int r=0;r<3;++r)data.levels[1][s][r][5]=50*signs[s];
        Require(Check(e,state,data).passed,"perfectly common angular noise failed to cancel");++tests;
    }
    {
        Engine e(args,config);AdaptiveState state;Init(e,state,root+"/ratio");auto data=Quiet();
        for(int s=0;s<8;++s){data.levels[0][s][0][5]=90;data.levels[0][s][1][5]=-90;
            for(int r=0;r<3;++r)data.levels[1][s][r][0]=5*signs[s];}
        e.ratios=false;Require(Check(e,state,data).passed,"common M11 noise did not cancel");
        e.ratios=true;Require(!Check(e,state,data).passed,"ratio contrast missed common-denominator noise");++tests;
    }
    {
        Engine e(args,config);AdaptiveState state;Init(e,state,root+"/invalid");auto data=Quiet();
        for(int s=0;s<8;++s)data.levels[0][s][1][0]=1+10*signs[s];
        Require(!Check(e,state,data).passed,"unresolved Fieller denominator passed");
        auto samples=Combine(data.levels);auto stats=Estimate(samples);
        for(int s=0;s<8;++s)data.levels[0][s][1][0]=-1;
        samples=Combine(data.levels);stats.mean[1][0]=-1;
        Require(!CheckAngularSampling(e,state,data,{0,1,2},samples,stats).passed,"nonpositive intermediate intensity accepted");
        samples.pop_back();bool threw=false;
        try{CheckAngularSampling(e,state,data,{0,1,2},samples,stats);}catch(const std::exception&){threw=true;}
        Require(threw,"fewer than eight replicas accepted");++tests;
    }
    {
        Engine e(args,config);AdaptiveState state;Init(e,state,root+"/scope");auto data=Quiet();
        for(int s=0;s<8;++s)data.levels[1][s][1][5]=5*signs[s];
        e.coherent=true;Require(!Check(e,state,data).enabled,"coherent mode changed");
        e.coherent=false;e.go=true;Require(!Check(e,state,data).enabled,"GO mode changed");
        e.go=false;e.local=false;Require(!Check(e,state,data).enabled,"explicit global orientation mode changed");++tests;
    }
    {
        Engine e(args,config);AdaptiveState state;Init(e,state,root+"/span-two");
        e.theta={160.,160.2,160.4,160.6,160.8};state.knots={0,1,2,3,4};
        state.counts=Budgets(2,std::vector<int>(5,512));e.ratios=false;
        AdaptiveData data;data.levels.assign(2,Samples(8,Matrix(5)));
        // Every short increment and span-one bend is below the allowance;
        // the broad correlated bump exceeds it at span two.
        std::array<double,8> unit{};for(int s=0;s<8;++s)unit[s]=signs[s];
        const double amplitude=.75*(e.eps/4)*e.eps*100/AngularNoise(unit);
        const double bump[5]={0,1,2,1,0};
        for(int s=0;s<8;++s)for(size_t r=0;r<5;++r){
            data.levels[0][s][r][0]=100;data.levels[0][s][r][5]=1;
            data.levels[1][s][r][5]=amplitude*bump[r]*signs[s];
        }
        const std::vector<size_t> ids={0,1,2,3,4};
        const auto samples=Combine(data.levels);const auto stats=Estimate(samples);
        const auto result=CheckAngularSampling(e,state,data,ids,samples,stats);
        Require(result.failed==0&&result.shapes==4&&result.failedShapes==1&&!result.passed,"broad shape noise escaped pair-only gate");
        Require(result.priorities.size()==1&&result.priorities[0].nodes==std::vector<size_t>({0,2,4}),"span-two cohort missing");
        Require(result.proposed[0]==state.counts[0]&&result.proposed[1]==std::vector<int>({1024,512,1024,512,1024}),"shape refinement grew unrelated rows");
        SaveAngularSampling(e,result,e.root);const auto saved=SavedAngularPriorities(e,state,e.root);
        Require(saved.size()==1&&saved[0].nodes==result.priorities[0].nodes,"shape-only resume lost the cohort");
        e.orientationBlock=4096;AdditionalWork work(e,e.seeds.data());auto proposed=result.proposed;
        BoundRoundWork(e,work,state.counts,proposed,{4,6},ids,stats,{}, {},nullptr,saved);
        Require(proposed==result.proposed,"unknown-cost limiter split the shape-only cohort");
        e.maxCount=512;const auto capped=CheckAngularSampling(e,state,data,ids,samples,stats);
        Require(!capped.passed&&capped.proposed==state.counts,"shape maxCount manufactured acceptance");
        ++tests;
    }
    {
        Engine e(args,config);AdaptiveState state;Init(e,state,root+"/uneven-linear");auto data=Quiet();
        e.theta={160.,160.01,170.};
        // Random affine angular offsets have zero curvature even on an
        // uneven grid; a large deterministic peak must also remain intact.
        for(int s=0;s<8;++s)for(int r=0;r<3;++r){
            data.levels[0][s][r][5]=(r==1?50.:0.);
            data.levels[1][s][r][5]=signs[s]*(10+e.theta[r]);
        }
        const auto result=Check(e,state,data);
        Require(result.failedShapes==0&&result.maximumShapeRatio<1e-9,"affine noise/physical peak constrained as curvature");
        ++tests;
    }
    {
        Engine e(args,config);AdaptiveState state;Init(e,state,root+"/custom-budget");auto data=Quiet();
        e.interpolation=.001;Close(Check(e,state,data).allowance[5],.001*.04,"strict explicit interpolation budget ignored");
        e.interpolation=.1;Close(Check(e,state,data).allowance[5],.01*.04,"loose interpolation consumed full component budget");
        ++tests;
    }
    {
        Engine e(args,config);AdaptiveState state;Init(e,state,root+"/migration");state.streak=2;e.independentPassed=true;
        const auto counts=state.counts;const auto previous=state.previous;const auto valid=state.valid;
        Atomic(root+"/migration/angular_sampling_policy.json","{\"policy\":\"legacy\"}");
        PrepareAngularSamplingPolicy(e,state,true);
        Require(state.streak==0&&!e.independentPassed&&counts==state.counts&&previous==state.previous&&valid==state.valid,"migration discarded data or retained obsolete pass");
        AdaptiveState restored=state;restored.streak=9;Require(RestoreAdaptive(e,restored)&&restored.counts==counts&&restored.streak==0,"checkpoint migration not restorable");
        state.streak=1;e.independentPassed=true;PrepareAngularSamplingPolicy(e,state,true);
        Require(state.streak==1&&e.independentPassed,"same-policy resume discarded acceptance history");++tests;
    }
    std::cout<<"PASS "<<tests<<" decision regressions (all 16 elements, ratios, covariance, caps, modes, resume)\n";return 0;
}catch(const std::exception& error){std::cerr<<error.what()<<'\n';return 1;}}
