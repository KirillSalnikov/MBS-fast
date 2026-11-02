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
    Require(result.enabled&&!result.passed&&result.failed==2&&result.priorities.size()==2,"real jagged supports were accepted");
    Close(result.priorities[0].error,.07789106231351053,"left paired noise changed");
    Close(result.priorities[1].error,.07982800520590279,"right paired noise changed");
    Require(data.levels==original&&state.counts==counts,"checking altered data or current allocation");
    bool grew=false;
    for(size_t l=0;l<levels;++l)for(size_t r=0;r<rows;++r){
        Require(result.proposed[l][r]>=counts[l][r]&&result.proposed[l][r]<=4*counts[l][r],"unbounded or shrinking real proposal");
        grew=grew||result.proposed[l][r]>counts[l][r];
    }
    Require(grew,"real pair did not schedule refinement");
    SaveAngularSampling(e,result,root);
    const auto saved=SavedAngularPriorities(e,state,root);
    Require(saved.size()==2&&saved[0].nodes==std::vector<size_t>({0,1})&&saved[1].nodes==std::vector<size_t>({1,2}),"resume split/lost measured pairs");
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
}
int main(int argc,char** argv){try{
    Require(argc==3,"expected fixture and private output directory");
    const std::string root=argv[2];Directory(root);RealReplay(argv[1],root+"/real562");
    ArgPP args;RunConfig config;size_t tests=0;
    for(int j=0;j<16;++j){
        Engine e(args,config);AdaptiveState state;Init(e,state,root+"/raw"+Number(j));e.ratios=false;
        auto data=Quiet();for(int s=0;s<8;++s)data.levels[1][s][1][j]=5*signs[s];
        const auto result=Check(e,state,data);
        Require(!result.passed&&result.failed==2&&result.proposed[0]==state.counts[0],"raw noise missed or quiet base grew");
        Require(result.proposed[1]==std::vector<int>({1024,1024,1024}),"equal-prefix pair did not grow together");++tests;
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
