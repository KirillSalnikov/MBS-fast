// Controller-only regressions. Optional measured fixtures contain the saved
// eight replicas at real support angles, not a synthetic optical reference.
#include "../src/FullAuto.cpp"
using namespace FullAuto;
namespace {
void Close(double a,double b,const char* message,double scale=1.){
    Require(std::isfinite(a)&&std::isfinite(b)&&std::fabs(a-b)<=2e-12*std::max(1.,std::fabs(scale)),message);
}
void NoExtraExtrema(const std::vector<double>& x,const Matrix& y,const std::vector<double>& q,const Matrix& out){
    for(size_t r=0;r<q.size();++r){const auto at=LocatePiecewiseLinear(x,q[r]);
        const auto& a=y[at.left];const auto& b=y[at.right];
        for(int j=0;j<16;++j){const double eps=2e-12*std::max(1.,std::max(std::fabs(a[j]),std::fabs(b[j])));
            Require(out[r][j]>=std::min(a[j],b[j])-eps&&out[r][j]<=std::max(a[j],b[j])+eps,"linear output overshot direct supports");
            if(a[0]>0&&b[0]>0){const double lo=std::min(a[j]/a[0],b[j]/b[0]),hi=std::max(a[j]/a[0],b[j]/b[0]);
                Require(out[r][0]>0&&out[r][j]/out[r][0]>=lo-2e-12&&out[r][j]/out[r][0]<=hi+2e-12,"normalized ratio overshot direct supports");}}
    }
}
void PureOperator(){
    const std::vector<double> x={0,.3,1,2},q={0,.1,.3,.7,1,1.25,1.75,2};
    Samples replicas(8,Matrix(x.size()));
    for(int s=0;s<8;++s)for(size_t r=0;r<x.size();++r){replicas[s][r][0]=100+r+2*s;
        for(int j=1;j<16;++j)replicas[s][r][j]=(r%2?-1.:1.)*(j+.1*r+.01*s*(r+1));}
    Samples output(8);
    for(int s=0;s<8;++s){
        output[s]=POOutputInterpolate(false,false,x,replicas[s],q);
        Require(POOutputInterpolate(false,false,x,replicas[s],x)==replicas[s],"support values changed");
        NoExtraExtrema(x,replicas[s],q,output[s]);
    }
    const auto mean=Mean(replicas),outputMean=Mean(output),meanOutput=POOutputInterpolate(false,false,x,mean,q);
    for(size_t r=0;r<q.size();++r){
        const auto weights=POOutputWeights(false,false,x,q[r]);Mueller predicted{};double sum=0;
        for(const auto& w:weights){Require(w.second>=0&&w.second<=1,"negative piecewise-linear weight");sum+=w.second;
            for(int j=0;j<16;++j)predicted[j]+=w.second*mean[w.first][j];}
        Close(sum,1,"linear weights do not sum to one");
        for(int j=0;j<16;++j){Close(predicted[j],meanOutput[r][j],"sparse weights disagree with output",meanOutput[r][0]);
            Close(outputMean[r][j],meanOutput[r][j],"interpolation did not commute with eight-replica mean",meanOutput[r][0]);}
    }
    Require(POOutputInterpolate(false,true,x,mean,q)==LocalCubicInterpolate(x,mean,q),"coherent estimator changed");
    Require(POOutputInterpolate(true,false,x,mean,q)==LocalCubicInterpolate(x,mean,q),"GO fallback estimator changed");
    Require(POOutputWeights(false,true,x,.7)==LocalCubicWeights(x,.7),"coherent covariance weights changed");
    Require(POOutputWeights(true,false,x,.7)==LocalCubicWeights(x,.7),"GO covariance weights changed");
    bool refused=false;try{PiecewiseLinearInterpolate(x,mean,{-1e-5});}catch(const std::runtime_error&){refused=true;}
    Require(refused,"linear extrapolation was allowed");
    Require(PiecewiseLinearInterpolate(x,mean,{-1e-13,2+1e-13})==Matrix({mean.front(),mean.back()}),"endpoint tolerance extrapolated");
}
void MidpointResidual(){
    ArgPP args;RunConfig config;Engine e(args,config);e.go=false;e.coherent=false;e.ratios=false;e.theta={0,.5,1};e.interpolation=.01;
    Samples samples(8,Matrix(3));
    for(int s=0;s<8;++s)for(int j=0;j<16;++j){
        const double a=(j?j:100)+.2*s,b=a+1.5+.1*j;
        samples[s][0][j]=a;samples[s][2][j]=b;
        samples[s][1][j]=(a+b)*.5+.005*(j+1)*(1.+(s<4?-1.:1.));
    }
    const auto result=InterpolateBudget(e,{0,2},{0,1,2},samples);
    for(int j=0;j<16;++j){double mean=0,variance=0,reference=0;
        for(int s=0;s<8;++s){mean+=((samples[s][0][j]+samples[s][2][j])*.5-samples[s][1][j])/8.;reference+=samples[s][1][0]/8.;}
        for(int s=0;s<8;++s)variance+=std::pow((samples[s][0][j]+samples[s][2][j])*.5-samples[s][1][j]-mean,2)/7.;
        Close(result.matrixErrors[1][j],(std::fabs(mean)+t95*std::sqrt(variance/8))/reference,"midpoint residual uses another operator");
        Close(result.matrixErrors[0][j],0,"support residual nonzero");Close(result.matrixErrors[2][j],0,"support residual nonzero");
    }
    // Independent angular control: high-count anchor interpolation plus the
    // paired low-count residual must use exactly this same linear operator.
    const std::vector<double> x={0,1},q={.5};
    for(int s=0;s<8;++s){Matrix anchors={samples[s][0],samples[s][2]},small={samples[7-s][0],samples[7-s][2]};
        const auto a=POOutputInterpolate(false,false,x,anchors,q),b=POOutputInterpolate(false,false,x,small,q);
        for(int j=0;j<16;++j)Close(a[0][j]+samples[7-s][1][j]-b[0][j],
            .5*(anchors[0][j]+anchors[1][j])+samples[7-s][1][j]-.5*(small[0][j]+small[1][j]),"angular control residual mismatch");}
}
void Migration(const std::string& root){
    ArgPP args;RunConfig config;Engine e(args,config);e.root=root;e.signature="linear-policy-regression";e.go=false;e.coherent=false;e.theta={0,.5,1};e.phi=8;
    Directory(root);AdaptiveState state;state.knots={0,2};state.counts=Budgets(1,std::vector<int>(3,512));state.valid.assign(3,true);state.coarse.assign(3,true);state.previous.resize(3);
    state.streak=3;e.independentPassed=true;
    Atomic(root+"/intermediate_current.json","{\"numerical_round_passed\":true}");Atomic(root+"/native_fixture.dat","saved-native-block");
    Require(!PublishedInterpolationCompatible(e),"legacy snapshot accepted as linear output");
    Require(PrepareOutputInterpolationPolicy(e,state,true),"legacy interpolation policy migration missed");
    Require(state.streak==0&&!e.independentPassed,"policy migration kept old acceptance");
    Require(Read(root+"/native_fixture.dat")=="saved-native-block","migration changed optical cache");
    Require(!PrepareOutputInterpolationPolicy(e,state,true),"same policy repeatedly migrates");
    Require(!PublishedInterpolationCompatible(e),"interrupted migration reused an old published decision");
    Atomic(root+"/intermediate_current.json","{\"interpolation_policy\":"+Quoted(POOutputInterpolationPolicy(false,false))+",\"numerical_round_passed\":true}");
    Require(PublishedInterpolationCompatible(e),"new published interpolation policy rejected");
    e.coherent=true;Require(PublishedInterpolationCompatible(e),"coherent legacy checkpoint behavior changed");
    // Same optical sources, different output operator: no reused audit and no
    // overwrite of a previously published snapshot with another estimator.
    AdaptiveData stage;stage.sources="same immutable optical sources\n";state.round=7;
    Statistics stats;stats.mean.resize(3);stats.interval.resize(3);stats.ratioInterval.resize(3);
    for(auto& row:stats.mean)row[0]=100;
    const std::vector<size_t> ids={0,1,2};const std::vector<double> errors(3,0);
    const auto legacy=root+"/round_snapshots/round7_"+Hash(stage.sources);Directory(legacy);Atomic(legacy+"/sentinel","legacy-preserved");
    const auto coherentFingerprint=AuditInputFingerprint(e,state,stage,ids);
    const auto coherentFolder=SaveRoundDiagnostics(e,state,stage,stats,ids,errors,e.theta);
    e.coherent=false;
    Require(AuditInputFingerprint(e,state,stage,ids)!=coherentFingerprint,"direct audit identity omitted output operator");
    const auto incoherentFolder=SaveRoundDiagnostics(e,state,stage,stats,ids,errors,e.theta);
    Require(coherentFolder!=incoherentFolder&&incoherentFolder!=legacy,"output policy reused another snapshot folder");
    Require(Read(legacy+"/sentinel")=="legacy-preserved","old published snapshot overwritten");
    state.streak=0;e.stablePasses=3;UpdatePassStreak(e,state,true,stage,incoherentFolder);
    Require(state.streak==1,"first new-policy pass was not counted once");
    UpdatePassStreak(e,state,true,stage,incoherentFolder);Require(state.streak==1,"same-policy pass counted twice");
    e.coherent=true;UpdatePassStreak(e,state,true,stage,coherentFolder);Require(state.streak==2,"changed operator reused old pass receipt");
}
void PublishedHints(const std::string& root){
    ArgPP args;RunConfig config;Engine e(args,config);e.root=root;e.signature="published-policy-regression";
    e.go=false;e.coherent=false;e.local=true;e.theta={0,.5,1};e.phi=8;e.interpolation=.01;
    Directory(root);AdaptiveState state;state.knots={0,2};state.counts=Budgets(1,std::vector<int>(3,512));
    state.valid.assign(3,true);state.coarse.assign(3,true);state.previous.resize(3);
    Statistics stats;stats.mean.resize(3);stats.interval.resize(3);stats.ratioInterval.resize(3);for(auto& row:stats.mean)row[0]=100;
    AdaptiveData stage;stage.sources="published hint fixture\n";const std::vector<size_t> ids={0,1,2};
    const auto folder=SaveRoundDiagnostics(e,state,stage,stats,ids,{0,0,0},e.theta);
    state.counts[0].assign(3,1024);e.finalChecks=false;e.stablePasses=1;state.streak=1;
    auto publish=[&](const std::string& policies,bool passed,double residual){
        const auto text="{\"snapshot_folder\":"+Quoted(folder)+policies+",\"numerical_round_passed\":"+(passed?"true":"false")+",\"interpolation_residual95\":"+Number(residual)+"}";
        Atomic(root+"/intermediate_current.json",text);return text;
    };
    const auto angular=",\"angular_sampling_policy\":"+Quoted(AngularSamplingPolicy());
    const auto interpolation=",\"interpolation_policy\":"+Quoted(POOutputInterpolationPolicy(false,false));
    const auto old=",\"angular_sampling_policy\":\"incoherent-component-envelope-noise-v2\"";
    auto check=[&](bool expectedPass,bool expectedResidual){
        Budgets counts;Statistics priorities;std::string saved;bool passed=true,residual=true;
        Require(PublishedRoundBudget(e,state,counts,priorities,saved,passed),"valid published round allocation lost");
        Require(passed==expectedPass,"published round pass policy guard failed");
        Require(counts[0]==std::vector<int>(3,512),"published allocation hints were discarded");
        Require(PublishedSchedulingBudget(e,state,counts,priorities,saved,passed,residual),"valid scheduling allocation lost");
        Require(passed==expectedPass&&residual==expectedResidual,"old published pass/residual hint bypassed policy guard");
    };
    auto text=publish(old,true,.02);Require(!PublishedAccuracyPoliciesCurrent(e,text),"old incoherent accuracy policy accepted");check(false,false);
    RevisitPublishedPass(e,state);Require(state.counts[0]==std::vector<int>(3,1024),"old published pass bypassed re-evaluation");
    text=publish(angular,true,.02);Require(!PublishedAccuracyPoliciesCurrent(e,text),"missing interpolation policy accepted");check(false,false);
    text=publish(interpolation+old,true,.02);Require(!PublishedAccuracyPoliciesCurrent(e,text),"old angular policy accepted with new interpolation");check(false,false);
    text=publish(interpolation+angular,true,.001);Require(PublishedAccuracyPoliciesCurrent(e,text),"current policies rejected");check(true,false);
    RevisitPublishedPass(e,state);Require(state.counts[0]==std::vector<int>(3,512),"compatible published allocation not restored");
    text=publish(interpolation+angular,false,.02);check(false,true);
    e.local=false;Require(PublishedAccuracyPoliciesCurrent(e,publish(interpolation,false,.02)),"disabled angular policy unnecessarily required");
    e.local=true;e.coherent=true;text=publish(old,true,.02);
    Require(PublishedAccuracyPoliciesCurrent(e,text),"coherent legacy hints behavior changed");check(true,true);
    const auto schemaFolder=root+"/schema";Directory(schemaFolder);
    e.coherent=false;WriteOutputGuide(e,schemaFolder,3,3,2,true);
    const auto linearSchema=Read(schemaFolder+"/matrix_output_schema.json");
    Require(linearSchema.find(POOutputInterpolationPolicy(false,false))!=std::string::npos&&linearSchema.find("fixed piecewise-linear interpolation")!=std::string::npos,"linear output schema mislabeled");
    e.coherent=true;WriteOutputGuide(e,schemaFolder,3,3,2,true);
    const auto coherentSchema=Read(schemaFolder+"/matrix_output_schema.json");
    Require(coherentSchema.find("local cubic interpolation")!=std::string::npos&&coherentSchema.find("fixed piecewise-linear interpolation")==std::string::npos,"coherent output schema mislabeled");
}
void RealFixture(const std::string& path){
    std::ifstream in(path);std::string magic;size_t knots,queries,replicas;
    Require(bool(in>>magic>>knots>>queries>>replicas)&&magic=="FULLAUTO_LINEAR_REAL_V1"&&replicas==8&&knots>=2,"invalid real interpolation fixture");
    std::vector<double> x(knots),q(queries);for(double& v:x)Require(bool(in>>v),"fixture support missing");for(double& v:q)Require(bool(in>>v),"fixture query missing");
    Samples data(8,Matrix(knots)),output(8);for(auto& replica:data)for(auto& row:replica)for(double& v:row)Require(bool(in>>v)&&std::isfinite(v),"fixture matrix missing");
    std::string extra;Require(!(in>>extra),"extra fixture fields");
    for(int s=0;s<8;++s){output[s]=POOutputInterpolate(false,false,x,data[s],q);Require(POOutputInterpolate(false,false,x,data[s],x)==data[s],"real knots changed");}
    const auto means=Mean(data),actual=Mean(output),expected=POOutputInterpolate(false,false,x,means,q);
    NoExtraExtrema(x,means,q,expected);double worst=0;
    for(size_t r=0;r<q.size();++r)for(int j=0;j<16;++j){Close(actual[r][j],expected[r][j],"real mean failed linearity",expected[r][0]);worst=std::max(worst,std::fabs(actual[r][j]-expected[r][j])/expected[r][0]);}
    std::cout<<"REAL "<<path<<" knots="<<knots<<" queries="<<queries<<" mean_difference_scaled_M11="<<Number(worst)<<" no_added_extrema_all16_and_ratios=PASS\n";
}
}
int main(int argc,char** argv){
    try{Require(argc>=2,"pass private output directory and optional real fixture files");PureOperator();MidpointResidual();Migration(argv[1]);PublishedHints(std::string(argv[1])+"/published_hints");
        for(int i=2;i<argc;++i)RealFixture(argv[i]);
        std::cout<<"PASS linear all16/ratios, exact knots, replica mean, paired midpoint/control residual, coherent/GO fallback, interruption-safe policy migration, published allocation/pass/residual guards and output schema\n";
        return 0;
    }catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}
}
