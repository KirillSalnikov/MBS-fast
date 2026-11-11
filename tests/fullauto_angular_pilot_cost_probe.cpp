#include "../src/FullAutoAngularPilotCost.h"
#include "../src/FullAuto.cpp"
using namespace FullAuto;
namespace {
void Close(double a,double b,const char* message){
    Require(std::isfinite(a)&&std::isfinite(b)&&std::fabs(a-b)<=1e-10*std::max(1.,std::fabs(b)),message);
}
void CheckUniform(Engine& e,const Samples& samples,const std::vector<double>& theta,
                  const FullAutoAngularPilotCost::Result& pilot){
    AdaptiveState state;state.counts=Budgets(1,std::vector<int>(theta.size(),512));
    std::vector<size_t> ids;for(size_t r=0;r<theta.size();++r){state.knots.push_back(r);ids.push_back(r);}
    e.theta=theta;e.eps=.04;e.interpolation=.01;e.coherent=false;e.go=false;e.local=true;
    e.maxCount=1<<24;e.fixedPhi=e.phi=8;
    AdaptiveData data;data.levels.push_back(samples);const auto original=data.levels;const auto counts=state.counts;
    const auto stats=Estimate(samples);const auto gate=CheckAngularSampling(e,state,data,ids,samples,stats);
    Require(pilot.resolved&&pilot.finite,"finite real pilot marked unresolved");
    Close(pilot.maximumRatio,gate.maximumRatio,"pilot score differs from actual angular acceptance gate");
    Close(pilot.maximumNoise95,gate.maximum,"pilot noise differs from actual angular acceptance gate");
    Require(pilot.pairs==gate.pairs&&pilot.shapes==gate.shapes,"pilot stencil coverage differs from actual gate");
    Require((pilot.maximumRatio<=1)==gate.passed,"pilot/gate decision mismatch on uniform supports");
    Require(data.levels==original&&state.counts==counts,"pilot/gate altered current sample allocation");
    Require(e.fixedPhi==8&&e.phi==8,"pilot cost modified an explicitly fixed phi");
}
}
int main(int argc,char** argv){try{
    Require(argc==2,"expected measured pilot fixture");std::ifstream file(argv[1]);size_t cases,seeds,rows;
    Require(bool(file>>cases>>seeds>>rows)&&cases==7&&seeds==8&&rows==33,"bad measured pilot fixture header");
    double bestOld=1e300,bestNew=1e300;int oldPhi=0,newPhi=0;Samples first;std::vector<double> firstTheta;
    for(size_t n=0;n<cases;++n){
        int phi,count,element,ratio,span;size_t left,center,right;
        double seconds,oldScore,expectedRatio,expectedScore,expectedInterval,expectedRatioInterval;
        Require(bool(file>>phi>>count>>seconds>>oldScore>>expectedRatio>>expectedScore>>element>>ratio>>span
            >>left>>center>>right>>expectedInterval>>expectedRatioInterval),"missing measured pilot metadata");
        std::vector<double> theta(rows);for(double& x:theta)Require(bool(file>>x),"missing pilot theta");
        Samples samples(seeds,Matrix(rows));for(auto& sample:samples)for(auto& row:sample)for(double& v:row)
            Require(bool(file>>v)&&std::isfinite(v),"missing measured pilot Mueller value");
        const auto before=samples;const auto stats=Estimate(samples);
        Close(stats.maxInterval,expectedInterval,"saved calibration interval not reproduced");
        Close(stats.maxRatioInterval,expectedRatioInterval,"saved calibration Fieller interval not reproduced");
        const auto pilot=FullAutoAngularPilotCost::Evaluate(samples,stats,theta,.04,.01,true,t95);
        Close(pilot.maximumRatio,expectedRatio,"measured Python angular replay mismatch");
        Require(pilot.worstElement==element&&pilot.worstNormalized==bool(ratio)&&pilot.worstSpan==span
            &&pilot.worstLeft==left&&pilot.worstCenter==center&&pilot.worstRight==right,"worst measured contrast metadata differs");
        const double statistical=.0190234375;
        const double oldMultiplier=std::max(1.,std::pow(std::max(stats.maxInterval,stats.maxRatioInterval)/statistical,2));
        const double setup=oldScore-seconds*oldMultiplier;
        const double score=seconds*std::max(oldMultiplier,pilot.maximumRatio*pilot.maximumRatio)+setup;
        Close(score,expectedScore,"measured angular work score mismatch");
        if(oldScore<bestOld){bestOld=oldScore;oldPhi=phi;}if(score<bestNew){bestNew=score;newPhi=phi;}
        ArgPP args;RunConfig config;Engine e(args,config);e.ratios=true;CheckUniform(e,samples,theta,pilot);
        e.ratios=false;const auto noRatios=FullAutoAngularPilotCost::Evaluate(samples,stats,theta,.04,.01,false,t95);
        CheckUniform(e,samples,theta,noRatios);
        Require(samples==before,"pilot changed measured samples");
        if(first.empty()){first=samples;firstTheta=theta;}
        std::cout<<"PASS measured phi="<<phi<<" old_score="<<oldScore<<" angular_ratio="<<pilot.maximumRatio<<" new_score="<<score<<'\n';
    }
    std::string extra;Require(!(file>>extra),"trailing pilot fixture data");
    Require(oldPhi==8&&newPhi==32,"measured ranking did not change from phi8 to phi32");
    {
        const std::vector<double> theta={150.,150.1,151.2,155.3,160.};
        Samples quiet(8,Matrix(theta.size()));
        for(auto& sample:quiet)for(size_t r=0;r<theta.size();++r){
            sample[r][0]=100+10*r;sample[r][5]=r%2?25:-15;sample[r][14]=2*r;
        }
        auto stats=Estimate(quiet);auto pilot=FullAutoAngularPilotCost::Evaluate(quiet,stats,theta,.04,.01,true,t95);
        Require(pilot.resolved&&pilot.maximumRatio<1,"resolved deterministic curvature penalized");
        ArgPP args;RunConfig config;Engine e(args,config);e.ratios=true;CheckUniform(e,quiet,theta,pilot);
        for(size_t seed=0;seed<8;++seed)quiet[seed][2][5]+=(seed<4?-1.:1.);
        stats=Estimate(quiet);pilot=FullAutoAngularPilotCost::Evaluate(quiet,stats,theta,.04,.01,true,t95);
        Require(pilot.resolved&&pilot.maximumRatio>1,"irregular-grid bend noise missed");
        CheckUniform(e,quiet,theta,pilot);
        for(auto& sample:quiet)sample.resize(2);
        const std::vector<double> endpoints={150.,160.};stats=Estimate(quiet);
        pilot=FullAutoAngularPilotCost::Evaluate(quiet,stats,endpoints,.04,.01,true,t95);
        Require(pilot.pairs==1&&pilot.shapes==0,"two-row pilot acquired nonexistent shape stencils");
        CheckUniform(e,quiet,endpoints,pilot);
        auto duplicate=endpoints;duplicate[1]=duplicate[0];
        Require(!FullAutoAngularPilotCost::Evaluate(quiet,stats,duplicate,.04,.01,true,t95).resolved,"duplicate pilot angles accepted");
    }
    auto stats=Estimate(first);
    stats.ratioInterval[2][5]=1e100;
    auto result=FullAutoAngularPilotCost::Evaluate(first,stats,firstTheta,.04,.01,true,t95);
    Require(!result.resolved&&result.finite&&result.maximumRatio==1e50,"unresolved Fieller interval became a cheap finite estimate");
    Require(FullAutoAngularPilotCost::Evaluate(first,stats,firstTheta,.04,.01,false,t95).resolved,"disabled ratios affected raw cost");
    stats=Estimate(first);stats.mean[0][0]=0;
    Require(!FullAutoAngularPilotCost::Evaluate(first,stats,firstTheta,.04,.01,true,t95).resolved,"zero pilot intensity accepted");
    stats=Estimate(first);first[0][0][0]=std::numeric_limits<double>::quiet_NaN();
    result=FullAutoAngularPilotCost::Evaluate(first,stats,firstTheta,.04,.01,true,t95);
    Require(!result.resolved&&!result.finite,"NaN pilot data accepted");
    first[0][0][0]=1;first.pop_back();
    Require(!FullAutoAngularPilotCost::Evaluate(first,stats,firstTheta,.04,.01,true,t95).resolved,"wrong replica count accepted");
    std::cout<<"PASS real ranking8->32, actual gate equivalence, fixed-phi immutability, irregular/2-row grids, deterministic curvature, ratios-off and invalid-input guards\n";
    return 0;
}catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}}
