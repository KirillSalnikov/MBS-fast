#include "../src/FullAuto.cpp"
using namespace FullAuto;
int main(int argc,char** argv){try{
    Require(argc==2||argc==3,"output and optional real fixture required");ArgPP args;RunConfig config;Engine e(args,config);
    Require(!DepthPilotFailureResolved(.010558647935480323,.0140348652046049,.01,true),"noisy real terminal pilot rejected prematurely");
    Require(DepthPilotFailureResolved(.010558647935480323,.0140348652046049,.01,false),"nonterminal depth scheduling changed");
    Require(DepthPilotFailureResolved(.02,.001,.01,true),"resolved excessive depth residual accepted");
    e.root=argv[1];Directory(e.root);Directory(e.root+"/jobs");e.signature="positive-probe";
    e.theta={150,155,160};e.seeds.assign(production,production+8);e.controls=false;
    e.local=true;e.maxCount=1024;e.initial=256;e.minimum=256;e.orientationBlock=256;
    e.incremental=true;e.gpus={""};e.phi=4;e.target=8;e.eps=.04;e.statistical=.03;e.interpolation=.01;
    AdaptiveState state;state.knots={0,1,2};state.counts=Budgets(2,std::vector<int>(3,256));
    state.valid.assign(3,true);state.coarse.assign(3,true);state.previous.resize(3);state.streak=1;
    AdaptiveData stage;stage.levels.assign(2,Samples(8,Matrix(3)));
    stage.costs.assign(2,std::vector<double>(3,.001));stage.analytic.assign(3,0);
    for(auto& sample:stage.levels[0])for(auto& row:sample)row[0]=10;
    for(auto& sample:stage.levels[1])sample[1][0]=-12;
    const auto samples=Combine(stage.levels),original=samples;
    Require(Mean(samples)[1][0]==-2,"fixture must be signed");
    Require(RecoverNonpositiveProduction(e,state,{6,8},stage,{0,1,2},samples,"cached_promotion")==1,"negative mean aborted instead of scheduling");
    Require(samples==original,"recovery clipped the estimator");
    Require(state.counts[0][0]==256&&state.counts[0][2]==256,"good angles recalculated");
    Require(state.counts[0][1]==512&&state.counts[1][1]==512,"all signed levels must grow");
    Require(state.streak==0&&!state.valid[1]&&!state.coarse[1],"negative row accepted");
    for(auto& level:state.counts)level[1]=1024;
    Require(RecoverNonpositiveProduction(e,state,{6,8},stage,{0,1,2},samples,"production")==2,"sampling limit ignored");
    for(auto& sample:stage.levels[1])sample[1][0]=0;
    Require(RecoverNonpositiveProduction(e,state,{6,8},stage,{0,1,2},Combine(stage.levels),"production")==0,"ordinary positive estimator changed");
    if(argc==3){
        std::ifstream file(argv[2]);size_t rows,levels,seeds;
        Require(bool(file>>rows>>levels>>seeds)&&rows==3&&levels==3&&seeds==8,"bad real fixture dimensions");
        stage.levels.assign(levels,Samples(8,Matrix(rows)));stage.costs.assign(levels,std::vector<double>(rows,.001));
        for(auto& level:stage.levels)for(auto& sample:level)for(auto& row:sample)for(double& v:row)Require(bool(file>>v),"truncated real fixture");
        const auto real=Combine(stage.levels);const auto originalReal=real;
        Require(std::fabs(Mean(real)[0][0]+430.376316226775)<1e-8,"real cached-promotion failure not reproduced");
        state.counts=Budgets(levels,std::vector<int>(rows,256));state.round=1;
        Require(RecoverNonpositiveProduction(e,state,{14,16,18},stage,{0,1,2},real,"cached_promotion")==1,"real failure not recovered");
        Require(real==originalReal,"real Mueller data clipped");
        std::cout<<"PASS: replayed signed cached promotion from L1000, theta160-170\n";
    }
    std::cout<<"PASS: signed estimates retained, bad angles grow, finite budget respected\n";return 0;
}catch(const std::exception& error){std::cerr<<error.what()<<'\n';return 1;}}
