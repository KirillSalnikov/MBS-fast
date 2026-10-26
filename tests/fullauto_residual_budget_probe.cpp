#define main residual_lazy_test_main
#include "fullauto_lazy_restore_probe.cpp"
#undef main
void SeparateBudget(const std::string& root){
    ArgPP args;RunConfig config;Engine e(args,config);SetupMemory(e,root,"/bin/true",false);
    e.theta={170,172,174,176,178,179,180};e.initial=524288;e.minimum=256;e.maxCount=67108864;e.local=true;e.phi=32;
    AdaptiveState state;state.knots={0,2,4,6};state.counts=Budgets(3,std::vector<int>(7,256));
    state.counts[0].assign(7,524288);state.counts[0][2]=33554432;
    state.counts[1][2]=131072;state.previous=Matrix(7);state.valid.assign(7,false);state.coarse.assign(7,false);
    const auto production=state.counts;PrepareResidualBudget(e,state);
    Require(state.counts==production,"residual initialization changed production");
    Require(PeakBudgets(state.interpolationCounts)==std::vector<int>({65536,256,256}),"large production prefix infected residual");
    const auto saved=state.interpolationCounts;state.counts[0][2]*=2;PrepareResidualBudget(e,state);
    Require(state.interpolationCounts==saved,"production growth changed frozen residual allocation");state.counts=production;
    const std::vector<size_t> ids={0,1,2,3,4,5,6};
    InterpolationResult fail;InitInterpolationNoise(fail,7);fail.bad={3};auto& p=fail.growth[3];
    p.level=1;p.upper=.02;p.mean=.001;p.noise=.019;p.counts={{0,512},{2,512},{3,512},{4,512},{6,512}};p.observedCounts.assign(5,256);
    Require(InterpolationSamplingProposal(state,fail)==production,"residual proposal leaked to production");
    RefineInterpolationFailures(e,state,ids,fail,true);
    Require(state.counts==production&&state.knots==std::vector<size_t>({0,2,4,6}),"noise-only residual changed production/grid");
    for(const auto& n:p.counts)Require(state.interpolationCounts[1][n.first]==512,"residual cohort growth incomplete");
    e.thetaPhi[3]=64;SaveAdaptive(e,state);
    AdaptiveState restored;restored.counts=Budgets(3,std::vector<int>(7,256));restored.previous=Matrix(7);restored.valid.assign(7,false);restored.coarse.assign(7,false);
    Require(RestoreAdaptive(e,restored,false)&&restored.counts==state.counts&&restored.interpolationCounts==state.interpolationCounts&&e.PhiForTheta(3)==64,"v5 lost residual allocation or phi");
    const auto text=Read(root+"/adaptive_checkpoint.tsv");Atomic(root+"/adaptive_checkpoint.tsv",text+"corrupt\n");bool rejected=false;
    try{RestoreAdaptive(e,restored,false);}catch(const std::exception&){rejected=true;}Require(rejected,"corrupt residual checkpoint accepted");Atomic(root+"/adaptive_checkpoint.tsv",text);
    // Old states keep their schema and can start a separate allocation later.
    state.interpolationCounts.clear();SaveAdaptive(e,state);Require(RestoreAdaptive(e,restored,false)&&restored.interpolationCounts.empty(),"legacy v4 migration failed");
    e.thetaPhi.clear();SaveAdaptive(e,state);Require(RestoreAdaptive(e,restored,false)&&restored.interpolationCounts.empty(),"legacy v3 migration failed");
    PrepareResidualBudget(e,state);fail.growth[3].biasDominated=true;fail.growth[3].level=-1;fail.growth[3].counts.clear();
    RefineInterpolationFailures(e,state,ids,fail,true);Require(std::binary_search(state.knots.begin(),state.knots.end(),3),"real angular bias did not refine support");
    Require(state.counts==production,"existing production counts lost on angular refinement");
}
void NativeResidual(const std::string& root,const std::string& binary,bool controls){
    ArgPP args;RunConfig config;Engine e(args,config);SetupMemory(e,root,binary,controls);
    e.workQueuePolicy="off";e.workerPolicy="process";e.orientationBlock=64;e.initial=512;e.minimum=64;e.local=true;e.depthSharing=false;
    e.theta={170,172,174,176,178,179,180};e.weightTheta=e.theta;e.weights.assign(7,{{.25,.125}});
    e.polarWeights.resize(7);for(auto& row:e.polarWeights)for(auto& control:row)for(double& value:control)value=.125;
    WriteGrid(root+"/requested_theta.csv",e.theta);
    AdaptiveState state;state.knots={0,2,4,6};state.counts=Budgets(2,std::vector<int>(7,64));state.counts[0].assign(7,512);state.counts[0][2]=33554432;
    PrepareResidualBudget(e,state);const auto production=state.counts;
    AdaptiveState check=state;check.counts=state.interpolationCounts;
    const auto result=CommonPrefixInterpolationImpl(e,check,{4,6},{0,1,2,3,4,5,6},e.seeds.data());
    e.forbidNativeWork=true;
    AdaptiveState legacy=check;legacy.interpolationCounts.clear();const auto reference=CommonPrefixInterpolationImpl(e,legacy,{4,6},{0,1,2,3,4,5,6},e.seeds.data());
    Require(result.matrixErrors==reference.matrixErrors&&result.ratioErrors==reference.ratioErrors&&result.bad==reference.bad&&result.maximum==reference.maximum,"independent allocation changed all16/ratio acceptance math");
    Require(state.counts==production,"residual verification modified production output estimator");
    e.forbidNativeWork=false;e.thetaPhi[1]=e.phi*2;const auto rules=e.thetaPhi;const int savedPhi=e.phi;
    const auto mixed=CommonPrefixInterpolation(e,check,{4,6},{0,1,2,3,4,5,6},e.seeds.data());
    Require(e.phi==savedPhi&&e.thetaPhi==rules,"common residual phi leaked into production");
    e.forbidNativeWork=true;
    {ResidualPhiScope scope(e,true);const auto uniform=CommonPrefixInterpolationImpl(e,legacy,{4,6},{0,1,2,3,4,5,6},e.seeds.data());
        Require(mixed.matrixErrors==uniform.matrixErrors&&mixed.ratioErrors==uniform.ratioErrors&&mixed.maximum==uniform.maximum,
            "common-phi residual differs from native uniform reference");}
    std::cout<<"PASS native residual controls="<<controls<<" exact all16/ratios, eight seeds, source validation, production preserved\n";
}
int main(int argc,char** argv){try{
    Require(argc==3,"OUTPUT BINARY");SeparateBudget(std::string(argv[1])+"/scheduling");
    for(bool controls:{false,true})NativeResidual(std::string(argv[1])+"/native"+Number(controls),argv[2],controls);
    std::cout<<"PASS separate residual growth, preserved large production prefixes, angular bias refinement, v3/v4/v5 resume and corruption rejection\n";return 0;
}catch(const std::exception& error){std::cerr<<error.what()<<'\n';return 1;}}
