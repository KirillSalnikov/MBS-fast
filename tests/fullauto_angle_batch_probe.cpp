#include "../src/FullAuto.cpp"
using namespace FullAuto;
int main(int argc,char** argv){try{
    Require(argc==4,"OUTPUT BINARY NATIVE_COMMAND_FILE");
    ArgPP args;RunConfig config;config.useGpu=true;Engine e(args,config);
    e.root=argv[1];Directory(e.root);Directory(e.root+"/jobs");e.exe=argv[2];e.binaryHash=FileHash(e.exe);e.signature="angle-batch-native-v1";
    std::istringstream input(Read(argv[3]));std::string line;std::getline(input,line);std::vector<std::string> command;
    while(std::getline(input,line))command.push_back(line);
    auto value=[&](const std::string& key){const auto i=std::find(command.begin(),command.end(),key);Require(i!=command.end()&&std::next(i)!=command.end(),"missing "+key);return *std::next(i);};
    e.phi=std::stoi(value("--phi-points"));e.workerThreads=std::stoi(value("--threads"));e.target=std::stoi(value("--max-reflections"));
    e.coherent=std::find(command.begin(),command.end(),"--incoherent")==command.end();
    e.controls=std::find(command.begin(),command.end(),"--analytic-facet-average")!=command.end();
    e.mirror=std::find(command.begin(),command.end(),"--mirror-gamma")!=command.end();
    e.sharedMeans=true;e.gpus={"0"};e.hostBudgetMb=8192;e.workerPolicy="persistent";e.resourcePolicy="auto";
    e.kernel.pipeline=true;e.kernel.block=256;e.eps=.02;e.ratios=true;e.orientationBlock=4096;e.maxCount=1<<22;e.dispatchThreads=6;
    e.depthSharing=e.coherent&&!getenv("TASK_MBS_NO_DEPTH_SHARING");e.depthVerified=true;e.seeds.assign(production,production+8);
    {std::istringstream grid(Read(value("--theta-grid-file")));double x;while(grid>>x)e.theta.push_back(x);}
    Require(e.theta.size()>16,"need a real grid wider than 16 angles");
    if(e.controls){
        e.analyticMeanThreads=6;
        e.requestedMeanCache=e.root+"/pilot_mean.cache";Atomic(e.requestedMeanCache,Read(value("--analytic-mean-cache")));
        Atomic(e.root+"/requested_theta.csv",Read(value("--analytic-mean-reference-grid")));
        std::istringstream weights(Read(value("--analytic-control-weights")));std::getline(weights,line);
        double x,a,b;while(weights>>x>>a>>b){e.weightTheta.push_back(x);e.weights.push_back({{a,b}});}
    }else WriteGrid(e.root+"/requested_theta.csv",e.theta);
    e.base=command;
    auto remove=[&](const std::string& name,size_t n){for(size_t i=0;i<e.base.size();)if(e.base[i]==name)e.base.erase(e.base.begin()+i,e.base.begin()+i+1+n);else ++i;};
    for(const auto& name:{"--max-reflections","--threads","--theta-grid-file","--phi-points","--sobol-offset","--analytic-mean-threads","--analytic-shadow-control","--analytic-mean-cache","--analytic-mean-reference-grid","--analytic-control-weights","--po-prepare-threads","--orientation-chunk","--output"})remove(name,1);
    remove("--sobol-seed",2);
    for(const auto& name:{"--close","--haar-alpha","--profile-phases","--analytic-facet-average","--orientation-pipeline","--mirror-gamma"})remove(name,0);
    // Exercise the corrected estimator independently of the direct off and
    // balanced native benchmarks.
    if(getenv("TASK_MBS_CORRECTION"))e.balancedCorrection=true;
    AdaptiveState state;std::vector<size_t> ids;
    const std::vector<int> depths={e.target-4,e.target-2,e.target};
    state.counts.assign(3,std::vector<int>(e.theta.size(),1024));
    for(size_t i=0;i<e.theta.size();++i){
        ids.push_back(i);if(i%4==0)state.knots.push_back(i);
        if(i>=e.theta.size()/2)state.counts[0][i]=2048;
    }
    if(state.knots.back()!=e.theta.size()-1)state.knots.push_back(e.theta.size()-1);
#ifdef MBS_GPU_FP32
    const auto frozen=state.counts;const auto seeds=e.seeds;const int phi=e.phi;const auto recipe=e.ExecutionHash();
    const bool accepted=VerifyInterpolationAngleBatches(e,state,depths,ids);
    Require(state.counts==frozen&&e.seeds==seeds&&e.phi==phi,"pilot changed production allocation or phi");
    const auto history=Read(e.root+"/fullauto_history.jsonl");
    Require(VerifyInterpolationAngleBatches(e,state,depths,ids)==accepted,"cached merge decision changed");
    Require(Read(e.root+"/fullauto_history.jsonl")==history,"cached pilot repeated optical work");
    Require(!e.collecting&&!e.shortPacketVerifying&&!e.phiCutoffPilot&&e.angleBatchPilotRecipe.empty(),"pilot leaked dispatch state");
    Require(e.ExecutionHash()==recipe,"pilot changed production cache identity");
    if(e.balancedCorrection)Require(history.find("cutoff_correction_start")!=std::string::npos,"pilot skipped the balanced correction");
    std::cout<<"PASS native FP32 merge pilot and cached decision; accepted="<<accepted<<'\n';
#endif
    return 0;
}catch(const std::exception& error){std::cerr<<error.what()<<'\n';return 1;}}
