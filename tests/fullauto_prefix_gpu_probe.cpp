// This probe explicitly tests a native build with prefix-output support.
#define MBS_NATIVE_ORIENTATION_PREFIX_OUTPUT 1
#include "../src/FullAuto.cpp"
using namespace FullAuto;
void SetupPrefixGpu(Engine& e,const std::string& root,const std::string& binary){
    e.root=root;e.signature="prefix-gpu-gate-regression-v1";Directory(root);Directory(root+"/jobs");e.exe=binary;e.binaryHash=FileHash(binary);
    e.controls=false;e.depthSharing=true;e.gpus={"0"};e.seeds.assign(production,production+8);e.workerThreads=2;e.hostBudgetMb=1024;
    e.phi=32;e.target=8;e.theta={0,170,175,180};e.initial=512;e.orientationBlock=512;e.workerPolicy="persistent";e.resourcePolicy="off";e.eps=.02;
    e.kernel.pipeline=true;e.base={"--method","po","--backend","gpu","--particle","1","2","1","--ri","1.3116","0","--wavelength-um",".532","--symmetry","1","1","--cutoff-profile","off"};e.selectedCutoff="off";
    e.balancedCorrection=true;WriteGrid(root+"/requested_theta.csv",e.theta);
}
int main(int argc,char** argv){try{
    Require(argc==3,"OUTPUT BINARY");ArgPP args;RunConfig config;config.useGpu=true;Engine e(args,config);SetupPrefixGpu(e,argv[1],argv[2]);e.WarmCache();
    auto first=e.DirectStage(512,e.seeds.data(),1,8,e.theta,true,false,0,6);
    Require(e.cutoffPrefixChecks.size()==1,"missing cached numerical gate");const bool accepted=e.cutoffPrefixChecks.begin()->second;
    auto next=e.DirectStage(512,e.seeds.data(),1,8,e.theta,true,false,512,6);
    Require(e.cutoffPrefixChecks.size()==1,"recalibrated unchanged execution profile");
    for(const auto* group:{&first,&next})Require(Exists((*group)[0].folder+"/cutoff_prefix_recipe.txt")==accepted,"gate decision not applied to production");
    e.eps=1e-13;auto strict=e.DirectStage(512,e.seeds.data(),1,8,e.theta,true,false,1024,6);
    Require(e.cutoffPrefixChecks.size()==2,"tighter tolerance reused loose validation");
    Require(!Exists(strict[0].folder+"/cutoff_prefix_recipe.txt"),"strict FP32 gate failed to restore original estimator");
    for(auto& service:e.services)service.reset();
    Engine restored(args,config);SetupPrefixGpu(restored,argv[1],argv[2]);restored.WarmCache();restored.forbidNativeWork=true;
    auto resumed=restored.DirectStage(512,e.seeds.data(),1,8,e.theta,true,false,0,6);
    Require(resumed[0].matrix==first[0].matrix,"completed prefix transaction changed on restart");
    std::cout<<"PASS GPU prefix runtime gate, cached decision, strict-tolerance fallback and checkpoint; accepted="<<accepted<<'\n';return 0;
}catch(const std::exception& ex){std::cerr<<ex.what()<<'\n';return 1;}}
