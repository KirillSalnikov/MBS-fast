// Build with the same CUDA/precision defines and objects as the native binary.
#include "../src/FullAuto.cpp"
using namespace FullAuto;
int main(int argc,char** argv){try{
    Require(argc==3,"OUTPUT NATIVE_BINARY");ArgPP args;RunConfig config;config.useGpu=true;
    Engine e(args,config);e.root=argv[1];e.exe=argv[2];e.binaryHash=FileHash(e.exe);e.signature="incoherent-tracing-pilot";
    Directory(e.root);Directory(e.root+"/jobs");e.coherent=false;e.controls=false;e.gpus={"0"};e.gpuCandidates=e.gpus;
    e.seeds.assign(production,production+8);e.initial=512;e.maxCount=512;e.phi=32;e.target=12;e.eps=.04;
    e.workerThreads=2;e.workerChunk=256;e.orientationBlock=512;e.hostBudgetMb=2048;e.workerPolicy="persistent";
    e.tracePolicy="auto";e.resourcePolicy="off";e.theta={163.75909223159732,167.13352685050799,170};
    e.base={"--method","po","--backend","cuda","--incoherent","--particle","1","631","174.8",
        "--ri","1.3116","0","--wavelength-um",".532","--symmetry","2","6","--cutoff-profile","off"};
    e.VerifyDepthSharing(e.theta);Require(e.depthSharing,"GPU incoherent shared-depth pilot rejected");
    e.TuneTracing(e.theta);const auto report=Read(e.root+"/trace_backend_calibration.json");
    Require(report.find("\"count_per_seed\":512")!=std::string::npos&&report.find("\"backend\":\"resident\"")!=std::string::npos
        &&report.find("\"backend\":\"hybrid\"")!=std::string::npos,"small incoherent startup skipped tracing calibration");
    const auto before=e.base;e.TuneTracing(e.theta);Require(e.base==before,"explicit saved tracing choice retuned");
    std::cout<<report<<'\n';return 0;
}catch(const std::exception& error){std::cerr<<error.what()<<'\n';return 1;}}
