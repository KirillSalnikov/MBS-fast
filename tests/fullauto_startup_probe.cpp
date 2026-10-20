#define main lazy_restore_test_main
#include "fullauto_lazy_restore_probe.cpp"
#undef main
int main(int argc,char** argv){try{
    Require(argc==3,"OUTPUT BINARY");
    for(int fixed:{0,3}){
        ArgPP args;RegisterCliOptions(args);const char* av[]={"probe","--threads","1","--orientation-chunk","64","--fullauto-orientation-block","64"};args.Parse(7,av);
        RunConfig config;Engine e(args,config);SetupMemory(e,std::string(argv[1])+"/fixed"+Number(fixed),argv[2],true);
        e.autoDepth=true;e.depthOrder="sampling-first";e.depthBudget=.005;e.pilot=128;e.kernelCount=128;e.minPhi=2;e.maxPhi=8;e.fixedPhi=fixed;
        e.workQueuePolicy="off";e.workerPolicy="persistent";e.depthPilot=128;e.maxCount=128;e.initial=128;
        PrepareDepthStartup(e);
        Require(e.controls&&e.pilot==128&&e.kernelCount==128&&e.startupKernel,"startup did not restore production settings");
        Require(fixed?e.phi==fixed:(e.phi==4||e.phi==8),"startup phi outside bounded/explicit candidates");
        const auto history=Read(e.root+"/fullauto_history.jsonl");
        Require(history.find("\"event\":\"phi_candidate\"")==std::string::npos&&history.find("phi_pilot_refinement")==std::string::npos,"startup ran broad or growing phi search");
        size_t candidates=0,at=0;while((at=history.find("\"event\":\"startup_phi_candidate\"",at))!=std::string::npos){++candidates;++at;}
        Require(candidates==size_t(fixed?0:2),"incorrect number of bounded startup candidates");
        for(const auto& entry:e.cache){Require(entry.second.count<=64,"startup pilot grew beyond bound");Require(!entry.second.hasControls,"preliminary pilot fitted controls");}
        e.resourcePolicy="auto";e.rotationPolicy="off";e.cutoffPolicy="off";e.workerChunk=64;
        e.PrepareTuning(e.theta);
        Require(Read(e.root+"/resource_calibration.json").find("\"deferred_until_final_phi\":true")!=std::string::npos,"resources benchmarked before final phi");
        e.phi=fixed?fixed:8;e.weightTheta=e.theta;e.weights.assign(e.theta.size(),{{0.25,0.125}});
        const auto coefficients=e.weights;const int phi=e.phi;e.TuneFinalResources(e.theta);
        const auto resources=Read(e.root+"/resource_calibration.json");
        Require(resources.find("\"phi\":"+Number(phi))!=std::string::npos&&e.phi==phi,"resource calibration used preliminary phi");
        Require(e.controls&&e.weights==coefficients&&resources.find("\"controls\":true")!=std::string::npos,"resource tuning changed final controls");
        Require(e.workerThreads==1&&e.workerChunk==64&&e.orientationBlock==64,"explicit resource controls were retuned");
        std::cout<<"PASS fixed_phi="<<fixed<<": bounded startup, full search deferred, controls restored, resources use final phi, explicit thread/chunk/block retained\n";
    }
    return 0;
}catch(const std::exception& error){std::cerr<<error.what()<<'\n';return 1;}}
