#include "../src/FullAuto.cpp"
using namespace FullAuto;
int main(int argc,char** argv){try{
    Require(argc==2,"OUTPUT");
    ArgPP args;RunConfig config;config.useGpu=true;Engine e(args,config);
    e.root=argv[1];Directory(e.root);
    e.resourcePolicy="auto";e.workerThreads=8;e.prepareThreadCap=12;e.phi=32;
    const auto recipe=e.ExecutionHash();
    for(int i=0;i<4;++i){Job j;j.identity="wait"+Number(i);j.folder=j.identity;j.executionHash=recipe;
        j.count=512;j.depth=14;j.phi=32;j.hasControls=e.controls;j.mirrored=e.mirror;j.preparationThreads=8;
        j.inputAngles.resize(13);j.traceSeconds=.13;j.diffractionSeconds=.14;j.processingSeconds=.20;j.preparedWaitSeconds=.045;
        e.CacheJob(j);}
    Require(e.AutoPrepareThreads(512,14,0,13)==12,"measured pipeline wait did not trigger useful trial");
    auto& p=e.preparationProfiles.at(e.PreparationKey(recipe,512,14,0,32,e.mirror,e.controls,13));
    for(int i=0;i<4;++i)p.ObserveCost(12,13,.30/512);
    Require(e.AutoPrepareThreads(512,14,0,13)==0,"slower candidate retained");
    for(int i=0;i<8;++i)p.ObserveCost(12,13,.18/512);
    for(int i=0;i<32;++i)p.ObserveWait(.001,.18);
    Require(e.AutoPrepareThreads(512,14,0,13)==12,"successful candidate lost after eliminating wait");
    p.costs.clear();for(int i=0;i<4;++i)p.ObserveCost(8,13,.20/512);
    Require(e.AutoPrepareThreads(512,14,0,13)==0,"diffraction-limited case needlessly explored");
    e.dispatchThreads=10;Require(e.AutoPrepareThreads(512,14,0,13)==10,"explicit dispatch lost priority");
    std::cout<<"PASS preparation: measured-wait trial, slower fallback, stable faster choice, no GPU-bound trial, explicit override\n";
    return 0;
}catch(const std::exception& x){std::cerr<<x.what()<<'\n';return 1;}}
