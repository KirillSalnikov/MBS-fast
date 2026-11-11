#include "../src/FullAuto.cpp"
using namespace FullAuto;
int main(int argc,char** argv){try{
    if(argc==2&&std::string(argv[1])=="--promotion-lifecycle"){
        ArgPP args;RunConfig config;std::unique_ptr<Engine> controller(new Engine(args,config));
        controller->exe="/isolated/test-worker";controller->hostBudgetMb=512;controller->gpus={"0","1"};
        const auto first=controller->WorkerIdentity("0"),second=controller->WorkerIdentity("1");
        controller->gpuChunkRetryCaps[first]=128;controller->gpuChunkRetryCaps[second]=64;
        const auto expected=controller->gpuChunkRetryCaps;
        auto workingGpuChunkCaps=std::move(controller->gpuChunkRetryCaps);
        controller.reset(new Engine(args,config));
        controller->exe="/isolated/test-worker";controller->hostBudgetMb=512;controller->gpus={"0","1"};
        controller->gpuChunkRetryCaps=std::move(workingGpuChunkCaps);
        Require(controller->gpuChunkRetryCaps==expected,"controller recreation lost per-worker chunk caps");
        Require(controller->gpuChunkRetryCaps.at(controller->WorkerIdentity("0"))==128,"same worker identity lost cap");
        std::vector<std::string> command={"--sobol-seed","512","23","--orientation-chunk","32"};
        Require(GpuChunkRetry::LimitCommand(command,0,controller->gpuChunkRetryCaps.at(first))==32,"promotion raised explicit small chunk");
        controller->hostBudgetMb=1024;
        Require(!controller->gpuChunkRetryCaps.count(controller->WorkerIdentity("0")),"different resource identity inherited cap");
        std::cout<<"PASS actual Engine destruction/recreation, two worker caps, explicit lower upper and changed resource identity\n";return 0;
    }
    Require(argc==5,"ROOT WORKER MEMORY CHUNK");const std::string root=argv[1];Directory(root);
    ArgPP args;RunConfig config;Engine e(args,config);e.root=root;e.signature="gpu-chunk-retry-fixture";e.exe=argv[2];e.binaryHash=FileHash(e.exe);
    e.controls=false;e.gpus={""};e.workerThreads=1;e.analyticMeanThreads=1;e.hostBudgetMb=512;e.phi=1;e.target=4;e.theta={150,160,170};
    e.workerPolicy="persistent";e.resourcePolicy="off";e.workerChunk=std::atoi(argv[4]);e.runStarted=Seconds();
    e.memoryMode=std::atoi(argv[3])!=0;
    if(e.memoryMode){e.memoryStore.reset(new FullAutoMemoryStore::Store(root,1024*1024,root));e.memoryScope.reset(new FullAutoMemoryStore::Scope(*e.memoryStore));e.memoryStore->Directory(root);}
    Directory(root+"/jobs");
    e.base={"--method","po","--backend","cpu","--particle","1",".2",".1","--ri","1.3116","0","--wavelength-um",".532","--symmetry","1","1","--cutoff-profile","off","--orientation-prefix-output","16"};
    unsigned seed=11;const auto jobs=e.DirectStage(32,&seed,1,4,e.theta,true,false,7,2);
    const pid_t successfulPid=e.services[0]?e.services[0]->pid:-1;
    // A later dispatch on the same worker identity must inherit the cap.
    seed=23;const auto next=e.DirectStage(32,&seed,1,4,e.theta,true,false,19,2);
    Require(successfulPid>0&&e.services[0]&&e.services[0]->pid==successfulPid,"successful persistent worker was not reused");
    std::ofstream processReceipt(root+"/successful_worker.pid");processReceipt<<successfulPid<<'\n';processReceipt.close();
    for(const auto& job:{jobs[0],next[0]}){
        const auto command=Read(job.folder+"/command.args");
        std::ofstream saved(root+(job.seed==11?"/first.command":"/next.command"));saved<<command;
        if(getenv("TASK_EXPECT_GPU_BACKOFF"))Require(job.dispatchChunk>0&&job.dispatchChunk<=std::min(128,e.workerChunk?e.workerChunk:128),"retry cap not retained across dispatches");
    }
    std::ofstream matrix(root+"/matrix.txt");matrix<<std::setprecision(17);for(const auto& row:jobs[0].matrix){for(double v:row)matrix<<v<<' ';matrix<<'\n';}
    if(e.memoryMode){
        // Preserve diagnostic state outside the transient store for assertions.
        e.memoryStore->Snapshot(root+"/test-only-memory-checkpoint.bin","retry fixture");
        std::ofstream listing(root+"/virtual_paths.txt");for(const auto& path:e.memoryStore->List(root,true))listing<<path<<'\n';
    }
    std::cout<<"PASS two native dispatches, complete validated results\n";return 0;
}catch(const std::exception& error){std::cerr<<error.what()<<'\n';return 1;}}
