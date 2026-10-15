// Included inside Engine. Processes outlive ExecutePrepared calls, but every
// request has its own result directory, PID ownership record and commit marker.
bool PersistentEnabled()const{return workerPolicy=="persistent"||(workerPolicy=="auto"&&config.useGpu);}
void CompleteNative(Job& j) {
    j=Load(j);Require(j.matrix.size()==j.inputAngles.size(),"worker returned wrong theta count");
    for(size_t r=0;r<j.inputAngles.size();++r)Require(std::fabs(j.angles[r]-j.inputAngles[r])<1e-10,"worker returned wrong theta coordinates");
    std::string pairHash;if(j.pairedDepth){const std::string name=j.folder+"/result/result_depth"+std::to_string(j.pairedDepth);std::ostringstream manifest;
        manifest<<name+".dat"<<' '<<FileHash(name+".dat")<<'\n';if(controls)manifest<<name+"_analytic_facets.tsv"<<' '<<FileHash(name+"_analytic_facets.tsv")<<'\n';
        Atomic(j.folder+"/paired_depth.tsv",manifest.str());pairHash=" "+Hash(manifest.str());}
    std::ostringstream integrals;
    for(int depth:{0,j.pairedDepth}){
        if(depth==0||j.pairedDepth){const auto path=j.folder+"/result/result"+(depth?"_depth"+Number(depth):"")+"_integrals.tsv";
            integrals<<path<<' '<<FileHash(path)<<'\n';}
        if(!j.pairedDepth)break;
    }
    if(j.pairedDepth){const auto path=j.folder+"/result/result_depth"+Number(j.pairedDepth)+"_integrals_revision.txt";
        if(Exists(path))integrals<<path<<' '<<FileHash(path)<<'\n';}
    Atomic(j.folder+"/integral_outputs.fnv",integrals.str());
    Atomic(j.folder+"/complete.txt",j.identity+" "+Number(j.seconds)+" "+FileHash(j.data)+" "+(controls?FileHash(j.summary):"none")
        +(pairHash.empty()?" none":pairHash)+" "+Hash(integrals.str())+"\n");CacheJob(j);ObservePreparation(j);
    Log("worker_complete","{\"folder\":"+Quoted(j.folder)+",\"seconds\":"+Number(j.seconds)+",\"trace_prepare_seconds\":"+Number(j.traceSeconds)
        +",\"diffraction_and_controls_seconds\":"+Number(j.diffractionSeconds)+",\"processing_wall_seconds\":"+Number(j.NativeWall())
        +",\"processing_wall_measured\":"+(j.processingSeconds>=0?"true":"false")+",\"outside_processing_seconds\":"+Number(j.OutsideProcessing())
        +",\"analytic_setup_seconds\":"+Number(j.controlSetup)+",\"pipeline_wait_for_prepared_seconds\":"+Number(j.preparedWaitSeconds)
        +",\"pipeline_producer_wait_seconds\":"+Number(j.producerWaitSeconds)+",\"pipeline_wait_measured\":"+(j.preparedWaitSeconds>=0?"true":"false")+",\"return_code\":0}");
}
void ExecutePersistent(std::vector<Job>& result,const std::vector<std::vector<std::string>>& commands,const std::vector<size_t>& todo){
    if(todo.empty())return;
    const auto& workerGpus=WorkerGpus();
    const size_t slots=std::min(workerGpus.size(),size_t(8));
    // Resource tuning may change the GPU list; never keep a second context on
    // a card owned by another slot or retain unused tuning workers.
    for(size_t slot=0;slot<services.size();++slot)if(services[slot]&&
        (slot>=slots||services[slot]->identity!=exe+"\n"+workerGpus[slot]+"\n"+Number(hostBudgetMb)))services[slot].reset();
    std::vector<int> slotJob(slots,-1);std::vector<double> begun(result.size());
    size_t next=0,completed=0;
    auto launch=[&](size_t slot,size_t index){
        Job& j=result[index];begun[index]=Seconds();
        std::string command=exe+"\n";for(const auto& v:commands[index])command+=v+"\n";
        Atomic(j.folder+"/command.args",command);
        if(!services[slot]){
            int sockets[2];Require(socketpair(AF_UNIX,SOCK_STREAM|SOCK_CLOEXEC,0,sockets)==0,"cannot create worker channel");
            pid_t pid=fork();if(pid<0){close(sockets[0]);close(sockets[1]);Require(false,"cannot fork persistent worker");}
            if(pid==0){signal(SIGTERM,SIG_DFL);signal(SIGINT,SIG_DFL);close(sockets[0]);fcntl(sockets[1],F_SETFD,0);
                if(config.useGpu){setenv("CUDA_VISIBLE_DEVICES",workerGpus[slot].c_str(),1);unsetenv("MBS_GPU_DEVICE");}
                if(!getenv("MBS_HOST_MEM_BUDGET_MB"))setenv("MBS_HOST_MEM_BUDGET_MB",Number(hostBudgetMb).c_str(),1);
                const std::string fd=std::to_string(sockets[1]);execl(exe.c_str(),exe.c_str(),"--fullauto-worker-server",fd.c_str(),static_cast<char*>(nullptr));_exit(127);}
            close(sockets[1]);services[slot].reset(new WorkerService(pid,sockets[0],exe+"\n"+workerGpus[slot]+"\n"+Number(hostBudgetMb)));
            activeWorkers[slot]=pid;
        }
        auto& service=*services[slot];++service.requests;j.pid=service.pid;slotJob[slot]=int(index);
        const std::string boot=Read("/proc/sys/kernel/random/boot_id");
        Atomic(j.folder+"/worker.pid",std::to_string(j.pid)+" "+ProcessStart(j.pid)+" "+boot);
        std::vector<std::string> request{j.folder+"/stdout.log",Number(j.block),j.warp};
        request.insert(request.end(),commands[index].begin(),commands[index].end());
        Require(WorkerProtocol::Send(service.fd,request),"persistent worker channel closed; see "+j.folder+"/stdout.log");
        Log("worker_start","{\"folder\":"+Quoted(j.folder)+",\"command_file\":"+Quoted(j.folder+"/command.args")+",\"pid\":"+Number(j.pid)
            +",\"gpu\":"+Quoted(workerGpus[slot])+",\"count\":"+Number(j.count)+",\"offset\":"+Number(j.offset)+",\"seed\":"+Number(j.seed)
            +",\"depth\":"+Number(j.depth)+",\"paired_depth\":"+Number(j.pairedDepth)+",\"phi\":"+Number(j.phi)+",\"theta_rows\":"+Number(j.inputAngles.size())
            +",\"block\":"+Number(j.block)+",\"warp\":"+Quoted(j.warp)+",\"pipeline\":"+(j.pipeline?"true":"false")+",\"threads\":"+Number(workerThreads)
            +",\"chunk\":"+Number(workerChunk)+",\"prepare_threads\":"+Number(j.preparationThreads?j.preparationThreads:workerThreads)+",\"persistent\":true,\"request_serial\":"+Number(service.requests)+"}");
    };
    try {
        for(size_t slot=0;slot<slots&&next<todo.size();++slot)launch(slot,todo[next++]);
        while(completed<todo.size()){
            std::vector<pollfd> events(slots);for(size_t slot=0;slot<slots;++slot){events[slot].fd=slotJob[slot]>=0?services[slot]->fd:-1;events[slot].events=POLLIN;}
            int rc;do{rc=poll(events.data(),events.size(),-1);}while(rc<0&&errno==EINTR);Require(rc>0,"cannot wait for persistent workers");
            for(size_t slot=0;slot<slots;++slot)if(events[slot].revents){
                const size_t index=size_t(slotJob[slot]);Job& j=result[index];int status=-1;
                Require(WorkerProtocol::Result(services[slot]->fd,status)&&status==0,"persistent worker failed; see "+j.folder+"/stdout.log");
                j.pid=-1;j.seconds=Seconds()-begun[index];slotJob[slot]=-1;std::remove((j.folder+"/worker.pid").c_str());
                // Feed the GPU before CPU validation and durable checkpoint I/O.
                if(next<todo.size())launch(slot,todo[next++]);
                CompleteNative(j);++completed;
            }
        }
    }catch(...){for(auto& service:services)service.reset();throw;}
}
