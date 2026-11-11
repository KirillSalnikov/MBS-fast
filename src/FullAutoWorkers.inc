// Included inside Engine. Processes outlive ExecutePrepared calls, but every
// request has its own result directory, PID ownership record and commit marker.
bool PersistentEnabled()const{return memoryMode||workerPolicy=="persistent"||(workerPolicy=="auto"&&config.useGpu);}
// Older checkpoint-compatible workers report the system RAM guard as status 1
// instead of the dedicated per-process budget status. Match only its exact
// diagnostic; numerical, CUDA and malformed-input failures remain fatal.
bool LegacyHostReserveFailure(const Job& job,int status)const{
    if(status!=1)return false;
    if(MemoryPath(job.folder)){if(!Exists(job.folder+"/stdout.log"))return false;
        return Read(job.folder+"/stdout.log").find("ERROR: calculation failed: host RAM reserve exhausted before preparing next orientation chunk;")!=std::string::npos;}
    std::ifstream file(job.folder+"/stdout.log",std::ios::binary);if(!file)return false;
    file.seekg(0,std::ios::end);const auto size=file.tellg();if(size<0)return false;
    file.seekg(std::max(std::streamoff(0),std::streamoff(size)-8192));
    const std::string tail((std::istreambuf_iterator<char>(file)),std::istreambuf_iterator<char>());
    return tail.find("ERROR: calculation failed: host RAM reserve exhausted before preparing next orientation chunk;")!=std::string::npos;
}
void WaitAfterHostReserveFailure(const Job& job){
    const double started=Seconds();bool waited=false;
    auto memory=HostMemoryReserve::Read();payloads.Clear();
    while(HostMemoryReserve::Low(memory)){
        const std::string detail="{\"folder\":"+Quoted(job.folder)+",\"available_kib\":"+Number(memory.availableKb)
            +",\"reserve_kib\":"+Number(HostMemoryReserve::ReserveKb(memory.totalKb))
            +",\"wait_seconds\":"+Number(Seconds()-started)+",\"maximum_wait_seconds\":60,\"same_sampling\":true,\"completed_blocks_preserved\":true}";
        if(!waited||Seconds()-memoryCheckAt>=5){Record(root+"/resource_wait.json",detail);Log("host_memory_wait",detail);memoryCheckAt=Seconds();}
        waited=true;WriteExecutionProgress();
        Require(Seconds()-started<60,"host RAM reserve remains exhausted after 60 seconds; completed checkpoints preserved; free RAM before resuming");
        usleep(250000);memory=HostMemoryReserve::Read();
    }
    memoryCheckAt=-1;
    Log("host_memory_retry_ready","{\"folder\":"+Quoted(job.folder)+",\"wait_seconds\":"+Number(Seconds()-started)
        +",\"available_kib\":"+Number(memory.availableKb)+",\"reserve_unchanged\":true}");
    if(waited)Record(root+"/resource_wait.json","{\"waiting\":false,\"folder\":"+Quoted(job.folder)+",\"same_sampling\":true}");
}
void CompleteNative(Job& j) {
    j=Load(j);Require(j.matrix.size()==j.inputAngles.size(),"worker returned wrong theta count");
    for(size_t r=0;r<j.inputAngles.size();++r)Require(std::fabs(j.angles[r]-j.inputAngles[r])<1e-10,"worker returned wrong theta coordinates");
    // Packed native jobs stay private: manifests are atomically visible but
    // their durability comes from the verified pack record and its frontier.
    std::unique_ptr<AtomicFileBatch> privateFiles;
    if(packJobs)privateFiles.reset(new AtomicFileBatch(j.folder));
    PackedJobStore::ExpectedMembers expected;
    auto remember=[&](const std::string& path,const std::string& digest){
        if(!packJobs)return;
        Require(path.compare(0,j.folder.size()+1,j.folder+"/")==0,"native receipt member outside job");
        const auto name=path.substr(j.folder.size()+1);const auto inserted=expected.emplace(name,digest);
        if(!inserted.second&&inserted.first->second!=digest)
            Require(false,"conflicting native member digest "+path
                +" (receipt="+inserted.first->second+", member="+digest+")");};
    auto rememberManifest=[&](const std::string& text){if(!packJobs)return;std::istringstream rows(text);std::string line;
        while(std::getline(rows,line)){const auto split=line.find_last_of(' ');Require(split!=std::string::npos,"invalid native integrity manifest");
            remember(line.substr(0,split),line.substr(split+1));}};
    if(packJobs)remember(j.folder+"/command.args",FileHash(j.folder+"/command.args"));
    std::string pairHash;if(j.pairedDepth){const std::string name=j.folder+"/result/result_depth"+std::to_string(j.pairedDepth);std::ostringstream manifest;
        manifest<<name+".dat"<<' '<<FileHash(name+".dat")<<'\n';if(controls)manifest<<name+"_analytic_facets.tsv"<<' '<<FileHash(name+"_analytic_facets.tsv")<<'\n';
        Atomic(j.folder+"/paired_depth.tsv",manifest.str());pairHash=" "+Hash(manifest.str());
        rememberManifest(manifest.str());remember(j.folder+"/paired_depth.tsv",Hash(manifest.str()));}
    std::ostringstream integrals;
    for(int depth:{0,j.pairedDepth}){
        if(depth==0||j.pairedDepth){const auto path=j.folder+"/result/result"+(depth?"_depth"+Number(depth):"")+"_integrals.tsv";
            integrals<<path<<' '<<FileHash(path)<<'\n';}
        if(!j.pairedDepth)break;
    }
    if(j.pairedDepth){const auto path=j.folder+"/result/result_depth"+Number(j.pairedDepth)+"_integrals_revision.txt";
        if(Exists(path))integrals<<path<<' '<<FileHash(path)<<'\n';}
    AppendNativePrefixIntegrity(j,integrals);
    Atomic(j.folder+"/integral_outputs.fnv",integrals.str());
    rememberManifest(integrals.str());remember(j.folder+"/integral_outputs.fnv",Hash(integrals.str()));
    const auto dataHash=FileHash(j.data),summaryHash=controls?FileHash(j.summary):"none";
    remember(j.data,dataHash);if(controls)remember(j.summary,summaryHash);
    // Leave the private staging scope without flushing loose manifests. The
    // typed oversize fallback below supplies a full raw durability fence.
    privateFiles.reset();
    PublishNativeComplete(j,j.identity+" "+Number(j.seconds)+" "+dataHash+" "+summaryHash
        +(pairHash.empty()?" none":pairHash)+" "+Hash(integrals.str())+"\n",expected);CacheJob(j);ObservePreparation(j);
    Log("worker_complete","{\"folder\":"+Quoted(j.folder)+",\"seconds\":"+Number(j.seconds)+",\"trace_prepare_seconds\":"+Number(j.traceSeconds)
        +",\"diffraction_and_controls_seconds\":"+Number(j.diffractionSeconds)+",\"processing_wall_seconds\":"+Number(j.NativeWall())
        +",\"processing_wall_measured\":"+(j.processingSeconds>=0?"true":"false")+",\"outside_processing_seconds\":"+Number(j.OutsideProcessing())
        +",\"analytic_setup_seconds\":"+Number(j.controlSetup)+",\"pipeline_wait_for_prepared_seconds\":"+Number(j.preparedWaitSeconds)
        +",\"pipeline_producer_wait_seconds\":"+Number(j.producerWaitSeconds)+",\"pipeline_wait_measured\":"+(j.preparedWaitSeconds>=0?"true":"false")+",\"return_code\":0}");
}
void ExecutePersistent(std::vector<Job>& result,std::vector<std::vector<std::string>>& commands,const std::vector<size_t>& todo,const std::function<void(size_t)>& onComplete={}){
    if(todo.empty())return;
    // Native completion markers remain individually durable. Only the
    // reconstructible discovery journal is flushed as a bounded batch.
    IndexAppendBatch indexWrites(*this);
    DispatchProgressScope progress(*this,todo.size());
    const auto& workerGpus=WorkerGpus();
    const size_t slots=std::min(workerGpus.size(),size_t(8));
    // Resource tuning may change the GPU list; never keep a second context on
    // a changed resource budget or retain unused tuning workers.
    for(size_t slot=0;slot<services.size();++slot)if(services[slot]&&
        (slot>=slots||services[slot]->identity!=WorkerIdentity(workerGpus[slot])))services[slot].reset();
    std::vector<int> slotJob(slots,-1),memoryRetries(result.size(),0);std::vector<bool> gpuRetried(result.size(),false);
    std::vector<double> begun(result.size()),firstBegun(result.size(),-1);
    size_t next=0;double checkpointSeconds=0,callbackSeconds=0;
    auto launch=[&](size_t slot,size_t index){
        CheckHostMemory();
        Job& j=result[index];ApplyAutomaticPreparation(j,commands[index]);
        // Apply after any dispatch tuning: a retry's host-safe upper bound
        // cannot be raised by a later automatic preparation decision.
        const auto safeCap=gpuChunkRetryCaps.find(WorkerIdentity(workerGpus[slot]));
        if(safeCap!=gpuChunkRetryCaps.end())j.dispatchChunk=GpuChunkRetry::LimitCommand(commands[index],j.dispatchChunk,safeCap->second);
        begun[index]=Seconds();
        if(firstBegun[index]<0)firstBegun[index]=begun[index];
        // Publish ownership before sending work, so a controller-only crash
        // remains detectable. Flush the launch files AFTER the worker starts:
        // raw jobs retain the durability fence while tracing overlaps it.
        // Packed jobs defer it: process-crash ownership remains visible in
        // pagecache, a machine crash kills the worker, and only a durable
        // verified pack frontier can publish completion.
        AtomicFileBatch launchFiles(j.folder);
        std::string command=exe+"\n";for(const auto& v:commands[index])command+=v+"\n";
        Atomic(j.folder+"/command.args",command);
        if(!services[slot]){
            int sockets[2];Require(socketpair(AF_UNIX,SOCK_STREAM|SOCK_CLOEXEC,0,sockets)==0,"cannot create worker channel");
            pid_t pid=fork();if(pid<0){close(sockets[0]);close(sockets[1]);Require(false,"cannot fork persistent worker");}
            if(pid==0){signal(SIGTERM,SIG_DFL);signal(SIGINT,SIG_DFL);close(sockets[0]);fcntl(sockets[1],F_SETFD,0);
                ConfigureWorkerEnvironment(workerGpus[slot]);
                const std::string fd=std::to_string(sockets[1]);execl(exe.c_str(),exe.c_str(),"--fullauto-worker-server",fd.c_str(),static_cast<char*>(nullptr));_exit(127);}
            close(sockets[1]);services[slot].reset(new WorkerService(pid,sockets[0],WorkerIdentity(workerGpus[slot]),memoryMode?root:""));
            activeWorkers[slot]=pid;
        }
        auto& service=*services[slot];++service.requests;j.pid=service.pid;slotJob[slot]=int(index);
        const std::string boot=Read("/proc/sys/kernel/random/boot_id");
        Atomic(j.folder+"/worker.pid",std::to_string(j.pid)+" "+ProcessStart(j.pid)+" "+boot);
        std::vector<std::string> request{j.folder+"/stdout.log",Number(j.block),j.warp};
        request.insert(request.end(),commands[index].begin(),commands[index].end());
        // The store already canonicalizes /./ and repeated slashes. Send the
        // same logical names to the private worker stage; keep the original
        // CLI spelling in the fingerprint and replay command.
        if(memoryMode)for(auto& value:request)if(memoryStore->Owns(value))value=FullAutoMemoryStore::Canonical(value);
        const bool sent=memoryMode?WorkerProtocol::SendMemoryStreamRequest(service.fd,root,MemoryWorkerOutputBudget(j),request,MemoryWorkerStreamInputs(commands[index]),
            [&](const WorkerProtocol::StreamFile& file,const WorkerProtocol::StreamWriter& writer){memoryStore->ReadStream(root+"/"+file.path,writer);}):WorkerProtocol::Send(service.fd,request);
        Require(sent,"persistent worker channel closed; see "+j.folder+"/stdout.log");
        if(!packJobs)launchFiles.Commit();
        Log("worker_start","{\"folder\":"+Quoted(j.folder)+",\"command_file\":"+Quoted(j.folder+"/command.args")+",\"pid\":"+Number(j.pid)
            +",\"gpu\":"+Quoted(workerGpus[slot])+WorkerResourceFields(workerGpus[slot])+",\"count\":"+Number(j.count)+",\"offset\":"+Number(j.offset)+",\"seed\":"+Number(j.seed)
            +",\"depth\":"+Number(j.depth)+",\"paired_depth\":"+Number(j.pairedDepth)+",\"phi\":"+Number(j.phi)+",\"theta_rows\":"+Number(j.inputAngles.size())
            +",\"block\":"+Number(j.block)+",\"warp\":"+Quoted(j.warp)+",\"pipeline\":"+(j.pipeline?"true":"false")+",\"threads\":"+Number(workerThreads)
            +",\"chunk\":"+Number(j.dispatchChunk?j.dispatchChunk:workerChunk)+",\"analytic_mean_threads\":"+Number(dispatchThreads?dispatchThreads:analyticMeanThreads)+",\"prepare_threads\":"+Number(j.preparationThreads?j.preparationThreads:workerThreads)+",\"persistent\":true,\"request_serial\":"+Number(service.requests)+"}");
        NativeStarted(j,workerGpus[slot]);
    };
    try {
        for(size_t slot=0;slot<slots&&next<todo.size();++slot)launch(slot,todo[next++]);
        const auto dispatch=WorkerProtocol::DispatchCompletions(todo.size(),slots,[&](bool block){
            std::vector<pollfd> events(slots);for(size_t slot=0;slot<slots;++slot){events[slot].fd=slotJob[slot]>=0?services[slot]->fd:-1;events[slot].events=POLLIN;}
            WriteExecutionProgress();
            int rc;do{rc=poll(events.data(),events.size(),block?1000:0);}while(rc<0&&errno==EINTR);Require(rc>=0,"cannot wait for persistent workers");
            std::vector<size_t> ready;for(size_t slot=0;slot<slots;++slot)if(events[slot].revents)ready.push_back(slot);return ready;
        },[&](size_t slot,size_t& index){
                index=size_t(slotJob[slot]);Job& j=result[index];int status=-1;
                const bool received=memoryMode?WorkerProtocol::ResultMemoryStream(services[slot]->fd,status,
                    [&](const WorkerProtocol::StreamFile& file,const WorkerProtocol::StreamReader& reader){AcceptMemoryWorkerStream(file,reader,j);}):WorkerProtocol::Result(services[slot]->fd,status);
                const bool hostReserve=received&&LegacyHostReserveFailure(j,status);
                if(received&&status==WorkerProtocol::GpuChunkRetry&&memoryRetries[index]++==0){
                    const double failedSeconds=Seconds()-begun[index];
                    NativeReturned(j,true);
                    const pid_t failedPid=services[slot]->pid;
                    services[slot].reset();slotJob[slot]=-1;j.pid=-1;
                    RemovePath(j.folder+"/worker.pid");
                    // A later controller resume may retry this same job again.
                    // Pick one unused stem for all three diagnostic members;
                    // never overwrite a previous attempt or collide with its
                    // nonempty result directory (including virtual storage).
                    std::string savedStem;unsigned attempt=0;
                    do{Require(attempt<1000000,"too many preserved GPU chunk retries");
                        savedStem=j.folder+"/gpu-chunk-retry_"+Number(failedPid)+"_"+Number(attempt++);
                    }while(Exists(savedStem+".stdout.log")||Exists(savedStem+".command.args")||Exists(savedStem+".result.failed"));
                    const std::string savedLog=savedStem+".stdout.log";
                    Require(RenamePath(j.folder+"/stdout.log",savedLog)==0,"cannot preserve GPU chunk retry log");
                    Require(RenamePath(j.folder+"/command.args",savedStem+".command.args")==0,"cannot preserve GPU chunk retry command");
                    // No output from the failed request may satisfy the new
                    // request. Keep it separate for diagnosis, without any
                    // complete marker or cache publication.
                    if(Exists(j.folder+"/result"))Require(RenamePath(j.folder+"/result",savedStem+".result.failed")==0,"cannot preserve GPU chunk retry partial output");
                    gpuChunkRetryCaps[WorkerIdentity(workerGpus[slot])]=GpuChunkRetry::SafeUpper;
                    j.dispatchChunk=GpuChunkRetry::LimitCommand(commands[index],j.dispatchChunk,workerChunk);
                    gpuRetried[index]=true;
                    Log("worker_gpu_chunk_retry","{\"folder\":"+Quoted(j.folder)+",\"old_pid\":"+Number(failedPid)
                        +",\"retry\":1,\"chunk_cap\":"+Number(j.dispatchChunk)+",\"failed_attempt_seconds\":"+Number(failedSeconds)
                        +",\"total_attempt_seconds\":"+Number(Seconds()-firstBegun[index])+",\"backoff_scope\":\"current_controller_worker_identity\",\"worker_identity\":"+Quoted(Hash(WorkerIdentity(workerGpus[slot])))
                        +",\"same_sampling\":true,\"completed_blocks_preserved\":true,\"failed_outputs_quarantined\":true,\"log\":"+Quoted(savedLog)+"}");
                    launch(slot,index);return false;
                }
                if(received&&(status==WorkerProtocol::MemoryPressure||hostReserve)&&memoryRetries[index]++==0){
                    NativeReturned(j,true);
                    const pid_t failedPid=services[slot]->pid;
                    services[slot].reset();slotJob[slot]=-1;j.pid=-1;
                    RemovePath(j.folder+"/worker.pid");
                    const std::string savedLog=j.folder+"/stdout.memory-retry.log";
                    Require(RenamePath(j.folder+"/stdout.log",savedLog)==0,"cannot preserve memory failure log");
                    Log("worker_memory_retry","{\"folder\":"+Quoted(j.folder)+",\"old_pid\":"+Number(failedPid)
                        +",\"retry\":1,\"system_ram_guard\":"+(hostReserve?"true":"false")+",\"same_sampling\":true,\"completed_blocks_preserved\":true,\"log\":"+Quoted(savedLog)+"}");
                    if(hostReserve)WaitAfterHostReserveFailure(j);
                    launch(slot,index);return false;
                }
                if(memoryMode&&(!received||status!=0)&&Exists(j.folder+"/stdout.log")){Atomic(root+"/last_worker_error.log",Read(j.folder+"/stdout.log"));ExportMemoryPublic(true);}
                Require(received&&status==0,"persistent worker failed; see "+(memoryMode?root+"/last_worker_error.log":j.folder+"/stdout.log"));
                if(gpuRetried[index])Log("worker_gpu_chunk_retry_complete","{\"folder\":"+Quoted(j.folder)
                    +",\"successful_attempt_seconds\":"+Number(Seconds()-begun[index])+",\"total_attempt_seconds\":"+Number(Seconds()-firstBegun[index])
                    +",\"chunk_cap\":"+Number(j.dispatchChunk)+",\"same_sampling\":true}");
                NativeReturned(j);j.pid=-1;j.seconds=Seconds()-begun[index];slotJob[slot]=-1;RemovePath(j.folder+"/worker.pid");
                if(memoryMode&&workerPolicy=="process")services[slot].reset();
                // Feed the GPU before CPU validation and durable checkpoint I/O.
                if(next<todo.size())launch(slot,todo[next++]);
                return true;
        },[&](size_t index){
            const double checkpointBegin=Seconds();CompleteNative(result[index]);const double committedSeconds=Seconds()-checkpointBegin;
            checkpointSeconds+=committedSeconds;NativeCommitted(result[index],committedSeconds);
            MaybeMemoryCheckpoint();
            if(onComplete){const double callbackBegin=Seconds();onComplete(index);callbackSeconds+=Seconds()-callbackBegin;}
        });
        indexWrites.Commit();
        Log("worker_dispatch_summary","{\"requests\":"+Number(todo.size())+",\"slots\":"+Number(slots)
            +",\"ready_batches\":"+Number(dispatch.readyBatches)+",\"ready_batch_peak\":"+Number(dispatch.readyBatchPeak)
            +",\"completed_pending_checkpoint_peak\":"+Number(dispatch.pendingPeak)+",\"completed_pending_checkpoint_limit\":"+Number(2*slots-1)
            +",\"checkpoint_seconds\":"+Number(checkpointSeconds)+",\"callback_seconds\":"+Number(callbackSeconds)
            +",\"refill_before_checkpoint\":true}");
    }catch(...){for(auto& service:services)service.reset();throw;}
}
