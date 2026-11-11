// Included inside Engine. Observation only: never used by allocation or by
// an accuracy/acceptance gate. A native batch is not a whole adaptive round.
struct LiveNativeRequest {
    std::string folder,gpu;int count,offset,depth,phi;unsigned seed;size_t rows;
    double thetaMin,thetaMax,started;
};
std::map<pid_t,LiveNativeRequest> liveNative;
std::string executionPhase="native_execution";
int executionRound=-1;
double executionPhaseStarted=Seconds(),executionProgressWritten=-1,lastNativeCommit=-1;
uint64_t nativeSubmitted=0,nativeReturned=0,nativeCommitted=0,nativeRetries=0,nativeOrientations=0;
double nativeWorkerSeconds=0,nativeCheckpointSeconds=0;
double nativeTraceSeconds=0,nativeConsumerSeconds=0,nativeProcessingSeconds=0;
double nativeGpuDiffractionSeconds=0,nativeAnalyticControlSeconds=0,nativeAnalyticSetupSeconds=0;
double nativePreparedWaitSeconds=0,nativeProducerWaitSeconds=0;
uint64_t nativeProfiledBlocks=0,nativeGpuProfiledBlocks=0,nativePipelineChunks=0;
struct DispatchProgress {size_t requested=0,committed=0;};
std::vector<DispatchProgress> dispatchProgress;
struct CompletedNativeBlock {std::string folder;unsigned seed;int offset,count,depth,phi;size_t rows;double lo,hi,completed;};
std::vector<CompletedNativeBlock> recentNativeBlocks;
uint64_t phaseNativeCommitted=0,phaseAngleOrientationSamples=0;

void WriteExecutionProgress(bool force=false){
    const double now=Seconds();
    if(!force&&executionProgressWritten>=0&&now-executionProgressWritten<(memoryMode?2:10))return;
    const double utc=std::chrono::duration<double>(std::chrono::system_clock::now().time_since_epoch()).count();
    std::ostringstream out;out<<std::setprecision(17)
        <<"{\"schema\":1,\"controller_pid\":"<<getpid()<<",\"unix_seconds\":"<<utc
        <<",\"round\":"<<executionRound<<",\"phase\":"<<Quoted(executionPhase)
        <<",\"phase_elapsed_seconds\":"<<now-executionPhaseStarted
        <<",\"phase_native_blocks_committed\":"<<phaseNativeCommitted
        <<",\"phase_angle_orientation_samples_committed\":"<<phaseAngleOrientationSamples
        <<",\"native_requests_submitted\":"<<nativeSubmitted<<",\"native_requests_returned\":"<<nativeReturned
        <<",\"native_requests_committed\":"<<nativeCommitted<<",\"native_requests_retried\":"<<nativeRetries
        <<",\"committed_orientation_samples_all_seeds\":"<<nativeOrientations
        <<",\"summed_worker_seconds\":"<<nativeWorkerSeconds<<",\"checkpoint_seconds\":"<<nativeCheckpointSeconds
        <<",\"native_phase_totals\":{\"profiled_blocks\":"<<nativeProfiledBlocks
        <<",\"gpu_split_profiled_blocks\":"<<nativeGpuProfiledBlocks<<",\"pipeline_chunks\":"<<nativePipelineChunks
        <<",\"trace_prepare\":"<<nativeTraceSeconds<<",\"diffraction_and_controls\":"<<nativeConsumerSeconds
        <<",\"gpu_diffraction_including_packing\":"<<nativeGpuDiffractionSeconds
        <<",\"analytic_controls\":"<<nativeAnalyticControlSeconds<<",\"analytic_setup\":"<<nativeAnalyticSetupSeconds
        <<",\"processing_wall\":"<<nativeProcessingSeconds
        <<",\"pipeline_wait_for_prepared\":"<<nativePreparedWaitSeconds
        <<",\"pipeline_producer_wait\":"<<nativeProducerWaitSeconds
        <<",\"scope\":\"completed native requests in this invocation; stage sums overlap; GPU stages include host overhead, not device busy time\"}"
        <<",\"completed_job_storage\":"<<Quoted(memoryMode?"memory_with_temporary_spill":packJobs?"packed":"individual_files")
        <<",\"native_blocks_durable\":"<<(memoryMode?"false":"true")
        <<",\"recovery_snapshot_interval_seconds\":"<<(memoryMode?checkpointInterval:0)
        <<",\"packed_inodes_released_this_invocation\":"<<packedInodesReleased
        <<",\"oversized_jobs_left_unpacked_this_invocation\":"<<unpackedOversizeJobs.size()
        <<",\"last_native_commit_age_seconds\":"<<(lastNativeCommit>=0?Number(now-lastNativeCommit):"null")
        <<",\"batch_requested\":"<<(dispatchProgress.empty()?0:dispatchProgress.back().requested)
        <<",\"batch_committed\":"<<(dispatchProgress.empty()?0:dispatchProgress.back().committed)
        <<",\"batch_active\":"<<(dispatchProgress.empty()?"false":"true")
        <<",\"counter_scope\":\"this controller invocation; child corrections count separately\""
        <<",\"batch_scope\":\"current native dispatch only; not whole-round completion\""
        <<",\"accuracy_source\":\"intermediate_current.json\",\"accuracy_recomputed\":false,\"active_workers\":[";
    bool first=true;for(const auto& item:liveNative){const auto& w=item.second;
        out<<(first?"":",")<<"{\"pid\":"<<item.first<<",\"gpu\":"<<Quoted(w.gpu)
            <<",\"folder\":"<<Quoted(w.folder)<<",\"count\":"<<w.count<<",\"offset\":"<<w.offset
            <<",\"seed\":"<<w.seed<<",\"depth\":"<<w.depth<<",\"phi\":"<<w.phi
            <<",\"theta_rows\":"<<w.rows<<",\"theta_min\":"<<w.thetaMin<<",\"theta_max\":"<<w.thetaMax
            <<",\"elapsed_seconds\":"<<now-w.started<<"}";first=false;}
    out<<"],\"recent_completed_blocks\":[";first=true;
    for(const auto& b:recentNativeBlocks){
        out<<(first?"":",")<<"{\"folder\":"<<Quoted(b.folder)<<",\"seed\":"<<b.seed<<",\"offset\":"<<b.offset
            <<",\"count\":"<<b.count<<",\"depth\":"<<b.depth<<",\"phi\":"<<b.phi<<",\"theta_rows\":"<<b.rows
            <<",\"theta_min\":"<<b.lo<<",\"theta_max\":"<<b.hi<<",\"committed_age_seconds\":"<<now-b.completed<<"}";first=false;
    }
    out<<"],\"completed_block_scope\":"<<Quoted(memoryMode
        ?"verified working values; recovery requires a published periodic snapshot; angular bounds may contain gaps; not a new accuracy estimate"
        :"durably saved native blocks; folder is a logical path when packed; angular bounds may contain gaps; not a new accuracy estimate")<<"}";
    Record(root+"/execution_progress.json",out.str());executionProgressWritten=now;
}
void SetExecutionPhase(const std::string& phase,int round=-1){
    executionPhase=phase;if(round>=0)executionRound=round;executionPhaseStarted=Seconds();
    phaseNativeCommitted=0;phaseAngleOrientationSamples=0;recentNativeBlocks.clear();
    WriteExecutionProgress(true);
}
void NativeStarted(const Job& job,const std::string& gpu){
    Require(job.pid>0&&!job.inputAngles.empty(),"invalid live native request");
    LiveNativeRequest info;info.folder=job.folder;info.gpu=gpu;info.count=job.count;info.offset=job.offset;
    info.depth=job.depth;info.phi=job.phi;info.seed=job.seed;info.rows=job.inputAngles.size();
    info.thetaMin=job.inputAngles.front();info.thetaMax=job.inputAngles.back();info.started=Seconds();
    liveNative[job.pid]=info;++nativeSubmitted;WriteExecutionProgress();
}
void NativeReturned(const Job& job,bool retry=false){
    const auto found=liveNative.find(job.pid);
    Require(found!=liveNative.end()&&found->second.folder==job.folder,"live worker identity changed");
    liveNative.erase(found);if(retry)++nativeRetries;else ++nativeReturned;
}
void NativeCommitted(const Job& job,double checkpointSeconds){
    ++nativeCommitted;nativeOrientations+=uint64_t(job.count);nativeWorkerSeconds+=job.seconds;
    if(job.processingSeconds>=0){++nativeProfiledBlocks;nativeTraceSeconds+=job.traceSeconds;
        nativeConsumerSeconds+=job.diffractionSeconds;nativeProcessingSeconds+=job.processingSeconds;
        nativePreparedWaitSeconds+=std::max(0.,job.preparedWaitSeconds);
        nativeProducerWaitSeconds+=std::max(0.,job.producerWaitSeconds);nativePipelineChunks+=job.pipelineChunks;}
    if(job.gpuDiffractionSeconds>=0&&job.analyticControlSeconds>=0){++nativeGpuProfiledBlocks;
        nativeGpuDiffractionSeconds+=job.gpuDiffractionSeconds;nativeAnalyticControlSeconds+=job.analyticControlSeconds;}
    nativeAnalyticSetupSeconds+=std::max(0.,job.controlSetup);
    nativeCheckpointSeconds+=checkpointSeconds;lastNativeCommit=Seconds();
    Require(!dispatchProgress.empty(),"native commit outside dispatch progress");
    ++dispatchProgress.back().committed;++phaseNativeCommitted;
    phaseAngleOrientationSamples+=uint64_t(job.count)*job.inputAngles.size()*job.phi;
    if(recentNativeBlocks.size()==16)recentNativeBlocks.erase(recentNativeBlocks.begin());
    recentNativeBlocks.push_back({job.folder,job.seed,job.offset,job.count,job.depth,job.phi,
        job.inputAngles.size(),job.inputAngles.front(),job.inputAngles.back(),lastNativeCommit});
    WriteExecutionProgress();
}
struct DispatchProgressScope {
    Engine& engine;
    DispatchProgressScope(Engine& e,size_t count):engine(e){DispatchProgress p;p.requested=count;engine.dispatchProgress.push_back(p);}
    ~DispatchProgressScope(){engine.dispatchProgress.pop_back();}
};
