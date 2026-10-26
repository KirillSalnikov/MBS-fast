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
struct DispatchProgress {size_t requested=0,committed=0;};
std::vector<DispatchProgress> dispatchProgress;

void WriteExecutionProgress(bool force=false){
    const double now=Seconds();
    if(!force&&executionProgressWritten>=0&&now-executionProgressWritten<10)return;
    const double utc=std::chrono::duration<double>(std::chrono::system_clock::now().time_since_epoch()).count();
    std::ostringstream out;out<<std::setprecision(17)
        <<"{\"schema\":1,\"controller_pid\":"<<getpid()<<",\"unix_seconds\":"<<utc
        <<",\"round\":"<<executionRound<<",\"phase\":"<<Quoted(executionPhase)
        <<",\"phase_elapsed_seconds\":"<<now-executionPhaseStarted
        <<",\"native_requests_submitted\":"<<nativeSubmitted<<",\"native_requests_returned\":"<<nativeReturned
        <<",\"native_requests_committed\":"<<nativeCommitted<<",\"native_requests_retried\":"<<nativeRetries
        <<",\"committed_orientation_samples_all_seeds\":"<<nativeOrientations
        <<",\"summed_worker_seconds\":"<<nativeWorkerSeconds<<",\"checkpoint_seconds\":"<<nativeCheckpointSeconds
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
    out<<"]}";Record(root+"/execution_progress.json",out.str());executionProgressWritten=now;
}
void SetExecutionPhase(const std::string& phase,int round=-1){
    executionPhase=phase;if(round>=0)executionRound=round;executionPhaseStarted=Seconds();
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
    nativeCheckpointSeconds+=checkpointSeconds;lastNativeCommit=Seconds();
    Require(!dispatchProgress.empty(),"native commit outside dispatch progress");
    ++dispatchProgress.back().committed;WriteExecutionProgress();
}
struct DispatchProgressScope {
    Engine& engine;
    DispatchProgressScope(Engine& e,size_t count):engine(e){DispatchProgress p;p.requested=count;engine.dispatchProgress.push_back(p);}
    ~DispatchProgressScope(){engine.dispatchProgress.pop_back();}
};
