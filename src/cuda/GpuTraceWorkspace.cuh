// The RawCache mutex serializes every lease. A single bounded workspace per
// process; device changes release the previous reservation on its own device.
// No optical result is reused, and every mutable counter is reset by Trace.
std::unique_ptr<Buffers>& SavedWorkspace(){static std::unique_ptr<Buffers> saved;return saved;}
struct WorkspaceLease {
    std::unique_ptr<Buffers> buffer;size_t limit=0;bool success=false;
    explicit WorkspaceLease(int device){
        if(const char* value=std::getenv("MBS_GPU_TRACE_CONVEX_WORKSPACE_MIB")){
            char* end=nullptr;const long mib=std::strtol(value,&end,10);
            if(end==value||*end||mib<0||mib>256)
                throw std::runtime_error("MBS_GPU_TRACE_CONVEX_WORKSPACE_MIB must be 0..256");
            limit=size_t(mib)*1024*1024;
        }
        auto& saved=SavedWorkspace();
        if(saved&&(!limit||saved->device!=device))saved.reset();
        if(saved){saved->Trim(limit);buffer=std::move(saved);}
        else buffer.reset(new Buffers(device));
        buffer->reused=buffer->allocated=0;
    }
    ~WorkspaceLease(){
        // A successful direct consumer has already downloaded its result;
        // raw-download mode also synchronizes before returning. Failed work
        // always destroys its workspace and follows the existing CPU retry.
        if(!success||!limit)return;
        buffer->Trim(limit);
        size_t free=0,total=0;
        if(cudaMemGetInfo(&free,&total)!=cudaSuccess||free<512ULL*1024*1024)return;
        std::fprintf(stderr,"GPU convex workspace: reused_buffers=%d allocations=%d retained_MiB=%.3f limit_MiB=%.3f\n",
            buffer->reused,buffer->allocated,buffer->Bytes()/1048576.,limit/1048576.);
        SavedWorkspace()=std::move(buffer);
    }
};
