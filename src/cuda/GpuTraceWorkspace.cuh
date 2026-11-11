// RawCache's mutex serializes every lease and boundary release. The same single
// owner slot is either a healthy cache or a permanently quarantined owner.
// No optical result is reused; Trace resets every mutable counter.
std::unique_ptr<Buffers>& SavedWorkspace(){static std::unique_ptr<Buffers> saved;return saved;}
std::atomic<bool>& WorkspaceQuarantined(){static std::atomic<bool> failed(false);return failed;}
void RequireHealthyWorkspace(){
    if(WorkspaceQuarantined().load(std::memory_order_acquire))
        throw std::runtime_error("GPU trace workspace cleanup failed; retained owner quarantined; restart this calculation process");
}
void QuarantineWorkspace(std::unique_ptr<Buffers>& owner) noexcept {
    // A lease moved the only saved owner out at construction. Boundary cleanup
    // already owns that saved slot, so neither path can create a second owner.
    auto& saved=SavedWorkspace();
    if(&owner!=&saved)saved=std::move(owner);
    WorkspaceQuarantined().store(true,std::memory_order_release);
}
bool RetireWorkspace(std::unique_ptr<Buffers>& owner) noexcept {
    if(!owner)return true;
    if(!owner->healthy||!owner->ReleaseAllOwned()){QuarantineWorkspace(owner);return false;}
    owner.reset();return true;
}
struct WorkspaceLease {
    std::unique_ptr<Buffers> buffer;size_t limit=0,requested=0;bool success=false,pooled=false;
    explicit WorkspaceLease(int device){
        RequireHealthyWorkspace();
        pooled=Capacity::BufferPoolEnabled();
        requested=pooled?size_t(2ULL*1024*1024*1024):0;
        if(const char* value=std::getenv("MBS_GPU_TRACE_CONVEX_WORKSPACE_MIB")){
            char* end=nullptr;const long mib=std::strtol(value,&end,10);
            if(end==value||*end||mib<0||mib>(pooled?2048:256))
                throw std::runtime_error("MBS_GPU_TRACE_CONVEX_WORKSPACE_MIB must be 0..256 (0..2048 with buffer pool)");
            requested=size_t(mib)*1024*1024;
        }
        // An explicit zero wins even when the new opt-in is enabled.
        if(!requested)pooled=false;
        auto& saved=SavedWorkspace();
        if(saved&&(!requested||saved->device!=device)&&!RetireWorkspace(saved))RequireHealthyWorkspace();
        limit=pooled?PoolLimit(device,saved.get()):requested;
        if(saved){
            saved->Trim(limit,pooled);
            if(!saved->healthy){QuarantineWorkspace(saved);RequireHealthyWorkspace();}
            buffer=std::move(saved);
        }else buffer.reset(new Buffers(device));
        buffer->reused=buffer->allocated=0;buffer->reuseLarger=pooled;
    }
    size_t PoolLimit(int device,const Buffers* owned,size_t* available=nullptr)const{
        size_t free=0,total=0;int current=-1;
        const bool queried=cudaGetDevice(&current)==cudaSuccess&&current==device
            &&cudaMemGetInfo(&free,&total)==cudaSuccess;
        if(available)*available=queried?free:0;
        return Capacity::RetentionLimit(queried,free,owned?owned->Bytes():0,requested);
    }
    bool KeepScratchForConsumer()const{
        if(!pooled||!buffer->healthy)return false;
        size_t free=0;const size_t cap=PoolLimit(buffer->device,buffer.get(),&free);
        return cap&&buffer->Bytes()<=cap&&free>=1024ULL*1024*1024;
    }
    ~WorkspaceLease() noexcept {
        if(!success||!buffer->healthy){RetireWorkspace(buffer);return;}
        if(pooled)limit=PoolLimit(buffer->device,buffer.get());
        if(!limit){RetireWorkspace(buffer);return;}
        buffer->Trim(limit,pooled);
        if(!buffer->healthy){QuarantineWorkspace(buffer);return;}
        size_t free=0,total=0;
        if(cudaMemGetInfo(&free,&total)!=cudaSuccess||free<512ULL*1024*1024){RetireWorkspace(buffer);return;}
        std::fprintf(stderr,"GPU convex workspace: reused_buffers=%d allocations=%d retained_MiB=%.3f limit_MiB=%.3f\n",
            buffer->reused,buffer->allocated,buffer->Bytes()/1048576.,limit/1048576.);
        if(pooled)std::fprintf(stderr,"GPU convex buffer pool retained: physical_MiB=%.3f cap_MiB=%.3f scratch_MiB=%.3f\n",
            buffer->Bytes()/1048576.,limit/1048576.,
            (buffer->sizes[Buffers::Stack]+buffer->sizes[Buffers::Queue0]+buffer->sizes[Buffers::Queue1])/1048576.);
        SavedWorkspace()=std::move(buffer);
    }
};
