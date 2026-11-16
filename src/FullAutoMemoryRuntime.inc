// Included inside Engine. Working values and optional durable snapshots have
// separate lifetimes. No numerical selection rule depends on this policy.
struct PhysicalStorageScope {
    FullAutoMemoryStore::Store* previous;
    PhysicalStorageScope():previous(FullAutoMemoryStore::Active()){FullAutoMemoryStore::Active()=nullptr;}
    ~PhysicalStorageScope(){FullAutoMemoryStore::Active()=previous;}
};
void ConfigureMemoryStorage(){
    const bool marked=Exists(root+"/memory_storage_v1.txt");
    const bool legacy=Exists(root+"/configuration.txt")&&!marked;
    const auto policy=args.IsCatched("fullauto_storage")?args.GetStringValue("fullauto_storage"):std::string();
    memoryMode=policy=="memory"||(policy.empty()&&!legacy&&!refinementWindow);
    Require(!(marked&&!memoryMode),"a memory calculation cannot resume through legacy storage");
    Require(!(memoryMode&&legacy),"memory storage requires a new output directory; existing legacy checkpoints are preserved");
    Require(!(memoryMode&&refinementWindow),"memory refinement import is not supported; use legacy storage for an existing refinement window");
    if(!memoryMode){Require(!args.IsCatched("fullauto_checkpoint_seconds")&&!args.IsCatched("fullauto_memory_mib"),"memory checkpoint/budget options cannot apply to an existing legacy output");return;}
    WorkerMemoryStage::CleanupDeadWorkers(root);
    CleanupMemoryAbandonedPublications();
    checkpointInterval=Int("fullauto_checkpoint_seconds",args.IsCatched("checkpoint")?300:0);
    const int storageShares=config.useGpu?std::max(1,VisibleGpuDeviceCount()):std::max(1,int(gpus.size()));
    const int defaultMiB=FullAutoStorageBudget::DefaultMiB(QueryRuntimeMemorySnapshot().ramAvailableKb,storageShares);
    const int mib=Int("fullauto_memory_mib",defaultMiB);
    Require(mib>0&&uint64_t(mib)<=uint64_t(SIZE_MAX)/(1024*1024),"invalid memory storage budget");
    const size_t storeBytes=size_t(mib)*1024*1024;
    // Move a large path catalogue out before it crowds hot matrix payloads out.
    // This is a split of the existing limit, not additional RAM.
    memoryStore.reset(new FullAutoMemoryStore::Store(root,storeBytes,root,storeBytes>=16*1024*1024?storeBytes/4:0));
    std::vector<std::pair<std::string,std::string>> ownedInputs;
    if(!marked){for(const auto& arg:base)if(memoryStore->Owns(arg)){
        struct stat st{};if(stat(arg.c_str(),&st)==0&&S_ISREG(st.st_mode)){
            ownedInputs.push_back({arg,Read(arg)});
            // Particle readers implicitly consult this optional sidecar.
            // Moving only the explicit CLI file into tmpfs would lose it.
            const auto sidecar=arg+".visibility";
            if(stat(sidecar.c_str(),&st)==0&&S_ISREG(st.st_mode))ownedInputs.push_back({sidecar,Read(sidecar)});
        }
    }}
    if(marked){
        Require(Exists(root+"/fullauto_checkpoint.bin"),"this memory run has no recovery checkpoint; use a new output directory (--checkpoint enables future recovery)");
        Require(memoryStore->Restore(root+"/fullauto_checkpoint.bin")==signature+"\n","memory checkpoint inputs/binary changed");
    }
    Atomic(root+"/memory_storage_v1.txt","memory_values_v1\n");
    memoryScope.reset(new FullAutoMemoryStore::Scope(*memoryStore));
    memoryStore->Directory(root);memoryStore->Directory(root+"/jobs");
    for(const auto& input:ownedInputs)memoryStore->Write(input.first,input.second);
    memoryCheckpointAt=Seconds();
    Record(root+"/fullauto_storage.json",MemoryStorageReport());
}
std::string MemoryStorageReport()const{
    if(!memoryMode)return "{\"storage\":\"legacy\"}";
    const auto s=memoryStore->Stats();
    return "{\"storage\":\"memory\",\"schema\":1,\"checkpoint_interval_seconds\":"+Number(checkpointInterval)
        +",\"resume_enabled\":"+(checkpointInterval>0?"true":"false")+",\"memory_limit_bytes\":"+Number(s.limitBytes)
        +",\"metadata_catalog_trigger_bytes\":"+Number(s.metadataLimitBytes)
        +",\"resident_bytes\":"+Number(s.residentBytes)+",\"metadata_bytes\":"+Number(s.metadataBytes)
        +",\"spill_bytes\":"+Number(s.spillBytes)+",\"spill_read_bytes\":"+Number(s.spillReadBytes)
        +",\"snapshot_spill_read_bytes\":"+Number(s.snapshotReadBytes)
        +",\"working_spill_read_bytes\":"+Number(s.spillReadBytes-s.snapshotReadBytes)
        +",\"spill_write_bytes\":"+Number(s.spillWriteBytes)+",\"virtual_files\":"+Number(s.files)
        +",\"temporary_catalog_active\":"+(s.catalogActive?"true":"false")+",\"temporary_catalog_bytes\":"+Number(s.catalogBytes)
        +",\"temporary_catalog_cache_bytes\":"+Number(s.catalogCacheBytes)
        +",\"checkpoint_count_this_invocation\":"+Number(memoryCheckpoints)
        +",\"abandoned_temporary_files_removed\":"+Number(memoryAbandonedTemps)
        +",\"worker_transport\":\"streamed_socket_files\",\"transport_chunk_bytes\":65536"
        +",\"worker_staging\":\"tmpfs for small jobs; disposable disk files for large jobs or low RAM\""
        +",\"durable_native_blocks\":false,\"snapshot\":"+Quoted(root+"/fullauto_checkpoint.bin")+"}";
}
bool MemoryPublicFile(const std::string& name)const{
    if(name=="adaptive_checkpoint.tsv"||name=="go_checkpoint.tsv"||name=="completed_jobs.journal"
        ||name=="job_index_format.txt"||name=="job_pack_format.txt"||name=="fullauto_checkpoint.bin")return false;
    if(name.compare(0,6,"theta_")==0&&name.size()>10&&name.compare(name.size()-4,4,".csv")==0){
        const auto digest=name.substr(6,name.size()-10);
        if(digest.find_first_not_of("0123456789abcdef")==std::string::npos)return false;
    }
    for(const char* prefix:{"weights_","mean_","means_","angular_sources_"})
        if(name.compare(0,std::strlen(prefix),prefix)==0)return false;
    return name.size()<4||name.compare(name.size()-4,4,".tmp")!=0;
}
void ExportMemoryObservation(const std::string& path,const std::string& bytes){
    // A live observation is replaceable, not a recovery transaction. Publish
    // one small, complete JSON without forcing SSD/metadata durability every
    // two seconds. Recovery still depends only on the explicit snapshot.
    Require(path==root+"/execution_progress.json","invalid memory observation path");
    const auto temporary=path+".observation.tmp";
    std::ofstream file(temporary,std::ios::binary);Require(bool(file),"cannot write live execution progress");
    file<<bytes;file.close();Require(bool(file),"cannot finish live execution progress");
    Require(rename(temporary.c_str(),path.c_str())==0,"cannot publish live execution progress");
}
std::vector<std::string> MemorySavedGeometryPaths()const{
    if(!args.IsCatched("save_geometry"))return {};
    const auto requested=args.GetStringValue("save_geometry");
    if(!memoryStore->Owns(requested))return {};
    const auto path=FullAutoMemoryStore::Canonical(requested);
    return {path,path+".visibility"};
}
#ifndef MBS_MEMORY_PUBLICATION_EVENT
#define MBS_MEMORY_PUBLICATION_EVENT(event,path) ((void)0)
#endif
std::string MemorySnapshotFolder(const std::string& json)const{
    const std::string key="\"snapshot_folder\":";size_t i=json.find(key);if(i==std::string::npos)return "";
    i=json.find('"',i+key.size());if(i==std::string::npos)return "";std::string result;
    for(++i;i<json.size();++i){char c=json[i];if(c=='"')return result;
        if(c=='\\'){if(++i>=json.size())break;c=json[i];
            if(c=='n')c='\n';else if(c=='r')c='\r';else if(c=='t')c='\t';
            else if(c!='\\'&&c!='"'&&c!='/')return "";}
        result+=c;
    }return "";
}
void CleanupMemoryAbandonedPublications(){
    // ConfigureMemoryStorage runs after the exclusive controller lock. These
    // temporary checkpoint files cannot be sources for recovery.
    const std::string prefix="fullauto_checkpoint.bin.tmp.";
    if(DIR* directory=opendir(root.c_str())){
        while(auto* entry=readdir(directory)){const std::string name=entry->d_name;
            if(name.size()!=prefix.size()+6||name.compare(0,prefix.size(),prefix)!=0
                ||name.substr(prefix.size()).find_first_not_of("0123456789abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ")!=std::string::npos)continue;
            const auto path=root+"/"+name;struct stat st{};
            if(lstat(path.c_str(),&st)==0&&S_ISREG(st.st_mode)&&st.st_uid==geteuid()&&(st.st_mode&0777)==0600&&st.st_nlink==1)
                if(unlink(path.c_str())==0)++memoryAbandonedTemps;
        }closedir(directory);
    }
    const auto parent=root+"/public_snapshots";DIR* directory=opendir(parent.c_str());if(!directory)return;
    while(auto* entry=readdir(directory)){const std::string name=entry->d_name;
        if(name.compare(0,9,".pending-")!=0)continue;const auto dash=name.find('-',9);
        if(dash==std::string::npos||dash==9||name.size()!=dash+7)continue;
        const auto digits=name.substr(9,dash-9),suffix=name.substr(dash+1);
        if(digits.find_first_not_of("0123456789")!=std::string::npos
            ||suffix.find_first_not_of("0123456789abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ")!=std::string::npos)continue;
        long long number=0;try{number=std::stoll(digits);}catch(...){continue;}
        if(number<=0||number>std::numeric_limits<pid_t>::max())continue;const auto pid=pid_t(number);
        errno=0;if(kill(pid,0)==0||errno!=ESRCH)continue;
        const auto path=parent+"/"+name;struct stat st{};
        if(lstat(path.c_str(),&st)!=0||!S_ISDIR(st.st_mode)||st.st_uid!=geteuid()||(st.st_mode&0777)!=0700)continue;
        std::ifstream marker(path+"/.mbs-public-pending-v1");if(!marker)continue;
        const std::string bytes((std::istreambuf_iterator<char>(marker)),{});
        if(bytes!="mbs_public_pending_v1\n"+root+"\n"+Number(pid)+"\n")continue;
        errno=0;if(kill(pid,0)==0||errno!=ESRCH)continue;
        WorkerMemoryStage::RemoveTree(path);if(lstat(path.c_str(),&st)!=0&&errno==ENOENT)++memoryAbandonedTemps;
    }closedir(directory);
}
std::string PublishMemoryMatrixGeneration(const std::string& source){
    std::vector<std::string> files;std::ostringstream manifest;
    manifest<<"mbs_public_matrix_v1\n";
    for(const auto& path:memoryStore->List(source,true))if(!memoryStore->IsDirectory(path)){
        files.push_back(path);manifest<<path.substr(source.size()+1)<<'\t'<<FileHash(path)<<'\n';}
    Require(!files.empty(),"empty public matrix generation");
    const auto parent=root+"/public_snapshots",destination=parent+"/"+Hash(manifest.str());
    if(memoryPublishedGenerations.count(destination))return destination;
    bool present=false;{PhysicalStorageScope disk;Directory(parent);present=Exists(destination);
        if(present){Require(Read(destination+"/.mbs-public-v1")==manifest.str(),"public matrix generation manifest changed");
            for(const auto& path:files){const auto published=destination+path.substr(source.size());
                Require(FileHash(published)==Hash(memoryStore->Read(path)),"public matrix generation payload changed");}}
    }
    if(!present){
        const auto pattern=parent+"/.pending-"+Number(getpid())+"-XXXXXX";
        std::vector<char> name(pattern.begin(),pattern.end());name.push_back('\0');
        char* created=mkdtemp(name.data());Require(created,"cannot stage public matrix generation");const std::string stage=created;
        try{
            {PhysicalStorageScope disk;Atomic(stage+"/.mbs-public-pending-v1","mbs_public_pending_v1\n"+root+"\n"+Number(getpid())+"\n");
                AtomicFileBatch::FlushPath(stage,true);}
            for(const auto& path:files){const auto bytes=Read(path),target=stage+path.substr(source.size());
                {PhysicalStorageScope disk;Directory(target.substr(0,target.find_last_of('/')));Atomic(target,bytes);}
                MBS_MEMORY_PUBLICATION_EVENT("member",target);
            }
            {PhysicalStorageScope disk;Atomic(stage+"/.mbs-public-v1",manifest.str());
                AtomicFileBatch::FlushPath(stage,true);
                Require(rename(stage.c_str(),destination.c_str())==0,"cannot publish matrix generation");
                AtomicFileBatch::FlushPath(parent,true);}
        }catch(...){WorkerMemoryStage::RemoveTree(stage);throw;}
    }
    memoryPublishedGenerations.insert(destination);return destination;
}
void PruneMemoryMatrixGenerations(const std::string& current,const std::string& previous){
    if(current.empty()||current==previous)return;
    const auto parent=root+"/public_snapshots";DIR* directory=opendir(parent.c_str());if(!directory)return;
    std::vector<std::string> obsolete;
    while(auto* entry=readdir(directory)){const std::string name=entry->d_name;
        if(name.empty()||name.size()>16||name.find_first_not_of("0123456789abcdef")!=std::string::npos)continue;
        const auto path=parent+"/"+name;if(path==current||path==previous)continue;struct stat st{};
        if(lstat(path.c_str(),&st)!=0||!S_ISDIR(st.st_mode)||st.st_uid!=geteuid()||(st.st_mode&0777)!=0700)continue;
        std::ifstream marker(path+"/.mbs-public-v1");if(!marker)continue;
        const std::string bytes((std::istreambuf_iterator<char>(marker)),{});
        if(bytes.compare(0,21,"mbs_public_matrix_v1\n")!=0||Hash(bytes)!=name)continue;
        obsolete.push_back(path);
    }closedir(directory);
    for(const auto& path:obsolete){WorkerMemoryStage::RemoveTree(path);memoryPublishedGenerations.erase(path);}
}
void ExportMemoryPublic(bool force=false){
    if(!memoryMode)return;
    const double now=Seconds();if(!force&&memoryPublicAt>=0&&now-memoryPublicAt<10)return;memoryPublicAt=now;
    memoryStore->Write(root+"/fullauto_storage.json",MemoryStorageReport());
    std::vector<std::pair<std::string,std::string>> exports;
    std::string intermediate,sourceFolder,publishedFolder,previousFolder;
    if(Exists(root+"/intermediate_current.json")){
        intermediate=Read(root+"/intermediate_current.json");
        sourceFolder=MemorySnapshotFolder(intermediate);
        if(publicRoot.empty()&&!sourceFolder.empty()&&memoryStore->Owns(sourceFolder)&&Exists(sourceFolder)){
            std::ifstream old(root+"/intermediate_current.json");
            previousFolder=MemorySnapshotFolder(std::string((std::istreambuf_iterator<char>(old)),{}));
            publishedFolder=PublishMemoryMatrixGeneration(sourceFolder);
            const auto from=Quoted(sourceFolder).substr(1,Quoted(sourceFolder).size()-2),to=Quoted(publishedFolder).substr(1,Quoted(publishedFolder).size()-2);
            for(size_t pos=0;(pos=intermediate.find(from,pos))!=std::string::npos;){intermediate.replace(pos,from.size(),to);pos+=to.size();}
        }
        exports.push_back({root+"/intermediate_current.json",root+"/intermediate_current.json"});
    }
    for(const auto& entry:memoryStore->ChildFiles(root)){
        const auto& path=entry.first;
        const auto name=path.substr(root.size()+1);if(name=="intermediate_current.json"||!MemoryPublicFile(name))continue;
        const auto old=memoryExportedGenerations.find(path);
        if(old!=memoryExportedGenerations.end()&&old->second==entry.second)continue;
        exports.push_back({path,path});
    }
    // Explicit geometry exports may be nested below the output directory.
    // Only these exact user-selected paths are allowed outside jobs/caches.
    for(const auto& path:MemorySavedGeometryPaths())if(memoryStore->Exists(path))exports.push_back({path,path});
    // A small current result is deliberately kept even without recovery.
    // Internal jobs/grids/source manifests never appear as physical files.
    for(const auto& file:exports){
        const uint64_t generation=memoryStore->Generation(file.second);
        const auto old=memoryExportedGenerations.find(file.first);
        if(old!=memoryExportedGenerations.end()&&old->second==generation
            &&(file.first!=root+"/intermediate_current.json"||publishedFolder.empty()||publishedFolder==previousFolder))continue;
        const auto bytes=file.first==root+"/intermediate_current.json"?intermediate:Read(file.second);
        {PhysicalStorageScope disk;Directory(file.first.substr(0,file.first.find_last_of('/')));Atomic(file.first,bytes);
            if(file.first==root+"/intermediate_current.json")AtomicFileBatch::FlushPath(root,true);}
        memoryExportedGenerations[file.first]=generation;
    }
    // The pointer is committed only after every member of its immutable
    // generation is durable. Keep the preceding generation for concurrent
    // readers; never overwrite a directory an old pointer still describes.
    PruneMemoryMatrixGenerations(publishedFolder,previousFolder);
}
void MaybeMemoryCheckpoint(bool force=false){
    if(!memoryMode||checkpointInterval<=0)return;
    if(!force&&Seconds()-memoryCheckpointAt<checkpointInterval)return;
    // Called only on the controller thread between complete transactions.
    // A dispatched worker may still be running; its unfinished request has no
    // completion marker and will be recomputed after a restart.
    if(activeAtomicFileBatch)return;
    const double checkpointStarted=Seconds();
    memoryStore->Snapshot(root+"/fullauto_checkpoint.bin",signature+"\n");
    memoryCheckpointAt=Seconds();++memoryCheckpoints;
    Log("memory_checkpoint_timing","{\"seconds\":"+Number(memoryCheckpointAt-checkpointStarted)+",\"generation\":"+Number(memoryStore->Epoch())+"}");
    Log("memory_checkpoint","{\"file\":"+Quoted(root+"/fullauto_checkpoint.bin")+",\"generation\":"+Number(memoryStore->Epoch())+"}");
    ExportMemoryPublic(true);
}
void FinishMemoryRun(){MaybeMemoryCheckpoint(true);ExportMemoryPublic(true);ExportCompactResults(true);}
std::vector<std::string> MemoryWorkerInputPaths(const std::vector<std::string>& command)const{
    std::vector<std::string> files;std::set<std::string> found;
    for(const auto& argument:command){
        for(const auto& input:{argument,argument+".visibility"}){
            if(!memoryStore->Owns(input)||!memoryStore->Exists(input)||memoryStore->IsDirectory(input))continue;
            const auto path=FullAutoMemoryStore::Canonical(input);
            if(!found.insert(path).second)continue;
            files.push_back(path);
        }
    }
    return files;
}
// Kept for small protocol/invariant fixtures. Production uses descriptors and
// bounded streaming, including for an expanding shared analytic cache.
std::vector<WorkerProtocol::File> MemoryWorkerInputs(const std::vector<std::string>& command)const{
    std::vector<WorkerProtocol::File> files;
    for(const auto& path:MemoryWorkerInputPaths(command))files.push_back({path.substr(root.size()+1),memoryStore->Read(path)});
    return files;
}
std::vector<WorkerProtocol::StreamFile> MemoryWorkerStreamInputs(const std::vector<std::string>& command)const{
    std::vector<WorkerProtocol::StreamFile> files;
    for(const auto& path:MemoryWorkerInputPaths(command))
        files.emplace_back(path.substr(root.size()+1),memoryStore->Size(path),memoryStore->Digest(path));
    return files;
}
uint64_t MemoryWorkerOutputBudget(const Job& job)const{
    // Conservative allowance for upper/lower and half-prefix DATs, control
    // summaries, scalar diagnostics and logs. It selects temporary staging,
    // never the numerical theta/orientation grid. Large jobs use ordinary
    // reclaimable disk cache instead of consuming unbounded tmpfs RAM.
    const uint64_t fixed=8ULL*1024*1024,bytesPerRow=6144;
    Require(job.inputAngles.size()<=(WorkerProtocol::MaximumStreamBytes-fixed)/bytesPerRow,"native angular output exceeds the 8 GiB streaming safety limit");
    return fixed+uint64_t(job.inputAngles.size())*bytesPerRow;
}
void AcceptMemoryWorkerStream(const WorkerProtocol::StreamFile& file,const WorkerProtocol::StreamReader& reader,const Job& job){
    const auto path=root+"/"+file.path;
    const bool native=path.compare(0,job.folder.size()+1,job.folder+"/")==0;
    const bool mean=(path==MeanCache("shared")||file.path.compare(0,5,"mean_")==0||file.path.compare(0,6,"means_")==0);
    const auto geometry=MemorySavedGeometryPaths();
    Require(native||mean||std::find(geometry.begin(),geometry.end(),path)!=geometry.end(),"worker returned a stream outside its result/cache namespace");
    RequireMutableMemoryPath(path);
    memoryStore->WriteStream(path,file.size,file.digest,reader);
}
void AcceptMemoryWorkerFiles(const std::vector<WorkerProtocol::File>& files,const Job& job){
    for(const auto& file:files){
        const auto path=root+"/"+file.path;
        const bool native=path.compare(0,job.folder.size()+1,job.folder+"/")==0;
        const bool mean=(path==MeanCache("shared")||file.path.compare(0,5,"mean_")==0||file.path.compare(0,6,"means_")==0);
        const auto geometry=MemorySavedGeometryPaths();
        Require(native||mean||std::find(geometry.begin(),geometry.end(),path)!=geometry.end(),"worker returned a file outside its result/cache namespace");
        Atomic(path,file.bytes);
    }
}
