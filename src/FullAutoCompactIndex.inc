// Included inside Engine. The journal is an acceleration index, never the
// authority for an optical result. Each referenced complete marker and source
// graph is verified before a restored payload can enter an estimator.
bool compactIndex=false,restoringIndex=false,indexJournalRestored=false;
bool packJobs=false;size_t packedInodesReleased=0;
std::set<std::string> unpackedOversizeJobs;
int indexImportFd=-1;size_t indexImportPendingBytes=0;
size_t indexJournalFlushes=0;
void BeginIndexImport(){
    if(memoryMode)return;
    Require(indexImportFd<0,"nested job index import");
    if(compactIndex){indexImportFd=open((root+"/completed_jobs.journal").c_str(),O_WRONLY|O_CREAT|O_APPEND,0644);
        Require(indexImportFd>=0,"cannot open imported job index");indexImportPendingBytes=0;}
}
void FinishIndexImport(){
    if(indexImportFd<0)return;
    const int fd=indexImportFd;indexImportFd=-1;indexImportPendingBytes=0;
    const bool okay=fsync(fd)==0;close(fd);Require(okay,"cannot flush imported job index");++indexJournalFlushes;
}
void CancelIndexImport(){
    // The index is only an accelerator. Its incomplete final append is repaired
    // on resume; authoritative completed jobs remain untouched.
    if(indexImportFd>=0){close(indexImportFd);indexImportFd=-1;indexImportPendingBytes=0;}
}
// The journal accelerates discovery only: each job's completion marker is
// already durable. Flush an eight-seed publication once, retaining the same
// checksummed records and the existing torn-tail/fallback recovery protocol.
struct IndexAppendBatch {
    Engine& engine;bool owner,done=false;
    explicit IndexAppendBatch(Engine& e):engine(e),owner(e.indexImportFd<0){if(owner)engine.BeginIndexImport();}
    ~IndexAppendBatch(){if(owner&&!done)engine.CancelIndexImport();}
    void Commit(){if(owner)engine.FinishIndexImport();done=true;}
};
std::string SamplingText(const Job& j)const{
    std::ostringstream out;
    out<<signature<<' '<<j.identity<<' '<<j.count<<' '<<j.offset<<' '<<j.seed<<' '<<j.depth<<' '<<j.phi<<' '<<j.block
       <<' '<<j.pipeline<<' '<<j.mirrored<<' '<<j.paired<<' '<<j.assembled<<' '<<j.hasControls<<' '<<j.pairedDepth
       <<' '<<j.warp<<' '<<j.weightHash<<' '<<j.executionHash<<'\n';
    for(double angle:j.inputAngles)out<<Number(angle)<<'\n';
    return out.str();
}
bool ParseSampling(const std::string& text,const std::string& folder,Job& j,bool currentControls=true)const{
    std::istringstream input(text);std::string configuration;
    Require(bool(input>>configuration>>j.identity>>j.count>>j.offset>>j.seed>>j.depth>>j.phi>>j.block
        >>j.pipeline>>j.mirrored>>j.paired>>j.assembled>>j.hasControls>>j.pairedDepth>>j.warp>>j.weightHash>>j.executionHash),"invalid sampling index");
    if(configuration!=signature||(currentControls&&j.hasControls!=controls))return false;
    double angle;while(input>>angle){Require(std::isfinite(angle),"nonfinite indexed angle");j.inputAngles.push_back(angle);}
    Require(input.eof()&&j.inputAngles.size()>=2&&j.count>0&&j.offset>=0,"invalid indexed grid");
    j.folder=folder;j.data=folder+"/result/result.dat";j.summary=folder+"/result/result_analytic_facets.tsv";
    return true;
}
void ConfigureIndex(bool fresh){
    if(memoryMode){compactIndex=true;packJobs=false;return;}
    const char* packed=getenv("MBS_FULLAUTO_PACK_JOBS");
    if(packed&&*packed)Require(std::string(packed)=="0"||std::string(packed)=="1","MBS_FULLAUTO_PACK_JOBS must be 0 or 1");
    const auto packedPolicy=root+"/job_pack_format.txt";
    if(Exists(packedPolicy)){Require(Read(packedPolicy)=="sharded_jobs_v1\n","unsupported packed job format");packJobs=true;}
    if(packed&&*packed)packJobs=std::string(packed)=="1";
    if(packJobs)Atomic(packedPolicy,"sharded_jobs_v1\n");
    Record(root+"/job_storage_policy.json","{\"format\":\"sharded_jobs_v1\",\"packing_new_jobs\":"+std::string(packJobs?"true":"false")
        +",\"maximum_archive_files\":256,\"logical_files_and_hashes_preserved\":true,\"physical_source_removed_after_verified_durable_pack\":true"
        +",\"commit_boundary\":\"two alternating checksummed 4 KiB frontier slots; data fsync before frontier fsync\""
        +",\"maximum_record_bytes\":"+Number(PackedJobStore::MaximumRecordBytes)+",\"oversize_action\":\"retain original job and log checkpoint_pack_oversize_retained\""
        +",\"unpack_tool\":\"fullauto_pack_jobs OUTPUT --unpack\",\"legacy_worker_output_policy_scope\":\"native duplicate text report suppression only\"}");
    const auto path=root+"/job_index_format.txt";
    if(Exists(path)){Require(Read(path)=="compact_v1\n","unsupported job index format");compactIndex=true;}
    else{
        // Existing unscoped campaigns also need the discovery journal. Without
        // this migration every restart rereads millions of small metadata
        // files, even when the controller supports compact indices.
        const auto identity=root+"/configuration.txt";
        if(!fresh&&Exists(identity))Require(Read(identity)==signature+"\n","cannot upgrade an incompatible legacy index");
        Atomic(path,"compact_v1\n");compactIndex=true;
        if(!fresh)Record(root+"/legacy_index_upgrade.json","{\"format\":\"compact_v1\",\"scope\":\"current output metadata only\",\"parent_unchanged\":true,\"sampling_unchanged\":true,\"payload_validation\":\"before use, including dependencies\"}");
    }
}
std::string DiagnosticText(const Job& j)const{
    std::ostringstream d;d<<std::setprecision(17)
        <<j.seconds<<' '<<j.controlRefinement<<' '<<j.controlSetup<<' '<<j.traceSeconds<<' '<<j.diffractionSeconds<<' '
        <<j.processingSeconds<<' '<<j.preparedWaitSeconds<<' '<<j.producerWaitSeconds<<' '<<j.assemblySeconds<<' '<<j.preparationThreads<<'\n';
    return d.str();
}
void ParseDiagnostics(const std::string& text,Job& j)const{
    std::istringstream d(text);std::string extra;
    Require(bool(d>>j.seconds>>j.controlRefinement>>j.controlSetup>>j.traceSeconds>>j.diffractionSeconds
        >>j.processingSeconds>>j.preparedWaitSeconds>>j.producerWaitSeconds>>j.assemblySeconds>>j.preparationThreads)
        &&!(d>>extra)&&j.seconds>0&&std::isfinite(j.seconds),"invalid journal diagnostics");
    for(double value:{j.controlRefinement,j.controlSetup,j.traceSeconds,j.diffractionSeconds,j.processingSeconds,
        j.preparedWaitSeconds,j.producerWaitSeconds,j.assemblySeconds})Require(std::isfinite(value),"nonfinite journal diagnostic");
}
void ReadCompactDiagnostics(const std::string& marker,const std::string& index,Job& j)const{
    const auto start=marker.find('\n',marker.find('\n')+1)+1+index.size();
    Require(start<marker.size(),"missing compact diagnostics");
    const auto end=marker.find('\n',start);Require(end!=std::string::npos,"incomplete compact diagnostics");
    std::istringstream header(marker.substr(start,end-start));std::string magic,hash,extra;
    Require(bool(header>>magic>>hash)&&!(header>>extra)&&magic=="MBS_DIAGNOSTICS_V1","invalid compact diagnostics header");
    const auto diagnostics=marker.substr(end+1);Require(Hash(diagnostics)==hash,"compact diagnostics changed");
    ParseDiagnostics(diagnostics,j);
}
std::string CompletionText(const Job& j,const std::string& metadata)const{
    if(!compactIndex)return metadata;
    const auto index=SamplingText(j),diagnostics=DiagnosticText(j);
    return metadata+"MBS_SAMPLING_V1 "+Hash(index)+" "+Number(index.size())+"\n"+index
        +"MBS_DIAGNOSTICS_V1 "+Hash(diagnostics)+"\n"+diagnostics;
}
void PublishComplete(const Job& j,const std::string& metadata){
    Atomic(j.folder+"/complete.txt",CompletionText(j,metadata));
}
// Oversize is the only packed-publication fallback. Deferred launch/manifests
// and native output must all be durable before a loose marker is published.
// This walks only the single current private job, never prior checkpoints.
void FlushNativePrivateTree(const std::string& path)const{
    struct stat before{},after{},named{};
    Require(lstat(path.c_str(),&before)==0,"cannot inspect private native checkpoint "+path);
    const bool directory=S_ISDIR(before.st_mode);
    Require(directory||S_ISREG(before.st_mode),"unsafe private native checkpoint member "+path);
    PackedJobStore::Descriptor fd(open(path.c_str(),O_RDONLY|O_NOFOLLOW|(directory?O_DIRECTORY:0)));
    Require(fd.fd>=0&&fstat(fd.fd,&after)==0&&PackedJobStore::Same(before,after),"private native checkpoint changed before flush "+path);
    if(directory){
        DIR* dir=fdopendir(dup(fd.fd));Require(dir,"cannot list private native checkpoint "+path);
        std::vector<std::string> names;errno=0;
        while(auto* entry=readdir(dir))if(std::string(entry->d_name)!="."&&std::string(entry->d_name)!="..")names.push_back(entry->d_name);
        const int error=errno;closedir(dir);Require(error==0,"cannot enumerate private native checkpoint "+path);
        std::sort(names.begin(),names.end());for(const auto& name:names)FlushNativePrivateTree(path+"/"+name);
    }
    PackedJobStore::Flush(fd.fd);
    Require(fstat(fd.fd,&after)==0&&lstat(path.c_str(),&named)==0&&PackedJobStore::Same(before,after)&&PackedJobStore::Same(before,named),"private native checkpoint changed during flush "+path);
}
void PublishPreparedComplete(const Job& j,const std::string& metadata,const PackedJobStore::ExpectedMembers& expected){
    if(!packJobs){PublishComplete(j,metadata);return;}
    const auto marker=CompletionText(j,metadata);
    try{packedInodesReleased+=PackedJobStore::PackPrepared(j.folder,marker,expected);return;}
    catch(const PackedJobStore::Oversize& error){
        unpackedOversizeJobs.insert(j.identity);
        Log("checkpoint_pack_oversize_retained","{\"folder\":"+Quoted(j.folder)+",\"message\":"+Quoted(error.what())+",\"original_preserved\":true}");
    }
    auto verifyExpected=[&](){for(const auto& member:expected){
        Require(PackedJobStore::Relative(member.first)&&member.first!="complete.txt","unsafe native expected member");
        Require(FileHash(j.folder+"/"+member.first)==member.second,"native receipt changed before durable fallback");}};
    verifyExpected();FlushNativePrivateTree(j.folder);
    PackedJobStore::FlushDirectory(j.folder.substr(0,j.folder.find_last_of('/')));verifyExpected();
    Require(!Exists(j.folder+"/complete.txt"),"native fallback would overwrite completed job");
    const auto temporary=j.folder+"/complete.txt.tmp";
    PackedJobStore::Descriptor file(open(temporary.c_str(),O_WRONLY|O_CREAT|O_TRUNC|O_NOFOLLOW,0644));
    Require(file.fd>=0,"cannot create native fallback completion");PackedJobStore::Write(file.fd,marker);PackedJobStore::Flush(file.fd);
    Require(rename(temporary.c_str(),(j.folder+"/complete.txt").c_str())==0,"cannot publish native fallback completion");
    PackedJobStore::FlushDirectory(j.folder);
}

// Derived data is private until publication. Preserve the legacy raw batching
// exactly; noncompact jobs still publish loose files before Load writes their
// separate sampling index. Compact packed jobs can commit their private bytes
// directly to the same durable V1 frontier used by native blocks.
struct DerivedFileBatch {
    bool prepared;std::unique_ptr<AtomicFileBatch> files;
    DerivedFileBatch(const std::string& folder,bool prepare,bool legacyBatch=true):prepared(prepare){
        if(prepared||legacyBatch)files.reset(new AtomicFileBatch(folder));
    }
    void Finish(){if(files&&!prepared)files->Commit();files.reset();}
};
bool PrepareDerivedPublication()const{return packJobs&&compactIndex;}
void PublishNativeComplete(const Job& j,const std::string& metadata,const PackedJobStore::ExpectedMembers& expected){
    PublishPreparedComplete(j,metadata,expected);
}
void PublishDerivedComplete(const Job& j,const std::string& metadata,DerivedFileBatch& batch){
    if(!batch.prepared){batch.Finish();PublishComplete(j,metadata);return;}
    Require(PrepareDerivedPublication()&&batch.files&&batch.files.get()==activeAtomicFileBatch,
            "derived publication lost private staging scope");
    PackedJobStore::ExpectedMembers expected;
    auto remember=[&](const std::string& path,const std::string& digest){
        Require(path.compare(0,j.folder.size()+1,j.folder+"/")==0,"derived receipt member outside job");
        const auto name=path.substr(j.folder.size()+1);const auto inserted=expected.emplace(name,digest);
        if(!inserted.second&&inserted.first->second!=digest)
            Require(false,"derived receipt differs from private member "+path
                +" (receipt="+inserted.first->second+", member="+digest+")");};
    // These are the hashes already serialized by the caller, not a fresh
    // receipt after a mutation. Binding them catches changes before Collect.
    std::istringstream receipt(metadata);std::string identity,seconds,dataHash,summaryHash,sourcesHash,integrityHash,extra;
    Require(bool(receipt>>identity>>seconds>>dataHash>>summaryHash>>sourcesHash)&&identity==j.identity,
            "invalid derived completion receipt");
    remember(j.data,dataHash);if(summaryHash!="none")remember(j.summary,summaryHash);
    remember(j.folder+"/angular_sources.tsv",sourcesHash);
    if(receipt>>integrityHash){Require(!(receipt>>extra),"extra derived completion fields");
        remember(j.folder+"/integral_outputs.fnv",integrityHash);}
    // Bind every staged local product, including command/diagnostic text and
    // hard-linked provenance. External source graphs retain their old hashes
    // and their complete dependency validation before use.
    for(const auto& path:batch.files->files)remember(path,FileHash(path));
    auto rememberManifest=[&](const std::string& name){
        const auto found=expected.find(name);Require(found!=expected.end(),"unbound derived integrity manifest");
        const auto text=ReadVerifiedText(j.folder+"/"+name,found->second);std::istringstream rows(text);std::string line;
        while(std::getline(rows,line)){const auto split=line.find_last_of(' ');
            Require(split!=std::string::npos,"invalid derived integrity manifest");
            remember(line.substr(0,split),line.substr(split+1));}
    };
    if(!integrityHash.empty())rememberManifest("integral_outputs.fnv");
    if(expected.count("paired_depth.tsv"))rememberManifest("paired_depth.tsv");
    batch.Finish();PublishPreparedComplete(j,metadata,expected);
}

void AppendIndex(const Job& j){
    if(memoryMode)return;
    if(!compactIndex||restoringIndex)return;
    const auto local=root+"/jobs/",parent=refinementParent+"/jobs/";
    const bool inLocal=j.folder.compare(0,local.size(),local)==0;
    Require(inLocal||(!refinementParent.empty()&&j.folder.compare(0,parent.size(),parent)==0),"job index source lies outside this continuation");
    const auto name=j.folder.substr(inLocal?local.size():parent.size());
    Require(!name.empty()&&name.find('/')==std::string::npos&&name!="."&&name!="..","invalid job index source name");
    std::ostringstream body;body<<name<<'\n'<<(j.indexMarkerHash.empty()?FileHash(j.folder+"/complete.txt"):j.indexMarkerHash)<<'\n'
        <<DiagnosticText(j)<<SamplingText(j);
    const auto data=body.str(),record=Number(data.size())+" "+Hash(data)+"\n"+data+"\n";
    const auto path=root+"/completed_jobs.journal";const bool batched=indexImportFd>=0;
    int fd=batched?indexImportFd:open(path.c_str(),O_WRONLY|O_CREAT|O_APPEND,0644);
    Require(fd>=0,"cannot append job index journal");size_t done=0;
    while(done<record.size()){ssize_t n=write(fd,record.data()+done,record.size()-done);if(n<0&&errno==EINTR)continue;
        if(n<=0){if(!batched)close(fd);Require(false,"job index journal write failed");}done+=size_t(n);}
    if(batched){indexImportPendingBytes+=record.size();
        if(indexImportPendingBytes>=1024*1024){Require(fsync(fd)==0,"imported job index flush failed");++indexJournalFlushes;indexImportPendingBytes=0;}}
    else{const bool okay=fsync(fd)==0;close(fd);Require(okay,"job index journal flush failed");++indexJournalFlushes;}
}
size_t RestoreIndexJournal(){
    if(memoryMode){indexJournalRestored=true;return 0;}
    if(!compactIndex||indexJournalRestored)return 0;
    const auto path=root+"/completed_jobs.journal";std::ifstream in(path.c_str(),std::ios::binary);
    if(!in){Require(!Exists(path),"cannot read job index journal");indexJournalRestored=true;return 0;}
    size_t restored=0;std::streamoff validEnd=0;bool torn=false;std::string header;
    while(std::getline(in,header)){
        if(in.eof()){torn=true;break;}
        size_t length=0;std::string hash,extra;std::istringstream h(header);
        Require(bool(h>>length>>hash)&&!(h>>extra)&&length>0&&length<=64*1024*1024,"invalid job index journal record");
        std::string body(length,'\0');in.read(&body[0],length);char newline=0;
        if(size_t(in.gcount())!=length||!in.get(newline)){torn=true;break;}
        Require(newline=='\n'&&Hash(body)==hash,"job index journal checksum mismatch");
        validEnd=in.tellg();std::istringstream row(body);std::string name,marker,diagnostics;
        Require(bool(std::getline(row,name)&&std::getline(row,marker)&&std::getline(row,diagnostics))
            &&!name.empty()&&name.find('/')==std::string::npos&&name!="."&&name!="..","invalid journal job path");
        std::ostringstream index;index<<row.rdbuf();Job j;
        if(ParseSampling(index.str(),CachedJobFolder(name),j,false)){
            ParseDiagnostics(diagnostics,j);
            j.indexMarkerHash=marker;CacheJob(j);++restored;
            if(restored%1024==0)CheckHostMemory();
            if(restored%16384==0)Record(root+"/cache_restore_status.json",
                "{\"status\":\"reading_journal\",\"controller_pid\":"+Number(getpid())+",\"indexed_jobs\":"+Number(cache.size())+"}");
        }
    }
    Require(!in.bad(),"job index journal read failed");
    in.close();
    // Under the exclusive controller lock, discard only an incomplete final
    // append. Its completed directory remains discoverable in the fallback scan.
    if(torn){int fd=open(path.c_str(),O_WRONLY);Require(fd>=0,"cannot repair torn index tail");
        const bool okay=ftruncate(fd,validEnd)==0&&fsync(fd)==0;close(fd);Require(okay,"cannot truncate torn index tail");
        Log("cache_index_tail_repaired","{\"bytes_preserved\":"+Number(validEnd)+"}");}
    indexJournalRestored=true;return restored;
}
