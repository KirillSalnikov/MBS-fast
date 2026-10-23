// Included inside Engine. Derived views often use identical native source
// lists on many angular stencils. Preserve the ordinary, independently
// readable angular_sources.tsv path, but share its immutable inode. The
// source contents and their directory are durable before any job links them.
std::set<std::string> durableSourceManifests;
size_t sharedSourceLinks=0,sharedSourceBytes=0,sharedSourceCreates=0;
void WriteSharedSources(const std::string& path,const std::string& text){
    // Small lists do not repay the extra shared-directory entry. Native jobs
    // and files outside an unfinished derived transaction use the old path.
    if(text.size()<4096||!activeAtomicFileBatch||!activeAtomicFileBatch->Accepts(path)){
        Atomic(path,text);return;
    }
    const auto directory=root+"/shared_source_manifests";
    const auto shared=directory+"/"+Hash(text)+".tsv";
    const bool alreadyDurable=durableSourceManifests.count(shared);
    if(!alreadyDurable){
        Directory(directory);
        if(Exists(shared))Require(Read(shared)==text,"shared source manifest changed or hash collision");
        else{Atomic(shared,text);++sharedSourceCreates;}
        // A matching file left by an interrupted first publication must cross
        // this fence too. No previous-process memo is trusted after restart.
        AtomicFileBatch::FlushPath(shared,false);
        AtomicFileBatch::FlushPath(directory,true);
        AtomicFileBatch::FlushPath(root,true);
        if(durableSourceManifests.size()>=4096)durableSourceManifests.clear();
        durableSourceManifests.insert(shared);
    }else Require(FileHash(shared)==Hash(text),"shared source manifest changed");
    struct stat source={},target={};
    Require(stat(shared.c_str(),&source)==0,"shared source manifest disappeared");
    if(stat(path.c_str(),&target)==0&&source.st_dev==target.st_dev&&source.st_ino==target.st_ino){
        activeAtomicFileBatch->AddedDurableLink(path);return;
    }
    const auto temporary=path+".link.tmp";
    Require(unlink(temporary.c_str())==0||errno==ENOENT,"cannot remove unfinished source link");
    if(link(shared.c_str(),temporary.c_str())!=0){
        const int error=errno;
        // Filesystems without hard links retain the fully durable old layout.
        Require(error==EXDEV||error==EPERM||error==EOPNOTSUPP||error==EMLINK,
                "cannot link shared source manifest "+path);
        Atomic(path,text);return;
    }
    Require(rename(temporary.c_str(),path.c_str())==0,"cannot publish shared source link");
    activeAtomicFileBatch->AddedDurableLink(path);
    ++sharedSourceLinks;sharedSourceBytes+=text.size();
}
