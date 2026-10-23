// Included inside Engine. Derived views often use identical native source
// lists on many angular stencils. Preserve the ordinary, independently
// readable angular_sources.tsv path, but share its immutable inode only when
// the list recurs. A unique list keeps the ordinary layout and write cost.
// Only completed local jobs can supply a link: their file and directory have
// already crossed the derived transaction's durability fence.
std::map<std::string,std::string> sharedSourcePaths;
size_t sharedSourceLinks=0,sharedSourceBytes=0,sharedSourceCreates=0;
void WriteSharedSources(const std::string& path,const std::string& text){
    // Small lists do not repay the extra lookup. Native jobs
    // and files outside an unfinished derived transaction use the old path.
    if(text.size()<4096||!activeAtomicFileBatch||!activeAtomicFileBatch->Accepts(path)){
        Atomic(path,text);return;
    }
    const auto key=Hash(text);const auto found=sharedSourcePaths.find(key);
    const std::string shared=found==sharedSourcePaths.end()?"":found->second;
    if(shared.empty()||shared==path||!Exists(shared)
            ||!Exists(shared.substr(0,shared.find_last_of('/'))+"/complete.txt")){
        // An interrupted first publication is not a durable source. Retry it
        // through the ordinary batch, or use this new job as the first copy.
        Atomic(path,text);++sharedSourceCreates;
        if(sharedSourcePaths.size()>=4096)sharedSourcePaths.clear();
        sharedSourcePaths[key]=path;return;
    }
    Require(Read(shared)==text,"shared source manifest changed or hash collision");
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
