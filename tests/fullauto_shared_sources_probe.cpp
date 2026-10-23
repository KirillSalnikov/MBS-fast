#define main original_lazy_main
#include "fullauto_lazy_restore_probe.cpp"
#undef main
#include <sys/syscall.h>

static std::string failFlushPath;
extern "C" int fsync(int fd){
    if(!failFlushPath.empty()){
        char path[4096];const auto name="/proc/self/fd/"+std::to_string(fd);
        const auto length=readlink(name.c_str(),path,sizeof(path));
        if(length>=0&&std::string(path,size_t(length))==failFlushPath){errno=EIO;return -1;}
    }
    return int(syscall(SYS_fsync,fd));
}

void SharedSources(const std::string& root){
    ArgPP args;RunConfig config;Engine engine(args,config);engine.root=root;
    Directory(root+"/jobs/leaf/result");Atomic(root+"/jobs/leaf/result/result.dat","leaf\n");
    Atomic(root+"/jobs/leaf/complete.txt","leaf 1 "+Hash("leaf\n")+" none\n");
    std::string sources;for(int i=0;i<100;++i)sources+=root+"/jobs/leaf\t"+FileHash(root+"/jobs/leaf/complete.txt")+"\n";
    Require(sources.size()>4096,"shared source fixture too short");
    auto write=[&](const std::string& folder,bool commit){
        Directory(folder+"/result");AtomicFileBatch files(folder);
        Atomic(folder+"/result/result.dat","derived\n");engine.WriteSharedSources(folder+"/angular_sources.tsv",sources);
        Require(files.durableLinks.count(folder+"/angular_sources.tsv"),"durable manifest link was not registered");
        if(!commit)return;
        files.Commit();Atomic(folder+"/complete.txt","derived 1 "+Hash("derived\n")+" none "+Hash(sources)+"\n");
        std::set<std::string> checked;ValidateStoredJob(folder,checked);
    };
    const auto first=root+"/jobs/first",second=root+"/jobs/second";
    write(first,true);write(second,true);
    struct stat a={},b={};Require(stat((first+"/angular_sources.tsv").c_str(),&a)==0&&stat((second+"/angular_sources.tsv").c_str(),&b)==0,"shared source stat failed");
    Require(a.st_dev==b.st_dev&&a.st_ino==b.st_ino&&engine.sharedSourceCreates==1,"identical source lists were rewritten");
    const auto interrupted=root+"/jobs/interrupted";
    const auto child=fork();Require(child>=0,"cannot fork interruption test");
    if(child==0){write(interrupted,false);_exit(86);}
    int status=0;Require(waitpid(child,&status,0)==child&&WIFEXITED(status)&&WEXITSTATUS(status)==86,"interruption fixture failed");
    Require(!Exists(interrupted+"/complete.txt"),"unfinished shared view is committed");
    engine.durableSourceManifests.clear();write(interrupted,true);
    for(bool sourceFence:{true,false}){
        const auto failed=root+"/jobs/failed"+Number(sourceFence);
        engine.durableSourceManifests.clear();
        failFlushPath=sourceFence?root+"/shared_source_manifests/"+Hash(sources)+".tsv":failed;
        bool failedSafely=false;try{write(failed,true);}catch(const std::exception&){failedSafely=true;}
        failFlushPath.clear();
        Require(failedSafely&&!Exists(failed+"/complete.txt")&&!activeAtomicFileBatch,"failed durable write published completion");
        write(failed,true);
    }
    // Ordinary files remain self-contained even after the shared store entry
    // is removed: existing hard links retain the exact source list.
    const auto shared=root+"/shared_source_manifests/"+Hash(sources)+".tsv";
    Require(unlink(shared.c_str())==0,"cannot unlink shared-store fixture");
    std::set<std::string> checked;ValidateStoredJob(first,checked);
    // In-place corruption through either alias must still invalidate a job.
    {std::ofstream corrupt((second+"/angular_sources.tsv").c_str(),std::ios::app);corrupt<<"changed\n";}
    bool rejected=false;try{std::set<std::string> seen;ValidateStoredJob(first,seen);}catch(const std::exception&){rejected=true;}
    Require(rejected,"modified shared source content was accepted");
    std::cout<<"PASS source manifests: one inode, exact content, interrupted publication recovery, source/directory flush failures, independent paths, corruption rejection\n";
}

void JournalBatch(const std::string& root,const std::string& binary){
    ArgPP args;RunConfig config;Engine engine(args,config);SetupMemory(engine,root,binary,false);
    engine.ConfigureIndex(true);engine.initial=128;engine.orientationBlock=64;
    for(int offset=0;offset<128;offset+=64)engine.DirectStage(64,engine.seeds.data(),8,6,engine.theta,true,false,offset);
    const auto before=Read(root+"/completed_jobs.journal");
    const auto flushes=engine.indexJournalFlushes;engine.forbidNativeWork=true;
    const auto result=engine.PrefixStage(128,engine.seeds.data(),8,6,engine.theta);
    Require(engine.indexJournalFlushes==flushes+1&&engine.indexImportFd<0,"eight-seed index was not committed once");
    const auto matrix=engine.Data(result);
    // Model a crash that loses all eight acceleration records while retaining
    // their independently durable completion markers. No optics may repeat.
    Atomic(root+"/completed_jobs.journal",before);
    Engine resumed(args,config);SetupMemory(resumed,root,binary,false);resumed.ConfigureIndex(false);
    resumed.initial=128;resumed.orientationBlock=64;resumed.forbidNativeWork=true;resumed.WarmCache();
    Require(resumed.cache.size()==engine.cache.size(),"missing journal batch lost completed results");
    const auto restored=resumed.PrefixStage(128,resumed.seeds.data(),8,6,resumed.theta);
    Require(resumed.Data(restored)==matrix,"journal batch recovery changed all16");
    std::cout<<"PASS journal: one flush for eight seeds, missing full batch recovered without native work, all16 exact\n";
}

int main(int argc,char** argv){try{
    Require(argc==3,"OUTPUT BINARY");SharedSources(std::string(argv[1])+"/sources");JournalBatch(std::string(argv[1])+"/journal",argv[2]);return 0;
}catch(const std::exception& error){std::cerr<<error.what()<<'\n';return 1;}}
