#include "../src/CheckedFileCache.h"
#include <fstream>
#include <functional>
#include <iostream>
#include <sys/wait.h>
using namespace PackedJobStore;
static void Put(const std::string& path,const std::string& data){std::ofstream stream(path,std::ios::binary);stream<<data;stream.close();Check(bool(stream),"fixture write");}
static std::string Make(const std::string& root,const std::string& name,const std::string& value){const auto folder=root+"/jobs/"+name;Check(mkdir(folder.c_str(),0755)==0,"fixture directory");Check(mkdir((folder+"/result").c_str(),0755)==0,"fixture result directory");Put(folder+"/complete.txt","original completion\n");Put(folder+"/result/result.dat",value);return folder;}
static bool Fails(const std::function<void()>& call){try{call();return false;}catch(const std::exception&){return true;}}
int main(){try{
    char temporary[]="/tmp/mbs-packed-store-XXXXXX";const char* made=mkdtemp(temporary);Check(made,"temporary directory");const std::string root=made;mkdir((root+"/jobs").c_str(),0755);
    const auto a=Make(root,"jabc","all16 numeric bytes\n");Check(Pack(a)==4,"inode release count");struct stat st{};Check(stat(a.c_str(),&st)!=0&&Exists(a)&&Exists(a+"/result"),"packed directory fallback");
    Check(CheckedFileCache::Read(a+"/result/result.dat")=="all16 numeric bytes\n","lossless member");
    CheckedFileCache cache;const auto hash=cache.Hash(a+"/result/result.dat");const auto opens=cache.opens;
    std::string colliding;for(unsigned k=1;;++k){const auto name="j"+Hex(k);if(name!="jabc"&&Shard(root,name)==Shard(root,"jabc")){colliding=name;break;}}
    const auto b=Make(root,colliding,"second job\n");Pack(b);Check(cache.Hash(a+"/result/result.dat")==hash&&cache.opens==opens,"own append invalidates old hashes");
    // An interrupted cleanup leaves a strict subset of identical files.
    mkdir(a.c_str(),0755);Put(a+"/complete.txt","original completion\n");Check(Pack(a)==2,"partial cleanup recovery");
    // Incomplete final append never loses earlier records. A writer repairs
    // only that tail and republishes a new complete job under the writer lock.
    const auto archive=Shard(root,"jabc");{Descriptor fd(open(archive.c_str(),O_RDWR));struct stat size{};Check(fstat(fd.fd,&size)==0,"fixture archive stat");const auto frontier=ReadFrontier(fd.fd,uint64_t(size.st_size));
        // Full visible length, unreadable new record, and a torn new watermark:
        // the last good watermark still limits the committed prefix.
        Check(lseek(fd.fd,0,SEEK_END)>=0,"fixture append seek");Write(fd.fd,std::string(32768,'X'));
        Check(lseek(fd.fd,off_t((frontier.sequence%2)*FrontierSlotBytes),SEEK_SET)>=0,"fixture frontier seek");Write(fd.fd,Fixed("MBS_JOB_FRONTIER_V1 torn",FrontierSlotBytes));Flush(fd.fd);}
    Check(CheckedFileCache::Read(a+"/result/result.dat")=="all16 numeric bytes\n","torn tail destroyed committed job");
    std::string next;for(unsigned k=1000;;++k){const auto name="j"+Hex(k);if(Shard(root,name)==archive&&name!=colliding){next=name;break;}}
    Pack(Make(root,next,"third job\n"));Check(Jobs(root).size()==3,"torn tail repair lost records");
    // A damaged trailer INSIDE the durable frontier must never be discarded
    // as an unfinished append, even when it is the final committed record.
    uint64_t tail=0;{Descriptor fd(open(archive.c_str(),O_RDWR));struct stat size{};fstat(fd.fd,&size);tail=ReadFrontier(fd.fd,uint64_t(size.st_size)).end-TrailerBytes;const char damaged='X';Check(pwrite(fd.fd,&damaged,1,off_t(tail))==1,"damage committed trailer");Flush(fd.fd);}
    Check(Fails([&]{Jobs(root);}),"committed last-record corruption accepted as a tail");
    {Descriptor fd(open(archive.c_str(),O_RDWR));const char repaired='M';Check(pwrite(fd.fd,&repaired,1,off_t(tail))==1,"restore committed trailer");Flush(fd.fd);}
    // A single damaged frontier after original cleanup must preserve the
    // LATEST job, not expose an older prefix that a writer might truncate.
    {Descriptor fd(open(archive.c_str(),O_RDWR));Check(lseek(fd.fd,0,SEEK_SET)==0,"fixture frontier corruption seek");Write(fd.fd,std::string(FrontierSlotBytes,'Z'));Flush(fd.fd);}
    Check(Jobs(root).size()==3&&CheckedFileCache::Read(root+"/jobs/"+next+"/result/result.dat")=="third job\n","one damaged frontier hid the latest committed job");
    // Discovery works without a completed_jobs.journal, in a fresh process.
    const auto child=fork();Check(child>=0,"fork");if(child==0){Global().indices.Clear();Global().members.Clear();try{Check(Jobs(root).size()==3,"cold discovery");Check(CheckedFileCache::Read(b+"/result/result.dat")=="second job\n","cold read");_exit(0);}catch(...){_exit(1);}}
    int status;waitpid(child,&status,0);Check(WIFEXITED(status)&&WEXITSTATUS(status)==0,"cold reopen");
    Check(Extract(b)==2&&Extract(b)==0,"lossless export/idempotence");Check(CheckedFileCache::Read(b+"/result/result.dat")=="second job\n","export changed bytes");Pack(b);
    const auto alias=root+"_alias";Check(symlink(root.c_str(),alias.c_str())==0,"fixture root alias");
    Check(Shard(alias,"jabc")==Shard(root,"jabc"),"root aliases use different append indices");
    CheckedFileCache aliasCache;const auto aliasHash=aliasCache.Hash(alias+"/jobs/"+colliding+"/result/result.dat");const auto aliasOpens=aliasCache.opens;
    std::string viaAlias;for(unsigned k=2000;;++k){const auto name="j"+Hex(k);if(Shard(root,name)==archive&&name!=next&&name!=colliding){viaAlias=name;break;}}
    Pack(Make(alias,viaAlias,"alias append\n"));Check(aliasCache.Hash(root+"/jobs/"+colliding+"/result/result.dat")==aliasHash&&aliasCache.opens==aliasOpens,"alias append invalidates old hashes");
    Member victim;Check(Locate(a+"/result/result.dat",victim),"locate corruption target");{Descriptor fd(open(archive.c_str(),O_WRONLY));const char bad='X';Check(pwrite(fd.fd,&bad,1,off_t(victim.offset))==1,"damage fixture");Flush(fd.fd);}
    Check(Fails([&]{cache.Hash(a+"/result/result.dat");}),"external payload corruption accepted by cached hash");
    // Refuse cleanup of an incomplete or symlinked source tree.
    const auto incomplete=root+"/jobs/j123";mkdir(incomplete.c_str(),0755);Put(incomplete+"/stdout.log","partial");Check(Fails([&]{Pack(incomplete);}),"incomplete packed");Check(stat((incomplete+"/stdout.log").c_str(),&st)==0,"incomplete source removed");
    const auto linked=Make(root,"j456","value");Check(symlink("/etc/passwd",(linked+"/escape").c_str())==0,"fixture symlink");Check(Fails([&]{Pack(linked);}),"symlink accepted");
    const auto huge=Make(root,"j789","value");{Descriptor file(open((huge+"/large").c_str(),O_WRONLY|O_CREAT,0644));Check(file.fd>=0&&ftruncate(file.fd,off_t(MaximumRecordBytes+1))==0,"oversize fixture");}
    bool oversize=false;try{Pack(huge);}catch(const Oversize&){oversize=true;}Check(oversize&&stat((huge+"/complete.txt").c_str(),&st)==0,"oversize job was not preserved");
    // A bounded number of archive inodes, independent of completed job count.
    for(unsigned k=20000;k<21024;++k)Pack(Make(root,"j"+Hex(k),"matrix bytes\n"));
    size_t files=0;DIR* dir=opendir((root+"/job_packs").c_str());while(auto* item=readdir(dir))if(item->d_name[0]!='.')++files;closedir(dir);Check(files<=256,"unbounded archive inode count");
    Check(Global().indices.Bytes()<=Global().indices.Limit()&&Global().members.Bytes()<=Global().members.Limit(),"unbounded metadata cache");
    std::vector<CompactIndex> scale(256);for(uint64_t k=0;k<2500000;++k){const auto name="j"+Hex(k);Check(scale[Digest(name)&255].Insert(name,k*1000),"large index insert");}
    size_t indexBytes=0;for(const auto& shard:scale)indexBytes+=shard.Bytes();Check(indexBytes<80*1024*1024,"2.5m job index exceeds bounded resident budget");
    for(uint64_t k=0;k<2500000;k+=97){const auto name="j"+Hex(k);uint64_t offset=0;Check(scale[Digest(name)&255].Find(name,offset)&&offset==k*1000,"large index lookup");}
    std::string savedFrontiers;{Descriptor fd(open(archive.c_str(),O_RDWR));savedFrontiers=Bytes(fd.fd,0,DataBegin);Check(lseek(fd.fd,0,SEEK_SET)==0,"fixture seek");Write(fd.fd,std::string(DataBegin,'Z'));Flush(fd.fd);}
    Check(Fails([&]{Jobs(root);}),"two corrupted frontier slots accepted");
    {Descriptor fd(open(archive.c_str(),O_RDWR));Check(lseek(fd.fd,0,SEEK_SET)==0,"fixture repair seek");Write(fd.fd,savedFrontiers);Flush(fd.fd);}
    std::cout<<"PASS 2500000 job index: "<<indexBytes<<" bytes, all 256 shards fit in 128 MiB cache\n";
    std::cout<<"PASS packed store: lossless bytes, <=256 archive inodes for 1027 completed jobs, cold discovery, interrupted cleanup, torn append, corruption rejection, append-stable hash cache. Fixture "<<root<<'\n';return 0;
}catch(const std::exception& error){std::cerr<<error.what()<<'\n';return 1;}}
