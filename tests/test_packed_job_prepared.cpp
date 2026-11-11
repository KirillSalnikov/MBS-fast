#include <sys/stat.h>
#include <unistd.h>
#include <fcntl.h>
#include <sys/wait.h>
#include <cerrno>
#include <functional>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>
static void Event(const char*,const std::string&);
static int FaultFsync(int);
#define MBS_PACK_EVENT(event,folder) Event(event,folder)
#define fsync FaultFsync
#include "../src/PackedJobStore.h"
#undef fsync
#undef MBS_PACK_EVENT
using namespace PackedJobStore;
static std::string mode,point,stage;static bool checkAbsent=false,failed=false;static int flushCount=0,events=0;
static void Put(const std::string& p,const std::string& text){std::ofstream f(p,std::ios::binary);f<<text;f.close();Check(bool(f),"fixture write failed");}
static bool Raw(const std::string& p){struct stat st;return lstat(p.c_str(),&st)==0;}
static void Reset(){mode.clear();point.clear();stage.clear();failed=false;checkAbsent=false;flushCount=0;events=0;}
static int FaultFsync(int fd){++flushCount;if(mode=="eio"&&!failed&&stage==point){failed=true;errno=EIO;return -1;}return ::fsync(fd);}
static void Event(const char* label,const std::string& folder){stage=label;++events;if(checkAbsent)Check(!Raw(folder+"/complete.txt"),"ephemeral loose complete appeared");
 if(stage==point){if(mode=="crash")_exit(73);if(mode=="change"&&!failed){failed=true;Put(folder+"/result/result.dat","changed after collection\n");}
 if(mode=="corrupt"&&!failed){failed=true;std::string root,name,part;Split(folder,root,name,part);const auto archive=Shard(root,name);Descriptor fd(open(archive.c_str(),O_RDWR));struct stat st;fstat(fd.fd,&st);const char bad='X';Check(pwrite(fd.fd,&bad,1,st.st_size-TrailerBytes-1)==1,"inject record corruption");}}
}
static bool Reject(const std::function<void()>& f){try{f();return false;}catch(const std::exception&){return true;}}
static std::string Fresh(){char name[]="/tmp/mbs-virtual-pack-XXXXXX";const auto p=mkdtemp(name);Check(p,"mkdtemp");mkdir((std::string(p)+"/jobs").c_str(),0755);return p;}
static std::string Fixture(const std::string& root,const std::string& name){const auto p=root+"/jobs/"+name;mkdir(p.c_str(),0755);mkdir((p+"/result").c_str(),0755);Put(p+"/command.args","immutable recipe\n");Put(p+"/result/result.dat","all16 matrix bytes\n");Put(p+"/result/result_analytic_facets.tsv","control bytes\n");Put(p+"/integral_outputs.fnv","integral provenance\n");return p;}
static ExpectedMembers Expected(){return {{"command.args",Hex(Digest("immutable recipe\n"))},{"result/result.dat",Hex(Digest("all16 matrix bytes\n"))},{"result/result_analytic_facets.tsv",Hex(Digest("control bytes\n"))},{"integral_outputs.fnv",Hex(Digest("integral provenance\n"))}};}
static size_t Prepared(const std::string& folder,const std::string& marker){return PackPrepared(folder,marker,Expected());}
static std::string Collide(const std::string& root,const std::string& other){for(unsigned i=1;;++i){const auto n="j"+Hex(i);if(n!=other&&Shard(root,n)==Shard(root,other))return n;}}
static std::string Get(const std::string& path){Member m;Check(Locate(path,m),"member missing");return Read(m);}
static void V1Read(const std::string& folder,const std::string& marker){PackedJobStore::Member m;Check(PackedJobStore::Locate(folder+"/complete.txt",m),"V1 cannot find completion");Check(PackedJobStore::Read(m)==marker,"V1 marker changed");Check(PackedJobStore::Locate(folder+"/result/result.dat",m)&&PackedJobStore::Read(m)=="all16 matrix bytes\n","V1 payload changed");}
static std::pair<std::string,std::string> FixturePair(){const auto root=Fresh();const auto old=Fixture(root,"jabc");Put(old+"/complete.txt","prior complete\n");Pack(old);return {old,Fixture(root,Collide(root,"jabc"))};}
int main(){try{
 const std::string marker="new compact completion and sampling\n";int cases=0;
 {Reset();const auto pair=FixturePair();checkAbsent=true;flushCount=0;const auto removed=Prepared(pair.second,marker);Check(removed==6,"inode count includes virtual nonexistent file");Check(flushCount==4,"existing shard must use exactly data/frontier/mirror+cleanup flush");V1Read(pair.second,marker);Check(!Raw(pair.second)&&Prepared(pair.second,marker)==0,"complete retry not idempotent");Check(Reject([&]{Prepared(pair.second,"different marker\n");}),"different prepared completion accepted");V1Read(pair.first,"prior complete\n");++cases;}
 // Every crash boundary retains the prior checkpoint; restart either ignores
 // an unpublished tail or resumes cleanup from a verified committed record.
 for(const auto& at:{"after_collect","before_record_write","before_data_fsync","after_data_fsync","before_frontier_write","before_frontier_fsync","after_frontier_fsync","before_mirror_fsync","after_mirror_fsync","before_cleanup","after_cleanup_member"}){
  Reset();const auto pair=FixturePair();const pid_t child=fork();Check(child>=0,"fork");if(!child){checkAbsent=true;mode="crash";point=at;Prepared(pair.second,marker);_exit(99);}int status=0;waitpid(child,&status,0);Check(WIFEXITED(status)&&WEXITSTATUS(status)==73,"fault point not reached");
  Reset();Check(!Raw(pair.second+"/complete.txt"),"crash left loose complete");V1Read(pair.first,"prior complete\n");Prepared(pair.second,marker);V1Read(pair.second,marker);Check(!Raw(pair.second),"retry left loose source");++cases;
 }
 for(const auto& at:{"before_data_fsync","before_frontier_fsync","before_mirror_fsync"}){
  Reset();const auto pair=FixturePair();checkAbsent=true;mode="eio";point=at;Check(Reject([&]{Prepared(pair.second,marker);})&&failed,"fsync EIO hidden");Check(Raw(pair.second+"/result/result.dat")&&!Raw(pair.second+"/complete.txt"),"fsync error removed private source or exposed raw complete");Reset();V1Read(pair.first,"prior complete\n");Prepared(pair.second,marker);V1Read(pair.second,marker);++cases;
 }
 for(const auto& at:{"after_collect","before_frontier_write"}){Reset();const auto pair=FixturePair();checkAbsent=true;mode="change";point=at;Check(Reject([&]{Prepared(pair.second,marker);}),"source mutation accepted");Check(Raw(pair.second)&&!Exists(pair.second+"/complete.txt"),"source mutation published completion");V1Read(pair.first,"prior complete\n");++cases;}
 {Reset();const auto pair=FixturePair();mode="corrupt";point="after_data_fsync";checkAbsent=true;Check(Reject([&]{Prepared(pair.second,marker);}),"readback corruption accepted");Check(Raw(pair.second)&&!Exists(pair.second+"/complete.txt"),"corrupt record published");V1Read(pair.first,"prior complete\n");++cases;}
 {Reset();const auto root=Fresh();const auto folder=Fixture(root,"j12");Descriptor fd(open((folder+"/large").c_str(),O_CREAT|O_WRONLY,0644));Check(fd.fd>=0&&ftruncate(fd.fd,off_t(MaximumRecordBytes+1))==0,"sparse oversize fixture");bool rejected=false;try{Prepared(folder,marker);}catch(const Oversize&){rejected=true;}Check(rejected&&Raw(folder+"/result/result.dat")&&!Raw(folder+"/complete.txt")&&!Exists(folder+"/complete.txt"),"oversize changed publication/source");++cases;}
 {Reset();const auto root=Fresh();const auto folder=Fixture(root,"j13");Put(folder+"/complete.txt",marker);Check(Reject([&]{Prepared(folder,marker);}),"already raw complete accepted by prepared API");Pack(folder);V1Read(folder,marker);++cases;}
 {Reset();const auto pair=FixturePair();Prepared(pair.second,marker);std::string root,name,part;Split(pair.second,root,name,part);Descriptor fd(open(Shard(root,name).c_str(),O_RDWR));std::string bad(FrontierSlotBytes,'X');Check(pwrite(fd.fd,bad.data(),bad.size(),0)==ssize_t(bad.size()),"damage frontier");::fsync(fd.fd);V1Read(pair.second,marker);Check(Prepared(pair.second,marker)==0,"missing-dir retry cannot repair frontier");V1Read(pair.first,"prior complete\n");++cases;}
 {Reset();const auto pair=FixturePair();const auto expected=Expected();Put(pair.second+"/result/result.dat","BAD16 matrix bytes\n");Check(Reject([&]{PackPrepared(pair.second,marker,expected);}),"pre-Collect same-size mutation accepted");Check(!Exists(pair.second+"/complete.txt")&&Raw(pair.second),"stale receipt published or source deleted");++cases;}
 {Reset();const auto pair=FixturePair();auto expected=Expected();expected["missing.dat"]=Hex(Digest("missing"));Check(Reject([&]{PackPrepared(pair.second,marker,expected);}),"missing expected member accepted");Check(Reject([&]{PackPrepared(pair.second,marker,{});}),"empty expected set accepted");expected=Expected();expected["../escape"]=Hex(Digest("x"));Check(Reject([&]{PackPrepared(pair.second,marker,expected);}),"unsafe expected path accepted");++cases;}
 {Reset();const auto pair=FixturePair();Prepared(pair.second,marker);auto expected=Expected();expected["result/result.dat"]=Hex(Digest("bad archived expectation"));Check(Reject([&]{PackPrepared(pair.second,marker,expected);}),"missing-dir retry ignored expected member digest");V1Read(pair.second,marker);++cases;}
 std::cout<<"PASS PackPrepared: "<<cases<<" cases; no ephemeral complete; V1 member reader; 4 existing-shard fsyncs; crash retry; EIO propagation; source and readback corruption; oversize retention; idempotence and mirror repair\n";return 0;
 }catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}}
