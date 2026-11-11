#include <unistd.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <map>
#include <fstream>
#include <functional>
static int TracedFsync(int);
static void PackEvent(const char*,const std::string&);
#define fsync TracedFsync
#define MBS_PACK_EVENT(event,folder) PackEvent(event,folder)
#include "../src/FullAuto.cpp"
#undef fsync
#undef MBS_PACK_EVENT
using namespace FullAuto;
static bool tracking=false,packedExpected=false,inject=false,injectAfterRename=false;
static std::string currentFolder;static size_t coreFsyncs=0;static std::vector<std::string> flushed;
static std::string FdName(int fd){char buffer[4096];const auto name="/proc/self/fd/"+std::to_string(fd);const auto n=readlink(name.c_str(),buffer,sizeof(buffer)-1);if(n<0)return "";buffer[n]=0;return buffer;}
static bool Raw(const std::string& p){struct stat st{};return lstat(p.c_str(),&st)==0;}
static int TracedFsync(int fd){const auto path=FdName(fd);if(tracking){flushed.push_back(path);if(path.find("/jobs")!=std::string::npos||path.find("/job_packs/")!=std::string::npos)++coreFsyncs;
 if(inject&&((!injectAfterRename&&path==currentFolder+"/complete.txt.tmp")||(injectAfterRename&&path==currentFolder&&Raw(currentFolder+"/complete.txt")))){inject=false;errno=EIO;return -1;}}
 return ::fsync(fd);}
static void PackEvent(const char*,const std::string& folder){if(tracking&&packedExpected)Require(!Raw(folder+"/complete.txt"),"virtual native path exposed loose completion");}
static void Put(const std::string& path,const std::string& text){std::ofstream out(path);out<<text;out.close();Require(bool(out),"fixture write failed");}
static bool Reject(const std::function<void()>& f){try{f();return false;}catch(const std::exception&){return true;}}
static void Setup(Engine& e,const std::string& root,bool packed,bool controls){e.root=root;e.signature="native-pack-proof-v1";Directory(root+"/jobs");e.compactIndex=true;e.packJobs=packed;e.controls=controls;e.workerThreads=1;e.phi=4;e.kernel.block=64;e.kernel.warp="auto";e.gpus={""};e.workerPolicy="persistent";}
static Job Fixture(Engine& e,bool pair,bool prefix){Job j;j.identity="abc";j.folder=e.root+"/jobs/jabc";Directory(j.folder+"/result");j.data=j.folder+"/result/result.dat";j.summary=j.folder+"/result/result_analytic_facets.tsv";j.inputAngles={150,153.5};j.count=32;j.seed=11;j.depth=pair?10:6;j.pairedDepth=pair?8:0;j.phi=4;j.block=64;j.warp="auto";j.hasControls=e.controls;j.executionHash="123";j.weightHash="456";j.seconds=1;
 Matrix values(2);for(size_t r=0;r<2;++r)for(size_t c=0;c<16;++c)values[r][c]=double(20+r+c)/32;
 auto result=[&](const std::string& stem){WriteMatrix(stem+".dat",j.inputAngles,values);Put(stem+"_integrals.tsv","theta_integral\n1 2 3\n");if(e.controls){std::ostringstream cv;cv<<"theta raw reflection shadow sampleR sampleS residual hybrid refinement setup wr ws\n";for(double angle:j.inputAngles)cv<<angle<<" 1 2 3 4 5 6 7 .01 .02 1 1\n";Put(stem+"_analytic_facets.tsv",cv.str());}};
 result(j.folder+"/result/result");if(pair){result(j.folder+"/result/result_depth8");Put(j.folder+"/result/result_depth8_integrals_revision.txt","revision\n");}
 if(prefix){result(j.folder+"/result/result_prefix16");if(pair){result(j.folder+"/result/result_prefix16_depth8");Put(j.folder+"/result/result_prefix16_depth8_integrals_revision.txt","revision\n");}}
 Put(j.folder+"/stdout.log","Hard tree-limit hits: 0\nOrientations still incomplete: 0\n");
 Put(j.folder+"/command.args","worker\n"+std::string(prefix?"--orientation-prefix-output\n16\n":""));return j;}
static void WarmShard(const Job& j){std::string root,name,member;PackedJobStore::Split(j.folder,root,name,member);std::string old;
 for(int n=1;;++n){const auto candidate="j"+PackedJobStore::Hex(n);if(candidate!=name&&PackedJobStore::Shard(root,candidate)==PackedJobStore::Shard(root,name)){old=root+"/jobs/"+candidate;break;}}
 Directory(old);Put(old+"/complete.txt","older checkpoint\n");PackedJobStore::Pack(old);}
static void Launch(Engine& e,const Job& j,bool prefix){AtomicFileBatch batch(j.folder);Atomic(j.folder+"/command.args","worker\n--real-launch-fixture\n"+std::string(prefix?"--orientation-prefix-output\n16\n":""));Atomic(j.folder+"/worker.pid","123 456 boot\n");Require(Raw(j.folder+"/worker.pid")&&Raw(j.folder+"/command.args"),"ownership invisible before send");if(!e.packJobs)batch.Commit();std::remove((j.folder+"/worker.pid").c_str());}
int main(int argc,char**argv){try{Require(argc==2,"fresh scratch path required");const std::string root=argv[1];int cases=0;
 for(bool controls:{false,true})for(bool pair:{false,true})for(bool prefix:{false,true}){
  Matrix reference;std::vector<std::array<double,5>> referenceControls;
  for(int mode=0;mode<3;++mode){ArgPP args;RunConfig config;Engine e(args,config);Setup(e,root+"/c"+Number(cases)+"/m"+Number(mode),mode==1,controls);auto j=Fixture(e,pair,prefix);if(mode<2)WarmShard(j);
   currentFolder=j.folder;coreFsyncs=0;flushed.clear();tracking=true;packedExpected=mode==1;Launch(e,j,prefix);e.CompleteNative(j);if(mode==0)PackedJobStore::Pack(j.folder);tracking=false;packedExpected=false;
   Require(e.cache.count(j.identity)==1,"completion not indexed");std::set<std::string> verified;ValidateStoredJob(j.folder,verified);Job restored;Require(e.ReadCachedIndex(j.folder,restored),"published native sampling missing");e.payloads.Clear();restored=e.Materialize(restored);Require(j.matrix==restored.matrix&&j.controls==restored.controls,"resume changed all16/controls");
   if(mode==0){reference=j.matrix;referenceControls=j.controls;}else Require(reference==j.matrix&&referenceControls==j.controls,"packed native estimator differs");
   if(mode<2)Require(!Raw(j.folder),"packed job loose tree retained");else Require(Raw(j.folder+"/complete.txt"),"no-pack did not publish raw marker");
   const size_t expected=mode==0?(pair?14:12):mode==1?4:(pair?10:8);Require(coreFsyncs==expected,"unexpected native fsync count got "+Number(coreFsyncs)+" expected "+Number(expected));
   std::cout<<"native controls="<<controls<<" pair="<<pair<<" prefix="<<prefix<<" mode="<<mode<<" fsync="<<coreFsyncs<<" exact16=true\n";
  }++cases;
 }
 for(int error=0;error<3;++error){ArgPP args;RunConfig config;Engine e(args,config);Setup(e,root+"/oversize"+Number(error),true,true);auto j=Fixture(e,true,true);int fd=open((j.folder+"/large-output").c_str(),O_CREAT|O_WRONLY,0644);Require(fd>=0&&ftruncate(fd,PackedJobStore::MaximumRecordBytes+1)==0,"sparse file");close(fd);Launch(e,j,true);currentFolder=j.folder;flushed.clear();tracking=true;inject=error!=0;injectAfterRename=error==2;bool rejected=Reject([&]{e.CompleteNative(j);});tracking=false;Require(rejected==(error!=0),"fallback fsync EIO hidden");
  if(error==1)Require(!Raw(j.folder+"/complete.txt")&&!e.cache.count(j.identity),"failed raw completion published");
  if(error==0){Require(Raw(j.folder+"/complete.txt")&&e.unpackedOversizeJobs.count(j.identity),"oversize fallback missing");std::set<std::string> verified;ValidateStoredJob(j.folder,verified);Require(std::find(flushed.begin(),flushed.end(),j.data)!=flushed.end()&&std::find(flushed.begin(),flushed.end(),j.folder+"/command.args")!=flushed.end()&&std::find(flushed.begin(),flushed.end(),j.folder+"/large-output")!=flushed.end(),"raw fallback omitted private tree flush");}
  ++cases;
 }
 std::cout<<"PASS native publication "<<cases<<" cases; all16/controls exact; prefix/paired receipts; resume; no-pack unchanged; 12/14->4 core fsync; durable oversize and EIO\n";return 0;
 }catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}}
