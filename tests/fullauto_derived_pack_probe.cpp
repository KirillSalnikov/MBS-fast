#include <unistd.h>
#include <fcntl.h>
#include <sys/stat.h>
#include <functional>
#include <string>
static int ReviewFsync(int);
static void ReviewEvent(const char*,const std::string&);
#define fsync ReviewFsync
#define MBS_PACK_EVENT(event,folder) ReviewEvent(event,folder)
#include "../src/FullAuto.cpp"
#undef fsync
#undef MBS_PACK_EVENT
using namespace FullAuto;
static bool failFlush=false,observe=false,trackFlush=false;
static size_t coreFlushes=0;
static std::string FdName(int fd){char b[4096];const auto p="/proc/self/fd/"+std::to_string(fd);const auto n=readlink(p.c_str(),b,sizeof(b)-1);if(n<0)return "";b[n]=0;return b;}
static size_t observedEvents=0;
static bool Raw(const std::string& p){struct stat st{};return lstat(p.c_str(),&st)==0;}
static int ReviewFsync(int fd){if(trackFlush){const auto p=FdName(fd);if(p.find("/jobs")!=std::string::npos||p.find("/job_packs/")!=std::string::npos)++coreFlushes;}if(failFlush){failFlush=false;errno=EIO;return -1;}return ::fsync(fd);}
static void ReviewEvent(const char*,const std::string& folder){if(observe){++observedEvents;Require(!Raw(folder+"/complete.txt"),"prepared derived completion visible loose before pack commit");}}
static bool Reject(const std::function<void()>& fn){try{fn();return false;}catch(const std::exception&){return true;}}
static void Put(const std::string& path,const std::string& text){std::ofstream s(path);s<<text;s.close();Require(bool(s),"fixture write");}
static void Mutate(const std::string& path){struct stat st{};Require(stat(path.c_str(),&st)==0,"fixture stat");auto text=Read(path);Require(!text.empty(),"empty mutation target");text[0]=text[0]=='x'?'y':'x';Put(path,text);struct timespec times[2]={st.st_atim,st.st_mtim};Require(utimensat(AT_FDCWD,path.c_str(),times,0)==0,"restore mtime");Require(FileHash(path)!=Hash(std::string()),"fixture hash");}
static void Setup(Engine& e,const std::string& path,bool pack,bool compact){e.root=path;e.signature="derived-independent-review";Directory(path+"/jobs");e.packJobs=pack;e.compactIndex=compact;e.controls=true;e.workerThreads=1;e.gpus={""};e.phi=4;e.orientationBlock=64;}
struct Fixture {Job job;std::unique_ptr<Engine::DerivedFileBatch> batch;std::string metadata;Matrix matrix;};
static Fixture Make(Engine& e,bool legacyBatch=true){
    Fixture f;Job& j=f.job;j.identity="abc";j.folder=e.root+"/jobs/pabc";Directory(j.folder+"/result");j.data=j.folder+"/result/result.dat";j.summary=j.folder+"/result/result_analytic_facets.tsv";j.inputAngles={150,153.5};j.count=128;j.seed=11;j.depth=10;j.pairedDepth=8;j.phi=4;j.block=64;j.warp="auto";j.hasControls=true;j.executionHash="123";j.weightHash="456";j.seconds=1;
    f.batch.reset(new Engine::DerivedFileBatch(j.folder,e.PrepareDerivedPublication(),legacyBatch));f.matrix.assign(2,Mueller{});for(size_t r=0;r<2;++r)for(size_t k=0;k<16;++k)f.matrix[r][k]=(20+r+k)/32.;WriteMatrix(j.data,j.inputAngles,f.matrix);
    std::string cv="theta raw reflection shadow sampleR sampleS residual hybrid refinement setup wr ws\n150 1 2 3 4 5 6 7 .01 .02 1 1\n153.5 1 2 3 4 5 6 7 .01 .02 1 1\n";Atomic(j.summary,cv);
    const auto lower=j.folder+"/result/result_depth8";WriteMatrix(lower+".dat",j.inputAngles,f.matrix);Atomic(lower+"_analytic_facets.tsv",cv);Atomic(j.folder+"/result/result_integrals.tsv","Cext Csca\n1 1\n");Atomic(lower+"_integrals.tsv","Cext Csca\n1 1\n");Atomic(lower+"_integrals_revision.txt","paired-v1\n");
    const auto pair=lower+".dat "+FileHash(lower+".dat")+"\n"+lower+"_analytic_facets.tsv "+FileHash(lower+"_analytic_facets.tsv")+"\n";Atomic(j.folder+"/paired_depth.tsv",pair);Atomic(j.folder+"/cutoff_sources.tsv","signed-fixture\n");
    std::string integrity;for(const auto& name:{std::string("result/result_integrals.tsv"),std::string("result/result_depth8_integrals.tsv"),std::string("result/result_depth8_integrals_revision.txt"),std::string("paired_depth.tsv"),std::string("cutoff_sources.tsv")})integrity+=j.folder+"/"+name+" "+FileHash(j.folder+"/"+name)+"\n";
    Atomic(j.folder+"/integral_outputs.fnv",integrity);Atomic(j.folder+"/angular_sources.tsv","");Atomic(j.folder+"/command.args","worker\n--review\n");Atomic(j.folder+"/stdout.log","Hard tree-limit hits: 0\nOrientations still incomplete: 0\n");
    f.metadata=j.identity+" 1 "+FileHash(j.data)+" "+FileHash(j.summary)+" "+Hash("")+" "+Hash(integrity)+"\n";return f;
}
int main(int argc,char**argv){try{Require(argc==2,"fresh scratch root required");const std::string root=argv[1];int cases=0;
    for(const auto& target:{std::string("result/result.dat"),std::string("result/result_analytic_facets.tsv"),std::string("result/result_depth8.dat"),std::string("result/result_depth8_analytic_facets.tsv"),std::string("result/result_integrals.tsv"),std::string("result/result_depth8_integrals_revision.txt"),std::string("paired_depth.tsv"),std::string("cutoff_sources.tsv"),std::string("integral_outputs.fnv")}){
        ArgPP args;RunConfig config;Engine e(args,config);Setup(e,root+"/tamper"+Number(cases),true,true);auto f=Make(e);Mutate(f.job.folder+"/"+target);Require(Reject([&]{e.PublishDerivedComplete(f.job,f.metadata,*f.batch);}),"stale derived receipt accepted: "+target);Require(!Exists(f.job.folder+"/complete.txt")&&Raw(f.job.folder),"failed derived publication visible or deleted");++cases;
    }
    {ArgPP args;RunConfig config;Engine e(args,config);Setup(e,root+"/angular_tamper",true,true);auto f=Make(e);Put(f.job.folder+"/angular_sources.tsv","changed\n");Require(Reject([&]{e.PublishDerivedComplete(f.job,f.metadata,*f.batch);}),"changed angular provenance accepted");Require(!Exists(f.job.folder+"/complete.txt"),"changed angular provenance published");++cases;}
    for(int mode=0;mode<3;++mode){ArgPP args;RunConfig config;Engine e(args,config);const bool packed=mode!=0,compact=mode!=2;Setup(e,root+"/mode"+Number(mode),packed,compact);auto f=Make(e,mode==1);observe=packed&&compact;e.PublishDerivedComplete(f.job,f.metadata,*f.batch);observe=false;std::set<std::string> visited;ValidateStoredJob(f.job.folder,visited);auto loaded=e.Load(f.job);Require(loaded.matrix==f.matrix&&loaded.controls.size()==2,"all16/controls changed after publication");if(!compact)Require(Raw(f.job.folder+"/sampling.index"),"noncompact Load lost index path");e.CacheJob(loaded);Require(Exists(f.job.folder+"/complete.txt"),"completion disappeared");Require(Raw(f.job.folder+"/complete.txt")==!packed,"wrong raw/packed visibility");++cases;}
    {ArgPP args;RunConfig config;Engine e(args,config);Setup(e,root+"/oversize",true,true);auto f=Make(e);const auto name=f.job.folder+"/oversize_payload";const int fd=open(name.c_str(),O_CREAT|O_WRONLY,0644);Require(fd>=0&&ftruncate(fd,PackedJobStore::MaximumRecordBytes+1)==0,"sparse oversize setup");close(fd);e.PublishDerivedComplete(f.job,f.metadata,*f.batch);Require(Raw(f.job.folder+"/complete.txt")&&e.unpackedOversizeJobs.count(f.job.identity),"oversize did not retain durable loose completion");std::set<std::string> visited;ValidateStoredJob(f.job.folder,visited);++cases;}
    {ArgPP args;RunConfig config;Engine e(args,config);Setup(e,root+"/eio",true,true);auto f=Make(e);failFlush=true;Require(Reject([&]{e.PublishDerivedComplete(f.job,f.metadata,*f.batch);}),"EIO swallowed");Require(!Exists(f.job.folder+"/complete.txt")&&Raw(f.job.folder)&&e.unpackedOversizeJobs.empty(),"EIO treated as successful/oversize publication");++cases;}
    for(bool packed:{false,true}){ArgPP args;RunConfig config;Engine e(args,config);Setup(e,root+(packed?"/fsync_prepared":"/fsync_raw"),packed,true);const auto target=e.root+"/jobs/pabc";std::string older;for(int n=1;;++n){const auto name="j"+PackedJobStore::Hex(n);if(PackedJobStore::Shard(e.root,name)==PackedJobStore::Shard(e.root,"pabc")){older=e.root+"/jobs/"+name;break;}}Directory(older);Put(older+"/complete.txt","older completed record\n");PackedJobStore::Pack(older);coreFlushes=0;trackFlush=true;auto f=Make(e,true);e.PublishDerivedComplete(f.job,f.metadata,*f.batch);if(!packed)PackedJobStore::Pack(f.job.folder);trackFlush=false;Require(coreFlushes==(packed?4:22),"unexpected publication-only fsync count "+Number(coreFlushes));std::set<std::string> visited;ValidateStoredJob(f.job.folder,visited);Require(e.Load(f.job).matrix==f.matrix,"fsync comparison changed matrix");std::cout<<"publication paired+controls warm-shard "<<(packed?"prepared":"legacy-batch+Pack")<<" core_fsync="<<coreFlushes<<"\n";++cases;}
    Require(observedEvents>0,"prepared path not exercised");std::cout<<"PASS independent derived publication "<<cases<<" cases: stale-metadata same-size/restored-mtime rejection; all16/controls load; raw and noncompact compatibility; oversize fallback; EIO propagation; no loose complete before frontier\n";return 0;
}catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}}
