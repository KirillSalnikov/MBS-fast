#include <string>
#include <stdexcept>
static int invariantSpillFd=-1;
static int invariantPublicMembers=-1;
#define MBS_MEMORY_STORE_SPILL_OPEN(fd) (invariantSpillFd=(fd))
#define MBS_MEMORY_PUBLICATION_EVENT(event,path) do { if(invariantPublicMembers>0&&!--invariantPublicMembers)throw std::runtime_error("injected public publication failure"); } while(false)
#include "../src/FullAuto.cpp"
#undef MBS_MEMORY_STORE_SPILL_OPEN

using namespace FullAuto;
namespace {
template<class Action> void Reject(Action action,const std::string& label){
    bool rejected=false;try{action();}catch(const std::exception&){rejected=true;}
    Require(rejected,"memory invariant was bypassed: "+label);
}
void CheckSealed(const std::string& root,const std::string& folder){
    const auto data=folder+"/result/result.dat",marker=folder+"/complete.txt";
    const auto bytes=Read(data),completion=Read(marker);
    auto& store=*FullAutoMemoryStore::Active();const auto epoch=store.Epoch();
    Atomic(data,bytes);Atomic(marker,completion);
    Require(store.Epoch()==epoch,"identical sealed Atomic changed generation");
    Reject([&]{Atomic(data,bytes+"changed\n");},"matrix overwrite");
    Reject([&]{Atomic(marker,completion+"changed\n");},"marker overwrite");
    Reject([&]{Atomic(folder+"/new_file","late");},"new file beneath sealed job");
    for(const auto& path:{data,marker,folder}){
        Reject([&]{RemovePath(path);},"remove "+path);
        Reject([&]{RenamePath(path,root+"/other");},"rename "+path);
    }
    Atomic(root+"/unsealed","other");
    Reject([&]{RenamePath(root+"/unsealed",data);},"replace sealed target");
    const auto tail=folder.substr(root.size());
    Reject([&]{Atomic(root+"/"+tail+"/result/result.dat",bytes+"alias\n");},"double-slash alias");
    Reject([&]{Atomic(root+"/."+tail+"/result/result.dat",bytes+"alias\n");},"dot alias");
    Reject([&]{RenamePath(root+"/jobs",root+"/all-moved");},"rename sealed ancestor");
    Require(Read(data)==bytes&&Read(marker)==completion,"rejected mutation changed sealed bytes");
}
void Validate(const std::string& folder,bool derived){
    ValidationBatch batch("memory invariant fixture",true);std::set<std::string> visited;
    ValidateStoredJob(folder,visited);
    if(derived)Require(visited.size()==1&&visited.count(folder),"memory graph shortcut unexpectedly traversed history");
}
}
int main(int argc,char** argv){try{
    Require(argc==2,"fresh scratch root required");const std::string disk=argv[1];
    Require(!Exists(disk),"use a new invariant-test output directory");Directory(disk);
    const std::string root=disk+"/virtual",leaf=root+"/jobs/native",parent=root+"/jobs/derived",checkpoint=disk+"/checkpoint.bin";
    FullAutoMemoryStore::Store store(root,32768,disk);
    {
        FullAutoMemoryStore::Scope scope(store);Directory(leaf+"/result");
        // The first payload exceeds the whole RAM budget, so offset zero in
        // this test-only anonymous spill is known without exposing internals.
        std::string data;for(int i=0;i<3000;++i)data+="150 1 1 .1 0 0 .1 .8 0 0 0 0 .6 0 0 0 0 .4\n";
        invariantSpillFd=-1;Atomic(leaf+"/result/result.dat",data);
        const int sourceSpillFd=invariantSpillFd;Require(sourceSpillFd>=0&&store.Stats().spillBytes>=data.size(),"source did not spill");
        Atomic(leaf+"/complete.txt","native 1 "+Hash(data)+" none\n");Validate(leaf,false);CheckSealed(root,leaf);

        Directory(parent+"/result");const std::string derived="150 1 1 .1 0 0 .1 .8 0 0 0 0 .6 0 0 0 0 .4\n";
        const auto sources=leaf+"\t"+FileHash(leaf+"/complete.txt")+"\n";
        Atomic(parent+"/result/result.dat",derived);Atomic(parent+"/angular_sources.tsv",sources);
        Atomic(parent+"/complete.txt","derived 1 "+Hash(derived)+" none "+Hash(sources)+"\n");
        Validate(parent,true);CheckSealed(root,parent);
        Require(ReadVerifiedText(parent+"/angular_sources.tsv",Hash(sources))==sources,"own provenance changed");
        Reject([&]{ReadVerifiedText(parent+"/angular_sources.tsv",Hash(sources+"changed"));},"wrong own-provenance digest");

        store.Snapshot(checkpoint,"invariant checkpoint");
        {
            FullAutoMemoryStore::Store restored(root,32768,disk);
            Require(restored.Restore(checkpoint)=="invariant checkpoint","checkpoint metadata mismatch");
            FullAutoMemoryStore::Scope restoredScope(restored);
            Validate(leaf,false);Validate(parent,true);CheckSealed(root,leaf);CheckSealed(root,parent);
            Require(Read(leaf+"/result/result.dat")==data&&Read(parent+"/result/result.dat")==derived,"snapshot changed all stored values");
        }
        Require(FullAutoMemoryStore::Active()==&store,"nested restore lost active store");
        // A stored hash can be used for immutable provenance without rereading
        // the old payload. When its bytes are needed, corruption must fail.
        const auto before=FileHash(leaf+"/result/result.dat");
        Require(pwrite(sourceSpillFd,"X",1,0)==1,"could not inject isolated spill corruption");
        Require(FileHash(leaf+"/result/result.dat")==before,"digest metadata unexpectedly changed");
        Validate(parent,true);
        Reject([&]{Read(leaf+"/result/result.dat");},"corrupt spilled source read");
        Reject([&]{ReadVerifiedText(leaf+"/result/result.dat",before);},"corrupt verified source read");
        Reject([&]{store.Snapshot(checkpoint,"must not publish corruption");},"corrupt spill snapshot");
    }
    Require(FullAutoMemoryStore::Active()==nullptr,"memory scope leaked");
    {FullAutoMemoryStore::Store intact(root,32768,disk);Require(intact.Restore(checkpoint)=="invariant checkpoint","failed checkpoint replaced valid generation");}
    {
        // Observation counters must survive log throttling, while worker
        // particle staging must retain the reader's implicit sidecar.
        ArgPP args;RunConfig config;Engine e(args,config);e.root=disk+"/observation";Directory(e.root);
        e.memoryMode=true;e.memoryStore.reset(new FullAutoMemoryStore::Store(e.root,1024*1024,disk));
        e.memoryScope.reset(new FullAutoMemoryStore::Scope(*e.memoryStore));
        e.memoryStore->Directory(e.root);Atomic(e.root+"/particle.dat","geometry");Atomic(e.root+"/particle.dat.visibility","visibility");
        const auto files=e.MemoryWorkerInputs({"--particle-file",e.root+"/./particle.dat",e.root+"/particle.dat"});
        Require(files.size()==2&&files[0].path=="particle.dat"&&files[0].bytes=="geometry"
            &&files[1].path=="particle.dat.visibility"&&files[1].bytes=="visibility","worker input sidecar/canonical dedup lost");
        Job job;job.folder=e.root+"/jobs/test";job.count=16;job.phi=4;job.inputAngles={150,170};job.seconds=5;
        ReadPhaseTimings(job,"phase\tseconds\tclock\ntrace_prepare\t2\tstage_wall_sum\n"
            "diffraction_and_controls\t3\tstage_wall_sum\ngpu_diffraction_including_packing\t1\tstage_wall_sum\n"
            "analytic_controls\t2\tstage_wall_sum\nprocessing\t4\twall\n"
            "pipeline_wait_for_prepared\t0.2\tconsumer_wall_subset_of_processing\n"
            "pipeline_producer_wait\t0.3\tproducer_wall_subset_of_processing\npipeline_chunks\t7\tcount\n");
        e.dispatchProgress.push_back({});e.NativeCommitted(job,0.1);
        Require(e.nativeProfiledBlocks==1&&e.nativeGpuProfiledBlocks==1&&e.nativePipelineChunks==7
            &&e.nativeTraceSeconds==2&&e.nativeAnalyticControlSeconds==2&&e.nativeGpuDiffractionSeconds==1,
            "native phase split lost or double-counted");
        std::ifstream observer(e.root+"/execution_progress.json");
        const std::string published((std::istreambuf_iterator<char>(observer)),{});
        Require(published.find("\"native_phase_totals\":")!=std::string::npos
            &&published.find("\"native_requests_committed\":1")!=std::string::npos,"live memory progress was not physically published");
        struct stat st{};Require(stat((e.root+"/fullauto_checkpoint.bin").c_str(),&st)!=0,"observer created a recovery snapshot");

        auto raw=[](const std::string& path){std::ifstream file(path,std::ios::binary);Require(bool(file),"missing physical publication "+path);
            return std::string((std::istreambuf_iterator<char>(file)),{});};
        auto stage=[&](const std::string& tag){const auto source=e.root+"/round_snapshots/"+tag;
            Directory(source);Atomic(source+"/matrix.dat",tag+" matrix");Atomic(source+"/parameters.csv",tag+" parameters");
            Atomic(e.root+"/intermediate_current.json","{\"snapshot_folder\":"+Quoted(source)+",\"matrix\":"+Quoted(source+"/matrix.dat")
                +",\"accuracy_and_parameters\":"+Quoted(source+"/parameters.csv")+"}");};
        const auto pointer=e.root+"/intermediate_current.json";
        stage("a");e.ExportMemoryPublic(true);const auto oldPointer=raw(pointer),a=e.MemorySnapshotFolder(oldPointer);
        stage("b");invariantPublicMembers=1;
        Reject([&]{e.ExportMemoryPublic(true);},"injected multi-file public publication failure");invariantPublicMembers=-1;
        Require(raw(pointer)==oldPointer&&raw(a+"/matrix.dat")=="a matrix"&&raw(a+"/parameters.csv")=="a parameters",
            "partial public export mixed matrix generations");
        e.ExportMemoryPublic(true);const auto b=e.MemorySnapshotFolder(raw(pointer));
        Require(a!=b&&raw(b+"/matrix.dat")=="b matrix"&&raw(b+"/parameters.csv")=="b parameters"
            &&raw(a+"/matrix.dat")=="a matrix","public snapshot switch lost old/current consistency");
        stage("c");e.ExportMemoryPublic(true);const auto c=e.MemorySnapshotFolder(raw(pointer));
        Require(raw(c+"/matrix.dat")=="c matrix"&&raw(b+"/parameters.csv")=="b parameters"
            &&lstat(a.c_str(),&st)!=0,"public generations were not bounded to current and previous");
        auto put=[](const std::string& path,const std::string& value){std::ofstream file(path,std::ios::binary);file<<value;file.close();Require(bool(file),"fixture write failed");};
        const auto checkpointPath=e.root+"/fullauto_checkpoint.bin",temporary=checkpointPath+".tmp.aBc123";
        put(checkpointPath,"committed checkpoint");put(temporary,"abandoned snapshot");chmod(temporary.c_str(),0600);
        const auto parent=e.root+"/public_snapshots";
        const pid_t child=fork();Require(child>=0,"pending publication fixture fork failed");
        if(child==0){const auto path=parent+"/.pending-"+Number(getpid())+"-aBc123";
            if(mkdir(path.c_str(),0700)!=0)_exit(1);
            put(path+"/.mbs-public-pending-v1","mbs_public_pending_v1\n"+e.root+"\n"+Number(getpid())+"\n");put(path+"/partial.dat","partial");_exit(0);}
        int status=0;Require(waitpid(child,&status,0)==child&&WIFEXITED(status)&&WEXITSTATUS(status)==0,"pending publication fixture failed");
        const auto abandoned=parent+"/.pending-"+Number(child)+"-aBc123",live=parent+"/.pending-"+Number(getpid())+"-aBc123";
        Require(mkdir(live.c_str(),0700)==0,"live pending fixture mkdir failed");
        put(live+"/.mbs-public-pending-v1","mbs_public_pending_v1\n"+e.root+"\n"+Number(getpid())+"\n");
        e.CleanupMemoryAbandonedPublications();
        Require(lstat(abandoned.c_str(),&st)!=0&&lstat(temporary.c_str(),&st)!=0&&lstat(live.c_str(),&st)==0
            &&raw(checkpointPath)=="committed checkpoint"&&raw(c+"/matrix.dat")=="c matrix"
            &&e.memoryAbandonedTemps==2,"abandoned publication cleanup removed live/committed data or leaked a dead temporary");
        WorkerMemoryStage::RemoveTree(live);
    }
    {
        // An explicitly requested geometry export is a narrow exception to
        // the per-job namespace, never permission to accept another root file.
        const auto root=disk+"/geometry_whitelist",geometry=root+"/geometry/deep/saved.particle";
        ArgPP args;args.AddRule("save_geometry",1,true);const char* argv[]={"probe","--save_geometry",geometry.c_str()};args.Parse(3,argv);
        RunConfig config;Engine e(args,config);e.root=root;Directory(root);e.memoryMode=true;
        e.memoryStore.reset(new FullAutoMemoryStore::Store(root,1024*1024,disk));e.memoryScope.reset(new FullAutoMemoryStore::Scope(*e.memoryStore));
        Job job;job.folder=root+"/jobs/test";
        const std::string value="geometry bytes";const auto digest=FullAutoMemoryStore::Hash(value);
        auto stream=[&](const std::string& path,bool& read){size_t offset=0;e.AcceptMemoryWorkerStream({path,value.size(),digest},
            [&](char* data,size_t requested)->size_t{read=true;const auto n=std::min(requested,value.size()-offset);std::memcpy(data,value.data()+offset,n);offset+=n;return n;},job);};
        for(const auto& path:{"unexpected.dat","geometry/deep/saved.particle.extra","geometry/deep/other.particle.visibility","jobs/another/result.dat"}){
            bool read=false;Reject([&]{stream(path,read);},"unauthorized streamed root export");Require(!read,"unauthorized stream consumed payload before namespace validation");
            Reject([&]{e.AcceptMemoryWorkerFiles({{path,value}},job);},"unauthorized bundled root export");Require(!Exists(root+"/"+path),"unauthorized file reached working state");
        }
        e.AcceptMemoryWorkerFiles({{"geometry/deep/saved.particle",value}},job);
        bool read=false;stream("geometry/deep/saved.particle.visibility",read);
        Require(read&&Read(geometry)==value&&Read(geometry+".visibility")==value,"explicit geometry/visibility exception was rejected");
    }
    std::cout<<"PASS memory invariants: sealed jobs/path aliases; spill corruption/restore; streamed inputs/progress; atomic public generations and interrupted export; dead temporary cleanup preserves live/current data\n";
    return 0;
}catch(const std::exception& error){std::cerr<<error.what()<<'\n';return 1;}}
