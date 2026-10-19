#include "../src/FullAuto.cpp"
using namespace FullAuto;
void SetupMemory(Engine& e,const std::string& root,const std::string& binary,bool controls){
    e.root=root;e.signature="bounded-memory-regression-v1";Directory(root);Directory(root+"/jobs");e.exe=binary;e.binaryHash=FileHash(binary);
    e.controls=controls;e.sharedMeans=controls;e.gpus={""};e.seeds.assign(production,production+8);e.workerThreads=1;e.analyticMeanThreads=1;e.hostBudgetMb=512;
    e.phi=1;e.target=6;e.theta={0,5,170,175,180};e.initial=512;e.orientationBlock=64;e.workerPolicy="persistent";e.resourcePolicy="off";e.workQueuePolicy="on";
    e.base={"--method","po","--backend","cpu","--particle","1",".2",".1","--ri","1.3116","0","--wavelength-um",".532","--symmetry","1","1","--cutoff-profile","off"};
    WriteGrid(root+"/requested_theta.csv",e.theta);
}
int main(int argc,char** argv){try{
    if(argc==2){
        ArgPP args;RunConfig config;Engine e(args,config);const bool bounded=std::string(argv[1])=="--bounded";
        e.payloads.SetLimit((bounded?8ULL:512ULL)*1024*1024);
        for(int i=0;i<1024;++i){Job j;j.identity=Number(i);j.folder="fixture/"+j.identity;j.count=64;j.seed=11;j.depth=6;j.offset=i*64;
            j.inputAngles.resize(2048);j.angles.resize(2048);j.matrix.resize(2048);for(size_t r=0;r<2048;++r){j.inputAngles[r]=j.angles[r]=double(r)/12;for(double& v:j.matrix[r])v=double(i+r);}
            e.CacheJob(j);}
        Require(e.payloads.Bytes()<=e.payloads.Limit(),"cache exceeded limit");
        std::istringstream proc(Read("/proc/self/status"));std::string line;while(std::getline(proc,line))if(line.find("VmRSS:")==0)std::cout<<line<<'\n';
        std::cout<<"resident_cache_bytes="<<e.payloads.Bytes()<<" evictions="<<e.payloads.evictions<<'\n';return 0;
    }
    Require(argc==3,"OUTPUT BINARY");const std::string root=argv[1];Directory(root);
    {
        ArgPP args;RunConfig config;Engine guard(args,config);guard.root=root+"/guard";Directory(guard.root);
        Atomic(guard.root+"/adaptive_checkpoint.tsv","saved checkpoint\n");
        const char* previous=getenv("MBS_HOST_MEM_RESERVE_MB");const std::string saved=previous?previous:"";
        const auto memory=HostMemoryReserve::Read();setenv("MBS_HOST_MEM_RESERVE_MB",Number(memory.totalKb/1024+1).c_str(),1);
        bool stopped=false;try{guard.CheckHostMemory();}catch(const std::exception& e){stopped=std::string(e.what()).find("host RAM reserve exhausted")!=std::string::npos;}
        if(previous)setenv("MBS_HOST_MEM_RESERVE_MB",saved.c_str(),1);else unsetenv("MBS_HOST_MEM_RESERVE_MB");
        Require(stopped&&Exists(guard.root+"/resource_stop.json")&&Read(guard.root+"/adaptive_checkpoint.tsv")=="saved checkpoint\n","memory guard failed to preserve checkpoint");
    }
    BoundedCache<int> bounded(512);auto pinned=std::make_shared<int>(7);bounded.Put("a",pinned,sizeof(int));
    for(int i=0;i<100;++i)bounded.Put(Number(i),std::make_shared<int>(i),sizeof(int));
    Require(bounded.Bytes()<=bounded.Limit()&&*pinned==7&&bounded.evictions>0,"bounded eviction or pin lifetime failed");
    const std::string hashFile=root+"/hash-test";std::string bytes(150000,'x');bytes[65536]='\0';bytes[130000]=char(255);Atomic(hashFile,bytes);
    Require(FileHash(hashFile)==Hash(bytes),"streaming hash changed");bytes[65537]='y';Atomic(hashFile,bytes);Require(FileHash(hashFile)==Hash(bytes),"changed file digest was reused");
    for(bool controls:{false,true}){
        ArgPP args;RunConfig config;Engine e(args,config);SetupMemory(e,root+(controls?"/controls":"/plain"),argv[2],controls);
        e.Prefetch({Engine::WorkRequest{e.initial,0,6,0,e.theta}},e.seeds.data());
        const auto full=e.PrefixStage(e.initial,e.seeds.data(),8,6,e.theta);const auto reference=e.Data(full);
        for(const auto& item:e.cache)Require(item.second.matrix.empty()&&item.second.controls.empty()&&item.second.angles.empty(),"metadata retains matrix payload");
        e.payloads.SetLimit(1);e.forbidNativeWork=true;
        auto same=e.PrefixStage(e.initial,e.seeds.data(),8,6,e.theta);Require(e.Data(same)==reference,"evicted exact prefix changed matrix");
        const std::vector<double> subset={0,170,180};auto rows=e.PrefixStage(e.initial,e.seeds.data(),8,6,subset);
        auto direct=e.Stage(64,e.seeds.data(),8,6,subset,true,false,0);
        for(size_t s=0;s<8;++s){for(size_t r=0;r<3;++r)for(size_t k=0;k<16;++k)Require(std::fabs(rows[s].matrix[r][k]-reference[s][r*2][k])<=1e-12*std::fabs(reference[s][r*2][0]),"evicted prefix row changed");
            IntegralReader reader;auto a=reader.Resolve(full[s].folder,6),b=reader.Resolve(rows[s].folder,6);Require(a.area==b.area&&a.ext==b.ext,"integrals changed after eviction");}
        Require(e.payloads.Bytes()==0,"one-byte cache retained payload");
        Engine resumed(args,config);SetupMemory(resumed,e.root,argv[2],controls);resumed.payloads.SetLimit(1);resumed.WarmCache();resumed.forbidNativeWork=true;
        const size_t before=resumed.cache.size();resumed.Prefetch({Engine::WorkRequest{resumed.initial,0,6,0,resumed.theta}},resumed.seeds.data());
        Require(resumed.Data(resumed.PrefixStage(resumed.initial,resumed.seeds.data(),8,6,resumed.theta))==reference,"resume changed matrices");
        Require(resumed.cache.size()==before,"resume retraced saved orientations");
        auto victim=resumed.cache.begin()->second;while(victim.assembled){auto it=resumed.cache.find(victim.identity);++it;Require(it!=resumed.cache.end(),"no native victim");victim=it->second;}
        const auto original=Read(victim.data);Atomic(victim.data,original+"corrupt\n");bool rejected=false;try{resumed.Materialize(victim);}catch(const std::exception&){rejected=true;}Atomic(victim.data,original);Require(rejected,"corrupt evicted source accepted");
        std::cout<<"PASS "<<(controls?"controls":"plain")<<": all16 bitwise equal, integral provenance, zero-residency resume, corrupt-source rejection\n";
    }
    std::cout<<"PASS bounded cache pinning and streaming hash mutation detection\n";return 0;
}catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}}
