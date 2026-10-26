// Legacy parent -> compact child index, with no optical work or parent writes.
#define main window_fixture_main
#include "fullauto_window_probe.cpp"
#undef main

int main(int argc,char** argv){try{
    Require(argc==2,"OUTPUT");const std::string root=argv[1];
    ArgPP args;RunConfig config;Engine parent(args,config);
    SetupCacheEngine(parent,root+"/parent");parent.compactIndex=false;
    Require(unlink((parent.root+"/job_index_format.txt").c_str())==0,"cannot make legacy fixture");
    for(size_t s=0;s<8;++s){
        Job j;parent.SetSampling(j,512,parent.seeds[s],8,parent.theta,false);
        j.identity="legacy-"+Number(s);j.folder=parent.root+"/jobs/"+j.identity;
        j.data=j.folder+"/result/result.dat";j.summary=j.folder+"/result/result_analytic_facets.tsv";
        j.weightHash=Hash("");j.seconds=1;j.angles=parent.theta;j.matrix=Matrix(2);
        for(size_t r=0;r<2;++r){j.matrix[r][0]=100+s+r;j.matrix[r][5]=20+s;}
        Directory(j.folder+"/result");WriteMatrix(j.data,j.angles,j.matrix);
        Atomic(j.folder+"/stdout.log","Hard tree-limit hits: 0\nOrientations still incomplete: 0\n");
        parent.PublishComplete(j,j.identity+" 1 "+FileHash(j.data)+" none\n");
        Job clean=j;clean.matrix.clear();clean.angles.clear();clean=parent.Load(clean);parent.CacheJob(clean);
    }
    Require(!Exists(parent.root+"/completed_jobs.journal"),"fixture unexpectedly has journal");
    const auto frozen=FrozenParentFiles(parent.root);
    auto setup=[&](Engine& e){e.root=root+"/child";e.signature=parent.signature;e.exe=parent.exe;
        e.binaryHash=parent.binaryHash;e.controls=false;e.forbidNativeWork=true;e.refinementParent=parent.root;
        Directory(e.root);Directory(e.root+"/jobs");e.ConfigureIndex(false);};
    Engine child(args,config);setup(child);Require(child.compactIndex,"legacy child was not upgraded");
    child.WarmCache(1);Require(child.cache.size()==8,"legacy import lost jobs");
    const auto path=child.root+"/completed_jobs.journal";const auto journal=Read(path);
    Require(!journal.empty()&&child.indexImportFd<0,"import journal not durably closed");
    Require(FrozenParentFiles(parent.root)==frozen,"legacy import mutated parent");
    Engine resumed(args,config);setup(resumed);resumed.WarmCache(1);
    Require(resumed.cache.size()==8&&resumed.payloads.Bytes()==0,"journal restore eagerly read matrices");
    const auto progress=Read(resumed.root+"/cache_restore_status.json");
    Require(progress.find("\"journal_jobs\":8")!=std::string::npos
        &&progress.find("\"fallback_folders\":0")!=std::string::npos,"resume reread legacy per-job metadata");
    for(const auto& item:parent.cache){const auto a=parent.Materialize(item.second),b=resumed.Materialize(resumed.cache.at(item.first));
        Require(a.matrix==b.matrix&&a.inputAngles==b.inputAngles&&a.seconds==b.seconds,"journal changed estimator or scheduling");}
    Require(FrozenParentFiles(parent.root)==frozen,"journal resume mutated parent");
    // Interrupted append still falls back to authoritative legacy directories.
    size_t offset=0,last=0;while(offset<journal.size()){last=offset;auto end=journal.find('\n',offset);
        Require(end!=std::string::npos,"bad journal fixture");size_t size=0;std::istringstream h(journal.substr(offset,end-offset));h>>size;offset=end+1+size+1;}
    Atomic(path,journal.substr(0,last)+"100 ");Engine interrupted(args,config);setup(interrupted);interrupted.WarmCache(1);
    Require(interrupted.cache.size()==8,"torn imported index lost completed job");
    Engine corrupt(args,config);setup(corrupt);corrupt.WarmCache(1);
    const Job victim=corrupt.cache.begin()->second;const auto marker=victim.folder+"/complete.txt",saved=Read(marker);
    Atomic(marker,saved+"\n");bool rejected=false;
    try{corrupt.Materialize(victim);}catch(const std::exception&){rejected=true;}
    Atomic(marker,saved);Require(rejected,"journal accepted changed authoritative marker");
    // Unscoped legacy campaigns must migrate too. The first import leaves
    // every authoritative result untouched; the next restore uses the journal.
    const auto beforeMigration=FrozenParentFiles(parent.root);
    auto standalone=[&](Engine& e){e.root=parent.root;e.signature=parent.signature;e.exe=parent.exe;
        e.binaryHash=parent.binaryHash;e.controls=false;e.forbidNativeWork=true;e.ConfigureIndex(false);};
    Engine migrated(args,config);standalone(migrated);Require(migrated.compactIndex,"unscoped legacy root was not upgraded");migrated.WarmCache(1);
    Engine fast(args,config);standalone(fast);fast.WarmCache(1);
    const auto fastProgress=Read(fast.root+"/cache_restore_status.json");
    Require(fastProgress.find("\"journal_jobs\":8")!=std::string::npos&&fastProgress.find("\"fallback_folders\":0")!=std::string::npos,"unscoped resume reread old metadata");
    for(const auto& item:parent.cache){const auto a=parent.Materialize(item.second),b=fast.Materialize(fast.cache.at(item.first));Require(a.matrix==b.matrix,"unscoped journal changed all16");}
    const auto after=FrozenParentFiles(parent.root);for(const auto& item:beforeMigration){const auto found=after.find(item.first);Require(found!=after.end()&&found->second==item.second,"migration changed an existing result: "+item.first);}
    Engine wrong(args,config);wrong.root=root+"/wrong";wrong.signature="different";Directory(wrong.root);Atomic(wrong.root+"/configuration.txt","original\n");
    rejected=false;try{wrong.ConfigureIndex(false);}catch(const std::exception&){rejected=true;}
    Require(rejected&&!Exists(wrong.root+"/job_index_format.txt"),"incompatible legacy root was modified");
    std::cout<<"PASS legacy scoped index: persistent metadata, zero fallback reads on resume, unchanged all16 and parent files, torn tail recovered, changed source rejected\n";
    return 0;
}catch(const std::exception& error){std::cerr<<error.what()<<'\n';return 1;}}
