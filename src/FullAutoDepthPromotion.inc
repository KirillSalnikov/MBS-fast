// Restartable metadata transition after a failed final depth audit. Native
// worker data and scientific matrices are never removed or rewritten here.
struct DepthPromotionPlan {
    std::string phase,signature;int previous=0,next=0;
    std::vector<std::pair<std::string,std::string>> files;
};
const std::vector<std::string>& DepthPromotionFiles(){
    static const std::vector<std::string> names={"adaptive_checkpoint.tsv","go_checkpoint.tsv","execution_tuning.tsv","physics_validation.json","fullauto_status.json","calibration.json","initial_calibration.json","depth_selection.json","automatic_depth_policy.tsv"};
    return names;
}
std::string DepthPromotionText(const DepthPromotionPlan& plan){
    std::ostringstream out;out<<1<<' '<<plan.phase<<' '<<plan.signature<<' '<<plan.previous<<' '<<plan.next<<' '<<plan.files.size()<<'\n';
    for(const auto& file:plan.files)out<<file.first<<' '<<file.second<<'\n';
    return Hash(out.str())+"\n"+out.str();
}
DepthPromotionPlan ReadDepthPromotionPlan(Engine& e){
    const auto text=Read(e.root+"/depth_promotion_pending.tsv");const auto split=text.find('\n');
    Require(split!=std::string::npos&&Hash(text.substr(split+1))==text.substr(0,split),"depth promotion pending receipt changed");
    std::istringstream in(text.substr(split+1));DepthPromotionPlan plan;int version=0;size_t count=0;
    Require(bool(in>>version>>plan.phase>>plan.signature>>plan.previous>>plan.next>>count)
        &&version==1&&(plan.phase=="archive"||plan.phase=="publish")&&plan.signature==e.signature
        &&plan.previous>=e.depthStart&&plan.previous<e.depthMaximum&&plan.next==std::min(e.depthMaximum,plan.previous+2)
        &&count<=DepthPromotionFiles().size(),"invalid depth promotion pending receipt");
    std::set<std::string> names;std::string name,digest;
    for(size_t i=0;i<count;++i){
        Require(bool(in>>name>>digest)&&std::find(DepthPromotionFiles().begin(),DepthPromotionFiles().end(),name)!=DepthPromotionFiles().end()
            &&names.insert(name).second&&!digest.empty(),"invalid depth promotion member");
        plan.files.push_back({name,digest});
    }
    Require(!(in>>name),"trailing depth promotion data");return plan;
}
int DepthPromotionFloor(Engine& e){
    const auto path=e.root+"/depth_promotion_floor.tsv";if(!Exists(path))return e.depthStart;
    const auto text=Read(path);const auto split=text.find('\n');
    Require(split!=std::string::npos&&Hash(text.substr(split+1))==text.substr(0,split),"depth promotion floor changed");
    std::istringstream in(text.substr(split+1));std::string signature,extra;int floor=0;
    Require(bool(in>>signature>>floor)&&!(in>>extra)&&signature==e.signature&&floor>=e.depthStart&&floor<=e.depthMaximum,
        "invalid depth promotion floor");return floor;
}
void SyncDepthPromotionDirectory(const std::string& path){
    if(MemoryPath(path))return;
    const int fd=open(path.c_str(),O_RDONLY|O_DIRECTORY);Require(fd>=0,"cannot open depth promotion directory");
    int result;do{result=fsync(fd);}while(result<0&&errno==EINTR);close(fd);Require(result==0,"cannot flush depth promotion directory");
}
void BeginDepthPromotion(Engine& e,int previous,int next){
    Require(!Exists(e.root+"/depth_promotion_pending.tsv"),"unfinished depth promotion transaction");
    Require(previous>=DepthPromotionFloor(e)&&previous<e.depthMaximum&&next==std::min(e.depthMaximum,previous+2),"invalid depth promotion transition");
    Require(!Exists(e.root+"/final_outputs.fnv"),"cannot promote an already published final estimator");
    DepthPromotionPlan plan;plan.phase="archive";plan.signature=e.signature;plan.previous=previous;plan.next=next;
    const auto archive=e.root+"/depth_attempts/n"+Number(previous);Directory(archive);
    for(const auto& name:DepthPromotionFiles()){
        Require(!Exists(archive+"/"+name),"depth promotion archive already contains "+name);
        if(Exists(e.root+"/"+name))plan.files.push_back({name,FileHash(e.root+"/"+name)});
    }
    Atomic(e.root+"/depth_promotion_pending.tsv",DepthPromotionText(plan));
}
bool ResumeDepthPromotion(Engine& e){
    const auto pending=e.root+"/depth_promotion_pending.tsv";if(!Exists(pending))return false;
    auto plan=ReadDepthPromotionPlan(e);const auto archive=e.root+"/depth_attempts/n"+Number(plan.previous);Directory(archive);
    Require(DepthPromotionFloor(e)<=plan.next,"stale depth promotion would lower the confirmed floor");
    // Validate the entire bounded manifest before moving anything. A conflict
    // leaves every remaining source and every existing archived file intact.
    for(const auto& file:plan.files){
        const auto source=e.root+"/"+file.first,destination=archive+"/"+file.first;
        const bool from=Exists(source),to=Exists(destination);
        if(plan.phase=="archive"){
            Require(from!=to,"missing or conflicting depth promotion member "+file.first);
            Require(FileHash(from?source:destination)==file.second,"depth promotion member changed "+file.first);
        }else Require(to&&FileHash(destination)==file.second,"published depth promotion archive changed "+file.first);
    }
    if(plan.phase=="archive"){
        for(const auto& file:plan.files){const auto source=e.root+"/"+file.first,destination=archive+"/"+file.first;
            if(Exists(source))Require(RenamePath(source,destination)==0,"cannot archive previous depth attempt");}
        SyncDepthPromotionDirectory(archive);SyncDepthPromotionDirectory(e.root+"/depth_attempts");SyncDepthPromotionDirectory(e.root);
        plan.phase="publish";Atomic(pending,DepthPromotionText(plan));
    }
    Require(!Exists(e.root+"/adaptive_checkpoint.tsv")&&!Exists(e.root+"/go_checkpoint.tsv")&&!Exists(e.root+"/final_outputs.fnv"),
        "depth promotion publish conflicts with an existing estimator");
    const auto floor=e.signature+" "+Number(plan.next)+"\n";
    Atomic(e.root+"/depth_promotion_floor.tsv",Hash(floor)+"\n"+floor);
    SaveSelectedDepth(e,plan.next);
    e.Record(e.root+"/depth_selection.json","{\"automatic\":true,\"target_depth\":"+Number(plan.next)+",\"reason\":\"promoted after fresh final depth validation\",\"previous\":"+Number(plan.previous)+",\"maximum\":"+Number(e.depthMaximum)+",\"budget\":"+Number(e.depthBudget)+"}");
    e.Log("physical_depth_promoted","{\"previous\":"+Number(plan.previous)+",\"next\":"+Number(plan.next)+",\"reason\":\"fresh final depth audit failed\",\"archive\":"+Quoted(archive)+",\"restartable_transaction\":true}");
    Require(RemovePath(pending)==0,"cannot finish depth promotion transaction");SyncDepthPromotionDirectory(e.root);
    return true;
}
