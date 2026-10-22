// Projection/import checks on frozen real production checkpoints. Native jobs
// are represented only by a completion marker; this probe never runs optics.
#include "../src/FullAuto.cpp"
using namespace FullAuto;

void SetupWindow(Engine& e,const std::string& parent,int initial){
    e.root=parent;e.signature=Read(parent+"/configuration.txt");
    if(!e.signature.empty()&&e.signature.back()=='\n')e.signature.pop_back();
    e.theta=ReadGrid(parent+"/requested_theta.csv");e.initial=initial;e.minimum=256;
    e.maxCount=67108864;e.minPhi=1;e.maxPhi=2147483647;e.eps=.02;e.interpolation=.005;
    e.autoDepth=true;e.depthStart=6;e.depthMaximum=18;e.gpuCandidates={"0","1","2","3"};
    e.refinementWindow=true;e.refinementLo=170;e.refinementHi=176;e.refinementEps=.05;
    e.forbidNativeWork=true;e.exe="/bin/false";
}
AdaptiveState EmptyWindowState(const Engine& e){
    AdaptiveState state;state.counts=Budgets(3,std::vector<int>(e.theta.size(),e.minimum));
    state.counts[0].assign(e.theta.size(),e.initial);state.previous=Matrix(e.theta.size());
    state.valid.assign(e.theta.size(),false);state.coarse.assign(e.theta.size(),false);return state;
}
std::map<std::string,std::string> FrozenParentFiles(const std::string& root){
    std::map<std::string,std::string> result;DIR* dir=opendir(root.c_str());Require(dir,"parent fixture missing");
    while(auto* entry=readdir(dir)){
        if(entry->d_name[0]=='.')continue;const std::string path=root+"/"+entry->d_name;struct stat st;
        Require(lstat(path.c_str(),&st)==0,"parent stat failed");
        if(S_ISDIR(st.st_mode)){const auto nested=FrozenParentFiles(path);result.insert(nested.begin(),nested.end());}
        else if(S_ISREG(st.st_mode))result[path]=FileHash(path)+" "+Number(st.st_mtim.tv_sec)+" "+Number(st.st_mtim.tv_nsec);
        else Require(false,"parent fixture unexpectedly contains symlink");
    }closedir(dir);return result;
}
void SetupCacheEngine(Engine& e,const std::string& root){
    e.root=root;e.signature="read-only-parent-cache-fixture";e.exe="/bin/false";e.binaryHash=FileHash(e.exe);
    e.controls=false;e.incremental=true;e.angularCache=true;e.forbidNativeWork=true;e.orientationBlock=512;
    e.phi=4;e.gpus={""};e.workerThreads=1;e.theta={170,176};e.seeds.assign(production,production+8);
    Directory(root);Directory(root+"/jobs");e.ConfigureIndex(true);
}
void ReadOnlyParentCache(const std::string& root){
    ArgPP args;RunConfig config;Engine parent(args,config);SetupCacheEngine(parent,root+"/parent");
    for(size_t s=0;s<8;++s)for(int offset:{0,512}){
        Job j;parent.SetSampling(j,512,parent.seeds[s],8,parent.theta,false);j.offset=offset;
        j.identity="native-block-"+Number(s)+"-"+Number(offset);j.folder=parent.root+"/jobs/"+j.identity;
        j.data=j.folder+"/result/result.dat";j.summary=j.folder+"/result/result_analytic_facets.tsv";
        j.weightHash=Hash("");j.seconds=1;j.angles=parent.theta;j.matrix=Matrix(2);
        for(size_t row=0;row<2;++row){j.matrix[row][0]=100+s+row+double(offset)/512;j.matrix[row][5]=20+s;}
        Directory(j.folder+"/result");WriteMatrix(j.data,j.angles,j.matrix);
        Atomic(j.folder+"/stdout.log","Hard tree-limit hits: 0\nOrientations still incomplete: 0\nTotal duration: 1 s\n");
        parent.PublishComplete(j,j.identity+" 1 "+FileHash(j.data)+" none\n");parent.CacheJob(j);
    }
    // Last seed has all native blocks but no assembled prefix in the parent.
    const auto seven=parent.PrefixStage(1024,parent.seeds.data(),7,8,parent.theta);parent.journal.flush();
    Require(seven.size()==7,"prefix fixture incomplete");const auto frozen=FrozenParentFiles(parent.root);
    Engine child(args,config);SetupCacheEngine(child,root+"/child");child.refinementParent=parent.root;
    child.WarmCache(1);Require(child.cache.size()==23,"child did not restore all parent metadata");
    for(const auto& j:seven)Require(child.cache.at(j.identity).folder==j.folder,"parent prefix redirected into child");
    const auto eight=child.PrefixStage(1024,child.seeds.data(),8,8,child.theta);
    Require(eight.size()==8,"partial prefix continuation failed");
    for(size_t s=0;s<8;++s){
        for(size_t row=0;row<2;++row)Require(eight[s].matrix[row][0]==100+s+row+.5,"prefix weighted mean changed");
        if(s<7)Require(eight[s].folder==seven[s].folder,"cached parent prefix replaced");
        else Require(eight[s].folder.find(child.root+"/jobs/")==0,"missing prefix was not created in child");
    }
    Require(FrozenParentFiles(parent.root)==frozen,"partial prefix changed parent payload, index or mtime");
    // A child's compact journal may refer to immutable parent directories.
    // On resume, choose the child directory only for completed local results.
    Engine resumed(args,config);SetupCacheEngine(resumed,child.root);resumed.refinementParent=parent.root;
    resumed.WarmCache(1);Require(resumed.cache.size()==24,"child journal lost inherited cache metadata");
    Require(resumed.cache.at(eight[7].identity).folder==eight[7].folder,"child prefix lost local path on resume");
    for(size_t s=0;s<7;++s)Require(resumed.cache.at(eight[s].identity).folder==seven[s].folder,"parent journal path resolution changed");
    Require(FrozenParentFiles(parent.root)==frozen,"resume wrote to frozen parent");
    Atomic(root+"/result.json","{\"status\":\"passed\",\"native_blocks\":16,\"inherited_complete_prefixes\":7,\"new_child_prefixes\":1,\"native_work\":false,\"parent_sha_and_mtimes_unchanged\":true,\"compact_journal_resume\":true}\n");
}

int main(int argc,char** argv){try{
    Require(argc==3,"usage: probe parent initial");const std::string parent=argv[1];const int initial=std::stoi(argv[2]);
    ArgPP args;RunConfig config;Engine original(args,config);SetupWindow(original,parent,initial);
    original.target=RefinementSelectedDepth(original);Require(original.RestoreTuning(),"missing saved tuning");
    auto old=EmptyWindowState(original);Require(RestoreAdaptive(original,old,false),"missing saved parent state");
    Require(original.controls&&original.sharedMeans,"fixture must exercise analytic shared means");
    const auto checkpoint=Read(parent+"/adaptive_checkpoint.tsv");const auto means=Read(parent+"/means_shared.cache");
    const auto weightFile=Read(parent+"/polarization_weights.tsv");
    const auto oldAngles=original.theta;const auto oldPhi=original.thetaPhi;
    const auto weights=original.weights;const auto polar=original.polarWeights;const auto wx=original.weightTheta;
    Engine child(args,config);SetupWindow(child,parent,initial);PrepareRefinementWindow(child);
    Require(child.eps==.05&&child.interpolation==.0125&&child.depthBudget==.0125,"wrong child accuracy budget");
    Require(child.theta.front()==170&&child.theta.back()==176,"wrong child range");
    Require(child.weights==weights&&child.polarWeights==polar&&child.weightTheta==wx,"analytic coefficients changed");
    Require(Read(child.root+"/polarization_weights.tsv")==weightFile,"weight file changed");
    Require(Read(child.root+"/means_shared.cache")==means,"authoritative analytic means changed during import");
    auto projected=EmptyWindowState(child);Require(RestoreAdaptive(child,projected,false),"missing projected checkpoint");
    Require(projected.round==0&&projected.streak==0&&projected.auditSerial==old.auditSerial,"old acceptance/audit handling incorrect");
    Require(std::none_of(projected.coarse.begin(),projected.coarse.end(),[](bool v){return v;}),"old quiet-half acceptance retained");
    size_t kept=0,newEndpoints=0;std::set<int> rules;
    for(size_t row=0;row<child.theta.size();++row){
        auto it=std::lower_bound(oldAngles.begin(),oldAngles.end(),child.theta[row]);const size_t j=it-oldAngles.begin();
        if(it!=oldAngles.end()&&*it==child.theta[row]){
            for(size_t level=0;level<3;++level)Require(projected.counts[level][row]==old.counts[level][j],"existing orientation prefix changed");
            Require(projected.previous[row]==old.previous[j]&&projected.valid[row]==old.valid[j],"saved matrix validity changed");
            Require(child.PhiForTheta(row)==original.PhiForTheta(j),"existing phi rule changed");++kept;
        }else{
            Require(j>0&&j<oldAngles.size(),"new endpoint outside parent");
            for(size_t level=0;level<3;++level)Require(projected.counts[level][row]==std::max(old.counts[level][j-1],old.counts[level][j]),"new endpoint underestimated count");
            Require(child.PhiForTheta(row)==std::max(original.PhiForTheta(j-1),original.PhiForTheta(j)),"new endpoint underestimated phi");
            Require(!projected.valid[row],"new endpoint marked directly calculated");++newEndpoints;
        }
        rules.insert(child.PhiForTheta(row));
    }
    Require(kept+newEndpoints==child.theta.size()&&newEndpoints==1,"expected only new exact 176 endpoint");
    auto reference=oldAngles;reference.push_back(176);std::sort(reference.begin(),reference.end());reference.erase(std::unique(reference.begin(),reference.end()),reference.end());
    Require(ReadGrid(child.MeanReferenceGrid())==reference,"analytic reference lost full parent union");
    const auto childRoot=child.root;const auto receipt=Read(childRoot+"/refinement_window_policy.json");
    // Simulate committed progress in the isolated controller, then invoke the
    // same startup path again. It must never import round0 over this state.
    projected.round=7;projected.streak=1;projected.coarse[0]=true;projected.counts[0][0]=std::min(child.maxCount,2*projected.counts[0][0]);
    SaveAdaptive(child,projected);const auto progress=Read(childRoot+"/adaptive_checkpoint.tsv");
    Engine resumed(args,config);SetupWindow(resumed,parent,initial);PrepareRefinementWindow(resumed);
    Require(resumed.root==childRoot&&Read(childRoot+"/adaptive_checkpoint.tsv")==progress,"resume overwrote child progress");
    auto restored=EmptyWindowState(resumed);Require(RestoreAdaptive(resumed,restored,false),"child resume unavailable");
    Require(restored.round==7&&restored.streak==1&&restored.coarse[0]&&restored.counts==projected.counts,"child progress not retained");
    Require(Read(childRoot+"/refinement_window_policy.json")==receipt,"resume changed scope policy");
    Require(Read(parent+"/adaptive_checkpoint.tsv")==checkpoint&&Read(parent+"/means_shared.cache")==means&&Read(parent+"/polarization_weights.tsv")==weightFile,"parent mutated");
    ReadOnlyParentCache(parent+"/../readonly_cache");
    std::cout<<"{\"status\":\"passed\",\"parent_rows\":"<<oldAngles.size()<<",\"child_rows\":"<<child.theta.size()<<",\"kept_rows\":"<<kept
        <<",\"new_endpoint_rows\":"<<newEndpoints<<",\"phi_rules\":[";bool first=true;for(int phi:rules){std::cout<<(first?"":",")<<phi;first=false;}
    std::cout<<"],\"parent_round\":"<<old.round<<",\"child_resume_round\":"<<restored.round<<",\"epsilon\":0.05,\"interpolation_budget\":0.0125,\"shared_means_V3\":true,\"parent_unchanged\":true,\"native_work\":false}\n";
    return 0;
}catch(const std::exception& error){std::cerr<<error.what()<<'\n';return 1;}}
