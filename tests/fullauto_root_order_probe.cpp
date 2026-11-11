#include "../src/FullAuto.cpp"
#include <thread>
using namespace FullAuto;

namespace {
void RootPlanNode(const std::string& folder,const std::vector<std::string>& children={}){
    Directory(folder+"/result");const std::string matrix="150 1 1 0.1 0 0 0.1 0.8 0 0 0 0 0.6 0 0 0 0 0.4\n";
    Atomic(folder+"/result/result.dat",matrix);std::string sources;
    for(const auto& child:children)sources+=child+'\t'+FileHash(child+"/complete.txt")+'\n';
    if(!children.empty())Atomic(folder+"/angular_sources.tsv",sources);
    Atomic(folder+"/complete.txt","node 1 "+Hash(matrix)+" none"+(children.empty()?"":" "+Hash(sources))+"\n");
}
template<class Action> void RejectRootPlan(Action action,const std::string& label){
    bool failed=false;try{action();}catch(const std::exception&){failed=true;}Require(failed,label);
}
}
int main(int argc,char** argv){try{
    Require(argc==2,"NEW_OUTPUT");const std::string root=argv[1];Directory(root+"/jobs");
    std::vector<std::string> leaves;for(int i=1;i<=32;++i){leaves.push_back(root+"/jobs/n"+PackedJobStore::Hex(i));RootPlanNode(leaves.back());}
    const auto a=root+"/jobs/a1",b=root+"/jobs/a2",c=root+"/jobs/a3";RootPlanNode(a,leaves);RootPlanNode(b,{a});RootPlanNode(c,{b});
    const std::vector<std::string> roots={a,b,c};size_t baseline=0;RootValidationStats stats;
    {ValidationBatch batch("old root order",true);ValidateStoredRootsPhase(roots,4,64*1024*1024);baseline=batch.checked;}
    {ValidationBatch batch("planned root order",true);bool planned=false;
        stats=ValidateStoredRoots(roots,4,[&](const RootValidationStats& plan){planned=true;Require(plan.ordered&&plan.firstRoots==1&&plan.secondRoots==2,"wrong planned phases");});
        Require(planned&&stats.checked<baseline&&stats.reused>=2,"shared descendants were not reused");Require(activeValidationBatch==&batch,"root plan scope leaked");}

    // A round memo may already cover the outer root. Corrupting a requested
    // descendant's direct payload must still fail in the SECOND phase. None
    // of this invocation's staged memo entries/counters may then be published.
    const auto original=Read(a+"/result/result.dat");
    {ValidationBatch batch("second phase failure",true);std::set<std::string> visited;ValidateStoredJob(c,visited);
        const auto bytes=batch.compounds.Bytes(),checked=batch.checked,reused=batch.reused;
        Atomic(a+"/result/result.dat",original+"corrupt\n");
        RejectRootPlan([&]{ValidateStoredRoots(roots,4);},"covered root payload skipped");
        Require(batch.compounds.Bytes()==bytes&&batch.checked==checked&&batch.reused==reused,"failed transaction changed parent memo");
        Require(activeValidationBatch==&batch,"failed transaction scope leaked");}
    Atomic(a+"/result/result.dat",original);
    {ValidationBatch batch("existing round memo",true);std::set<std::string> visited;ValidateStoredJob(a,visited);
        const auto result=ValidateStoredRoots(roots,4);Require(result.reused>2,"fallback chain lost existing round memo");}
    {ValidationBatch batch("duplicate requested roots",true);const auto result=ValidateStoredRoots({a,a},4);
        Require(!result.ordered&&result.firstRoots==2&&std::string(result.reason)=="duplicate roots","duplicates did not use original traversal");}

    // More roots than workers must reuse earlier waves without publishing
    // partially checked graphs to the enclosing round. These are independent
    // requested roots with a common descendant, so root ordering alone cannot
    // eliminate the duplicate walk.
    std::vector<std::string> waveRoots;
    for(int i=0;i<13;++i){waveRoots.push_back(root+"/jobs/w"+PackedJobStore::Hex(i+1));RootPlanNode(waveRoots.back(),{a});}
    size_t independentChecks=0;
    {ValidationBatch batch("independent baseline",true);ValidateStoredRootsPhase(waveRoots,4,64*1024*1024);independentChecks=batch.checked;}
    {ValidationBatch batch("bounded waves",true);size_t progressCalls=0,lastRoots=0;
        const auto owner=std::this_thread::get_id();
        const auto result=ValidateStoredRoots(waveRoots,4,{},[&](const RootValidationStats& progress){
            Require(std::this_thread::get_id()==owner,"progress escaped controller thread");
            Require(batch.checked==0&&batch.reused==0&&batch.compounds.Bytes()==0,"progress published parent prematurely");
            Require(progress.validatedRoots>lastRoots&&progress.validatedRoots-lastRoots<=progress.waveWidth,"invalid wave progress");
            lastRoots=progress.validatedRoots;++progressCalls;
        },0); // Force the bounded planning-budget fallback, without huge files.
        Require(result.incompletePlan&&result.boundedWaves,"planning limit did not enable bounded waves");
        Require(result.checked<independentChecks&&result.validatedRoots==waveRoots.size(),"waves did not reuse common sources");
        Require(progressCalls==result.plannedWaves&&progressCalls==result.completedWaves,"missing wave progress");
        for(const auto& path:waveRoots)Require(bool(batch.compounds.Peek(path)),"earlier requested root was not retained");
    }
    {ValidationBatch batch("complete independent plan",true);const auto result=ValidateStoredRoots(waveRoots,4);
        Require(!result.incompletePlan&&!result.boundedWaves&&result.checked==independentChecks,
            "complete independent roots lost their original scheduling");}
    auto memoDigest=[](const ValidationBatch& batch){std::string text;
        batch.compounds.LeastRecentFirst([&](const std::string& key,std::shared_ptr<const std::string> value){text+=key+' '+*value+'\n';});
        return Hash(text);
    };
    const auto lastMatrix=Read(waveRoots.back()+"/result/result.dat");
    for(int threads:{1,4}){
        ValidationBatch batch("late wave failure",true);std::set<std::string> seen;ValidateStoredJob(a,seen);
        const auto digest=memoDigest(batch);const auto checked=batch.checked,reused=batch.reused;size_t progressCalls=0;
        RejectRootPlan([&]{ValidateStoredRoots(waveRoots,threads,{},[&](const RootValidationStats&){
            if(++progressCalls==1)Atomic(waveRoots.back()+"/result/result.dat",lastMatrix+"corrupt\n");
        },0);},"late corrupt root accepted");
        Require(progressCalls>1&&batch.checked==checked&&batch.reused==reused&&memoDigest(batch)==digest,
            "late wave failure changed parent memo, LRU or counters");
        Require(activeValidationBatch==&batch,"late wave failure leaked scope");
        Atomic(waveRoots.back()+"/result/result.dat",lastMatrix);
    }
    {ValidationBatch batch("progress failure",true);const auto digest=memoDigest(batch);
        RejectRootPlan([&]{ValidateStoredRoots(waveRoots,4,{},[&](const RootValidationStats&){
            throw std::runtime_error("progress sink failed");
        },0);},"progress failure ignored");
        Require(batch.checked==0&&batch.reused==0&&memoDigest(batch)==digest,"failed progress published parent memo");
    }

    // Partial-known DAG plus waves: root 12 is an unread root that is already
    // known to be a child of root 0. Unknown-first planning is only an ordering
    // hint; root 12 and the covered root 0 both require their own direct checks.
    {
        const auto partial=root+"/partial_known";std::vector<std::string> partialRoots;
        for(int i=0;i<13;++i)partialRoots.push_back(partial+"/jobs/r"+PackedJobStore::Hex(i));
        RootPlanNode(partialRoots[12]);RootPlanNode(partialRoots[0],{partialRoots[12]});
        for(int i=1;i<12;++i)RootPlanNode(partialRoots[i],{partialRoots[i-1]});
        const size_t budget=Read(partialRoots[0]+"/angular_sources.tsv").size()+Read(partialRoots[1]+"/angular_sources.tsv").size();
        const auto plan=PlanValidationRoots(partialRoots,[&](size_t i,size_t available,std::string& text){
            return ReadValidationPlanManifest(partialRoots[i]+"/angular_sources.tsv",available,text);
        },512,budget);
        Require(plan.incomplete&&plan.Ordered()&&plan.edges==2&&plan.first.size()==12&&plan.second==std::vector<size_t>({0}),"partial-known plan lost discovered links");
        Require(std::find(plan.first.begin(),plan.first.end(),12)!=plan.first.end(),"unread known child was omitted");
        auto selected=plan.first;selected.insert(selected.end(),plan.second.begin(),plan.second.end());std::sort(selected.begin(),selected.end());
        for(size_t i=0;i<partialRoots.size();++i)Require(selected.at(i)==i,"partial-known plan duplicated or skipped a root");
        Require(selected.size()==partialRoots.size(),"partial-known plan changed root count");
        auto partialMemoDigest=[](const ValidationBatch& batch){std::string text;
            batch.compounds.LeastRecentFirst([&](const std::string& key,std::shared_ptr<const std::string> value){text+=key+' '+*value+'\n';});return Hash(text);};
        for(int threads:{1,4}){
            {ValidationBatch batch("successful partial-known waves",true);size_t progressed=0;
                const auto result=ValidateStoredRoots(partialRoots,threads,{},[&](const RootValidationStats& progress){
                    Require(batch.checked==0&&batch.reused==0&&batch.compounds.Bytes()==0,"partial-known progress published parent");
                    Require(progress.validatedRoots>progressed,"partial-known progress did not advance");progressed=progress.validatedRoots;
                },budget);
                Require(result.incompletePlan&&result.ordered&&result.boundedWaves&&result.edges==2,"partial-known wave policy changed");
                Require(result.validatedRoots==partialRoots.size()&&progressed==partialRoots.size()&&result.completedWaves==result.plannedWaves,"partial-known waves omitted roots");}
            for(size_t bad:{size_t(0),size_t(12)}){
                ValidationBatch batch("late partial-known failure",true);std::set<std::string> seen;ValidateStoredJob(partialRoots[11],seen);
                const auto before=partialMemoDigest(batch),matrix=Read(partialRoots[bad]+"/result/result.dat");const auto checked=batch.checked,reused=batch.reused;
                size_t progressCalls=0,lastRoots=0;
                RejectRootPlan([&]{ValidateStoredRoots(partialRoots,threads,{},[&](const RootValidationStats& progress){
                    lastRoots=progress.validatedRoots;if(++progressCalls==1)Atomic(partialRoots[bad]+"/result/result.dat",matrix+"corrupt\n");
                },budget);},"partial-known covered/unread root payload escaped validation");
                Require(progressCalls>1&&batch.checked==checked&&batch.reused==reused&&partialMemoDigest(batch)==before,"late partial-known failure published memo/LRU/counters");
                if(bad==0)Require(lastRoots==12,"covered root did not reach the second phase");
                else Require(lastRoots<12,"unread known child was not checked in first phase");
                Require(activeValidationBatch==&batch,"partial-known exception leaked validation scope");Atomic(partialRoots[bad]+"/result/result.dat",matrix);
            }
        }
        // Even an existing parent memo for the entire chain cannot authorize
        // skipping the direct payload of ANY of the 13 explicitly requested roots.
        for(size_t bad=0;bad<partialRoots.size();++bad){
            ValidationBatch batch("all partial-known root payloads",true);std::set<std::string> seen;ValidateStoredJob(partialRoots[11],seen);
            const auto before=partialMemoDigest(batch),matrix=Read(partialRoots[bad]+"/result/result.dat");const auto checked=batch.checked,reused=batch.reused;
            Atomic(partialRoots[bad]+"/result/result.dat",matrix+"corrupt\n");
            RejectRootPlan([&]{ValidateStoredRoots(partialRoots,4,{},{},budget);},"partial-known requested root payload skipped");
            Require(batch.checked==checked&&batch.reused==reused&&partialMemoDigest(batch)==before,"partial-known failure changed caller memo");
            Atomic(partialRoots[bad]+"/result/result.dat",matrix);
        }
    }

    std::string text="old";Require(ReadValidationPlanManifest(root+"/absent",0,text)&&text.empty(),"missing optional manifest");
    const auto large=root+"/large_manifest";int fd=open(large.c_str(),O_CREAT|O_TRUNC|O_WRONLY,0600);
    Require(fd>=0&&ftruncate(fd,8*1024*1024+1)==0,"large fixture");close(fd);
    Require(!ReadValidationPlanManifest(large,8*1024*1024,text)&&text.empty(),"reader allocated oversized manifest");
    const auto pipe=root+"/manifest_fifo";Require(mkfifo(pipe.c_str(),0600)==0,"fifo fixture");
    Require(!ReadValidationPlanManifest(pipe,1024,text),"nonregular planning input accepted");

    // Bounded planning must also work through the immutable packed view and
    // an interrupted loose-copy cleanup; it never authenticates acceptance.
    PackedJobStore::Pack(a);
    {ValidationBatch batch("packed plan",true);const auto result=ValidateStoredRoots(roots,4);Require(result.ordered&&result.firstRoots==1,"packed adjacency lost");}
    PackedJobStore::Extract(a);Require(unlink((a+"/angular_sources.tsv").c_str())==0,"partial cleanup fixture");
    {ValidationBatch batch("partial packed plan",true);Require(ValidateStoredRoots(roots,4).ordered,"partial packed adjacency lost");}
    Require(Read(a+"/result/result.dat")==original,"root matrix changed");
    PackedJobStore::Member member;Require(PackedJobStore::Locate(a+"/angular_sources.tsv",member),"packed sources missing");
    fd=open(member.archive.c_str(),O_RDWR);const char bad='X';Require(fd>=0&&pwrite(fd,&bad,1,member.offset)==1,"corrupt fixture");close(fd);
    Require(!ReadValidationPlanManifest(a+"/angular_sources.tsv",8*1024*1024,text),"corrupt plan member read");
    {ValidationBatch batch("corrupt packed source",true);RejectRootPlan([&]{ValidateStoredRoots(roots,4);},"planner fallback hid source corruption");Require(batch.compounds.Bytes()==0,"failed fallback published memo");}
    Require(activeValidationBatch==nullptr,"final scope leaked");
    std::cout<<"PASS root validation ordering: checked="<<baseline<<"->"<<stats.checked
        <<"; every root/direct payload checked; existing round memo; bounded limited-plan waves; late-wave/progress rollback; partial-known links and unread known child direct checks; controller-thread progress; complete scheduling retained; bounded reader; packed/partial fallback; corrupt sources rejected\n";
    return 0;
}catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}}
