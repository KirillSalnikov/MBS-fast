#include "../src/FullAuto.cpp"
using namespace FullAuto;

namespace {
template<class Action> void RejectMarker(Action action,const std::string& label){
    bool rejected=false;try{action();}catch(const std::exception&){rejected=true;}
    Require(rejected,"source marker accepted: "+label);
}
void ValidateFresh(const std::string& folder){
    ValidationBatch batch("marker regression",true);std::set<std::string> visited;ValidateStoredJob(folder,visited);
}
void MarkerLeaf(const std::string& folder){
    Directory(folder+"/result");
    const std::string data="150 1 1 0.1 0 0 0.1 0.8 0 0 0 0 0.6 0 0 0 0 0.4\n";
    const std::string integrals="Cext\tprojected_area\n1.2\t1\n";
    Atomic(folder+"/result/result.dat",data);Atomic(folder+"/result/result_integrals.tsv",integrals);
    const auto manifest=folder+"/result/result_integrals.tsv "+Hash(integrals)+"\n";
    Atomic(folder+"/integral_outputs.fnv",manifest);
    Atomic(folder+"/complete.txt","native 1 "+Hash(data)+" none none "+Hash(manifest)+"\n");
}
void MarkerParent(const std::string& folder,const std::string& sources){
    Directory(folder+"/result");const std::string data="150 1 1 0.1 0 0 0.1 0.8 0 0 0 0 0.6 0 0 0 0 0.4\n";
    Atomic(folder+"/result/result.dat",data);Atomic(folder+"/angular_sources.tsv",sources);
    Atomic(folder+"/complete.txt","parent 1 "+Hash(data)+" none "+Hash(sources)+"\n");
}
void MarkerPaired(const std::string& folder,bool withIntegrals){
    MarkerLeaf(folder);const auto data=Read(folder+"/result/result.dat");
    const auto path=folder+"/result/result_depth6.dat";Atomic(path,data);
    const auto pair=path+" "+Hash(data)+"\n";Atomic(folder+"/paired_depth.tsv",pair);
    Atomic(folder+"/complete.txt","paired 1 "+Hash(data)+" none "+Hash(pair)
        +(withIntegrals?" "+Hash(Read(folder+"/integral_outputs.fnv")):"")+"\n");
}
void MarkerBalancedPaired(const std::string& folder,const std::string& sources){
    MarkerPaired(folder,true);const auto data=Read(folder+"/result/result.dat");
    const std::string cutoff="6 1 source\n";Atomic(folder+"/cutoff_sources.tsv",cutoff);
    const auto integrity=Read(folder+"/integral_outputs.fnv")
        +folder+"/paired_depth.tsv "+Hash(Read(folder+"/paired_depth.tsv"))+"\n"
        +folder+"/cutoff_sources.tsv "+Hash(cutoff)+"\n";
    Atomic(folder+"/angular_sources.tsv",sources);Atomic(folder+"/integral_outputs.fnv",integrity);
    Atomic(folder+"/complete.txt","balanced 1 "+Hash(data)+" none "+Hash(sources)+" "+Hash(integrity)+"\n");
}
}

int main(int argc,char** argv){try{
    Require(argc==2,"OUTPUT");const std::string root=argv[1],leaf=root+"/jobs/n1",parent=root+"/jobs/a1";
    Directory(root+"/jobs");MarkerLeaf(leaf);
    const auto marker=Read(leaf+"/complete.txt"),expected=Hash(marker),edge=leaf+"\t"+expected+"\n";
    MarkerParent(parent,edge);ValidateFresh(parent);const auto savedMatrix=Read(parent+"/result/result.dat");

    // Parent hashes are correct in these malformed-edge fixtures. Rejection
    // must come from the empty child digest, not a stale parent manifest hash.
    MarkerParent(parent,leaf+"\t\n");RejectMarker([&]{ValidateFresh(parent);},"empty first edge");
    MarkerParent(parent,edge+leaf+"\t\n");RejectMarker([&]{ValidateFresh(parent);},"empty duplicate edge");
    MarkerParent(parent,edge+leaf+"\t"+expected+"0\n");RejectMarker([&]{ValidateFresh(parent);},"wrong duplicate digest");
    MarkerParent(parent,edge+edge);ValidateFresh(parent);

    Atomic(leaf+"/complete.txt",marker+"changed\n");
    RejectMarker([&]{ValidateFresh(parent);},"changed first-visit marker");
    {std::set<std::string> visited{leaf};RejectMarker([&]{ValidateStoredJob(leaf,visited,expected);},"changed visited marker");}
    Atomic(leaf+"/complete.txt",marker);ValidateFresh(parent);

    {ValidationBatch batch("memo marker regression",true);std::set<std::string> first;ValidateStoredJob(parent,first);
        Require(bool(batch.compounds.Peek(parent)),"parent graph not memoized");
        Atomic(leaf+"/complete.txt",marker+"changed\n");
        RejectMarker([&]{std::set<std::string> next;ValidateStoredJob(parent,next);},"memo hit skipped direct marker");
        Atomic(leaf+"/complete.txt",marker);
        Atomic(parent+"/angular_sources.tsv",leaf+"\t\n");
        RejectMarker([&]{std::set<std::string> next;ValidateStoredJob(parent,next);},"memo hit skipped source manifest hash");
        MarkerParent(parent,edge+edge);
    }
    const auto data=Read(leaf+"/result/result.dat");Atomic(leaf+"/result/result.dat",data+"changed\n");
    RejectMarker([&]{ValidateFresh(parent);},"child payload hash");Atomic(leaf+"/result/result.dat",data);
    const auto integrals=Read(leaf+"/result/result_integrals.tsv");Atomic(leaf+"/result/result_integrals.tsv",integrals+"changed\n");
    RejectMarker([&]{ValidateFresh(parent);},"integral output hash");Atomic(leaf+"/result/result_integrals.tsv",integrals);

    // A real fifth-field digest requires its manifest even when the main DAT
    // and completion marker are intact, or the ancestry is already memoized.
    {const auto folder=root+"/jobs/missing-angular";MarkerParent(folder,edge);ValidateFresh(folder);
        Require(unlink((folder+"/angular_sources.tsv").c_str())==0,"cannot remove raw angular manifest");
        RejectMarker([&]{ValidateFresh(folder);},"missing raw angular manifest");
        Atomic(leaf+"/result/result.dat",data+"changed\n");
        RejectMarker([&]{ValidateFresh(folder);},"missing angular manifest hides corrupt dependency");
        Atomic(leaf+"/result/result.dat",data);MarkerParent(folder,edge);
        ValidationBatch batch("missing manifest memo",true);std::set<std::string> seen;ValidateStoredJob(folder,seen);
        Require(bool(batch.compounds.Peek(folder)),"missing manifest fixture not memoized");
        Require(unlink((folder+"/angular_sources.tsv").c_str())==0,"cannot remove memo angular manifest");
        RejectMarker([&]{std::set<std::string> next;ValidateStoredJob(folder,next);},"memo missing angular manifest");
    }
    // Four-field native checkpoints, explicit none, and paired five/six-field
    // native markers remain valid. Native pairs authenticate paired_depth.tsv.
    {const auto folder=root+"/jobs/legacy";MarkerLeaf(folder);
        for(const auto& suffix:std::vector<std::string>{""," none"," none "+Hash(Read(folder+"/integral_outputs.fnv"))}){
            Atomic(folder+"/complete.txt","native 1 "+Hash(Read(folder+"/result/result.dat"))+" none"+suffix+"\n");ValidateFresh(folder);}
        for(bool withIntegrals:{false,true}){const auto paired=root+"/jobs/paired"+Number(withIntegrals);MarkerPaired(paired,withIntegrals);ValidateFresh(paired);
            Require(unlink((paired+"/paired_depth.tsv").c_str())==0,"cannot remove paired manifest");
            RejectMarker([&]{ValidateFresh(paired);},"missing native paired manifest");}
    }
    // Balanced paired jobs have both manifests, but their fifth-field digest
    // belongs to angular_sources. A surviving paired file cannot replace it.
    {const auto folder=root+"/jobs/balanced-paired";MarkerBalancedPaired(folder,edge);ValidateFresh(folder);
        ValidationBatch batch("balanced missing manifest memo",true);std::set<std::string> seen;ValidateStoredJob(folder,seen);
        Require(unlink((folder+"/angular_sources.tsv").c_str())==0,"cannot remove balanced angular manifest");
        RejectMarker([&]{std::set<std::string> next;ValidateStoredJob(folder,next);},"balanced paired memo missing angular manifest");
        RejectMarker([&]{ValidateFresh(folder);},"balanced paired missing angular manifest");
    }
    // Existence includes a verified packed member. Raw cleanup must therefore
    // remain accepted for both source-manifest formats.
    {const auto angular=root+"/jobs/a90",paired=root+"/jobs/n90";
        MarkerParent(angular,edge);MarkerPaired(paired,true);
        for(const auto& item:std::vector<std::pair<std::string,std::string>>{{angular,"angular_sources.tsv"},{paired,"paired_depth.tsv"}}){
            PackedJobStore::Pack(item.first);ValidateFresh(item.first);PackedJobStore::Extract(item.first);
            Require(unlink((item.first+"/"+item.second).c_str())==0,"cannot remove loose packed manifest");ValidateFresh(item.first);
            PackedJobStore::Member member;Require(PackedJobStore::Locate(item.first+"/"+item.second,member),"packed manifest missing");
            const int fd=open(member.archive.c_str(),O_RDWR);const char corrupt='X';
            Require(fd>=0&&pwrite(fd,&corrupt,1,member.offset)==1,"cannot mutate test packed manifest");close(fd);
            RejectMarker([&]{ValidateFresh(item.first);},"corrupt packed source manifest");
        }
    }

    // Same validation entry point for serial and independent parallel roots.
    std::vector<std::string> roots;for(int i=1;i<=8;++i){const auto folder=root+"/jobs/a"+Number(i);MarkerParent(folder,edge+edge);roots.push_back(folder);}
    {ValidationBatch batch("parallel markers",true);ValidateStoredRoots(roots,4);}
    Require(Read(parent+"/result/result.dat")==savedMatrix,"matrix changed during validation");

    // Missing loose files still resolve from a verified pack, including an
    // interrupted cleanup with both a directory and an archived copy present.
    PackedJobStore::Pack(leaf);ValidateFresh(parent);PackedJobStore::Extract(leaf);
    Require(unlink((leaf+"/result/result.dat").c_str())==0,"cannot make partial packed fixture");ValidateFresh(parent);
    PackedJobStore::Member member;Require(PackedJobStore::Locate(leaf+"/result/result.dat",member),"packed payload missing");
    const int fd=open(member.archive.c_str(),O_RDWR);const char corrupt='X';
    Require(fd>=0&&pwrite(fd,&corrupt,1,member.offset)==1,"cannot mutate test pack");close(fd);
    RejectMarker([&]{ValidateFresh(parent);},"corrupt packed child payload");
    Require(activeValidationBatch==nullptr,"validation scope leaked");
    std::cout<<"PASS source marker validation: first/duplicate/memo edges; empty hashes; payload/integral hashes; missing angular/native-paired/balanced-paired manifests; legacy four/five/six fields; parallel roots; packed and partial cleanup; matrix unchanged\n";
    return 0;
}catch(const std::exception& error){std::cerr<<error.what()<<'\n';return 1;}}
