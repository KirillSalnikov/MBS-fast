#include "FullAutoTraceEvidence.h"
#include <cassert>
#include <cstdlib>
#include <iostream>
#include <sys/stat.h>
#include <unistd.h>
int main(){
    char path[]="/tmp/mbs-trace-evidence-XXXXXX";assert(mkdtemp(path));
    const std::string root=path,parent=root+"/assembled",child=root+"/native child",view=root+"/prefix";
    for(const auto& d:{parent,child,view})assert(mkdir(d.c_str(),0700)==0);
    std::ofstream(parent+"/stdout.log")<<"Balanced control-variate assembly\n";
    std::ofstream(parent+"/angular_sources.tsv")<<view<<"\tchecksum\n";
    std::ofstream(view+"/angular_sources.tsv")<<parent<<"\tcycle\n"<<child<<"\tchecksum\n";
    std::ofstream(child+"/stdout.log")<<"GPU convex resident direct: accepted=0 CPU_fallback=64\n";
    assert(!FullAutoTraceEvidence::Executed(parent));
    std::ofstream(child+"/stdout.log")<<"GPU convex resident direct: accepted=64 CPU_fallback=0\n";
    assert(FullAutoTraceEvidence::Executed(parent));
    assert(FullAutoTraceEvidence::Executed(child));
    assert(!FullAutoTraceEvidence::Executed(root+"/missing"));
    assert(mkdir((root+"/jobs").c_str(),0700)==0);
    const auto packedParent=root+"/jobs/c123",packedView=root+"/jobs/v123",packedChild=root+"/jobs/j123";
    for(const auto& d:{packedParent,packedView,packedChild}){
        assert(mkdir(d.c_str(),0700)==0);std::ofstream(d+"/complete.txt")<<"test committed block\n";
    }
    std::ofstream(packedParent+"/stdout.log")<<"Balanced control-variate assembly\n";
    std::ofstream(packedParent+"/angular_sources.tsv")<<packedView<<"\tchecksum\n";
    std::ofstream(packedView+"/angular_sources.tsv")<<packedParent<<"\tcycle\n"<<packedChild<<"\tchecksum\n";
    std::ofstream(packedChild+"/stdout.log")<<"GPU convex resident direct: accepted=64 CPU_fallback=0\n";
    for(const auto& d:{packedParent,packedView,packedChild})PackedJobStore::Pack(d);
    assert(FullAutoTraceEvidence::Executed(packedParent));
    assert(FullAutoTraceEvidence::Executed(packedChild));
    assert(!FullAutoTraceEvidence::Executed(root+"/jobs/j999"));
    for(const auto& d:{parent,child,view}){
        unlink((d+"/stdout.log").c_str());unlink((d+"/angular_sources.tsv").c_str());rmdir(d.c_str());}
    rmdir(root.c_str());std::cout<<"Resident execution evidence: PASS\n";
}
