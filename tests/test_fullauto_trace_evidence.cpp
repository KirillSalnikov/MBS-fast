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
    for(const auto& d:{parent,child,view}){
        unlink((d+"/stdout.log").c_str());unlink((d+"/angular_sources.tsv").c_str());rmdir(d.c_str());}
    rmdir(root.c_str());std::cout<<"Resident execution evidence: PASS\n";
}
