#pragma once
#include <fstream>
#include <set>
#include <sstream>
#include <string>

// Corrected jobs contain an assembly log, not the logs of their native
// children. Follow the recorded source graph to verify actual GPU execution.
namespace FullAutoTraceEvidence {
inline bool Executed(const std::string& folder,std::set<std::string>& seen){
    if(!seen.insert(folder).second)return false;
    std::ifstream log(folder+"/stdout.log");std::string line;
    const std::string marker="GPU convex resident direct: accepted=";
    while(std::getline(log,line)){
        const auto pos=line.find(marker);if(pos==std::string::npos)continue;
        std::istringstream value(line.substr(pos+marker.size()));long long count=0;
        if(value>>count&&count>0)return true;
    }
    std::ifstream sources(folder+"/angular_sources.tsv");
    while(std::getline(sources,line)){
        const auto tab=line.find('\t');if(tab==std::string::npos||!tab)continue;
        if(Executed(line.substr(0,tab),seen))return true;
    }
    return false;
}
inline bool Executed(const std::string& folder){std::set<std::string> seen;return Executed(folder,seen);}
}
