#pragma once
#include <algorithm>
#include <cstdlib>
#include <fstream>
#include <sstream>
#include <string>
#include <stdexcept>

namespace HostMemoryReserve {
// A recoverable per-process limit, distinct from invalid input or CUDA errors.
struct BudgetExhausted : std::runtime_error { using std::runtime_error::runtime_error; };
struct Snapshot {long long totalKb=0,availableKb=0;};
inline Snapshot Read(){
    Snapshot s;std::ifstream input("/proc/meminfo");std::string line;
    while(std::getline(input,line)){std::istringstream row(line);std::string key;long long value;if(!(row>>key>>value))continue;
        if(key=="MemTotal:")s.totalKb=value;
        if(key=="MemAvailable:")s.availableKb=value;
        }
    return s;
}
inline long long ReserveKb(long long totalKb){
    const char* setting=std::getenv("MBS_HOST_MEM_RESERVE_MB");char* end=nullptr;
    const long long requested=setting?std::strtoll(setting,&end,10):0;
    if(setting&&end!=setting&&!*end&&requested>0&&requested<=1048576)return requested*1024;
    // The safety reserve must not shrink as available RAM approaches zero.
    return std::min(4096LL*1024,std::max(256LL*1024,totalKb/16));
}
inline bool Low(const Snapshot& s){return s.totalKb>0&&s.availableKb<ReserveKb(s.totalKb)+16LL*1024;}
}
