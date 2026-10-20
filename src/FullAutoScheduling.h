#pragma once
#include <algorithm>
#include <cmath>
#include <cstddef>
#include <vector>

namespace FullAutoScheduling {
// An unbounded ratio interval is a sampling failure, never a work estimate.
inline bool Resolved(double error) { return std::isfinite(error)&&error>=0&&error<1e50; }
inline int PilotLimit(int count,int maximum,bool fixed) {
    if(fixed)return count;
    return std::max(count,std::min(maximum,count>4096/8?4096:count*8));
}
inline int NextPilot(int count,int limit) { return count>=limit?count:std::min(limit,count>limit/2?limit:2*count); }
inline std::vector<std::vector<double>> AngleBatches(const std::vector<double>& angles,size_t cap) {
    if(angles.empty())return {};
    if(!cap||angles.size()<=cap)return {angles};
    cap=std::max(size_t(2),cap);
    // Balance the last two packets: never emit a singleton or exceed the cap.
    const size_t packets=(angles.size()+cap-1)/cap;
    const size_t width=angles.size()/packets,extra=angles.size()%packets;
    std::vector<std::vector<double>> result;size_t begin=0;
    for(size_t i=0;i<packets;++i){const size_t end=begin+width+(i<extra);
        result.emplace_back(angles.begin()+begin,angles.begin()+end);begin=end;}
    return result;
}
inline int ShortPacketCandidate(size_t rows,bool fp32) {
    // Very narrow angular packets cannot hide long producer bursts. This is
    // only a candidate: the controller still checks two independent pilots
    // for all16/ratio agreement AND a measured reduction in processing time.
    return fp32||rows<=8?512:1024;
}
inline bool ShortPacket(size_t rows,int phi,int count,int chunk,bool explicitChunk,int candidateChunk=1024) {
    return !explicitChunk&&chunk!=candidateChunk&&count>=2048&&rows>0&&rows<=256&&phi>0&&rows<=size_t(1024/phi);
}
inline bool RetryPhi(int oldPeak,int peak,int attempts) {
    return oldPeak>0&&attempts<3&&peak/oldPeak>=4;
}
}
