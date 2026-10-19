#pragma once
#include <algorithm>
#include <cmath>
#include <cstddef>

namespace FullAutoScheduling {
// An unbounded ratio interval is a sampling failure, never a work estimate.
inline bool Resolved(double error) { return std::isfinite(error)&&error>=0&&error<1e50; }
inline int PilotLimit(int count,int maximum,bool fixed) {
    if(fixed)return count;
    return std::max(count,std::min(maximum,count>4096/8?4096:count*8));
}
inline int NextPilot(int count,int limit) { return count>=limit?count:std::min(limit,count>limit/2?limit:2*count); }
inline bool ShortPacket(size_t rows,int phi,int count,int chunk,bool explicitChunk,int candidateChunk=1024) {
    return !explicitChunk&&chunk!=candidateChunk&&count>=2048&&rows>0&&rows<=256&&phi>0&&rows<=size_t(1024/phi);
}
inline bool RetryPhi(int oldPeak,int peak,int attempts) {
    return oldPeak>0&&attempts<3&&peak/oldPeak>=4;
}
}
