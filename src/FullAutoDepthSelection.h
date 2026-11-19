#pragma once
#include <algorithm>
#include <cmath>
#include <vector>

namespace FullAutoDepthSelection {
inline const char* Policy(){return "common-reference-local-pilot-v4-max-fallback";}
inline int ProductionDepth(bool passed,int candidate,int maximum){return passed?candidate:maximum;}
inline bool FiniteBound(double v){return std::isfinite(v)&&v>=0&&v<1e50;}
// The worker supports at most 30 interactions. Keep explicit depth-max=28
// usable: its reference guard is 28/30 rather than attempting depth 32.
inline int Reference(int maximum){return std::min(28,maximum+2);}
inline int Anchor(int maximum){return std::min(30,maximum+4);}
inline int PilotFloor(int initial,int maximum){return std::min(maximum,std::max(initial,4096));}
inline int PilotLimit(int initial,int maximum){
    return std::min(maximum,std::max(4096,initial>maximum/256?maximum:initial*256));
}
inline int DoubleCount(int count,int limit){return count>=limit?limit:(count>limit/2?limit:2*count);}
enum Decision { Pass, Fail, Unresolved, ReferenceFail };
inline Decision Classify(double upper,double lower,double referenceUpper,double referenceLower,double budget){
    if(!(budget>0)||!std::isfinite(budget))return Unresolved;
    if(!FiniteBound(referenceUpper)||!FiniteBound(upper)||!FiniteBound(lower)||!FiniteBound(referenceLower))return Unresolved;
    if(referenceLower>budget)return ReferenceFail;
    if(referenceUpper>budget)return Unresolved;
    if(FiniteBound(upper)&&upper<=budget)return Pass;
    if(FiniteBound(lower)&&lower>budget)return Fail;
    return Unresolved;
}
inline std::vector<int> Candidates(int start,int maximum,bool shallowestFirst=false){
    std::vector<int> result;
    for(int n=start;;n=std::min(maximum,n>maximum/2?maximum:2*n)){
        result.push_back(n);if(n==maximum)break;
    }
    // No monotonicity assumption: inspect every skipped depth below the first
    // passing coarse candidate, not just the last coarse interval.
    for(int n=start;n<=maximum;n+=2)
        if(std::find(result.begin(),result.end(),n)==result.end())result.push_back(n);
    if(shallowestFirst)std::sort(result.begin(),result.end());
    return result;
}
}
