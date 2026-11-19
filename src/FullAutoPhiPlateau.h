#pragma once
#include <algorithm>
#include <cmath>
#include <map>
namespace PhiPlateau {
struct Observation {
    int count;
    double raw, ratio;
    bool resolved;
};
// Compare the measured uncertainties themselves. Equal orientation floors can
// conceal a material difference in the underlying validation measurements.
inline bool Same(double a,double b) {
    return std::isfinite(a)&&std::isfinite(b)&&a>=0&&b>=0
        &&std::abs(a-b)<=2e-4*std::max(a,b);
}
inline int DenseRule(int original,const std::map<int,Observation>& observations,bool ratios) {
    auto at=observations.find(original);
    if(at==observations.end()||!at->second.resolved)return original;
    const auto& ref=at->second;
    int best=original;
    while(best>1&&best%2==0) {
        const auto lower=observations.find(best/2);
        if(lower==observations.end())break;
        const auto& entry=lower->second;
        if(!entry.resolved||entry.count!=ref.count||!Same(ref.raw,entry.raw)
           ||(ratios&&!Same(ref.ratio,entry.ratio)))break;
        best=lower->first;
    }
    return best;
}
}
