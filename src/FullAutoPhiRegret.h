#pragma once
#include <algorithm>
#include <cmath>
#include <limits>
#include <map>

namespace FullAutoPhiRegret {
struct Cost {
    double row;
    double paired;
};
struct Choice {
    int phi;
    double relative;
};
// Each model supplies a positive cost in its own units. Choose the rule with
// the smallest worst relative slowdown to the best rule under either model.
inline Choice Select(const std::map<int,Cost>& costs,int original) {
    Choice selected{original,std::numeric_limits<double>::infinity()};
    if(!costs.count(original))return selected;
    double row=selected.relative,paired=selected.relative;
    for(const auto& entry:costs){
        const auto& t=entry.second;
        if(!(t.row>0 && t.paired>0 && std::isfinite(t.row) && std::isfinite(t.paired)))return selected;
        row=std::min(row,t.row);paired=std::min(paired,t.paired);
    }
    for(const auto& entry:costs){
        const auto& t=entry.second;
        const double ratio=std::max(t.row/row,t.paired/paired);
        if(ratio<selected.relative || (ratio==selected.relative && entry.first==original))
            selected={entry.first,ratio};
    }
    return selected;
}
}
