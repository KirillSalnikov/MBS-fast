#pragma once
#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <string>

// PreparedChunkController already reports a padded maximum over observed orientations.
namespace FullAutoPhiTileMemory {
inline double Peak(const std::string& log){
    const std::string prefix="Orientation average adaptive memory: measured ";
    const std::string field=", conservative ";
    double maximum=0;
    for(size_t start=log.find(prefix);start!=std::string::npos;start=log.find(prefix,start+prefix.size())){
        const size_t end=log.find('\n',start);
        const size_t at=log.find(field,start+prefix.size());
        if(at==std::string::npos || (end!=std::string::npos && at>=end))continue;
        const char* value=log.c_str()+at+field.size();
        char* suffix=nullptr;
        const double parsed=std::strtod(value,&suffix);
        if(value!=suffix && std::isfinite(parsed) && parsed>0 &&
           log.compare(static_cast<size_t>(suffix-log.c_str()),16," KiB/orientation")==0)
            maximum=std::max(maximum,parsed);
    }
    return maximum;
}
inline int Fit(long long cacheMiB,double safeKiB,int ceiling){
    if(!(safeKiB>0) || !std::isfinite(safeKiB))return ceiling;
    int n=1;
    while(n<=ceiling/2 && safeKiB*(2*n)<=1024.0*cacheMiB)n*=2;
    return n;
}
}
