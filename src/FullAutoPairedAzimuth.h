#pragma once
#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <vector>

namespace FullAutoPairedAzimuth {
struct Difference {
    double scale=0;
    size_t row=0, element=0;
    bool ratio=false, valid=false;
};

// Compare matched scrambled means. The normalized statistic uses the
// influence function of each ratio of means, including its own denominator.
// This schedules pilot work only; independent production acceptance is intact.
template<class Samples>
Difference Measure(const Samples& lower,const Samples& upper,bool ratios) {
    Difference result;
    if(lower.size()!=8||upper.size()!=8||lower[0].empty())return result;
    const size_t rows=lower[0].size();
    for(size_t s=0;s<8;++s)
        if(lower[s].size()!=rows||upper[s].size()!=rows)return result;
    for(size_t r=0;r<rows;++r){
        std::array<double,16> lo{},hi{};
        for(size_t s=0;s<8;++s)for(size_t j=0;j<16;++j){
            lo[j]+=lower[s][r][j]/8;
            hi[j]+=upper[s][r][j]/8;
        }
        if(!(lo[0]>0&&hi[0]>0))return Difference{};
        for(size_t j=0;j<16;++j){
            const double lr=lo[j]/lo[0],hr=hi[j]/hi[0];
            for(int norm=0;norm<=(ratios&&j!=0);++norm){
                const double mean=norm?lr-hr:(lo[j]-hi[j])/hi[0];
                double sum=0;
                for(size_t s=0;s<8;++s){
                    const auto& a=lower[s][r];const auto& b=upper[s][r];
                    const double residual=norm?(a[j]-lr*a[0])/lo[0]-(b[j]-hr*b[0])/hi[0]
                        :(a[j]-b[j])/hi[0]-mean;
                    sum+=residual*residual;
                }
                const double se=std::sqrt(sum/56.);
                if(!std::isfinite(mean)||!std::isfinite(se))return Difference{};
                const double width=std::max(std::abs(mean),2.364624251*se);
                if(width>result.scale){result.scale=width;result.row=r;result.element=j;result.ratio=norm;}
            }
        }
    }
    result.valid=true;
    return result;
}
}
