#pragma once
#include <algorithm>
#include <array>
#include <cmath>
#include <deque>
#include <vector>

namespace FullAutoPreparation {
inline size_t Bucket(size_t n) {
    size_t b=1;while(b<n && b<=n/2)b*=2;
    return b; // floor power of two, without overflow
}
inline double Median(std::vector<double> v) {
    if(v.empty())return 0;
    std::sort(v.begin(),v.end());const size_t n=v.size();
    return n%2?v[n/2]:.5*(v[n/2-1]+v[n/2]);
}
struct Profile {
    bool initialized=false;int lastSelected=-1;size_t restoredSamples=0,liveSamples=0;
    std::deque<std::array<double,2>> recent;
    void Restore(const std::vector<double>& trace,const std::vector<double>& diffraction) {
        initialized=true;restoredSamples=trace.size();
        if(trace.empty())return;
        // A bounded prior survives resume without freezing the scheduler to
        // historical GPU contention. Nine fresh jobs can replace this prior.
        const std::array<double,2> prior={{Median(trace),Median(diffraction)}};
        for(int i=0;i<8;++i)recent.push_back(prior);
    }
    void Observe(double trace,double diffraction) {
        if(!std::isfinite(trace)||!std::isfinite(diffraction)||trace<=0||diffraction<0)return;
        recent.push_back({{trace,diffraction}});++liveSamples;
        while(recent.size()>32)recent.pop_front();
    }
    std::array<double,2> Estimate()const {
        std::vector<double> t,d;for(const auto& p:recent){t.push_back(p[0]);d.push_back(p[1]);}
        return {{Median(t),Median(d)}};
    }
};
}
