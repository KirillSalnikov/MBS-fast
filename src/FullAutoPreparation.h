#pragma once
#include <algorithm>
#include <array>
#include <cmath>
#include <deque>
#include <vector>
#include <fstream>
#include <set>
#include <map>
#include <string>
#include <sched.h>

namespace FullAutoPreparation {
// Count physical cores in this controller's cpuset, including SMT siblings only once.
inline int AllowedPhysicalCores() {
    cpu_set_t allowed;CPU_ZERO(&allowed);
    if(sched_getaffinity(0,sizeof(allowed),&allowed)!=0)return 0;
    std::set<std::string> cores;
    for(int cpu=0;cpu<CPU_SETSIZE;++cpu)if(CPU_ISSET(cpu,&allowed)){
        const std::string base="/sys/devices/system/cpu/cpu"+std::to_string(cpu)+"/topology/";
        int package=-1,core=-1;std::ifstream p(base+"physical_package_id"),c(base+"core_id");
        if(p>>package && c>>core && package>=0 && core>=0)cores.insert(std::to_string(package)+":"+std::to_string(core));
        else cores.insert("cpu:"+std::to_string(cpu));
    }
    return int(cores.size());
}
inline int ThreadCap(int physical,int allowed,int visibleGpus) {
    return std::max(1,std::min(allowed>0?allowed:physical,physical/std::max(1,visibleGpus)));
}

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
    std::deque<double> waits;
    void ObserveWait(double seconds,double processing) {
        if(!std::isfinite(seconds)||!std::isfinite(processing)||seconds<0||processing<=0||seconds>processing)return;
        waits.push_back(seconds/processing);while(waits.size()>32)waits.pop_front();
    }
    double WaitFraction()const {return Median(std::vector<double>(waits.begin(),waits.end()));}
    // Same exact theta count, same recipe/count bucket; compare wall time per
    // orientation. Never select CPU resources merely to increase GPU load.
    std::map<std::pair<int,size_t>,std::deque<double>> costs;
    void ObserveCost(int threads,size_t rows,double secondsPerOrientation) {
        if(threads<1||!std::isfinite(secondsPerOrientation)||secondsPerOrientation<=0)return;
        auto& values=costs[{threads,rows}];values.push_back(secondsPerOrientation);
        while(values.size()>8)values.pop_front();
    }
    double Cost(int threads,size_t rows)const {
        auto it=costs.find({threads,rows});if(it==costs.end()||it->second.size()<4)return 0;
        return Median(std::vector<double>(it->second.begin(),it->second.end()));
    }
    int SelectThreads(int baseline,int candidate,size_t rows)const {
        const double base=Cost(baseline,rows),next=Cost(candidate,rows);
        if(!base)return baseline; // gather real baseline work before exploration
        if(!next)return candidate;
        // One bounded retry in 64 live observations lets changed contention
        // replace old evidence. This uses useful work, never extra orientations.
        if(liveSamples>0&&liveSamples%64==63)return candidate;
        return next<base*.97?candidate:baseline;
    }
    void Restore(const std::vector<double>& trace,const std::vector<double>& diffraction,const std::vector<double>& waiting={}) {
        initialized=true;restoredSamples=trace.size();
        if(trace.empty())return;
        // A bounded prior survives resume without freezing the scheduler to
        // historical GPU contention. Nine fresh jobs can replace this prior.
        const std::array<double,2> prior={{Median(trace),Median(diffraction)}};
        for(int i=0;i<8;++i)recent.push_back(prior);
        if(!waiting.empty())for(int i=0;i<8;++i)waits.push_back(Median(waiting));
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
