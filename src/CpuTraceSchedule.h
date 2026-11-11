#pragma once
#include "TraceTreeReserve.h"
#include <algorithm>
#include <chrono>
#include <cstdlib>
#include <ctime>
#include <limits>
#include <vector>

// Scheduling only: orientation indices, optical work and reduction order are
// owned by the caller. ScratchBytes is an admission estimate, not a bound on
// arbitrary beam-tree growth. The ordinary scheduler remains unchanged.
namespace CpuTraceSchedule {
inline bool Fine(const char* value) { return value && value[0]=='1' && value[1]=='\0'; }
inline int Grain(bool fine) { return fine ? 1 : 4; }
inline size_t ReserveCount(const char* value,bool convex=false) { return TraceTreeReserve::Count(value,convex); }

inline bool ScratchBytes(size_t beamBytes,size_t reserve,size_t base,size_t& bytes) {
    const size_t max=std::numeric_limits<size_t>::max();
    if(beamBytes && reserve>(max-base)/2/beamBytes)return false;
    bytes=base+2*reserve*beamBytes;return bytes>0;
}
inline int Threads(bool fine,int requested,int count,size_t budget,size_t scratch) {
    if(!fine)return std::max(1,requested);
    if(!scratch || budget<scratch)return 0;
    return static_cast<int>(std::min<size_t>(std::max(1,std::min(requested,count)),budget/scratch));
}
inline double ThreadCpuSeconds() {
#if defined(CLOCK_THREAD_CPUTIME_ID)
    timespec value;
    if(clock_gettime(CLOCK_THREAD_CPUTIME_ID,&value)==0)
        return double(value.tv_sec)+double(value.tv_nsec)*1e-9;
#endif
    return -1;
}
using Clock=std::chrono::steady_clock;
inline double Seconds(Clock::time_point begin,Clock::time_point end) {
    return std::chrono::duration<double>(end-begin).count();
}
struct ThreadProfile {
    bool present=false;
    size_t tasks=0;
    double setup=0,work=0,loop=0,barrier=0,teardown=0,cpu=-1;
};
// Counts cache hits and caught exceptions as assigned work too; no clock call
// or counter traffic when profiling is disabled.
class TaskProfile {
    ThreadProfile* profile;
    Clock::time_point started;
public:
    explicit TaskProfile(ThreadProfile* p):profile(p) {
        if(profile){started=Clock::now();++profile->tasks;}
    }
    ~TaskProfile(){if(profile)profile->work+=Seconds(started,Clock::now());}
};
struct Summary {
    size_t threads=0,active=0,tasks=0,minTasks=0,maxTasks=0;
    double setup=0,work=0,maxWork=0,loop=0,maxLoop=0,barrier=0,maxBarrier=0,teardown=0,cpu=0;
    bool cpuAvailable=true;
};
inline Summary Summarize(const std::vector<ThreadProfile>& profiles) {
    Summary s;s.minTasks=std::numeric_limits<size_t>::max();
    for(const auto& p:profiles){if(!p.present)continue;
        ++s.threads;s.active+=p.tasks>0;s.tasks+=p.tasks;
        s.minTasks=std::min(s.minTasks,p.tasks);s.maxTasks=std::max(s.maxTasks,p.tasks);
        s.setup+=p.setup;s.work+=p.work;s.maxWork=std::max(s.maxWork,p.work);
        s.loop+=p.loop;s.maxLoop=std::max(s.maxLoop,p.loop);
        s.barrier+=p.barrier;s.maxBarrier=std::max(s.maxBarrier,p.barrier);
        s.teardown+=p.teardown;if(p.cpu<0)s.cpuAvailable=false;else s.cpu+=p.cpu;
    }
    if(!s.threads)s.minTasks=0;
    return s;
}
}
