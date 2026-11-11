#pragma once
#include "RuntimeInfo.h"
#include "HostMemoryReserve.h"
#include <cmath>
#include <cstdlib>
#include <cerrno>

// Releasing every short request's reusable malloc arenas stalls the next CPU
// preparation and faults the same pages back in. Keep bounded slack between
// requests; memory pressure always takes precedence over the time cadence.
namespace PersistentWorkerMemory {
inline long long ProcessBudgetKb() {
    const char* text=std::getenv("MBS_HOST_MEM_BUDGET_MB");
    if(!text||!*text)return 0;
    char* end=nullptr;errno=0;const long long mb=std::strtoll(text,&end,10);
    return !errno&&end!=text&&!*end&&mb>0&&mb<=1048576?mb*1024:0;
}
struct TrimPolicy {
    double lastTrim=-1;
    long long retainedAfterTrimKb=0;
    unsigned requestsSinceTrim=0;
    const char* Decide(const RuntimeMemorySnapshot& memory,double now,
                       long long processBudgetKb,long long reserveKb) {
        ++requestsSinceTrim;
        // Unknown live counters must not silently disable the old cleanup.
        if(memory.processRssKb<=0||memory.ramTotalKb<=0||memory.ramAvailableKb<=0
           ||!std::isfinite(now))return "memory_unknown";
        if(memory.ramAvailableKb<=reserveKb+256LL*1024)return "system_reserve";
        if(processBudgetKb>0&&memory.processRssKb>=processBudgetKb-processBudgetKb/4)
            return "process_budget";
        if(memory.processRssKb<=256LL*1024)return nullptr;
        if(lastTrim<0||now<lastTrim)return "initial";
        if(memory.processRssKb>retainedAfterTrimKb+256LL*1024)return "rss_growth";
        if(now-lastTrim>=5)return "periodic";
        if(requestsSinceTrim>=64)return "request_limit";
        return nullptr;
    }
    void Trimmed(double now,long long retainedKb) {
        lastTrim=now;retainedAfterTrimKb=retainedKb;requestsSinceTrim=0;
    }
};
}
