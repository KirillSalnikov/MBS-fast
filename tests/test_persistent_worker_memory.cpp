#include "PersistentWorkerMemory.h"
#include <cassert>
#include <cstring>
#include <iostream>
#include <limits>
int main(){
    using PersistentWorkerMemory::TrimPolicy;
    RuntimeMemorySnapshot m;m.ramTotalKb=64LL*1024*1024;m.ramAvailableKb=16LL*1024*1024;
    m.processRssKb=800LL*1024;const long long budget=4LL*1024*1024,reserve=4LL*1024*1024;
    TrimPolicy p;assert(std::strcmp(p.Decide(m,10,budget,reserve),"initial")==0);
    p.Trimmed(10,600LL*1024);
    for(int i=0;i<20;++i)assert(!p.Decide(m,10+i*.1,budget,reserve));
    assert(std::strcmp(p.Decide(m,15,budget,reserve),"periodic")==0);
    p.Trimmed(15,600LL*1024);m.processRssKb=857LL*1024;
    assert(std::strcmp(p.Decide(m,15.01,budget,reserve),"rss_growth")==0);
    m.processRssKb=3LL*1024*1024;
    assert(std::strcmp(p.Decide(m,15.01,budget,reserve),"process_budget")==0);
    m.processRssKb=800LL*1024;m.ramAvailableKb=reserve+256LL*1024;
    assert(std::strcmp(p.Decide(m,15.01,budget,reserve),"system_reserve")==0);
    m.ramAvailableKb=16LL*1024*1024;p.Trimmed(20,600LL*1024);
    for(int i=0;i<63;++i)assert(!p.Decide(m,20+i*.01,budget,reserve));
    assert(std::strcmp(p.Decide(m,20.64,budget,reserve),"request_limit")==0);
    p.Trimmed(21,600LL*1024);m.processRssKb=128LL*1024;
    assert(!p.Decide(m,100,budget,reserve));
    // Low process RSS does not bypass pressure or unknown-memory cleanup.
    m.ramAvailableKb=reserve;
    assert(std::strcmp(p.Decide(m,100,budget,reserve),"system_reserve")==0);
    m.ramAvailableKb=0;
    assert(std::strcmp(p.Decide(m,100,budget,reserve),"memory_unknown")==0);
    m.ramAvailableKb=16LL*1024*1024;
    assert(std::strcmp(p.Decide(m,std::numeric_limits<double>::quiet_NaN(),budget,reserve),"memory_unknown")==0);
    m.processRssKb=800LL*1024;
    assert(std::strcmp(p.Decide(m,20,budget,reserve),"initial")==0);
    unsetenv("MBS_HOST_MEM_BUDGET_MB");assert(PersistentWorkerMemory::ProcessBudgetKb()==0);
    for(const char* value:{"-1","0","x","99999999999999999999999999","2048garbage"}){
        setenv("MBS_HOST_MEM_BUDGET_MB",value,1);assert(PersistentWorkerMemory::ProcessBudgetKb()==0);}
    setenv("MBS_HOST_MEM_BUDGET_MB","2048",1);assert(PersistentWorkerMemory::ProcessBudgetKb()==2LL*1024*1024);
    unsetenv("MBS_HOST_MEM_BUDGET_MB");
    std::cout<<"Persistent worker memory policy: PASS\n";
}
