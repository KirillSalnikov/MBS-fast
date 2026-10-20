// Deterministic event-clock regression: callbacks represent slow descriptor
// preparation/checkpoint I/O, and workers represent independent native jobs.
// Exercise the production dispatch policy without GPU scheduling noise.
#include "../src/WorkerProtocol.h"
#include <algorithm>
#include <cassert>
#include <cmath>
#include <iostream>
#include <limits>
#include <stdexcept>

struct Simulation {
    const size_t slots,total;
    const double work,callback;
    double now=0,idle=0;
    size_t next=0,committed=0,retries=0;
    bool injectRetry=false;
    std::vector<int> jobs;
    std::vector<double> due;
    std::vector<unsigned> completions,attempts;
    std::vector<std::string> events;
    Simulation(size_t s,size_t n,double w,double c,bool retry=false):
        slots(s),total(n),work(w),callback(c),injectRetry(retry),
        jobs(s,-1),due(s,0),completions(n,0),attempts(n,0){
        for(size_t slot=0;slot<slots&&next<total;++slot)launch(slot,next++);
    }
    void launch(size_t slot,size_t index){
        jobs[slot]=int(index);due[slot]=now+work;++attempts[index];
        events.push_back("launch:"+std::to_string(index));
    }
    std::vector<size_t> ready(bool block){
        double first=std::numeric_limits<double>::infinity();
        for(size_t slot=0;slot<slots;++slot)if(jobs[slot]>=0)first=std::min(first,due[slot]);
        if(block){assert(std::isfinite(first));now=std::max(now,first);}
        std::vector<size_t> result;for(size_t slot=0;slot<slots;++slot)if(jobs[slot]>=0&&due[slot]<=now)result.push_back(slot);
        return result;
    }
    bool reap(size_t slot,size_t& index){
        index=size_t(jobs[slot]);assert(due[slot]<=now);idle+=now-due[slot];jobs[slot]=-1;
        if(injectRetry&&index==1&&attempts[index]==1){++retries;launch(slot,index);return false;}
        if(next<total)launch(slot,next++);
        return true;
    }
    void commit(size_t index){
        events.push_back("commit:"+std::to_string(index));
        assert(++completions[index]==1);++committed;now+=callback;
    }
    void oldDispatch(){
        while(committed<total)for(size_t slot:ready(true)){size_t index=0;if(reap(slot,index))commit(index);}
    }
    WorkerProtocol::CompletionDispatchStats newDispatch(){
        return WorkerProtocol::DispatchCompletions(total,slots,
            [&](bool block){return ready(block);},
            [&](size_t slot,size_t& index){return reap(slot,index);},
            [&](size_t index){commit(index);});
    }
    void check()const{
        assert(committed==total);
        for(size_t i=0;i<total;++i){assert(completions[i]==1);assert(attempts[i]==(injectRetry&&i==1?2u:1u));}
        for(int job:jobs)assert(job==-1);
    }
};
int main(){
    for(size_t slots:{size_t(1),size_t(2),size_t(4),size_t(8)}){
        for(size_t count:{size_t(1),slots,slots*12})for(bool retry:{false,true}){
            Simulation run(slots,count,100,40,retry);auto stats=run.newDispatch();run.check();
            assert(stats.pendingPeak<=2*slots-1);assert(stats.readyBatchPeak<=slots);
            assert(run.retries==(retry&&count>1?1u:0u));
            if(count>=2*slots&&slots>1){
                const auto firstCommit=std::find_if(run.events.begin(),run.events.end(),[](const std::string& e){return e.find("commit:")==0;});
                // Every ready worker, including a one-shot retry, is fed before
                // the first durable checkpoint or expensive caller callback.
                assert(size_t(firstCommit-run.events.begin())==2*slots);
            }
        }
    }
    Simulation oldRun(4,48,100,40),newRun(4,48,100,40);
    oldRun.oldDispatch();auto stats=newRun.newDispatch();oldRun.check();newRun.check();
    assert(newRun.idle<oldRun.idle);
    std::cout<<"{\"old_event_clock_wall\":"<<oldRun.now<<",\"new_event_clock_wall\":"<<newRun.now
        <<",\"old_worker_idle\":"<<oldRun.idle<<",\"new_worker_idle\":"<<newRun.idle
        <<",\"pending_peak\":"<<stats.pendingPeak<<",\"pending_limit\":7,\"simulation\":true}\n";
}
