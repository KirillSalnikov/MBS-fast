// Transport-level wall benchmark. The child delays stand in for optical work;
// the parent delays stand in for measured checkpoint/lookahead preparation.
// It measures controller scheduling only, not an optical/GPU speedup.
#include "../src/WorkerProtocol.h"
#include <algorithm>
#include <cassert>
#include <chrono>
#include <iostream>
#include <poll.h>
#include <signal.h>
#include <sys/prctl.h>
#include <sys/wait.h>

double seconds(){return std::chrono::duration<double>(std::chrono::steady_clock::now().time_since_epoch()).count();}
struct Worker {
    int fd;pid_t pid;
    Worker(){
        int channel[2];assert(socketpair(AF_UNIX,SOCK_STREAM,0,channel)==0);pid=fork();assert(pid>=0);
        if(pid==0){close(channel[0]);prctl(PR_SET_PDEATHSIG,SIGTERM);
            std::vector<std::string> request;
            while(WorkerProtocol::Receive(channel[1],request)&&!request.empty()){
                usleep(unsigned(std::stoul(request[0]))*1000);
                if(!WorkerProtocol::Reply(channel[1],0))break;
            }
            close(channel[1]);_exit(0);
        }
        close(channel[1]);fd=channel[0];
    }
    ~Worker(){close(fd);kill(pid,SIGTERM);int status;while(waitpid(pid,&status,0)<0&&errno==EINTR){}}
};
double run(bool grouped){
    const size_t slots=4,total=64;Worker workers[slots];
    std::vector<int> active(slots,-1);std::vector<unsigned> committed(total,0);size_t next=0,finished=0;
    const unsigned work[]={12,25,15,40,5,18,10,60},callback[]={5,5,20,4,10,7,15,5};
    auto launch=[&](size_t slot,size_t index){active[slot]=int(index);assert(WorkerProtocol::Send(workers[slot].fd,{std::to_string(work[index%8])}));};
    auto ready=[&](bool block){
        pollfd events[slots];for(size_t slot=0;slot<slots;++slot)events[slot]={active[slot]>=0?workers[slot].fd:-1,POLLIN,0};
        int rc;do{rc=poll(events,slots,block?-1:0);}while(rc<0&&errno==EINTR);assert(rc>=0);
        std::vector<size_t> result;for(size_t slot=0;slot<slots;++slot)if(events[slot].revents)result.push_back(slot);return result;
    };
    auto reap=[&](size_t slot,size_t& index){
        index=size_t(active[slot]);int status=-1;assert(WorkerProtocol::Result(workers[slot].fd,status)&&status==0);active[slot]=-1;
        if(next<total)launch(slot,next++);return true;
    };
    auto commit=[&](size_t index){assert(++committed[index]==1);++finished;usleep(callback[index%8]*1000);};
    const double begun=seconds();for(size_t slot=0;slot<slots;++slot)launch(slot,next++);
    if(grouped){const auto stats=WorkerProtocol::DispatchCompletions(total,slots,ready,reap,commit);assert(stats.pendingPeak<2*slots);}
    else while(finished<total)for(size_t slot:ready(true)){size_t index=0;reap(slot,index);commit(index);}
    for(unsigned n:committed)assert(n==1);
    return seconds()-begun;
}
int main(){
    std::vector<double> oldTimes,newTimes;
    for(int trial=0;trial<5;++trial){
        if(trial%2){newTimes.push_back(run(true));oldTimes.push_back(run(false));}
        else{oldTimes.push_back(run(false));newTimes.push_back(run(true));}
    }
    auto median=[](std::vector<double> v){std::sort(v.begin(),v.end());return v[v.size()/2];};
    std::cout<<"{\"transport_benchmark\":true,\"optical_work\":false,\"trials\":5,\"workers\":4,\"requests\":64,\"old_seconds\":[";
    for(size_t i=0;i<oldTimes.size();++i)std::cout<<(i?",":"")<<oldTimes[i];
    std::cout<<"],\"new_seconds\":[";for(size_t i=0;i<newTimes.size();++i)std::cout<<(i?",":"")<<newTimes[i];
    std::cout<<"],\"median_speedup\":"<<median(oldTimes)/median(newTimes)<<"}\n";
}
