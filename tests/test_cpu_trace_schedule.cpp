#include "CpuTraceSchedule.h"
#include <omp.h>
#include <atomic>
#include <cassert>
#include <cmath>
#include <iostream>
#include <stdexcept>
#include <thread>

using namespace CpuTraceSchedule;

static void Run(int count,int requested,bool fine,bool synchronizeFirstWave) {
    const int grain=Grain(fine);
    std::vector<int> seen(count,0),owner(count,-1);
    std::vector<unsigned long long> values(count,0);
    std::vector<ThreadProfile> profiles(requested);
    std::atomic<int> arrived(0);
    int actual=0;
    #pragma omp parallel num_threads(requested)
    {
        const int tid=omp_get_thread_num(),team=omp_get_num_threads();
        #pragma omp single
        actual=team;
        auto& p=profiles[tid];p.present=true;
        const double cpuStart=ThreadCpuSeconds();
        const auto started=Clock::now();
        bool first=true;
        #pragma omp for schedule(dynamic,grain) nowait
        for(int i=0;i<count;++i){
            TaskProfile task(&p);
            if(first && synchronizeFirstWave){
                first=false;arrived.fetch_add(1);
                // A bounded first-wave rendezvous makes the scheduler's group
                // parallelism observable without a fragile elapsed-time test.
                const int target=std::min(team,(count+grain-1)/grain);
                const auto deadline=Clock::now()+std::chrono::seconds(10);
                while(arrived.load()<target){assert(Clock::now()<deadline);std::this_thread::yield();}
            }
            ++seen[i];owner[i]=tid;
            // Indexed cached and exception paths must still count one task.
            if(i%7==1){values[i]=123;continue;}
            try {
                if(i%11==3)throw std::runtime_error("synthetic caught trace failure");
                unsigned long long value=17;
                const int repetitions=2000+(i%5)*137;
                for(int k=0;k<repetitions;++k)value=value*6364136223846793005ULL+1;
                values[i]=value;
            }catch(const std::runtime_error&){values[i]=456;}
        }
        p.loop=Seconds(started,Clock::now());
        const auto wait=Clock::now();
        #pragma omp barrier
        p.barrier=Seconds(wait,Clock::now());
        const double cpuEnd=ThreadCpuSeconds();
        p.cpu=cpuStart>=0&&cpuEnd>=cpuStart?cpuEnd-cpuStart:-1;
    }
    const auto s=Summarize(profiles);
    assert(s.threads==size_t(actual)&&s.tasks==size_t(count));
    assert(s.active<=size_t(std::min(actual,(count+grain-1)/grain)));
    if(synchronizeFirstWave)assert(s.active==size_t(std::min(actual,(count+grain-1)/grain)));
    assert(s.work>=0&&s.maxWork>=0&&s.barrier>=0&&s.maxBarrier>=0);
    if(s.cpuAvailable)assert(std::isfinite(s.cpu)&&s.cpu>=0);
    for(int i=0;i<count;++i){
        assert(seen[i]==1&&owner[i]>=0&&owner[i]<actual);
        unsigned long long reference=17;
        for(int k=0;k<2000+(i%5)*137;++k)reference=reference*6364136223846793005ULL+1;
        if(i%7==1)reference=123;else if(i%11==3)reference=456;
        assert(values[i]==reference);
    }
    std::cout<<"count="<<count<<" grain="<<grain<<" team="<<s.threads<<" active="<<s.active<<" tasks="<<s.tasks<<" PASS\n";
}

int main(){
    omp_set_dynamic(0);
    assert(!Fine(nullptr)&&!Fine("")&&!Fine("0")&&!Fine("true")&&!Fine("10")&&Fine("1"));
    assert(ReserveCount(nullptr)==65536&&ReserveCount("0")==65536&&ReserveCount("-1")==65536);
    assert(ReserveCount("8x")==65536&&ReserveCount("8")==8);
    size_t scratch=0;assert(ScratchBytes(2248,65536,17*1024*1024,scratch));
    assert(scratch==2ULL*2248*65536+17*1024*1024);
    assert(Threads(false,16,16,0,0)==16);
    assert(Threads(true,16,16,3*scratch,scratch)==3);
    assert(Threads(true,16,2,32*scratch,scratch)==2);
    assert(Threads(true,16,16,scratch-1,scratch)==0);
    assert(Threads(true,16,16,16*scratch,scratch)==16);
    assert(!ScratchBytes(std::numeric_limits<size_t>::max(),2,0,scratch));
    assert(!ScratchBytes(1,std::numeric_limits<size_t>::max(),17,scratch));
    assert(!ScratchBytes(0,1,0,scratch));
    for(bool fine:{false,true}){
        Run(16,16,fine,true);
        for(int count:{0,1,2,3,4,7,17,97})Run(count,16,fine,false);
        Run(16,1,fine,false);
    }
    std::cout<<"PASS fine CPU orientation scheduling, exact-once indexed results, cached/error continuations, thread timing and bounded scratch admission\n";
}
