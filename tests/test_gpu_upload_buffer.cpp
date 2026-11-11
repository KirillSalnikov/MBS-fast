#include "cuda/GpuUploadBuffer.h"
#include <sys/wait.h>
#include <unistd.h>
#include <cassert>
#include <cstdint>
#include <map>
#include <cstdio>
struct Cell {
    static size_t alive;
    uint64_t value[16];
    Cell(){++alive;for(auto& x:value)x=0;}
    Cell(const Cell& other){++alive;for(size_t i=0;i<16;++i)value[i]=other.value[i];}
    ~Cell(){--alive;}
};
size_t Cell::alive=0;
struct MockRegistration {
    static std::map<void*,size_t> active;
    static bool failRegister,failUnregister;
    static size_t registrations,unregistrations,clears;
    static int Register(void* p,size_t n){++registrations;if(failRegister)return 1;assert(!active.count(p));active[p]=n;return 0;}
    static int Unregister(void* p){++unregistrations;assert(active.count(p));if(failUnregister)return 2;active.erase(p);return 0;}
    static const char* Error(int){return "injected registration error";}
    static void ClearError(){++clears;}
};
std::map<void*,size_t> MockRegistration::active;
bool MockRegistration::failRegister=false,MockRegistration::failUnregister=false;
size_t MockRegistration::registrations=0,MockRegistration::unregistrations=0,MockRegistration::clears=0;
using Buffer=MbsGpuUpload::Vector<Cell,MockRegistration>;
using Budget=MbsGpuUpload::Budget;
static void option(const char* name,const char* value){assert(setenv(name,value,1)==0);}
static void fill(Buffer& b){for(size_t i=0;i<b.size();++i)for(size_t k=0;k<16;++k)b[i].value[k]=i*37+k;}
static void prefix(Buffer& b,size_t count){for(size_t i=0;i<count;++i)for(size_t k=0;k<16;++k)assert(b[i].value[k]==i*37+k);}
static void clean(){assert(MockRegistration::active.empty()&&Cell::alive==0);}
int main(){
    const size_t mib=1024*1024,initial=mib/sizeof(Cell);
    option("MBS_GPU_PINNED_UPLOAD","1");option("MBS_GPU_TIMING","1");
    option("MBS_GPU_UPLOAD_RESERVE_GROWTH","0");
    {Budget budget(true);budget.limit=16*mib;
        {Buffer b(&budget);b.resize(initial);fill(b);auto* original=b.data();assert(budget.used==mib);
            const size_t calls=MockRegistration::unregistrations;b.resize(initial/2);assert(b.data()==original);assert(MockRegistration::unregistrations==calls);prefix(b,initial/2);
            b.resize(12000);assert(MockRegistration::unregistrations==calls+1);prefix(b,initial/2);b.data();assert(budget.used==b.capacity_bytes());
            assert(budget.timings.reserveCalls==0&&budget.timings.growthCalls==2);}
        assert(budget.used==0);}
    clean();
    option("MBS_GPU_UPLOAD_RESERVE_GROWTH","1");
    {Budget budget(true);budget.limit=8*mib;
        {Buffer b(&budget);b.resize(initial);fill(b);b.data();b.resize(12000);prefix(b,initial);
            assert(b.capacity_bytes()==48000*sizeof(Cell));auto* original=b.data();const auto calls=MockRegistration::unregistrations;
            b.resize(20000);assert(b.data()==original&&MockRegistration::unregistrations==calls);prefix(b,initial);
            assert(budget.used<=budget.limit&&budget.timings.reserveCalls==1);}
        assert(budget.used==0);}
    clean();
    {Budget budget(true);budget.limit=8*mib;
        {Buffer a(&budget),b(&budget),c(&budget);a.resize(2*initial);a.data();b.resize(initial);b.data();assert(budget.used==3*mib);
            a.resize(24000);a.data();assert(a.capacity_bytes()+b.capacity_bytes()<=budget.limit);assert(budget.used==8*mib);
            const auto calls=MockRegistration::registrations;c.resize(initial);assert(c.data());assert(MockRegistration::registrations==calls&&budget.used==8*mib);}
        assert(budget.used==0);}
    clean();
    {Budget budget(true);budget.limit=8*mib;
        {Buffer b(&budget);b.resize(initial);b.data();option("MBS_GPU_PINNED_UPLOAD","0");b.resize(12000);b.data();
            assert(budget.used==0&&budget.timings.reserveCalls==0&&b.capacity_bytes()<48000*sizeof(Cell));}
        option("MBS_GPU_PINNED_UPLOAD","1");}
    clean();
    {Budget budget(true);budget.limit=8*mib;
        {Buffer b(&budget);b.resize(initial);fill(b);MockRegistration::failRegister=true;const auto clears=MockRegistration::clears;
            assert(b.data());prefix(b,initial);assert(budget.used==0&&MockRegistration::clears==clears+1);
            MockRegistration::failRegister=false;const auto calls=MockRegistration::registrations;b.data();assert(MockRegistration::registrations==calls);}
        assert(budget.used==0);}
    clean();
    {Budget budget(true);budget.limit=8*mib;
        {Buffer b(&budget);b.resize(initial);fill(b);auto* original=b.data();MockRegistration::failUnregister=true;
            bool threw=false;try{b.resize(12000);}catch(const std::runtime_error&){threw=true;}
            assert(threw&&b.size()==initial&&b.capacity_bytes()==mib&&b.data()==original&&budget.used==mib);prefix(b,initial);
            option("MBS_GPU_PINNED_UPLOAD","0");threw=false;try{b.data();}catch(const std::runtime_error&){threw=true;}
            assert(threw&&budget.used==mib);option("MBS_GPU_PINNED_UPLOAD","1");MockRegistration::failUnregister=false;
            // Corrupt accounting must not wrap size_t or authorize freeing a pin.
            budget.used=0;const auto calls=MockRegistration::unregistrations;threw=false;try{b.resize(12000);}catch(const std::runtime_error&){threw=true;}
            assert(threw&&MockRegistration::unregistrations==calls&&b.size()==initial);budget.used=mib;}
        assert(budget.used==0);}
    clean();
    // The intentionally retained error allocation is owned until this child
    // exits; no test leak is carried into the other cases/process.
    const pid_t child=fork();assert(child>=0);
    if(!child){Budget budget(true);budget.limit=8*mib;Buffer* b=new Buffer(&budget);b->resize(initial);fill(*b);auto* pointer=b->data();
        MockRegistration::failUnregister=true;delete b;assert(Cell::alive==initial&&budget.used==mib&&MockRegistration::active.count(pointer));
        assert(pointer[initial-1].value[15]==(initial-1)*37+15);_exit(0);}
    int status=0;assert(waitpid(child,&status,0)==child&&WIFEXITED(status)&&WEXITSTATUS(status)==0);clean();
    assert(MbsGpuUpload::GrowthReserve(12,8,800)==48);
    assert(MbsGpuUpload::GrowthReserve(40,8,800)==100);
    assert(MbsGpuUpload::GrowthReserve(101,8,800)==101);
    assert(MbsGpuUpload::GrowthReserve(std::numeric_limits<size_t>::max(),8,std::numeric_limits<size_t>::max())==std::numeric_limits<size_t>::max());
    option("MBS_HOST_MEM_BUDGET_MB","64");Budget bounded(true);assert(bounded.limit==4*mib);
    std::puts("PASS GPU upload lifetime, opt-out, shared reserve bound, failure handling and staged-value preservation");
}
