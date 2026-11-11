#include "cuda/GpuResidentPoolPolicy.h"
#include <cassert>
#include <iostream>
#include <cstdlib>
#include <set>

struct Memory {
    std::set<void*> owned;bool failAllocate=false,failFree=false;
    bool Allocate(void**p,size_t n){if(failAllocate)return false;*p=std::malloc(n);if(!*p)return false;owned.insert(*p);return true;}
    bool Free(void*p){assert(owned.count(p));if(failFree)return false;owned.erase(p);std::free(p);return true;}
    ~Memory(){assert(owned.empty());}
};

int main(){
    using namespace GpuResidentPoolPolicy;
    const size_t MiB=1024*1024;Sizes empty={},small={},old={};Plan p;
    small[0]=720*MiB;small[1]=16*1024;
    assert(MakePlan(small,empty,12ULL*1024*MiB,.2,p)&&p.cache&&p.budget==1024*MiB);
    assert(p.retainedBytes==p.requestedBytes);
    old=p.capacities;assert(MakePlan(small,old,10ULL*1024*MiB,.2,p)&&p.cache&&p.capacities==old);
    // Reduced pressure cap trims old oversized slots before any replacement.
    small[0]=32*MiB;
    assert(MakePlan(small,old,1024*MiB,.2,p)&&p.cache&&p.capacities==small);
    assert(p.retainedBytes<=p.budget);
    // Legal work larger than the retention cap uses ephemeral storage.
    small[0]=2048*MiB;
    assert(MakePlan(small,old,12ULL*1024*MiB,.2,p)&&!p.cache&&p.retainedBytes==0);
    assert(MakePlan(small,empty,0,.2,p)&&!p.cache&&p.budget==0);
    small[0]=64*MiB;old={};old[13]=900*MiB;
    assert(MakePlan(small,old,1024*MiB,.2,p)&&p.cache&&p.capacities[13]==0);
    small={};small[0]=128*MiB;
    assert(MakePlan(small,empty,1024*MiB,.25,p)&&p.cache&&p.budget==p.requestedBytes);
    size_t value=0;assert(!Multiply(std::numeric_limits<size_t>::max(),2,value));
    assert(!Add(std::numeric_limits<size_t>::max(),1,value));
    small[0]=std::numeric_limits<size_t>::max();small[1]=1;
    assert(!MakePlan(small,empty,1024*MiB,.2,p));
    old[0]=std::numeric_limits<size_t>::max();small={};
    assert(!MakePlan(small,old,1,.2,p));
    for(double fraction:{0.,-1.,2.,std::numeric_limits<double>::infinity(),std::numeric_limits<double>::quiet_NaN()})
        assert(!MakePlan(empty,empty,1024*MiB,fraction,p));
    Memory memory;Storage storage;void* pointer=nullptr;
    assert(storage.Get(0,128,memory,pointer)&&pointer&&storage.slots[0].bytes==128);
    void* original=pointer;assert(storage.Get(0,64,memory,pointer)&&pointer==original&&storage.reuses==1);
    memory.failAllocate=true;assert(!storage.Get(1,128,memory,pointer)&&!pointer&&storage.slots[1].bytes==0);
    memory.failAllocate=false;assert(storage.Get(1,128,memory,pointer));
    memory.failFree=true;assert(!storage.Clear(memory)&&storage.poisoned&&storage.slots[0].pointer==original);
    assert(!storage.Get(0,64,memory,pointer)&&!pointer);
    memory.failFree=false;assert(storage.Clear(memory)&&!storage.poisoned&&memory.owned.empty());
    assert(storage.Get(0,128,memory,pointer));p=Plan();p.capacities[0]=64;
    assert(storage.Apply(p,memory)&&memory.owned.empty());
    assert(storage.Get(0,64,memory,pointer)&&storage.slots[0].bytes==64);
    assert(!storage.Get(Slots,64,memory,pointer)&&!pointer);
    assert(storage.Clear(memory));
    std::cout<<"PASS bounded resident pool policy: reuse, pressure, oversize, exact fit, unused slots and overflow\n";
}
