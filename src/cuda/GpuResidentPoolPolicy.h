#pragma once
#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <limits>

// Allocation-free accounting for the experimental resident preparation cache.
// Bytes describe capacities only; no optical result is cached.
namespace GpuResidentPoolPolicy {
static const size_t Slots=14;
typedef std::array<size_t,Slots> Sizes;
inline bool Add(size_t a,size_t b,size_t& result){
    if(b>std::numeric_limits<size_t>::max()-a)return false;
    result=a+b;return true;
}
inline bool Multiply(size_t a,size_t b,size_t& result){
    if(a&&b>std::numeric_limits<size_t>::max()/a)return false;
    result=a*b;return true;
}
inline bool Sum(const Sizes& values,size_t& result){
    result=0;for(size_t value:values)if(!Add(result,value,result))return false;return true;
}
struct Plan {
    Sizes capacities={};
    size_t budget=0,requestedBytes=0,retainedBytes=0;
    bool cache=false;
};
inline bool MakePlan(const Sizes& requested,const Sizes& old,size_t freeBytes,double fraction,Plan& result){
    if(!std::isfinite(fraction)||fraction<=0||fraction>1)return false;
    size_t cached=0,available=0;
    Plan plan;
    if(!Sum(old,cached)||!Add(freeBytes,cached,available)||!Sum(requested,plan.requestedBytes))return false;
    const size_t fixedCap=size_t(1024)*1024*1024;
    // Compute the fractional bound after limiting the scale to avoid rounding
    // a SIZE_MAX-sized floating value back outside size_t's range.
    const long double fractionBytes=static_cast<long double>(available)*fraction/2;
    const size_t share=fractionBytes>=fixedCap?fixedCap:static_cast<size_t>(fractionBytes);
    plan.budget=std::min(fixedCap,std::min(available/8,share));
    if(plan.requestedBytes>plan.budget){result=plan;return true;}
    plan.cache=true;
    for(size_t i=0;i<Slots;++i)plan.capacities[i]=std::max(old[i],requested[i]);
    if(!Sum(plan.capacities,plan.retainedBytes))return false;
    if(plan.retainedBytes>plan.budget){plan.capacities=requested;plan.retainedBytes=plan.requestedBytes;}
    result=plan;return true;
}

// Small ownership model shared by CUDA and mock-allocation tests. A failed
// free keeps its pointer/accounting and poisons reuse; it is never forgotten.
struct Storage {
    struct Slot {void* pointer=nullptr;size_t bytes=0;};
    std::array<Slot,Slots> slots;
    bool poisoned=false;
    size_t allocations=0,releases=0,reuses=0;
    Sizes Capacities()const{Sizes sizes={};for(size_t i=0;i<Slots;++i)sizes[i]=slots[i].bytes;return sizes;}
    template<class Memory>bool Release(size_t i,Memory& memory){
        auto& slot=slots[i];
        if(slot.pointer&&!memory.Free(slot.pointer)){poisoned=true;return false;}
        if(slot.pointer)++releases;
        slot.pointer=nullptr;slot.bytes=0;return true;
    }
    template<class Memory>bool Clear(Memory& memory){
        bool ok=true;for(size_t i=0;i<Slots;++i)if(!Release(i,memory))ok=false;
        if(ok)poisoned=false;
        return ok;
    }
    template<class Memory>bool Apply(const Plan& plan,Memory& memory){
        if(poisoned)return false;
        // Free every superseded allocation before growing any replacement.
        for(size_t i=0;i<Slots;++i)if(slots[i].bytes!=plan.capacities[i]&&!Release(i,memory))return false;
        return true;
    }
    template<class Memory>bool Get(size_t i,size_t bytes,Memory& memory,void*& pointer){
        pointer=nullptr;
        if(poisoned||i>=Slots)return false;
        auto& slot=slots[i];
        if(slot.pointer){if(slot.bytes<bytes)return false;++reuses;pointer=slot.pointer;return true;}
        if(!bytes)return false;
        void* fresh=nullptr;if(!memory.Allocate(&fresh,bytes))return false;
        slot.pointer=fresh;slot.bytes=bytes;++allocations;pointer=fresh;return true;
    }
};
}
