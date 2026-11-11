#pragma once
// Host-side upload storage. The policy adapter makes CUDA registration lifetime
// testable without a GPU. Physical beam data and ordering are not changed here.
#include <algorithm>
#include <chrono>
#include <cstdio>
#include <cstdlib>
#include <limits>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

namespace MbsGpuUpload {
struct Timings {
    double unregisterMs=0,resizeMs=0,registerMs=0;
    size_t unregisterCalls=0,resizeCalls=0,growthCalls=0,registerCalls=0,registerFailures=0,reserveCalls=0;
};
inline bool Enabled(const char* name){const char* p=std::getenv(name);return p&&p[0]=='1'&&p[1]=='\0';}
inline double NowMs(){return std::chrono::duration<double,std::milli>(std::chrono::steady_clock::now().time_since_epoch()).count();}
struct Budget {
    size_t used=0,limit=256ULL<<20;
    bool defaultPinned;
    Timings timings;
    explicit Budget(bool pinByDefault):defaultPinned(pinByDefault){
        if(const char* value=std::getenv("MBS_HOST_MEM_BUDGET_MB")){
            char* end=nullptr;const unsigned long long mib=std::strtoull(value,&end,10);
            if(value[0]!='-'&&end!=value&&*end=='\0'&&mib>0)limit=size_t(std::min(256ULL,mib/16))<<20;
        }
    }
    bool PinEnabled()const{
        const char* value=std::getenv("MBS_GPU_PINNED_UPLOAD");
        return value ? value[0]=='1'&&value[1]=='\0':defaultPinned;
    }
};
// A reserve is optional and only used after an already registered buffer has
// actually outgrown its capacity. Bounds are in elements to avoid overflow.
inline size_t GrowthReserve(size_t requested,size_t elementBytes,size_t availableBytes){
    if(!elementBytes)return requested;
    const size_t maximum=availableBytes/elementBytes;
    if(requested>maximum)return requested;
    return requested>maximum/4 ? maximum:requested*4;
}

template<class T,class Registration> class Vector {
    std::unique_ptr<std::vector<T>> values{new std::vector<T>};
    Budget* budget;
    void* registered=nullptr;
    size_t lockedBytes=0;
    bool registrationFailed=false;
    bool FailUnregister(bool throwOnError,const char* detail){
        if(throwOnError)throw std::runtime_error(std::string("CUDA upload unregister failed: ")+detail+"; registered host allocation retained");
        std::fprintf(stderr,"CUDA upload unregister failed: %s; registered host allocation retained\n",detail);return false;
    }
    bool Unregister(bool throwOnError=true){
        if(!registered)return true;
        if(lockedBytes>budget->used)return FailUnregister(throwOnError,"pin-budget accounting underflow");
        const bool timing=Enabled("MBS_GPU_TIMING");const double started=timing?NowMs():0;
        const int result=Registration::Unregister(registered);
        if(timing){budget->timings.unregisterMs+=NowMs()-started;++budget->timings.unregisterCalls;}
        if(result)return FailUnregister(throwOnError,Registration::Error(result));
        budget->used-=lockedBytes;registered=nullptr;lockedBytes=0;return true;
    }
public:
    explicit Vector(Budget* owner):budget(owner){}
    Vector(const Vector&)=delete;Vector& operator=(const Vector&)=delete;
    ~Vector()noexcept{
        // A failed CUDA unregistration does not authorize releasing its host
        // allocation. Leak its owner until process teardown rather than free
        // memory the driver may still access. This is only an error path.
        if(!Unregister(false))values.release();
    }
    void resize(size_t count){
        const bool grows=count>values->capacity(),wasRegistered=registered!=nullptr;
        if(grows)Unregister(); // Failure leaves the old pointer/size untouched.
        const bool timing=Enabled("MBS_GPU_TIMING");const double started=timing?NowMs():0;
        if(grows&&wasRegistered&&!registrationFailed&&budget->PinEnabled()
            &&Enabled("MBS_GPU_UPLOAD_RESERVE_GROWTH")&&budget->used<=budget->limit){
            const size_t target=GrowthReserve(count,sizeof(T),budget->limit-budget->used);
            if(target>count&&target<=values->max_size()){
                values->reserve(target);if(timing)++budget->timings.reserveCalls;
            }
        }
        values->resize(count);
        if(timing){budget->timings.resizeMs+=NowMs()-started;++budget->timings.resizeCalls;
            if(grows)++budget->timings.growthCalls;}
    }
    size_t size()const{return values->size();}
    bool empty()const{return values->empty();}
    size_t capacity_bytes()const{return values->capacity()*sizeof(T);}
    T& operator[](size_t i){return (*values)[i];}
    T* data(){
        const char* mode=std::getenv("MBS_GPU_PINNED_UPLOAD");const bool enabled=budget->PinEnabled();
        if(!enabled)Unregister();
        const size_t bytes=capacity_bytes(),minimum=(mode?1ULL:8ULL)<<20;
        if(!registered&&!registrationFailed&&enabled&&!values->empty()&&bytes>=minimum
            &&budget->used<=budget->limit&&bytes<=budget->limit-budget->used){
            const bool timing=Enabled("MBS_GPU_TIMING");const double started=timing?NowMs():0;
            const int result=Registration::Register(values->data(),bytes);
            if(timing){budget->timings.registerMs+=NowMs()-started;++budget->timings.registerCalls;
                if(result)++budget->timings.registerFailures;}
            if(!result){registered=values->data();lockedBytes=bytes;budget->used+=bytes;}
            else{registrationFailed=true;Registration::ClearError();}
        }
        return values->data();
    }
};
}
