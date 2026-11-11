#pragma once
#include <algorithm>
#include <cstddef>
#include <cstdint>
namespace FullAutoThetaPilotPlan {
// Only a startup geometry proposal, never a production accuracy certificate.
// Keep every requested theta; bound implicit orientation/azimuth work instead.
inline int Count(int requested,size_t rows,int phi,bool explicitCount,
                 uint64_t budget=UINT64_C(16)*1024*1024){
    if(explicitCount||requested<=32||!rows||phi<=0)return requested;
    const uint64_t per=uint64_t(rows)>UINT64_MAX/uint64_t(phi)?UINT64_MAX:uint64_t(rows)*uint64_t(phi);
    int count=requested;
    while(count>32&&uint64_t(count)>budget/per)count=std::max(32,count/2);
    return count;
}
}
