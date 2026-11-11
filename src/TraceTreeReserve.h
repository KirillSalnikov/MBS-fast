#pragma once
#include <cstdlib>
#include <cstddef>

namespace TraceTreeReserve {
// Initial capacity only. std::vector growth and physical trace limits are
// unchanged. Preserve the larger reserve for concave/aggregate geometries.
inline size_t Count(const char* value, bool convex = false) {
    if(value && *value) {
        char* end=nullptr;
        const long count=std::strtol(value,&end,10);
        if(end && *end=='\0' && count>0)return static_cast<size_t>(count);
    }
    return convex ? 4096u : 65536u;
}
}
