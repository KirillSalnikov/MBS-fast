#pragma once
#include <algorithm>

namespace FullAutoStorageBudget {
// Limit, not preallocation. Leave three quarters of currently available RAM
// for native workers, the OS and other applications. Divide that allowance
// by visible GPUs so independent controllers do not each claim the same share.
inline int DefaultMiB(long long availableKiB,int shares) {
    if(availableKiB<=0)return 512;
    const long long share=availableKiB/1024/4/std::max(1,shares);
    return int(std::max(512LL,std::min(16384LL,share)));
}
}
