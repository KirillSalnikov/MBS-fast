#pragma once
#include <algorithm>

// CPU preparation and GPU raw-beam storage serve different reuse patterns.
// This policy changes only the implicit GPU limit when the bounded trace
// workspace is enabled. Explicit CLI limits keep their existing RAM clamps.
namespace TraceCachePolicy {
struct Limits {
    long long hostMiB;
    long long gpuRawMiB;
};
inline Limits Resolve(bool persistent, bool explicitLimit, long long requestedMiB,
                      long long ramAvailableKiB, bool hasHostBudget,
                      long long hostBudgetMiB, bool traceBufferPool) {
    long long hostMiB = explicitLimit ? requestedMiB : (persistent ? 256 : 0);
    hostMiB = std::min(hostMiB, std::max(0LL, ramAvailableKiB / 1024 / 32));
    if (hasHostBudget)
        hostMiB = std::min(hostMiB, std::max(0LL, hostBudgetMiB / 16));
    return {hostMiB, traceBufferPool && !explicitLimit ? 0 : hostMiB};
}
}
