#include "TraceCachePolicy.h"
#include "FullAutoPhiTileMemory.h"
#include "cuda/GpuTraceCapacity.h"
#include <cassert>
#include <cstdlib>
#include <cstdio>
#include <limits>

static void expect(bool persistent, bool explicitLimit, long long requested,
                   long long ramKiB, bool budgetPresent, long long budgetMiB,
                   bool pool, long long expectedHost, long long expectedGpu) {
    const auto limits = TraceCachePolicy::Resolve(persistent, explicitLimit,
        requested, ramKiB, budgetPresent, budgetMiB, pool);
    assert(limits.hostMiB == expectedHost && limits.gpuRawMiB == expectedGpu);
}
int main() {
    const std::string memoryReport =
        "Orientation average adaptive memory: measured 1 MiB / 1 orientations, conservative 16000 KiB/orientation, budget 10 MiB\n"
        "Orientation average adaptive memory: measured 1 MiB / 1 orientations, conservative 32000 KiB/orientation, budget 10 MiB\n";
    assert(FullAutoPhiTileMemory::Peak(memoryReport) == 32000);
    assert(FullAutoPhiTileMemory::Peak("unrelated conservative 32000 KiB/orientation") == 0);
    assert(FullAutoPhiTileMemory::Fit(256, 32000, 16) == 8);
    assert(FullAutoPhiTileMemory::Fit(256, 16000, 16) == 16);
    const long long ampleRam = 64LL * 1024 * 1024;
    // Implicit persistent CPU preparation stays 256MiB. Only GPU raw storage
    // changes under the opt-in; a standalone native run still defaults to 0.
    expect(true, false, 0, ampleRam, false, 0, false, 256, 256);
    expect(true, false, 0, ampleRam, false, 0, true, 256, 0);
    expect(false, false, 0, ampleRam, false, 0, false, 0, 0);
    expect(false, false, 0, ampleRam, false, 0, true, 0, 0);
    for (bool persistent : {false, true}) {
        for (bool pool : {false, true}) {
            expect(persistent, true, 0, ampleRam, false, 0, pool, 0, 0);
            expect(persistent, true, 512, ampleRam, false, 0, pool, 512, 512);
            // Explicit settings retain both original host-memory clamps.
            expect(persistent, true, 512, 1024LL * 1024, false, 0, pool, 32, 32);
            expect(persistent, true, 512, ampleRam, true, 1024, pool, 64, 64);
            expect(persistent, true, 512, 1024LL * 1024, true, 256, pool, 16, 16);
            expect(persistent, true, 512, ampleRam, true, 0, pool, 0, 0);
            expect(persistent, true, 512, -1, false, 0, pool, 0, 0);
        }
    }
    expect(true, false, 9999, ampleRam, true, 1024, true, 64, 0);
    expect(true, false, 9999, ampleRam, true, 1024, false, 64, 64);
    expect(true, false, 0, 1024LL * 1024, false, 0, true, 32, 0);
    expect(true, false, 0, ampleRam, true, -1, true, 0, 0);
    expect(true, false, 0, std::numeric_limits<long long>::max(), false, 0, true, 256, 0);
    // Match the exact production opt-in parser. Other spellings do not enable
    // the pool and must not unexpectedly change legacy GPU cache defaults.
    for (const char* value : {"0", "true", "yes", "01", "1x", ""}) {
        setenv("MBS_GPU_TRACE_BUFFER_POOL", value, 1);
        const bool enabled = GpuConvexQueue::Capacity::BufferPoolEnabled();
        assert(!enabled);
        expect(true, false, 0, ampleRam, false, 0, enabled, 256, 256);
    }
    setenv("MBS_GPU_TRACE_BUFFER_POOL", "1", 1);
    assert(GpuConvexQueue::Capacity::BufferPoolEnabled());
    expect(true, false, 0, ampleRam, false, 0,
        GpuConvexQueue::Capacity::BufferPoolEnabled(), 256, 0);
    unsetenv("MBS_GPU_TRACE_BUFFER_POOL");
    assert(!GpuConvexQueue::Capacity::BufferPoolEnabled());
    std::puts("PASS CPU/GPU cache split, explicit priority, RAM clamps, persistent/native and exact opt-in flags");
}
