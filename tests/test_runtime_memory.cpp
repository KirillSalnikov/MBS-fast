#include "../src/RuntimeInfo.h"
#include <algorithm>
#include <chrono>
#include <cmath>
#include <iostream>
#include <stdexcept>
#include <sys/mman.h>
#include <unistd.h>

static void Require(bool ok, const char* message)
{
    if (!ok) throw std::runtime_error(message);
}

template<class Query> double TimeQueries(Query query, int count)
{
    volatile long long observed = 0;
    const auto begin = std::chrono::steady_clock::now();
    for (int i = 0; i < count; ++i) observed += query().processRssKb;
    Require(observed > 0, "missing process memory during benchmark");
    return std::chrono::duration<double>(
        std::chrono::steady_clock::now() - begin).count();
}

int main()
{
    try {
        const auto complete = QueryRuntimeResourceSnapshot();
        const auto before = QueryRuntimeMemorySnapshot();
        Require(before.ramTotalKb > 0 && before.ramTotalKb == complete.ramTotalKb,
                "memory-only query changed total RAM");
        Require(before.ramAvailableKb > 0 && before.ramAvailableKb <= before.ramTotalKb,
                "invalid current available RAM");
        Require(before.processRssKb > 0 && before.processPeakRssKb >= before.processRssKb,
                "invalid current RSS or peak RSS");
        Require(complete.physicalCores > 0 && complete.logicalProcessors > 0
                && !complete.cpuModel.empty(), "full resource report lost CPU fields");

        const size_t allocation = 64 * 1024 * 1024;
        void* block = mmap(nullptr, allocation, PROT_READ | PROT_WRITE,
                           MAP_PRIVATE | MAP_ANONYMOUS, -1, 0);
        Require(block != MAP_FAILED, "cannot allocate memory freshness fixture");
        auto bytes = static_cast<volatile unsigned char*>(block);
        const size_t page = size_t(sysconf(_SC_PAGESIZE));
        for (size_t i = 0; i < allocation; i += page) bytes[i] = 1;
        const auto allocated = QueryRuntimeMemorySnapshot();
        Require(munmap(block, allocation) == 0, "cannot release fixture mapping");
        const auto released = QueryRuntimeMemorySnapshot();
        // Linux RSS accounting can lag slightly; a 32 MiB margin distinguishes
        // a live query from cached counters without relying on exact page totals.
        Require(allocated.processRssKb - before.processRssKb > 32 * 1024,
                "RSS query did not observe allocation");
        Require(allocated.processRssKb - released.processRssKb > 32 * 1024,
                "RSS query did not observe release");
        Require(released.processPeakRssKb + 1024 >= allocated.processRssKb,
                "peak RSS lost high-water mark");

        std::vector<double> full, memory;
        const int count = 200;
        for (int trial = 0; trial < 7; ++trial) {
            if (trial % 2) {
                memory.push_back(TimeQueries(QueryRuntimeMemorySnapshot, count));
                full.push_back(TimeQueries(QueryRuntimeResourceSnapshot, count));
            } else {
                full.push_back(TimeQueries(QueryRuntimeResourceSnapshot, count));
                memory.push_back(TimeQueries(QueryRuntimeMemorySnapshot, count));
            }
        }
        std::sort(full.begin(), full.end());
        std::sort(memory.begin(), memory.end());
        std::cout << "{\"passed\":true,\"queries_per_trial\":" << count
                  << ",\"trials\":7,\"full_query_median_seconds\":" << full[3]
                  << ",\"memory_query_median_seconds\":" << memory[3]
                  << ",\"query_speedup\":" << full[3] / memory[3]
                  << ",\"rss_before_kib\":" << before.processRssKb
                  << ",\"rss_allocated_kib\":" << allocated.processRssKb
                  << ",\"rss_released_kib\":" << released.processRssKb
                  << ",\"cpu_logical_processors\":" << complete.logicalProcessors
                  << ",\"scope\":\"memory monitoring only; no native physics benchmark\"}\n";
    } catch (const std::exception& error) {
        std::cerr << error.what() << '\n';
        return 1;
    }
}
