#pragma once

#include <string>
#include <utility>
#include <vector>

// Live memory counters only. Short-lived native requests must not rescan CPU
// topology merely to enforce their host-memory budget or decide whether to trim.
struct RuntimeMemorySnapshot
{
    long long ramTotalKb = 0;
    long long ramAvailableKb = 0;
    long long processRssKb = 0;
    long long processPeakRssKb = 0;
};

struct RuntimeResourceSnapshot
{
    std::string hostname;
    std::string cpuModel;
    int logicalProcessors = 0;
    int physicalCores = 0;
    int openmpMaxThreads = 1;
    long long ramTotalKb = 0;
    long long ramAvailableKb = 0;
    long long processRssKb = 0;
    long long processPeakRssKb = 0;
};

struct RuntimeEnvironmentSnapshot
{
    std::vector<std::pair<std::string, std::string>> active;
    std::vector<std::string> experimental;
};

RuntimeMemorySnapshot QueryRuntimeMemorySnapshot();
RuntimeResourceSnapshot QueryRuntimeResourceSnapshot();
RuntimeEnvironmentSnapshot QueryRuntimeEnvironmentSnapshot();
std::string FormatRuntimeResourceReport(const std::string &stage,
                                        bool includeActiveGpu);
std::string FormatRuntimeEnvironmentReport(
    const RuntimeEnvironmentSnapshot &snapshot);
