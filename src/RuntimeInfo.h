#pragma once

#include <string>
#include <utility>
#include <vector>

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

RuntimeResourceSnapshot QueryRuntimeResourceSnapshot();
RuntimeEnvironmentSnapshot QueryRuntimeEnvironmentSnapshot();
std::string FormatRuntimeResourceReport(const std::string &stage,
                                        bool includeActiveGpu);
std::string FormatRuntimeEnvironmentReport(
    const RuntimeEnvironmentSnapshot &snapshot);
