#pragma once

#include <string>

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

RuntimeResourceSnapshot QueryRuntimeResourceSnapshot();
std::string FormatRuntimeResourceReport(const std::string &stage,
                                        bool includeActiveGpu);
