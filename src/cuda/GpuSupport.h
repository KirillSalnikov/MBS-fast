#pragma once

#include <string>

struct GpuDeviceInfo
{
    int deviceId = -1;
    std::string name;
    long long totalGlobalMem = 0;
    int computeMajor = 0;
    int computeMinor = 0;
    int runtimeVersion = 0;
    int driverVersion = 0;
    int buildArch = 0;
    int visibleDeviceCount = 0;
};

bool CheckGpuRuntime(GpuDeviceInfo &info, std::string &error);
std::string FormatGpuInfo(const GpuDeviceInfo &info);
