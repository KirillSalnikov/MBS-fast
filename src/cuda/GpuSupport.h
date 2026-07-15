#pragma once

#include <string>
#include <vector>

struct GpuVisibleDeviceInfo
{
    int deviceId = -1;
    std::string name;
    long long totalGlobalMem = 0;
    int computeMajor = 0;
    int computeMinor = 0;
};

struct GpuDeviceInfo
{
    int deviceId = -1;
    std::string name;
    long long totalGlobalMem = 0;
    long long freeGlobalMem = 0;
    int computeMajor = 0;
    int computeMinor = 0;
    int runtimeVersion = 0;
    int driverVersion = 0;
    int buildArch = 0;
    int visibleDeviceCount = 0;
    std::vector<GpuVisibleDeviceInfo> visibleDevices;
};

bool CheckGpuRuntime(GpuDeviceInfo &info, std::string &error);
bool QueryActiveGpuMemory(long long &freeBytes, long long &totalBytes,
                          std::string &error);
std::string FormatGpuInfo(const GpuDeviceInfo &info);
