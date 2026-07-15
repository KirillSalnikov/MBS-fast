#include "GpuSupport.h"

#include <cstdlib>
#include <sstream>

#ifdef USE_CUDA
#include <cuda_runtime.h>
#endif

#ifndef MBS_GPU_BUILD_ARCH
#define MBS_GPU_BUILD_ARCH 0
#endif

bool CheckGpuRuntime(GpuDeviceInfo &info, std::string &error)
{
#ifndef USE_CUDA
    error = "binary was built without CUDA support; rebuild with make USE_CUDA=1";
    return false;
#else
    auto getEnvInt = [](const char *name, int fallback) {
        const char *value = std::getenv(name);
        if (!value || !*value)
            return fallback;
        char *end = nullptr;
        long parsed = std::strtol(value, &end, 10);
        return (end && *end == '\0' && parsed >= 0) ? (int)parsed : fallback;
    };

    int count = 0;
    cudaError_t err = cudaGetDeviceCount(&count);
    if (err != cudaSuccess)
    {
        error = std::string("cudaGetDeviceCount failed: ") + cudaGetErrorString(err);
        return false;
    }
    if (count <= 0)
    {
        error = "no CUDA devices found";
        return false;
    }

    int localRank = getEnvInt("OMPI_COMM_WORLD_LOCAL_RANK",
                    getEnvInt("MPI_LOCALRANKID",
                    getEnvInt("SLURM_LOCALID", 0)));
    int deviceId = localRank % count;

    info.visibleDevices.clear();
    for (int index = 0; index < count; ++index)
    {
        cudaDeviceProp visibleProp;
        err = cudaGetDeviceProperties(&visibleProp, index);
        if (err != cudaSuccess)
        {
            error = std::string("cudaGetDeviceProperties failed for device ")
                + std::to_string(index) + ": " + cudaGetErrorString(err);
            return false;
        }
        GpuVisibleDeviceInfo visible;
        visible.deviceId = index;
        visible.name = visibleProp.name;
        visible.totalGlobalMem = (long long)visibleProp.totalGlobalMem;
        visible.computeMajor = visibleProp.major;
        visible.computeMinor = visibleProp.minor;
        info.visibleDevices.push_back(visible);
    }

    cudaDeviceProp prop;
    err = cudaGetDeviceProperties(&prop, deviceId);
    if (err != cudaSuccess)
    {
        error = std::string("cudaGetDeviceProperties failed: ") + cudaGetErrorString(err);
        return false;
    }

    size_t freeBytes = 0;
    size_t totalBytes = 0;
    err = cudaMemGetInfo(&freeBytes, &totalBytes);
    if (err != cudaSuccess)
    {
        error = std::string("cudaMemGetInfo failed: ") + cudaGetErrorString(err);
        return false;
    }
    err = cudaSetDevice(deviceId);
    if (err != cudaSuccess)
    {
        error = std::string("cudaSetDevice failed: ") + cudaGetErrorString(err);
        return false;
    }
    err = cudaFree(0);
    if (err != cudaSuccess)
    {
        error = std::string("CUDA context initialization failed: ") + cudaGetErrorString(err);
        return false;
    }

    int runtimeVersion = 0;
    int driverVersion = 0;
    cudaRuntimeGetVersion(&runtimeVersion);
    cudaDriverGetVersion(&driverVersion);

    info.deviceId = deviceId;
    info.name = prop.name;
    info.totalGlobalMem = (long long)prop.totalGlobalMem;
    info.freeGlobalMem = (long long)freeBytes;
    info.computeMajor = prop.major;
    info.computeMinor = prop.minor;
    info.runtimeVersion = runtimeVersion;
    info.driverVersion = driverVersion;
    info.buildArch = MBS_GPU_BUILD_ARCH;
    info.visibleDeviceCount = count;
    return true;
#endif
}

bool QueryActiveGpuMemory(long long &freeBytes, long long &totalBytes,
                          std::string &error)
{
#ifndef USE_CUDA
    freeBytes = 0;
    totalBytes = 0;
    error = "CUDA support is not compiled into this binary";
    return false;
#else
    size_t freeValue = 0;
    size_t totalValue = 0;
    const cudaError_t err = cudaMemGetInfo(&freeValue, &totalValue);
    if (err != cudaSuccess)
    {
        error = std::string("cudaMemGetInfo failed: ") + cudaGetErrorString(err);
        return false;
    }
    freeBytes = (long long)freeValue;
    totalBytes = (long long)totalValue;
    error.clear();
    return true;
#endif
}

std::string FormatGpuInfo(const GpuDeviceInfo &info)
{
    std::ostringstream out;
    out << "CUDA device " << info.deviceId << ": " << info.name
        << ", cc " << info.computeMajor << "." << info.computeMinor
        << ", memory " << (info.totalGlobalMem / (1024LL * 1024LL)) << " MiB"
        << ", free " << (info.freeGlobalMem / (1024LL * 1024LL)) << " MiB";
    if (info.visibleDeviceCount > 1)
        out << ", visible devices " << info.visibleDeviceCount
            << " (auto multi-GPU)";
    if (info.buildArch > 0)
        out << ", build sm_" << info.buildArch;
    if (info.runtimeVersion > 0 || info.driverVersion > 0)
        out << ", CUDA runtime " << info.runtimeVersion
            << ", driver " << info.driverVersion;
    if (!info.visibleDevices.empty())
    {
        out << "\nVisible CUDA devices:";
        for (const GpuVisibleDeviceInfo &device : info.visibleDevices)
        {
            out << "\n  [" << device.deviceId << "] " << device.name
                << ", cc " << device.computeMajor << "." << device.computeMinor
                << ", " << device.totalGlobalMem / (1024LL * 1024LL)
                << " MiB";
        }
    }
    return out.str();
}
