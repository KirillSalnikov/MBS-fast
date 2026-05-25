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

    cudaDeviceProp prop;
    err = cudaGetDeviceProperties(&prop, deviceId);
    if (err != cudaSuccess)
    {
        error = std::string("cudaGetDeviceProperties failed: ") + cudaGetErrorString(err);
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
    info.computeMajor = prop.major;
    info.computeMinor = prop.minor;
    info.runtimeVersion = runtimeVersion;
    info.driverVersion = driverVersion;
    info.buildArch = MBS_GPU_BUILD_ARCH;
    return true;
#endif
}

std::string FormatGpuInfo(const GpuDeviceInfo &info)
{
    std::ostringstream out;
    out << "CUDA device " << info.deviceId << ": " << info.name
        << ", cc " << info.computeMajor << "." << info.computeMinor
        << ", memory " << (info.totalGlobalMem / (1024LL * 1024LL)) << " MB";
    if (info.buildArch > 0)
        out << ", build sm_" << info.buildArch;
    if (info.runtimeVersion > 0 || info.driverVersion > 0)
        out << ", CUDA runtime " << info.runtimeVersion
            << ", driver " << info.driverVersion;
    return out.str();
}
