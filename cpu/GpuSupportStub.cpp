#include "../src/cuda/GpuSupport.h"

bool CheckGpuRuntime(GpuDeviceInfo &/*info*/, std::string &error)
{
    error = "CPU MPI/OpenMP binary was built without CUDA support; use gpu/bin/mbs_po_gpu_float_fast";
    return false;
}

std::string FormatGpuInfo(const GpuDeviceInfo &/*info*/)
{
    return "CUDA unavailable in CPU MPI/OpenMP binary";
}
