#include "../src/cuda/GpuSupport.h"
#include "../src/cuda/GpuAnalyticFacet.h"

bool EvaluateAnalyticFacetGpu(const AnalyticGpuModel&,
    const std::vector<std::array<AnalyticBackscatter::Vec,3>>&,
    const std::vector<double>&,const std::vector<double>&,int,
    std::vector<double>&,std::vector<double>&,
    std::vector<AnalyticFacetAverage::Components>*){return false;}

bool EvaluateAnalyticAzimuthGaussianGpu(
    const std::vector<std::vector<AnalyticAzimuthGaussian::Beam>>&,
    const std::vector<double>&,const std::vector<double>&,int,
    std::vector<double>&,std::vector<double>&,
    std::vector<AnalyticAzimuthGaussian::Values>*){return false;}

bool CheckGpuRuntime(GpuDeviceInfo &/*info*/, std::string &error)
{
    error = "CPU MPI/OpenMP binary was built without CUDA support; use gpu/bin/mbs_po_gpu_float_fast";
    return false;
}

bool QueryActiveGpuMemory(long long &freeBytes, long long &totalBytes,
                          std::string &error)
{
    freeBytes = 0;
    totalBytes = 0;
    error = "CPU binary has no active CUDA device";
    return false;
}

int VisibleGpuDeviceCount()
{
    return 0;
}

std::string FormatGpuInfo(const GpuDeviceInfo &/*info*/)
{
    return "CUDA unavailable in CPU MPI/OpenMP binary";
}
