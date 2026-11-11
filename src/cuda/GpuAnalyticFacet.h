#pragma once
#include "AnalyticFacetAverage.h"
#include "AnalyticAzimuthGaussian.h"
#include <array>
#include <cstddef>
#include <memory>
#include <string>

namespace AnalyticGpuWorkspace
{
// Capacities and requests are in doubles. Account for buffers that must grow
// as well as retained larger buffers; the latter alone is not a safe bound.
inline bool CanRetain(size_t budget,const std::array<size_t,6>&capacity,
    const std::array<size_t,6>&requested)
{
    size_t remaining=budget/sizeof(double);
    for(size_t i=0;i<capacity.size();++i){
        const size_t needed=capacity[i]>requested[i] ? capacity[i]:requested[i];
        if(needed>remaining)return false;
        remaining-=needed;
    }
    return true;
}
// Prefer an existing large scratch allocation to growing it with every
// orientation chunk. Bound the additional launch count, rather than assuming
// that arbitrarily small tiles are faster than allocation. The actual speed
// tradeoff is measured separately; this helper only enforces resource bounds.
inline size_t ReuseTile(size_t budget,const std::array<size_t,6>&capacity,
    size_t orientations,size_t planned,size_t sampleElements,size_t maximumTileFactor=4)
{
    if(!orientations || !planned || !maximumTileFactor)return planned;
    size_t tile=planned;
    for(size_t i=0;i<2;++i)if(capacity[i]/orientations<tile)tile=capacity[i]/orientations;
    for(size_t i=2;i<4;++i)if(capacity[i]<tile)tile=capacity[i];
    const size_t minimum=planned/maximumTileFactor+(planned%maximumTileFactor!=0);
    if(tile<minimum)return planned;
    const std::array<size_t,6>requested={{orientations*tile,orientations*tile,tile,tile,sampleElements,sampleElements}};
    return CanRetain(budget,capacity,requested) ? tile:planned;
}
}

struct AnalyticGpuPatch
{
    double normal[3],horizontal[3],vertical[3],area;
    int firstVertex,vertices;
};
struct AnalyticGpuModel
{
    std::vector<AnalyticGpuPatch> patches;
    std::vector<double> verticesXY;
    double indexReal,indexImag,wave,circleRadius;
    int shadowMode;
};
// FP64 radial preparation only; the small-z core and convergence checks stay
// on CPU. A bounded workspace is used per host thread/device. False requests
// an exact CPU fallback and never authorizes an approximate mean.
int AnalyticMeanGpuDevice();
bool EvaluateAnalyticMeanTailGpu(int device,int degree,double step,size_t first,size_t end,
    const std::vector<double>&factor,double z0,double offset,double fluxOffset,
    double*values,double*derivatives,std::string*error=nullptr);
bool EvaluateAnalyticFacetGpu(const AnalyticGpuModel &model,
    // Each frame contains transverse azimuth axes followed by the direction
    // toward the source, all in particle coordinates. Incidence is arbitrary;
    // frames must be orthonormal and right handed (see SourceFrame).
    const std::vector<std::array<AnalyticBackscatter::Vec,3>>&frames,
    const std::vector<double>&weights,const std::vector<double>&theta,int nPhi,
    std::vector<double>&reflectedCells,std::vector<double>&shadowCells,
    std::vector<AnalyticFacetAverage::Components>*orientationRows,
    std::string*error=nullptr);

// Experimental one-tile submit/collect. All methods, including destruction,
// belong to the same host thread/device. Prepare launches no optical kernels;
// unsupported tiling returns supported=false without changing any matrix.
// A submitted operation is drained by the destructor even on exceptions.
class AnalyticFacetGpuPending
{
public:
    struct State;
    AnalyticFacetGpuPending();
    ~AnalyticFacetGpuPending();
    AnalyticFacetGpuPending(const AnalyticFacetGpuPending&)=delete;
    AnalyticFacetGpuPending& operator=(const AnalyticFacetGpuPending&)=delete;
    bool Prepare(const AnalyticGpuModel &model,
        const std::vector<std::array<AnalyticBackscatter::Vec,3>>&frames,
        const std::vector<double>&weights,const std::vector<double>&theta,int nPhi,
        bool&supported,std::string*error=nullptr);
    bool Submit(std::string*error=nullptr);
    bool Collect(std::vector<double>&reflection,std::vector<double>&shadow,
        std::string*error=nullptr);
private:
    std::unique_ptr<State> state;
};
// Diagnostic environment options (not part of the physical model):
// MBS_CUDA_ANALYTIC_PROFILE=1 prints device-stage timings and transfer counts.
// MBS_CUDA_ANALYTIC_CPU_SAMPLES=1 uses the reference CPU azimuth reduction.
