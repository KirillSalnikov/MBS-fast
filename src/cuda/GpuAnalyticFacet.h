#pragma once
#include "AnalyticFacetAverage.h"
#include "AnalyticAzimuthGaussian.h"
#include <array>

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
bool EvaluateAnalyticFacetGpu(const AnalyticGpuModel &model,
    // Each frame contains transverse azimuth axes followed by the direction
    // toward the source, all in particle coordinates. Incidence is arbitrary;
    // frames must be orthonormal and right handed (see SourceFrame).
    const std::vector<std::array<AnalyticBackscatter::Vec,3>>&frames,
    const std::vector<double>&weights,const std::vector<double>&theta,int nPhi,
    std::vector<double>&reflectedCells,std::vector<double>&shadowCells,
    std::vector<AnalyticFacetAverage::Components>*orientationRows);
