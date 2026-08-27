#pragma once

#include <vector>

#include "Beam.h"
#include "Facet.h"
#include "geometry_lib.h"

struct GpuTraceExactHit
{
    bool evaluated = false;
    int candidateIndex = -1;
    int facetId = -1;
    Polygon intersection;
    bool partitionEvaluated = false;
    bool partitionOverflow = false;
    double partitionOverlapArea = 0.0;
    Polygon reachedIntersection;
    std::vector<Polygon> remaining;

    void Reset()
    {
        evaluated = false;
        candidateIndex = -1;
        facetId = -1;
        intersection.Clear();
        partitionEvaluated = false;
        partitionOverflow = false;
        partitionOverlapArea = 0.0;
        reachedIntersection.Clear();
        remaining.clear();
    }
};

struct GpuTraceBeamFacets
{
    const Beam *beam;
    IntArray *facetIds;
    std::vector<unsigned char> *mayIntersect;
    bool *sortedOnGpu;
    GpuTraceExactHit *exactHit = nullptr;
};

bool GpuTracePrefilterBeamFacets(const Beam &beam,
                                 const Facet *facets,
                                 const IntArray &facetIds,
                                 std::vector<unsigned char> &mayIntersect);

bool GpuTracePrepareBeamFacetBatch(const Facet *facets,
                                   double geometryScale,
                                   const std::vector<GpuTraceBeamFacets> &items);

// Bind the calling host worker to its stable CUDA tracing device.  The
// assignment is round-robin across the devices selected by MBS_GPU_MULTI and
// MBS_GPU_MULTI_MAX, matching the diffraction backend's device policy.
int GpuTraceBindWorkerDevice();
int GpuTraceWorkerDeviceCount();

// Reports whether the immediately preceding batch preparation on this host
// thread failed because CUDA could not reserve device memory (including
// kernel stack backing). Other unsupported cases must keep using the CPU
// fallback instead of retrying smaller batches.
bool GpuTraceLastFailureWasOutOfMemory();

void GpuTraceInvalidateFacetCache();
