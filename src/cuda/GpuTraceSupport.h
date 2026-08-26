#pragma once

#include <vector>

#include "Beam.h"
#include "Facet.h"
#include "geometry_lib.h"

struct GpuTraceBeamFacets
{
    const Beam *beam;
    IntArray *facetIds;
    std::vector<unsigned char> *mayIntersect;
    bool *sortedOnGpu;
};

bool GpuTracePrefilterBeamFacets(const Beam &beam,
                                 const Facet *facets,
                                 const IntArray &facetIds,
                                 std::vector<unsigned char> &mayIntersect);

bool GpuTracePrepareBeamFacetBatch(const Facet *facets,
                                   double geometryScale,
                                   const std::vector<GpuTraceBeamFacets> &items);

// Reports whether the immediately preceding batch preparation on this host
// thread failed because CUDA could not reserve device memory (including
// kernel stack backing). Other unsupported cases must keep using the CPU
// fallback instead of retrying smaller batches.
bool GpuTraceLastFailureWasOutOfMemory();

void GpuTraceInvalidateFacetCache();
