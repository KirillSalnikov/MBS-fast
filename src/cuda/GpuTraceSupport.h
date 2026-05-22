#pragma once

#include <vector>

#include "Beam.h"
#include "Facet.h"
#include "geometry_lib.h"

struct GpuTraceBeamFacets
{
    const Beam *beam;
    const IntArray *facetIds;
    std::vector<unsigned char> *mayIntersect;
};

bool GpuTracePrefilterBeamFacets(const Beam &beam,
                                 const Facet *facets,
                                 const IntArray &facetIds,
                                 std::vector<unsigned char> &mayIntersect);

bool GpuTracePrefilterBeamFacetBatch(const Facet *facets,
                                     const std::vector<GpuTraceBeamFacets> &items);

void GpuTraceInvalidateFacetCache();
