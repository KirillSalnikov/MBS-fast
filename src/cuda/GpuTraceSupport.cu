#include "GpuTraceSupport.h"

#ifdef USE_CUDA

#include <cuda_runtime.h>

#include <algorithm>
#include <cstdlib>
#include <cstdio>
#include <vector>

struct GpuTraceCandidate
{
    float beam[64][4];
    float facet[64][4];
    float dir[4];
    float normal[4];
    int beamN;
    int facetN;
};

struct GpuTraceBeamRecord
{
    double vertices[64][4];
    double dir[4];
    int nVertices;
    int location;
};

struct GpuTraceFacetRecord
{
    double vertices[64][4];
    double normalIn[4];
    double normalOut[4];
    double margin;
    int nVertices;
};

struct GpuTracePair
{
    int beam;
    int facet;
};

struct GpuTraceWorkspace
{
    GpuTraceCandidate *candidates = nullptr;
    GpuTraceBeamRecord *beams = nullptr;
    GpuTraceFacetRecord *facets = nullptr;
    GpuTracePair *pairs = nullptr;
    unsigned char *results = nullptr;
    size_t cap = 0;
    size_t beamCap = 0;
    size_t facetCap = 0;
    const Facet *facetOwner = nullptr;
    int copiedMaxFacetId = -1;
};

static thread_local GpuTraceWorkspace g_traceWorkspace;

static float gpu_trace_margin()
{
    const char *value = std::getenv("MBS_GPU_TRACE_MARGIN");
    if (!value || !*value)
        return 1.0f;
    char *end = nullptr;
    float parsed = std::strtof(value, &end);
    if (!end || *end != '\0' || parsed < 0.0f)
        return 1.0f;
    return parsed;
}

static bool gpu_trace_cache_facets()
{
    const char *value = std::getenv("MBS_GPU_TRACE_CACHE_FACETS");
    return !(value && *value == '0');
}

static bool trace_cuda_ok(cudaError_t err, const char *where)
{
    if (err == cudaSuccess)
        return true;
    std::fprintf(stderr, "CUDA trace prefilter error at %s: %s\n",
                 where, cudaGetErrorString(err));
    return false;
}

static bool ensure_trace_capacity(size_t count)
{
    if (g_traceWorkspace.cap >= count)
        return true;
    cudaFree(g_traceWorkspace.candidates);
    cudaFree(g_traceWorkspace.pairs);
    cudaFree(g_traceWorkspace.results);
    g_traceWorkspace.candidates = nullptr;
    g_traceWorkspace.pairs = nullptr;
    g_traceWorkspace.results = nullptr;
    g_traceWorkspace.cap = 0;
    if (count == 0)
        return true;
    if (cudaMalloc(&g_traceWorkspace.pairs,
                   count * sizeof(GpuTracePair)) != cudaSuccess)
        return false;
    if (cudaMalloc(&g_traceWorkspace.results,
                   count * sizeof(unsigned char)) != cudaSuccess)
        return false;
    g_traceWorkspace.cap = count;
    return true;
}

static bool ensure_trace_beam_capacity(size_t count)
{
    if (g_traceWorkspace.beamCap >= count)
        return true;
    cudaFree(g_traceWorkspace.beams);
    g_traceWorkspace.beams = nullptr;
    g_traceWorkspace.beamCap = 0;
    if (count == 0)
        return true;
    if (cudaMalloc(&g_traceWorkspace.beams,
                   count * sizeof(GpuTraceBeamRecord)) != cudaSuccess)
        return false;
    g_traceWorkspace.beamCap = count;
    return true;
}

static bool ensure_trace_facet_capacity(size_t count)
{
    if (g_traceWorkspace.facetCap >= count)
        return true;
    cudaFree(g_traceWorkspace.facets);
    g_traceWorkspace.facets = nullptr;
    g_traceWorkspace.facetCap = 0;
    g_traceWorkspace.facetOwner = nullptr;
    g_traceWorkspace.copiedMaxFacetId = -1;
    if (count == 0)
        return true;
    if (cudaMalloc(&g_traceWorkspace.facets,
                   count * sizeof(GpuTraceFacetRecord)) != cudaSuccess)
        return false;
    g_traceWorkspace.facetCap = count;
    return true;
}

static bool upload_trace_facets_if_needed(const Facet *facets, int maxFacetId)
{
    if (maxFacetId < 0)
        return false;
    if (!ensure_trace_facet_capacity((size_t)maxFacetId + 1))
        return false;
    if (gpu_trace_cache_facets()
        && g_traceWorkspace.facetOwner == facets
        && g_traceWorkspace.copiedMaxFacetId >= maxFacetId)
        return true;

    std::vector<GpuTraceFacetRecord> hostFacets((size_t)maxFacetId + 1);
    for (int facetId = 0; facetId <= maxFacetId; ++facetId)
    {
        const Facet &facet = facets[facetId];
        GpuTraceFacetRecord &record = hostFacets[(size_t)facetId];
        record.nVertices = facet.nVertices;
        if (record.nVertices > 64)
            record.nVertices = 64;
        for (int k = 0; k < 4; ++k)
        {
            record.normalIn[k] = facet.in_normal.coordinates[k];
            record.normalOut[k] = facet.ex_normal.coordinates[k];
        }
        record.margin = gpu_trace_margin();
        for (int v = 0; v < record.nVertices; ++v)
            for (int k = 0; k < 4; ++k)
                record.vertices[v][k] = facet.arr[v].coordinates[k];
    }

    if (!trace_cuda_ok(cudaMemcpy(g_traceWorkspace.facets, hostFacets.data(),
                                  hostFacets.size() * sizeof(GpuTraceFacetRecord),
                                  cudaMemcpyHostToDevice), "copy trace facets"))
        return false;
    g_traceWorkspace.facetOwner = facets;
    g_traceWorkspace.copiedMaxFacetId = maxFacetId;
    return true;
}

void GpuTraceInvalidateFacetCache()
{
    g_traceWorkspace.facetOwner = nullptr;
    g_traceWorkspace.copiedMaxFacetId = -1;
}

__device__ inline float dot3_dev(const float *a, const float *b)
{
    return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
}

__device__ inline void add_bounds(float u, float v,
                                  float &minU, float &maxU,
                                  float &minV, float &maxV)
{
    minU = fminf(minU, u);
    maxU = fmaxf(maxU, u);
    minV = fminf(minV, v);
    maxV = fmaxf(maxV, v);
}

__global__ void trace_prefilter_kernel(const GpuTraceCandidate *candidates,
                                       unsigned char *results,
                                       int count)
{
    int idx = (int)(blockIdx.x * blockDim.x + threadIdx.x);
    if (idx >= count) return;

    const GpuTraceCandidate &c = candidates[idx];
    float dp0 = dot3_dev(c.dir, c.normal);
    if (fabsf(dp0) < 1.7453284e-5f)
    {
        results[idx] = 0;
        return;
    }

    float ax = fabsf(c.normal[0]);
    float ay = fabsf(c.normal[1]);
    float az = fabsf(c.normal[2]);
    int drop = (ax > ay && ax > az) ? 0 : ((ay > az) ? 1 : 2);

    float bMinU = 1e30f, bMaxU = -1e30f;
    float bMinV = 1e30f, bMaxV = -1e30f;
    float fMinU = 1e30f, fMaxU = -1e30f;
    float fMinV = 1e30f, fMaxV = -1e30f;

    for (int i = 0; i < c.beamN; ++i)
    {
        float p[3] = {c.beam[i][0], c.beam[i][1], c.beam[i][2]};
        float t = (dot3_dev(p, c.normal) + c.normal[3]) / dp0;
        p[0] -= c.dir[0] * t;
        p[1] -= c.dir[1] * t;
        p[2] -= c.dir[2] * t;

        float u = drop == 0 ? p[1] : p[0];
        float v = drop == 2 ? p[1] : p[2];
        add_bounds(u, v, bMinU, bMaxU, bMinV, bMaxV);
    }

    for (int i = 0; i < c.facetN; ++i)
    {
        float u = drop == 0 ? c.facet[i][1] : c.facet[i][0];
        float v = drop == 2 ? c.facet[i][1] : c.facet[i][2];
        add_bounds(u, v, fMinU, fMaxU, fMinV, fMaxV);
    }

    const float margin = 0.1f;
    bool overlap = !(bMaxU < fMinU - margin || fMaxU < bMinU - margin ||
                     bMaxV < fMinV - margin || fMaxV < bMinV - margin);
    results[idx] = overlap ? 1 : 0;
}

__device__ inline double dot3_dev_d(const double *a, const double *b)
{
    return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
}

__device__ inline void add_bounds_d(double u, double v,
                                    double &minU, double &maxU,
                                    double &minV, double &maxV)
{
    minU = fmin(minU, u);
    maxU = fmax(maxU, u);
    minV = fmin(minV, v);
    maxV = fmax(maxV, v);
}

__global__ void trace_prefilter_pair_kernel(const GpuTraceBeamRecord *beams,
                                            const GpuTraceFacetRecord *facets,
                                            const GpuTracePair *pairs,
                                            unsigned char *results,
                                            int count)
{
    int idx = (int)(blockIdx.x * blockDim.x + threadIdx.x);
    if (idx >= count) return;

    const GpuTracePair pair = pairs[idx];
    const GpuTraceBeamRecord &b = beams[pair.beam];
    const GpuTraceFacetRecord &f = facets[pair.facet];
    const double *normal = b.location == 0 ? f.normalIn : f.normalOut;
    double dp0 = dot3_dev_d(b.dir, normal);
    if (fabs(dp0) < 1.7453284e-5)
    {
        results[idx] = 0;
        return;
    }

    double ax = fabs(normal[0]);
    double ay = fabs(normal[1]);
    double az = fabs(normal[2]);
    int drop = (ax > ay && ax > az) ? 0 : ((ay > az) ? 1 : 2);

    double bMinU = 1e300, bMaxU = -1e300;
    double bMinV = 1e300, bMaxV = -1e300;
    double fMinU = 1e300, fMaxU = -1e300;
    double fMinV = 1e300, fMaxV = -1e300;

    for (int i = 0; i < b.nVertices; ++i)
    {
        double p[3] = {b.vertices[i][0], b.vertices[i][1], b.vertices[i][2]};
        double t = (dot3_dev_d(p, normal) + normal[3]) / dp0;
        p[0] -= b.dir[0] * t;
        p[1] -= b.dir[1] * t;
        p[2] -= b.dir[2] * t;

        double u = drop == 0 ? p[1] : p[0];
        double v = drop == 2 ? p[1] : p[2];
        add_bounds_d(u, v, bMinU, bMaxU, bMinV, bMaxV);
    }

    for (int i = 0; i < f.nVertices; ++i)
    {
        double u = drop == 0 ? f.vertices[i][1] : f.vertices[i][0];
        double v = drop == 2 ? f.vertices[i][1] : f.vertices[i][2];
        add_bounds_d(u, v, fMinU, fMaxU, fMinV, fMaxV);
    }

    const double margin = f.margin;
    bool overlap = !(bMaxU < fMinU - margin || fMaxU < bMinU - margin ||
                     bMaxV < fMinV - margin || fMaxV < bMinV - margin);
    results[idx] = overlap ? 1 : 0;
}

static int gpu_trace_threshold()
{
    const char *value = std::getenv("MBS_GPU_TRACE_MIN_CANDIDATES");
    if (!value || !*value)
        return 8192;
    char *end = nullptr;
    long parsed = std::strtol(value, &end, 10);
    if (!end || *end != '\0' || parsed < 1)
        return 8192;
    return (int)parsed;
}

bool GpuTracePrefilterBeamFacets(const Beam &beam,
                                 const Facet *facets,
                                 const IntArray &facetIds,
                                 std::vector<unsigned char> &mayIntersect)
{
    GpuTraceBeamFacets item;
    item.beam = &beam;
    item.facetIds = &facetIds;
    item.mayIntersect = &mayIntersect;
    std::vector<GpuTraceBeamFacets> items(1, item);
    return GpuTracePrefilterBeamFacetBatch(facets, items);
}

bool GpuTracePrefilterBeamFacetBatch(const Facet *facets,
                                     const std::vector<GpuTraceBeamFacets> &items)
{
    size_t total = 0;
    for (size_t itemIdx = 0; itemIdx < items.size(); ++itemIdx)
    {
        GpuTraceBeamFacets item = items[itemIdx];
        item.mayIntersect->assign(item.facetIds->size, 1);
        if (item.beam->nVertices <= 0 || item.beam->nVertices > 64)
            continue;
        total += item.facetIds->size;
    }

    if (total < (size_t)gpu_trace_threshold())
        return false;
    if (!ensure_trace_capacity(total))
        return false;

    if (!ensure_trace_beam_capacity(items.size()))
        return false;

    int maxFacetId = -1;
    for (size_t itemIdx = 0; itemIdx < items.size(); ++itemIdx)
    {
        const IntArray &facetIds = *items[itemIdx].facetIds;
        for (size_t facetIdx = 0; facetIdx < facetIds.size; ++facetIdx)
            maxFacetId = std::max(maxFacetId, facetIds.arr[facetIdx]);
    }
    if (maxFacetId < 0)
        return false;
    if (!upload_trace_facets_if_needed(facets, maxFacetId))
        return false;

    std::vector<GpuTraceBeamRecord> hostBeams(items.size());
    std::vector<GpuTracePair> hostPairs(total);
    std::vector<GpuTraceBeamFacets> mappedItems;
    std::vector<size_t> mappedFacetOffsets;
    mappedItems.reserve(total);
    mappedFacetOffsets.reserve(total);

    for (size_t itemIdx = 0; itemIdx < items.size(); ++itemIdx)
    {
        const Beam &beam = *items[itemIdx].beam;
        GpuTraceBeamRecord &record = hostBeams[itemIdx];
        record.nVertices = beam.nVertices;
        record.location = beam.location == Location::In ? 0 : 1;
        record.dir[0] = beam.direction.coordinates[0];
        record.dir[1] = beam.direction.coordinates[1];
        record.dir[2] = beam.direction.coordinates[2];
        record.dir[3] = beam.direction.coordinates[3];
        for (int v = 0; v < beam.nVertices; ++v)
            for (int k = 0; k < 4; ++k)
                record.vertices[v][k] = beam.arr[v].coordinates[k];
    }

    size_t out = 0;
    for (size_t itemIdx = 0; itemIdx < items.size(); ++itemIdx)
    {
        GpuTraceBeamFacets item = items[itemIdx];
        const Beam &beam = *item.beam;
        if (beam.nVertices <= 0 || beam.nVertices > 64)
            continue;

        for (size_t facetIdx = 0; facetIdx < item.facetIds->size; ++facetIdx)
        {
            hostPairs[out].beam = (int)itemIdx;
            hostPairs[out].facet = item.facetIds->arr[facetIdx];
            mappedItems.push_back(item);
            mappedFacetOffsets.push_back(facetIdx);
            ++out;
        }
    }

    if (out == 0)
        return false;

    if (!trace_cuda_ok(cudaMemcpy(g_traceWorkspace.beams, hostBeams.data(),
                                  hostBeams.size() * sizeof(GpuTraceBeamRecord),
                                  cudaMemcpyHostToDevice), "copy trace beams"))
        return false;
    if (!trace_cuda_ok(cudaMemcpy(g_traceWorkspace.pairs, hostPairs.data(),
                                  out * sizeof(GpuTracePair),
                                  cudaMemcpyHostToDevice), "copy trace pairs"))
        return false;
    int block = 128;
    int grid = ((int)out + block - 1) / block;
    trace_prefilter_pair_kernel<<<grid, block>>>(g_traceWorkspace.beams,
                                                 g_traceWorkspace.facets,
                                                 g_traceWorkspace.pairs,
                                                 g_traceWorkspace.results,
                                                 (int)out);
    if (!trace_cuda_ok(cudaGetLastError(), "trace_prefilter_batch_kernel"))
        return false;

    std::vector<unsigned char> results(out);
    if (!trace_cuda_ok(cudaMemcpy(results.data(), g_traceWorkspace.results,
                                  out * sizeof(unsigned char),
                                  cudaMemcpyDeviceToHost), "copy batch results"))
        return false;

    for (size_t idx = 0; idx < out; ++idx)
        (*mappedItems[idx].mayIntersect)[mappedFacetOffsets[idx]] = results[idx];
    return true;
}

#else

bool GpuTracePrefilterBeamFacets(const Beam &/*beam*/,
                                 const Facet */*facets*/,
                                 const IntArray &facetIds,
                                 std::vector<unsigned char> &mayIntersect)
{
    mayIntersect.assign(facetIds.size, 1);
    return false;
}

bool GpuTracePrefilterBeamFacetBatch(const Facet */*facets*/,
                                     const std::vector<GpuTraceBeamFacets> &items)
{
    for (size_t i = 0; i < items.size(); ++i)
        items[i].mayIntersect->assign(items[i].facetIds->size, 1);
    return false;
}

void GpuTraceInvalidateFacetCache()
{
}

#endif
