#include "GpuTraceSupport.h"

#ifdef USE_CUDA

#include <cuda_runtime.h>

#include <algorithm>
#include <cfloat>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <cstdio>
#include <vector>

constexpr int GPU_TRACE_COMPACT_VERTICES = 8;
constexpr int GPU_TRACE_TRIANGLE_VERTICES = 3;
constexpr int GPU_TRACE_TRIANGLE_INTERSECTION_VERTICES = 6;
constexpr int GPU_TRACE_SORT_CANDIDATES = 24;
constexpr int GPU_TRACE_LARGE_SORT_CANDIDATES = 256;
constexpr int GPU_TRACE_LARGE_SORT_WORDS =
    GPU_TRACE_LARGE_SORT_CANDIDATES / 32;
constexpr uint32_t GPU_TRACE_SORTED_MARKER = 0x80000000u;
constexpr uint32_t GPU_TRACE_SKIP_MARKER = 0x40000000u;
constexpr uint32_t GPU_TRACE_FIRST_SKIP_MARKER = 0x00200000u;
constexpr uint32_t GPU_TRACE_COUNT_SHIFT = 22u;
constexpr uint32_t GPU_TRACE_COUNT_MASK = 0x7fc00000u;
constexpr uint32_t GPU_TRACE_FACET_ID_MASK = 0x001fffffu;

struct GpuTraceBeamRecord
{
    double dir[3];
    int vertexOffset;
    int nVertices;
    int location;
};

struct GpuTraceFacetRecord
{
    double vertices[GPU_TRACE_COMPACT_VERTICES][3];
    double surfaceNormal[3];
    double normalIn[4];
    double normalOut[4];
    double bounds[3][4];
    double margin;
    int nVertices;
};

struct GpuTraceItemRecord
{
    int offset;
    int count;
};

struct GpuTraceInputLayout
{
    size_t beams;
    size_t items;
    size_t itemIndices;
    size_t vertices;
    size_t candidates;
    size_t bytes;
};

struct GpuTraceWorkspace
{
    unsigned char *input = nullptr;
    GpuTraceFacetRecord *facets = nullptr;
    unsigned char *hostInput = nullptr;
    size_t inputCap = 0;
    size_t facetCap = 0;
    size_t hostInputCap = 0;
    const Facet *facetOwner = nullptr;
    int copiedMaxFacetId = -1;
    bool copiedFacetsAreTriangles = false;
    cudaEvent_t timingStart = nullptr;
    cudaEvent_t timingSmallDone = nullptr;
    cudaEvent_t timingLargeDone = nullptr;
    bool timingEventsReady = false;
    cudaStream_t stream = nullptr;
    bool streamReady = false;
};

static thread_local GpuTraceWorkspace g_traceWorkspace;
static thread_local bool g_traceLastFailureWasOutOfMemory = false;

static bool trace_cuda_ok(cudaError_t err, const char *where);

static double gpu_trace_margin_override()
{
    const char *value = std::getenv("MBS_GPU_TRACE_MARGIN");
    if (!value || !*value)
        return -1.0;
    char *end = nullptr;
    double parsed = std::strtod(value, &end);
    if (!end || *end != '\0' || parsed < 0.0)
        return -1.0;
    return parsed;
}

static bool gpu_trace_cache_facets()
{
    const char *value = std::getenv("MBS_GPU_TRACE_CACHE_FACETS");
    return !(value && *value == '0');
}

static bool gpu_trace_prefilter_first()
{
    const char *value = std::getenv("MBS_GPU_TRACE_PREFILTER_FIRST");
    return value && value[0] == '1' && value[1] == '\0';
}

static bool gpu_trace_full_sort()
{
    const char *value = std::getenv("MBS_GPU_TRACE_FULL_SORT");
    return !(value && value[0] == '0' && value[1] == '\0');
}

static bool gpu_trace_large_skip_flags()
{
    const char *value = std::getenv("MBS_GPU_TRACE_LARGE_SKIP_FLAGS");
    return !(value && value[0] == '0' && value[1] == '\0');
}

static bool gpu_trace_timing_enabled()
{
    const char *value = std::getenv("MBS_GPU_TRACE_TIMING");
    return value && value[0] == '1' && value[1] == '\0';
}

static bool gpu_trace_nonblocking_stream_enabled()
{
    const char *value = std::getenv("MBS_GPU_TRACE_NONBLOCKING_STREAM");
    return value && value[0] == '1' && value[1] == '\0';
}

static bool ensure_trace_stream()
{
    if (!gpu_trace_nonblocking_stream_enabled())
        return true;
    if (g_traceWorkspace.streamReady)
        return true;
    if (!trace_cuda_ok(
            cudaStreamCreateWithFlags(&g_traceWorkspace.stream,
                                      cudaStreamNonBlocking),
            "create nonblocking trace stream"))
        return false;
    g_traceWorkspace.streamReady = true;
    return true;
}

static cudaStream_t trace_stream()
{
    return g_traceWorkspace.streamReady ? g_traceWorkspace.stream : nullptr;
}

static bool ensure_trace_timing_events()
{
    if (g_traceWorkspace.timingEventsReady)
        return true;
    if (cudaEventCreate(&g_traceWorkspace.timingStart) != cudaSuccess
        || cudaEventCreate(&g_traceWorkspace.timingSmallDone) != cudaSuccess
        || cudaEventCreate(&g_traceWorkspace.timingLargeDone) != cudaSuccess)
        return false;
    g_traceWorkspace.timingEventsReady = true;
    return true;
}

static bool trace_cuda_ok(cudaError_t err, const char *where)
{
    if (err == cudaSuccess)
        return true;
    if (err == cudaErrorMemoryAllocation)
        g_traceLastFailureWasOutOfMemory = true;
    std::fprintf(stderr, "CUDA trace prefilter error at %s: %s\n",
                 where, cudaGetErrorString(err));
    return false;
}

static size_t trace_growth_capacity(size_t current, size_t required)
{
    size_t capacity = std::max<size_t>(current, 256);
    while (capacity < required)
    {
        const size_t grown = capacity + capacity / 2;
        if (grown <= capacity)
            return required;
        capacity = grown;
    }
    return capacity;
}

static size_t trace_align_offset(size_t offset, size_t alignment)
{
    return (offset + alignment - 1) & ~(alignment - 1);
}

static GpuTraceInputLayout trace_input_layout(size_t beamCount,
                                              size_t vertexCount,
                                              size_t candidateCount)
{
    GpuTraceInputLayout layout = {};
    layout.beams = 0;
    layout.items = trace_align_offset(
        beamCount * sizeof(GpuTraceBeamRecord),
        alignof(GpuTraceItemRecord));
    layout.itemIndices = trace_align_offset(
        layout.items + beamCount * sizeof(GpuTraceItemRecord),
        alignof(int));
    layout.vertices = trace_align_offset(
        layout.itemIndices + beamCount * sizeof(int),
        alignof(double));
    layout.candidates = trace_align_offset(
        layout.vertices + vertexCount * 3 * sizeof(double),
        alignof(uint32_t));
    layout.bytes = layout.candidates
                 + candidateCount * sizeof(uint32_t);
    return layout;
}

static bool ensure_trace_input_capacity(size_t bytes)
{
    if (g_traceWorkspace.inputCap >= bytes)
        return true;
    if (bytes == 0)
        return true;
    const size_t capacity = trace_growth_capacity(
        g_traceWorkspace.inputCap, bytes);
    unsigned char *input = nullptr;
    const cudaError_t status = cudaMalloc(&input, capacity);
    if (status != cudaSuccess)
    {
        if (status == cudaErrorMemoryAllocation)
            g_traceLastFailureWasOutOfMemory = true;
        return false;
    }
    cudaFree(g_traceWorkspace.input);
    g_traceWorkspace.input = input;
    g_traceWorkspace.inputCap = capacity;
    return true;
}

static bool ensure_trace_host_input_capacity(size_t bytes)
{
    if (g_traceWorkspace.hostInputCap >= bytes)
        return true;
    if (bytes == 0)
        return true;
    const size_t capacity = trace_growth_capacity(
        g_traceWorkspace.hostInputCap, bytes);
    unsigned char *input = nullptr;
    if (cudaMallocHost(reinterpret_cast<void **>(&input), capacity)
        != cudaSuccess)
        return false;
    if (g_traceWorkspace.hostInput != nullptr)
        cudaFreeHost(g_traceWorkspace.hostInput);
    g_traceWorkspace.hostInput = input;
    g_traceWorkspace.hostInputCap = capacity;
    return true;
}

static bool ensure_trace_facet_capacity(size_t count)
{
    if (g_traceWorkspace.facetCap >= count)
        return true;
    if (count == 0)
        return true;
    const size_t capacity = trace_growth_capacity(g_traceWorkspace.facetCap, count);
    GpuTraceFacetRecord *facets = nullptr;
    const cudaError_t status = cudaMalloc(
        &facets, capacity * sizeof(GpuTraceFacetRecord));
    if (status != cudaSuccess)
    {
        if (status == cudaErrorMemoryAllocation)
            g_traceLastFailureWasOutOfMemory = true;
        return false;
    }
    cudaFree(g_traceWorkspace.facets);
    g_traceWorkspace.facets = facets;
    g_traceWorkspace.facetCap = capacity;
    g_traceWorkspace.facetOwner = nullptr;
    g_traceWorkspace.copiedMaxFacetId = -1;
    g_traceWorkspace.copiedFacetsAreTriangles = false;
    return true;
}

static bool upload_trace_facets_if_needed(const Facet *facets, int maxFacetId,
                                          double geometryScale)
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
    bool allFacetsAreTriangles = true;
    const double marginOverride = gpu_trace_margin_override();
    for (int facetId = 0; facetId <= maxFacetId; ++facetId)
    {
        const Facet &facet = facets[facetId];
        allFacetsAreTriangles = allFacetsAreTriangles
                             && facet.nVertices
                                == GPU_TRACE_TRIANGLE_VERTICES;
        GpuTraceFacetRecord &record = hostFacets[(size_t)facetId];
        if (facet.nVertices < 3 || facet.nVertices > 64)
            return false;
        record.nVertices = facet.nVertices;
        for (int k = 0; k < 4; ++k)
        {
            record.normalIn[k] = facet.in_normal.coordinates[k];
            record.normalOut[k] = facet.ex_normal.coordinates[k];
        }
        const int copiedVertices = std::min(
            facet.nVertices, GPU_TRACE_COMPACT_VERTICES);
        for (int vertex = 0; vertex < copiedVertices; ++vertex)
            for (int coordinate = 0; coordinate < 3; ++coordinate)
                record.vertices[vertex][coordinate] =
                    facet.arr[vertex].coordinates[coordinate];
        const Point3f surfaceNormal = facet.Normal();
        for (int coordinate = 0; coordinate < 3; ++coordinate)
            record.surfaceNormal[coordinate] =
                surfaceNormal.coordinates[coordinate];
        for (int drop = 0; drop < 3; ++drop)
        {
            double minU = DBL_MAX, maxU = -DBL_MAX;
            double minV = DBL_MAX, maxV = -DBL_MAX;
            for (int vertex = 0; vertex < facet.nVertices; ++vertex)
            {
                const Point3f &point = facet.arr[vertex];
                const double u = drop == 0 ? point.cy : point.cx;
                const double v = drop == 2 ? point.cy : point.cz;
                minU = std::min(minU, u);
                maxU = std::max(maxU, u);
                minV = std::min(minV, v);
                maxV = std::max(maxV, v);
            }
            record.bounds[drop][0] = minU;
            record.bounds[drop][1] = maxU;
            record.bounds[drop][2] = minV;
            record.bounds[drop][3] = maxV;
        }
    }
    const double margin = marginOverride >= 0.0
        ? marginOverride
        : 0x1p-40
            * std::max(std::fabs(geometryScale), DBL_MIN);
    for (GpuTraceFacetRecord &record : hostFacets)
        record.margin = margin;

    if (!trace_cuda_ok(cudaMemcpy(g_traceWorkspace.facets, hostFacets.data(),
                                  hostFacets.size() * sizeof(GpuTraceFacetRecord),
                                  cudaMemcpyHostToDevice), "copy trace facets"))
        return false;
    g_traceWorkspace.facetOwner = facets;
    g_traceWorkspace.copiedMaxFacetId = maxFacetId;
    g_traceWorkspace.copiedFacetsAreTriangles = allFacetsAreTriangles;
    return true;
}

void GpuTraceInvalidateFacetCache()
{
    g_traceWorkspace.facetOwner = nullptr;
    g_traceWorkspace.copiedMaxFacetId = -1;
    g_traceWorkspace.copiedFacetsAreTriangles = false;
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

struct GpuTraceFacetDepth
{
    int facetId;
    double depth;
    double maximumDepth;
};

union GpuTraceDepthScratch
{
    double storage[2];
    float projectedBounds[4];
};

__device__ inline void trace_store_projected_bounds(
    GpuTraceFacetDepth &entry,
    float minU, float maxU, float minV, float maxV)
{
    GpuTraceDepthScratch scratch;
    scratch.projectedBounds[0] = minU;
    scratch.projectedBounds[1] = maxU;
    scratch.projectedBounds[2] = minV;
    scratch.projectedBounds[3] = maxV;
    entry.depth = scratch.storage[0];
    entry.maximumDepth = scratch.storage[1];
}

__device__ inline GpuTraceDepthScratch trace_load_projected_bounds(
    const GpuTraceFacetDepth &entry)
{
    GpuTraceDepthScratch scratch;
    scratch.storage[0] = entry.depth;
    scratch.storage[1] = entry.maximumDepth;
    return scratch;
}

__device__ inline bool trace_depth_less(const GpuTraceFacetDepth &left,
                                        const GpuTraceFacetDepth &right)
{
    return left.depth < right.depth
        || (left.depth == right.depth && left.facetId < right.facetId);
}

struct GpuTracePoint2
{
    double x;
    double y;
};

__device__ inline double projected_cross_dev(const GpuTracePoint2 &a,
                                             const GpuTracePoint2 &b,
                                             const GpuTracePoint2 &c)
{
    return (b.x - a.x)*(c.y - a.y) - (b.y - a.y)*(c.x - a.x);
}

__device__ double projected_area_dev(const GpuTracePoint2 *polygon, int size)
{
    if (size <= 0)
        return 0.0;
    double twiceArea = 0.0;
    for (int index = 0; index + 1 < size; ++index)
    {
        const GpuTracePoint2 &a = polygon[index];
        const GpuTracePoint2 &b = polygon[index + 1];
        twiceArea += a.x*b.y - a.y*b.x;
    }
    const GpuTracePoint2 &last = polygon[size - 1];
    const GpuTracePoint2 &first = polygon[0];
    twiceArea += last.x*first.y - last.y*first.x;
    return 0.5*twiceArea;
}

__device__ bool projected_overlap_point_dev(
    const GpuTraceFacetRecord &first,
    const GpuTraceFacetRecord &second,
    const double *axisU,
    const double *axisV,
    double areaTolerance,
    GpuTracePoint2 &point)
{
    constexpr int capacity = 2*GPU_TRACE_COMPACT_VERTICES;
    GpuTracePoint2 firstProjected[GPU_TRACE_COMPACT_VERTICES];
    GpuTracePoint2 secondProjected[GPU_TRACE_COMPACT_VERTICES];
    for (int vertex = 0; vertex < first.nVertices; ++vertex)
    {
        firstProjected[vertex].x = dot3_dev_d(first.vertices[vertex], axisU);
        firstProjected[vertex].y = dot3_dev_d(first.vertices[vertex], axisV);
    }
    for (int vertex = 0; vertex < second.nVertices; ++vertex)
    {
        secondProjected[vertex].x = dot3_dev_d(second.vertices[vertex], axisU);
        secondProjected[vertex].y = dot3_dev_d(second.vertices[vertex], axisV);
    }

    const double secondOrientation = projected_area_dev(
        secondProjected, second.nVertices) >= 0.0 ? 1.0 : -1.0;
    GpuTracePoint2 buffers[2][capacity];
    for (int vertex = 0; vertex < first.nVertices; ++vertex)
        buffers[0][vertex] = firstProjected[vertex];
    GpuTracePoint2 *result = buffers[0];
    int resultSize = first.nVertices;
    const double orientation = secondOrientation;
    for (int edge = 0; edge < second.nVertices && resultSize != 0; ++edge)
    {
        const GpuTracePoint2 &a = secondProjected[edge];
        const GpuTracePoint2 &b =
            secondProjected[(edge + 1) % second.nVertices];
        GpuTracePoint2 *input = result;
        const int inputSize = resultSize;
        result = result == buffers[0] ? buffers[1] : buffers[0];
        resultSize = 0;
        GpuTracePoint2 previous = input[inputSize - 1];
        double previousSide = orientation*projected_cross_dev(a, b, previous);
        bool previousInside = previousSide >= -areaTolerance;
        for (int index = 0; index < inputSize; ++index)
        {
            const GpuTracePoint2 current = input[index];
            const double currentSide =
                orientation*projected_cross_dev(a, b, current);
            const bool currentInside = currentSide >= -areaTolerance;
            if (currentInside != previousInside)
            {
                const double denominator = previousSide - currentSide;
                if (fabs(denominator) > DBL_MIN)
                {
                    if (resultSize >= capacity)
                        return false;
                    const double t = previousSide/denominator;
                    result[resultSize++] = {
                        previous.x + t*(current.x - previous.x),
                        previous.y + t*(current.y - previous.y)};
                }
            }
            if (currentInside)
            {
                if (resultSize >= capacity)
                    return false;
                result[resultSize++] = current;
            }
            previous = current;
            previousSide = currentSide;
            previousInside = currentInside;
        }
    }
    if (resultSize < 3
        || fabs(projected_area_dev(result, resultSize)) <= areaTolerance)
        return false;

    point = {0.0, 0.0};
    for (int vertex = 0; vertex < resultSize; ++vertex)
    {
        point.x += result[vertex].x;
        point.y += result[vertex].y;
    }
    point.x /= resultSize;
    point.y /= resultSize;
    return true;
}

__device__ __forceinline__ int projected_clip_edge_dev(
    const GpuTracePoint2 *input,
    int inputSize,
    const GpuTracePoint2 &a,
    const GpuTracePoint2 &b,
    double orientation,
    double areaTolerance,
    GpuTracePoint2 *output,
    int outputCapacity)
{
    if (inputSize <= 0)
        return 0;
    int outputSize = 0;
    GpuTracePoint2 previous = input[inputSize - 1];
    double previousSide =
        orientation*projected_cross_dev(a, b, previous);
    bool previousInside = previousSide >= -areaTolerance;
    for (int index = 0; index < inputSize; ++index)
    {
        const GpuTracePoint2 current = input[index];
        const double currentSide =
            orientation*projected_cross_dev(a, b, current);
        const bool currentInside = currentSide >= -areaTolerance;
        if (currentInside != previousInside)
        {
            const double denominator = previousSide - currentSide;
            if (fabs(denominator) > DBL_MIN)
            {
                if (outputSize >= outputCapacity)
                    return -1;
                const double t = previousSide/denominator;
                output[outputSize++] = {
                    previous.x + t*(current.x - previous.x),
                    previous.y + t*(current.y - previous.y)};
            }
        }
        if (currentInside)
        {
            if (outputSize >= outputCapacity)
                return -1;
            output[outputSize++] = current;
        }
        previous = current;
        previousSide = currentSide;
        previousInside = currentInside;
    }
    return outputSize;
}

struct GpuTracePolygonAccumulator
{
    GpuTracePoint2 sum;
    GpuTracePoint2 first;
    GpuTracePoint2 previous;
    double twiceArea;
    int count;
};

__device__ __forceinline__ void projected_accumulate_point_dev(
    const GpuTracePoint2 &value,
    GpuTracePolygonAccumulator &accumulator)
{
    if (accumulator.count == 0)
    {
        accumulator.first = value;
    }
    else
    {
        accumulator.twiceArea +=
            accumulator.previous.x*value.y
            - accumulator.previous.y*value.x;
    }
    accumulator.previous = value;
    accumulator.sum.x += value.x;
    accumulator.sum.y += value.y;
    ++accumulator.count;
}

__device__ __forceinline__ bool projected_triangle_overlap_clip_dev(
    const GpuTracePoint2 *firstProjected,
    const GpuTracePoint2 *secondProjected,
    double secondOrientation,
    double areaTolerance,
    GpuTracePoint2 &point)
{
    GpuTracePoint2 firstBuffer[4];
    GpuTracePoint2 secondBuffer[5];
    int firstSize = projected_clip_edge_dev(
        firstProjected, GPU_TRACE_TRIANGLE_VERTICES,
        secondProjected[0], secondProjected[1],
        secondOrientation, areaTolerance, firstBuffer, 4);
    if (firstSize <= 0)
        return false;
    int secondSize = projected_clip_edge_dev(
        firstBuffer, firstSize,
        secondProjected[1], secondProjected[2],
        secondOrientation, areaTolerance, secondBuffer, 5);
    if (secondSize <= 0)
        return false;

    const GpuTracePoint2 &a = secondProjected[2];
    const GpuTracePoint2 &b = secondProjected[0];
    GpuTracePolygonAccumulator accumulator = {
        {0.0, 0.0}, {0.0, 0.0}, {0.0, 0.0}, 0.0, 0};
    GpuTracePoint2 previous = secondBuffer[secondSize - 1];
    double previousSide =
        secondOrientation*projected_cross_dev(a, b, previous);
    bool previousInside = previousSide >= -areaTolerance;
    for (int index = 0; index < secondSize; ++index)
    {
        const GpuTracePoint2 current = secondBuffer[index];
        const double currentSide =
            secondOrientation*projected_cross_dev(a, b, current);
        const bool currentInside = currentSide >= -areaTolerance;
        if (currentInside != previousInside)
        {
            const double denominator = previousSide - currentSide;
            if (fabs(denominator) > DBL_MIN)
            {
                const double t = previousSide/denominator;
                const GpuTracePoint2 intersection = {
                    previous.x + t*(current.x - previous.x),
                    previous.y + t*(current.y - previous.y)};
                projected_accumulate_point_dev(intersection, accumulator);
            }
        }
        if (currentInside)
            projected_accumulate_point_dev(current, accumulator);
        previous = current;
        previousSide = currentSide;
        previousInside = currentInside;
    }
    if (accumulator.count < 3)
        return false;
    accumulator.twiceArea +=
        accumulator.previous.x*accumulator.first.y
        - accumulator.previous.y*accumulator.first.x;
    if (fabs(0.5*accumulator.twiceArea) <= areaTolerance)
        return false;
    point.x = accumulator.sum.x / accumulator.count;
    point.y = accumulator.sum.y / accumulator.count;
    return true;
}

template <bool TriangleFacets>
__device__ __forceinline__ bool projected_overlap_point_cached_first_dev(
    const GpuTracePoint2 *firstProjected,
    int firstVertices,
    const GpuTraceFacetRecord &second,
    const GpuTracePoint2 *cachedSecondProjected,
    const double *axisU,
    const double *axisV,
    double firstOrientation,
    double secondOrientation,
    double areaTolerance,
    GpuTracePoint2 &point)
{
    constexpr int projectedCapacity = TriangleFacets
        ? GPU_TRACE_TRIANGLE_VERTICES : GPU_TRACE_COMPACT_VERTICES;
    constexpr int capacity = TriangleFacets
        ? GPU_TRACE_TRIANGLE_INTERSECTION_VERTICES
        : 2*GPU_TRACE_COMPACT_VERTICES;
    GpuTracePoint2 localSecondProjected[projectedCapacity];
    const GpuTracePoint2 *secondProjected = cachedSecondProjected;
    const int clippedFirstVertices = TriangleFacets
        ? GPU_TRACE_TRIANGLE_VERTICES : firstVertices;
    const int clippedSecondVertices = TriangleFacets
        ? GPU_TRACE_TRIANGLE_VERTICES : second.nVertices;
    if (!TriangleFacets || cachedSecondProjected == nullptr)
    {
        for (int vertex = 0; vertex < clippedSecondVertices; ++vertex)
        {
            localSecondProjected[vertex].x =
                dot3_dev_d(second.vertices[vertex], axisU);
            localSecondProjected[vertex].y =
                dot3_dev_d(second.vertices[vertex], axisV);
        }
        secondProjected = localSecondProjected;
    }
    if (TriangleFacets)
    {
        return projected_triangle_overlap_clip_dev(
            firstProjected, secondProjected,
            secondOrientation, areaTolerance, point);
    }
    GpuTracePoint2 buffers[2][capacity];
    for (int vertex = 0; vertex < clippedFirstVertices; ++vertex)
        buffers[0][vertex] = firstProjected[vertex];
    GpuTracePoint2 *result = buffers[0];
    int resultSize = clippedFirstVertices;
    for (int edge = 0; edge < clippedSecondVertices && resultSize != 0;
         ++edge)
    {
        const GpuTracePoint2 &a = secondProjected[edge];
        const int nextEdge = edge + 1 == clippedSecondVertices ? 0 : edge + 1;
        const GpuTracePoint2 &b = secondProjected[nextEdge];
        GpuTracePoint2 *input = result;
        const int inputSize = resultSize;
        result = result == buffers[0] ? buffers[1] : buffers[0];
        resultSize = 0;
        GpuTracePoint2 previous = input[inputSize - 1];
        double previousSide =
            secondOrientation*projected_cross_dev(a, b, previous);
        bool previousInside = previousSide >= -areaTolerance;
        for (int index = 0; index < inputSize; ++index)
        {
            const GpuTracePoint2 current = input[index];
            const double currentSide =
                secondOrientation*projected_cross_dev(a, b, current);
            const bool currentInside = currentSide >= -areaTolerance;
            if (currentInside != previousInside)
            {
                const double denominator = previousSide - currentSide;
                if (fabs(denominator) > DBL_MIN)
                {
                    if (resultSize >= capacity)
                        return false;
                    const double t = previousSide/denominator;
                    result[resultSize++] = {
                        previous.x + t*(current.x - previous.x),
                        previous.y + t*(current.y - previous.y)};
                }
            }
            if (currentInside)
            {
                if (resultSize >= capacity)
                    return false;
                result[resultSize++] = current;
            }
            previous = current;
            previousSide = currentSide;
            previousInside = currentInside;
        }
    }
    if (resultSize < 3
        || fabs(projected_area_dev(result, resultSize)) <= areaTolerance)
        return false;

    point = {0.0, 0.0};
    for (int vertex = 0; vertex < resultSize; ++vertex)
    {
        point.x += result[vertex].x;
        point.y += result[vertex].y;
    }
    point.x /= resultSize;
    point.y /= resultSize;
    return true;
}

__device__ bool normalize3_dev(double *vector)
{
    const double length = sqrt(dot3_dev_d(vector, vector));
    if (!(length > DBL_MIN))
        return false;
    vector[0] /= length;
    vector[1] /= length;
    vector[2] /= length;
    return true;
}

__device__ inline void cross3_dev(const double *a, const double *b,
                                  double *result)
{
    result[0] = a[1]*b[2] - a[2]*b[1];
    result[1] = a[2]*b[0] - a[0]*b[2];
    result[2] = a[0]*b[1] - a[1]*b[0];
}

__device__ bool trace_sort_facets_dev(const GpuTraceBeamRecord &beam,
                                      const GpuTraceFacetRecord *facets,
                                      uint32_t *candidateIds,
                                      int count,
                                      double geometryScale)
{
    if (count < 2)
        return true;
    if (count > GPU_TRACE_SORT_CANDIDATES
        || beam.nVertices <= 0
        || beam.nVertices > GPU_TRACE_COMPACT_VERTICES)
        return false;

    GpuTraceFacetDepth ordered[GPU_TRACE_SORT_CANDIDATES];
    for (int index = 0; index < count; ++index)
    {
        const int facetId = static_cast<int>(candidateIds[index]);
        const GpuTraceFacetRecord &facet = facets[facetId];
        if (facet.nVertices < 3
            || facet.nVertices > GPU_TRACE_COMPACT_VERTICES)
            return false;
        double minimum = dot3_dev_d(facet.vertices[0], beam.dir);
        double maximum = minimum;
        for (int vertex = 1; vertex < facet.nVertices; ++vertex)
        {
            const double depth = dot3_dev_d(facet.vertices[vertex], beam.dir);
            minimum = fmin(minimum, depth);
            maximum = fmax(maximum, depth);
        }
        ordered[index] = {facetId, minimum, maximum};
    }

    for (int index = 1; index < count; ++index)
    {
        const GpuTraceFacetDepth value = ordered[index];
        int position = index;
        while (position > 0
               && (value.depth < ordered[position - 1].depth
                   || (value.depth == ordered[position - 1].depth
                       && value.facetId < ordered[position - 1].facetId)))
        {
            ordered[position] = ordered[position - 1];
            --position;
        }
        ordered[position] = value;
    }

    const double scale = fmax(fabs(geometryScale), DBL_MIN);
    const double depthTolerance = 0x1p-40*scale;
    for (int begin = 0; begin < count;)
    {
        int end = begin + 1;
        while (end < count
               && ordered[end].depth - ordered[begin].depth <= depthTolerance)
            ++end;
        for (int index = begin + 1; index < end; ++index)
        {
            const GpuTraceFacetDepth value = ordered[index];
            int position = index;
            while (position > begin
                   && value.facetId < ordered[position - 1].facetId)
            {
                ordered[position] = ordered[position - 1];
                --position;
            }
            ordered[position] = value;
        }
        begin = end;
    }

    double direction[3] = {beam.dir[0], beam.dir[1], beam.dir[2]};
    if (!normalize3_dev(direction))
    {
        for (int index = 0; index < count; ++index)
            candidateIds[index] =
                static_cast<uint32_t>(ordered[index].facetId);
        return true;
    }
    double reference[3];
    const double ax = fabs(direction[0]);
    const double ay = fabs(direction[1]);
    const double az = fabs(direction[2]);
    if (ax <= ay && ax <= az)
    {
        reference[0] = 1.0; reference[1] = 0.0; reference[2] = 0.0;
    }
    else if (ay <= az)
    {
        reference[0] = 0.0; reference[1] = 1.0; reference[2] = 0.0;
    }
    else
    {
        reference[0] = 0.0; reference[1] = 0.0; reference[2] = 1.0;
    }
    double axisU[3], axisV[3];
    cross3_dev(direction, reference, axisU);
    if (!normalize3_dev(axisU))
        return false;
    cross3_dev(direction, axisU, axisV);
    if (!normalize3_dev(axisV))
        return false;

    unsigned int precedes[GPU_TRACE_SORT_CANDIDATES] = {};
    int indegree[GPU_TRACE_SORT_CANDIDATES] = {};
    const double areaTolerance = 0x1p-40*scale*scale;
    for (int first = 0; first < count; ++first)
    {
        for (int second = first + 1; second < count; ++second)
        {
            if (ordered[second].depth
                > ordered[first].maximumDepth + depthTolerance)
                break;
            const GpuTraceFacetRecord &firstFacet =
                facets[ordered[first].facetId];
            const GpuTraceFacetRecord &secondFacet =
                facets[ordered[second].facetId];
            GpuTracePoint2 overlapPoint;
            if (!projected_overlap_point_dev(
                    firstFacet, secondFacet, axisU, axisV,
                    areaTolerance, overlapPoint))
                continue;

            double rayOrigin[3] = {
                axisU[0]*overlapPoint.x + axisV[0]*overlapPoint.y,
                axisU[1]*overlapPoint.x + axisV[1]*overlapPoint.y,
                axisU[2]*overlapPoint.x + axisV[2]*overlapPoint.y};
            const double *firstNormal = firstFacet.surfaceNormal;
            const double *secondNormal = secondFacet.surfaceNormal;
            const double firstDenominator = dot3_dev_d(direction, firstNormal);
            const double secondDenominator = dot3_dev_d(direction, secondNormal);
            if (fabs(firstDenominator) <= 0x1p-40
                || fabs(secondDenominator) <= 0x1p-40)
                continue;
            double firstDelta[3], secondDelta[3];
            for (int coordinate = 0; coordinate < 3; ++coordinate)
            {
                firstDelta[coordinate] =
                    firstFacet.vertices[0][coordinate] - rayOrigin[coordinate];
                secondDelta[coordinate] =
                    secondFacet.vertices[0][coordinate] - rayOrigin[coordinate];
            }
            const double firstDepth =
                dot3_dev_d(firstDelta, firstNormal)/firstDenominator;
            const double secondDepth =
                dot3_dev_d(secondDelta, secondNormal)/secondDenominator;
            if (fabs(firstDepth - secondDepth) <= depthTolerance)
                continue;
            const int before = firstDepth < secondDepth ? first : second;
            const int after = firstDepth < secondDepth ? second : first;
            const unsigned int edge = 1u << after;
            if (!(precedes[before] & edge))
            {
                precedes[before] |= edge;
                ++indegree[after];
            }
        }
    }

    unsigned int emitted = 0;
    int output = 0;
    while (output < count)
    {
        int selected = count;
        for (int candidate = 0; candidate < count; ++candidate)
        {
            if (!(emitted & (1u << candidate)) && indegree[candidate] == 0)
            {
                selected = candidate;
                break;
            }
        }
        if (selected == count)
        {
            for (int candidate = 0; candidate < count; ++candidate)
                if (!(emitted & (1u << candidate)))
                    candidateIds[output++] =
                        static_cast<uint32_t>(ordered[candidate].facetId);
            break;
        }
        emitted |= 1u << selected;
        candidateIds[output++] =
            static_cast<uint32_t>(ordered[selected].facetId);
        for (int after = 0; after < count; ++after)
            if (precedes[selected] & (1u << after))
                --indegree[after];
    }
    return true;
}

template <bool TriangleFacets>
__device__ void trace_sort_facets_warp_dev(
    const GpuTraceBeamRecord &beam,
    const GpuTraceFacetRecord *facets,
    uint32_t *candidateIds,
    int count,
    double geometryScale,
    GpuTraceFacetDepth *ordered,
    double *sortAxes,
    unsigned int *precedes,
    double *planeDenominators,
    float (*projectedBounds)[4],
    GpuTracePoint2 *cachedProjected,
    int &state)
{
    const int lane = static_cast<int>(threadIdx.x);
    if (lane == 0)
    {
        state = count <= GPU_TRACE_SORT_CANDIDATES
             && beam.nVertices > 0
             && beam.nVertices <= GPU_TRACE_COMPACT_VERTICES
            ? 1 : 0;
    }
    __syncthreads();
    if (state == 0 || count < 2)
        return;

    for (int index = lane; index < count; index += static_cast<int>(blockDim.x))
    {
        const int facetId = static_cast<int>(candidateIds[index]);
        const GpuTraceFacetRecord &facet = facets[facetId];
        if (facet.nVertices < 3
            || facet.nVertices > GPU_TRACE_COMPACT_VERTICES)
        {
            atomicExch(&state, 0);
            continue;
        }
        double minimum = dot3_dev_d(facet.vertices[0], beam.dir);
        double maximum = minimum;
        for (int vertex = 1; vertex < facet.nVertices; ++vertex)
        {
            const double depth = dot3_dev_d(facet.vertices[vertex], beam.dir);
            minimum = fmin(minimum, depth);
            maximum = fmax(maximum, depth);
        }
        ordered[index] = {facetId, minimum, maximum};
    }
    __syncthreads();
    if (state == 0)
        return;

    const double scale = fmax(fabs(geometryScale), DBL_MIN);
    const double depthTolerance = 0x1p-40*scale;
    int sortWidth = 1;
    while (sortWidth < count)
        sortWidth <<= 1;
    GpuTraceFacetDepth sortValue = {0x7fffffff, DBL_MAX, DBL_MAX};
    if (lane < count)
        sortValue = ordered[lane];
    if (lane < sortWidth)
    {
        const unsigned int mask = sortWidth == 32
            ? 0xffffffffu : ((1u << sortWidth) - 1u);
        for (int sequence = 2; sequence <= sortWidth; sequence <<= 1)
        {
            for (int stride = sequence >> 1; stride > 0; stride >>= 1)
            {
                GpuTraceFacetDepth other;
                other.facetId = __shfl_xor_sync(
                    mask, sortValue.facetId, stride, sortWidth);
                other.depth = __shfl_xor_sync(
                    mask, sortValue.depth, stride, sortWidth);
                other.maximumDepth = __shfl_xor_sync(
                    mask, sortValue.maximumDepth, stride, sortWidth);
                const bool ascending = (lane & sequence) == 0;
                const bool lowerLane = (lane & stride) == 0;
                const bool takeMinimum = ascending == lowerLane;
                const bool otherLess = trace_depth_less(other, sortValue);
                if ((takeMinimum && otherLess)
                    || (!takeMinimum && !otherLess))
                {
                    sortValue = other;
                }
            }
        }
        if (lane < count)
            ordered[lane] = sortValue;
    }
    __syncthreads();

    if (lane == 0)
    {
        for (int begin = 0; begin < count;)
        {
            int end = begin + 1;
            while (end < count
                   && ordered[end].depth - ordered[begin].depth <= depthTolerance)
                ++end;
            for (int index = begin + 1; index < end; ++index)
            {
                const GpuTraceFacetDepth value = ordered[index];
                int position = index;
                while (position > begin
                       && value.facetId < ordered[position - 1].facetId)
                {
                    ordered[position] = ordered[position - 1];
                    --position;
                }
                ordered[position] = value;
            }
            begin = end;
        }

        bool overlappingIntervals = false;
        for (int first = 0; first + 1 < count; ++first)
        {
            if (ordered[first + 1].depth
                <= ordered[first].maximumDepth + depthTolerance)
            {
                overlappingIntervals = true;
                break;
            }
        }
        if (!overlappingIntervals)
        {
            for (int index = 0; index < count; ++index)
                candidateIds[index] =
                    static_cast<uint32_t>(ordered[index].facetId);
            state = 2;
        }
        else
        {
            double *direction = sortAxes;
            direction[0] = beam.dir[0];
            direction[1] = beam.dir[1];
            direction[2] = beam.dir[2];
            if (!normalize3_dev(direction))
            {
                for (int index = 0; index < count; ++index)
                    candidateIds[index] =
                        static_cast<uint32_t>(ordered[index].facetId);
                state = 2;
            }
            else
            {
                double reference[3];
                const double ax = fabs(direction[0]);
                const double ay = fabs(direction[1]);
                const double az = fabs(direction[2]);
                if (ax <= ay && ax <= az)
                {
                    reference[0] = 1.0;
                    reference[1] = 0.0;
                    reference[2] = 0.0;
                }
                else if (ay <= az)
                {
                    reference[0] = 0.0;
                    reference[1] = 1.0;
                    reference[2] = 0.0;
                }
                else
                {
                    reference[0] = 0.0;
                    reference[1] = 0.0;
                    reference[2] = 1.0;
                }
                double *axisU = sortAxes + 3;
                double *axisV = sortAxes + 6;
                cross3_dev(direction, reference, axisU);
                if (!normalize3_dev(axisU))
                {
                    state = 0;
                }
                else
                {
                    cross3_dev(direction, axisU, axisV);
                    if (!normalize3_dev(axisV))
                        state = 0;
                }
            }
        }
    }
    __syncthreads();
    if (state != 1)
        return;

    for (int index = lane; index < count; index += static_cast<int>(blockDim.x))
        precedes[index] = 0;
    __syncthreads();

    const double *direction = sortAxes;
    const double *axisU = sortAxes + 3;
    const double *axisV = sortAxes + 6;
    const double areaTolerance = 0x1p-40*scale*scale;
    if (TriangleFacets)
    {
        for (int candidate = lane; candidate < count;
             candidate += static_cast<int>(blockDim.x))
        {
            const GpuTraceFacetRecord &facet =
                facets[ordered[candidate].facetId];
            double minU = DBL_MAX, maxU = -DBL_MAX;
            double minV = DBL_MAX, maxV = -DBL_MAX;
            for (int vertex = 0; vertex < GPU_TRACE_TRIANGLE_VERTICES;
                 ++vertex)
            {
                const double u =
                    dot3_dev_d(facet.vertices[vertex], axisU);
                const double v =
                    dot3_dev_d(facet.vertices[vertex], axisV);
                cachedProjected[candidate*GPU_TRACE_TRIANGLE_VERTICES
                                + vertex] = {u, v};
                minU = fmin(minU, u);
                maxU = fmax(maxU, u);
                minV = fmin(minV, v);
                maxV = fmax(maxV, v);
            }
            projectedBounds[candidate][0] = nextafterf(
                static_cast<float>(minU - depthTolerance), -FLT_MAX);
            projectedBounds[candidate][1] = nextafterf(
                static_cast<float>(maxU + depthTolerance), FLT_MAX);
            projectedBounds[candidate][2] = nextafterf(
                static_cast<float>(minV - depthTolerance), -FLT_MAX);
            projectedBounds[candidate][3] = nextafterf(
                static_cast<float>(maxV + depthTolerance), FLT_MAX);
            planeDenominators[candidate] =
                dot3_dev_d(direction, facet.surfaceNormal);
        }
        __syncthreads();
    }
    const int pairCount = count*(count - 1)/2;
    for (int pairIndex = lane; pairIndex < pairCount;
         pairIndex += static_cast<int>(blockDim.x))
    {
        int first = 0;
        int rowSize = count - 1;
        int rowOffset = pairIndex;
        while (rowOffset >= rowSize)
        {
            rowOffset -= rowSize;
            ++first;
            --rowSize;
        }
        const int second = first + 1 + rowOffset;
        if (ordered[second].depth
            > ordered[first].maximumDepth + depthTolerance)
            continue;

        const GpuTraceFacetRecord &firstFacet =
            facets[ordered[first].facetId];
        const GpuTraceFacetRecord &secondFacet =
            facets[ordered[second].facetId];
        if (TriangleFacets
            && (projectedBounds[first][1] < projectedBounds[second][0]
                || projectedBounds[second][1]
                    < projectedBounds[first][0]
                || projectedBounds[first][3] < projectedBounds[second][2]
                || projectedBounds[second][3]
                    < projectedBounds[first][2]))
        {
            continue;
        }
        const double firstDenominator = TriangleFacets
            ? planeDenominators[first]
            : dot3_dev_d(direction, firstFacet.surfaceNormal);
        const double secondDenominator = TriangleFacets
            ? planeDenominators[second]
            : dot3_dev_d(direction, secondFacet.surfaceNormal);
        if (fabs(firstDenominator) <= 0x1p-40
            || fabs(secondDenominator) <= 0x1p-40)
            continue;
        GpuTracePoint2 overlapPoint;
        if (TriangleFacets)
        {
            if (!projected_overlap_point_cached_first_dev<true>(
                    cachedProjected
                        + first*GPU_TRACE_TRIANGLE_VERTICES,
                    GPU_TRACE_TRIANGLE_VERTICES, secondFacet,
                    cachedProjected
                        + second*GPU_TRACE_TRIANGLE_VERTICES,
                    axisU, axisV,
                    firstDenominator >= 0.0 ? 1.0 : -1.0,
                    secondDenominator >= 0.0 ? 1.0 : -1.0,
                    areaTolerance, overlapPoint))
                continue;
        }
        else if (!projected_overlap_point_dev(
                     firstFacet, secondFacet, axisU, axisV,
                     areaTolerance, overlapPoint))
        {
            continue;
        }

        double rayOrigin[3] = {
            axisU[0]*overlapPoint.x + axisV[0]*overlapPoint.y,
            axisU[1]*overlapPoint.x + axisV[1]*overlapPoint.y,
            axisU[2]*overlapPoint.x + axisV[2]*overlapPoint.y};
        const double *firstNormal = firstFacet.surfaceNormal;
        const double *secondNormal = secondFacet.surfaceNormal;
        double firstDelta[3], secondDelta[3];
        for (int coordinate = 0; coordinate < 3; ++coordinate)
        {
            firstDelta[coordinate] =
                firstFacet.vertices[0][coordinate] - rayOrigin[coordinate];
            secondDelta[coordinate] =
                secondFacet.vertices[0][coordinate] - rayOrigin[coordinate];
        }
        const double firstDepth =
            dot3_dev_d(firstDelta, firstNormal)/firstDenominator;
        const double secondDepth =
            dot3_dev_d(secondDelta, secondNormal)/secondDenominator;
        if (fabs(firstDepth - secondDepth) <= depthTolerance)
            continue;
        const int before = firstDepth < secondDepth ? first : second;
        const int after = firstDepth < secondDepth ? second : first;
        atomicOr(precedes + before, 1u << after);
    }
    __syncthreads();

    if (lane == 0)
    {
        int indegree[GPU_TRACE_SORT_CANDIDATES] = {};
        for (int before = 0; before < count; ++before)
            for (int after = 0; after < count; ++after)
                if (precedes[before] & (1u << after))
                    ++indegree[after];

        unsigned int emitted = 0;
        int output = 0;
        while (output < count)
        {
            int selected = count;
            for (int candidate = 0; candidate < count; ++candidate)
            {
                if (!(emitted & (1u << candidate))
                    && indegree[candidate] == 0)
                {
                    selected = candidate;
                    break;
                }
            }
            if (selected == count)
            {
                for (int candidate = 0; candidate < count; ++candidate)
                    if (!(emitted & (1u << candidate)))
                        candidateIds[output++] =
                            static_cast<uint32_t>(ordered[candidate].facetId);
                break;
            }
            emitted |= 1u << selected;
            candidateIds[output++] =
                static_cast<uint32_t>(ordered[selected].facetId);
            for (int after = 0; after < count; ++after)
                if (precedes[selected] & (1u << after))
                    --indegree[after];
        }
    }
    __syncthreads();
}

template <bool TriangleFacets>
__device__ void trace_sort_facets_block_dev(
    const GpuTraceBeamRecord &beam,
    const GpuTraceFacetRecord *facets,
    uint32_t *candidateIds,
    int count,
    double geometryScale,
    GpuTraceFacetDepth *ordered,
    double *sortAxes,
    unsigned int *precedes,
    int *indegree,
    int *intervalEnds,
    double *planeDenominators,
    GpuTracePoint2 *cachedProjected,
    int &state)
{
    const int lane = static_cast<int>(threadIdx.x);
    if (lane == 0)
    {
        state = count > GPU_TRACE_SORT_CANDIDATES
             && count <= GPU_TRACE_LARGE_SORT_CANDIDATES
            ? 1 : 0;
    }
    __syncthreads();
    if (state == 0)
        return;

    int sortWidth = 1;
    while (sortWidth < count)
        sortWidth <<= 1;
    for (int index = lane; index < sortWidth;
         index += static_cast<int>(blockDim.x))
    {
        if (index < count)
        {
            const int facetId = static_cast<int>(candidateIds[index]);
            const GpuTraceFacetRecord &facet = facets[facetId];
            if (facet.nVertices < 3
                || facet.nVertices > GPU_TRACE_COMPACT_VERTICES)
            {
                atomicExch(&state, 0);
                ordered[index] = {0x7fffffff, DBL_MAX, DBL_MAX};
                continue;
            }
            double minimum = dot3_dev_d(facet.vertices[0], beam.dir);
            double maximum = minimum;
            for (int vertex = 1; vertex < facet.nVertices; ++vertex)
            {
                const double depth = dot3_dev_d(
                    facet.vertices[vertex], beam.dir);
                minimum = fmin(minimum, depth);
                maximum = fmax(maximum, depth);
            }
            ordered[index] = {facetId, minimum, maximum};
        }
        else
        {
            ordered[index] = {0x7fffffff, DBL_MAX, DBL_MAX};
        }
    }
    __syncthreads();
    if (state == 0)
        return;

    for (int sequence = 2; sequence <= sortWidth; sequence <<= 1)
    {
        for (int stride = sequence >> 1; stride > 0; stride >>= 1)
        {
            for (int index = lane; index < sortWidth;
                 index += static_cast<int>(blockDim.x))
            {
                const int otherIndex = index ^ stride;
                if (otherIndex <= index)
                    continue;
                const GpuTraceFacetDepth left = ordered[index];
                const GpuTraceFacetDepth right = ordered[otherIndex];
                const bool ascending = (index & sequence) == 0;
                const bool swap = ascending
                    ? trace_depth_less(right, left)
                    : trace_depth_less(left, right);
                if (swap)
                {
                    ordered[index] = right;
                    ordered[otherIndex] = left;
                }
            }
            __syncthreads();
        }
    }

    const double scale = fmax(fabs(geometryScale), DBL_MIN);
    const double depthTolerance = 0x1p-40*scale;
    if (lane == 0)
    {
        for (int begin = 0; begin < count;)
        {
            int end = begin + 1;
            while (end < count
                   && ordered[end].depth - ordered[begin].depth
                      <= depthTolerance)
                ++end;
            for (int index = begin + 1; index < end; ++index)
            {
                const GpuTraceFacetDepth value = ordered[index];
                int position = index;
                while (position > begin
                       && value.facetId < ordered[position - 1].facetId)
                {
                    ordered[position] = ordered[position - 1];
                    --position;
                }
                ordered[position] = value;
            }

            begin = end;
        }

        bool overlappingIntervals = false;
        for (int first = 0; first + 1 < count; ++first)
        {
            if (ordered[first + 1].depth
                <= ordered[first].maximumDepth + depthTolerance)
            {
                overlappingIntervals = true;
                break;
            }
        }
        if (!overlappingIntervals)
        {
            for (int index = 0; index < count; ++index)
                candidateIds[index] =
                    static_cast<uint32_t>(ordered[index].facetId);
            state = 2;
        }
        else
        {
            double *direction = sortAxes;
            direction[0] = beam.dir[0];
            direction[1] = beam.dir[1];
            direction[2] = beam.dir[2];
            if (!normalize3_dev(direction))
            {
                for (int index = 0; index < count; ++index)
                    candidateIds[index] =
                        static_cast<uint32_t>(ordered[index].facetId);
                state = 2;
            }
            else
            {
                double reference[3];
                const double ax = fabs(direction[0]);
                const double ay = fabs(direction[1]);
                const double az = fabs(direction[2]);
                if (ax <= ay && ax <= az)
                {
                    reference[0] = 1.0;
                    reference[1] = 0.0;
                    reference[2] = 0.0;
                }
                else if (ay <= az)
                {
                    reference[0] = 0.0;
                    reference[1] = 1.0;
                    reference[2] = 0.0;
                }
                else
                {
                    reference[0] = 0.0;
                    reference[1] = 0.0;
                    reference[2] = 1.0;
                }
                double *axisU = sortAxes + 3;
                double *axisV = sortAxes + 6;
                cross3_dev(direction, reference, axisU);
                if (!normalize3_dev(axisU))
                {
                    state = 0;
                }
                else
                {
                    cross3_dev(direction, axisU, axisV);
                    if (!normalize3_dev(axisV))
                        state = 0;
                }
            }
        }
    }
    __syncthreads();
    if (state != 1)
        return;

    const int wordCount = (count + 31) >> 5;
    const int bitsetEntries = count*wordCount;
    for (int entry = lane; entry < bitsetEntries;
         entry += static_cast<int>(blockDim.x))
        precedes[entry] = 0;
    for (int candidate = lane; candidate < count;
         candidate += static_cast<int>(blockDim.x))
        indegree[candidate] = 0;
    __syncthreads();

    const double *direction = sortAxes;
    const double *axisU = sortAxes + 3;
    const double *axisV = sortAxes + 6;
    const double areaTolerance = 0x1p-40*scale*scale;
    for (int first = lane; first < count;
         first += static_cast<int>(blockDim.x))
    {
        int end = first + 1;
        while (end < count
               && ordered[end].depth
                  <= ordered[first].maximumDepth + depthTolerance)
            ++end;
        intervalEnds[first] = end;
    }
    __syncthreads();

    for (int candidate = lane; candidate < count;
         candidate += static_cast<int>(blockDim.x))
    {
        const GpuTraceFacetRecord &facet =
            facets[ordered[candidate].facetId];
        double minU = DBL_MAX, maxU = -DBL_MAX;
        double minV = DBL_MAX, maxV = -DBL_MAX;
        for (int vertex = 0; vertex < facet.nVertices; ++vertex)
        {
            const double u = dot3_dev_d(facet.vertices[vertex], axisU);
            const double v = dot3_dev_d(facet.vertices[vertex], axisV);
            minU = fmin(minU, u);
            maxU = fmax(maxU, u);
            minV = fmin(minV, v);
            maxV = fmax(maxV, v);
        }
        trace_store_projected_bounds(
            ordered[candidate],
            nextafterf(static_cast<float>(minU - depthTolerance),
                       -FLT_MAX),
            nextafterf(static_cast<float>(maxU + depthTolerance),
                       FLT_MAX),
            nextafterf(static_cast<float>(minV - depthTolerance),
                       -FLT_MAX),
            nextafterf(static_cast<float>(maxV + depthTolerance),
                       FLT_MAX));
        planeDenominators[candidate] =
            dot3_dev_d(direction, facet.surfaceNormal);
    }
    __syncthreads();

    constexpr int groupSize = 8;
    const int group = lane / groupSize;
    const int groupLane = lane & (groupSize - 1);
    const int groupCount = static_cast<int>(blockDim.x) / groupSize;
    const int groupOffset = (lane & 31) & ~(groupSize - 1);
    const unsigned int groupMask =
        ((1u << groupSize) - 1u) << groupOffset;
    constexpr int cachedVertices = TriangleFacets
        ? GPU_TRACE_TRIANGLE_VERTICES : GPU_TRACE_COMPACT_VERTICES;
    for (int first = group; first < count; first += groupCount)
    {
        const GpuTraceFacetRecord &firstFacet =
            facets[ordered[first].facetId];
        GpuTracePoint2 *firstProjected =
            cachedProjected + group*cachedVertices;
        for (int vertex = groupLane; vertex < firstFacet.nVertices;
             vertex += groupSize)
        {
            firstProjected[vertex].x =
                dot3_dev_d(firstFacet.vertices[vertex], axisU);
            firstProjected[vertex].y =
                dot3_dev_d(firstFacet.vertices[vertex], axisV);
        }
        __syncwarp(groupMask);
        const GpuTraceDepthScratch firstBounds =
            trace_load_projected_bounds(ordered[first]);
        const double firstDenominator = planeDenominators[first];
        if (fabs(firstDenominator) <= 0x1p-40)
            continue;
        for (int second = first + 1 + groupLane;
             second < intervalEnds[first]; second += groupSize)
        {
            const GpuTraceDepthScratch secondBounds =
                trace_load_projected_bounds(ordered[second]);
            if (firstBounds.projectedBounds[1]
                    < secondBounds.projectedBounds[0]
                || secondBounds.projectedBounds[1]
                    < firstBounds.projectedBounds[0]
                || firstBounds.projectedBounds[3]
                    < secondBounds.projectedBounds[2]
                || secondBounds.projectedBounds[3]
                    < firstBounds.projectedBounds[2])
            {
                continue;
            }
            const GpuTraceFacetRecord &secondFacet =
                facets[ordered[second].facetId];
            const double secondDenominator = planeDenominators[second];
            if (fabs(secondDenominator) <= 0x1p-40)
                continue;
            GpuTracePoint2 overlapPoint;
            if (!projected_overlap_point_cached_first_dev<TriangleFacets>(
                    firstProjected, firstFacet.nVertices,
                    secondFacet,
                    nullptr,
                    axisU, axisV,
                    firstDenominator >= 0.0 ? 1.0 : -1.0,
                    secondDenominator >= 0.0 ? 1.0 : -1.0,
                    areaTolerance, overlapPoint))
                continue;

            double rayOrigin[3] = {
                axisU[0]*overlapPoint.x + axisV[0]*overlapPoint.y,
                axisU[1]*overlapPoint.x + axisV[1]*overlapPoint.y,
                axisU[2]*overlapPoint.x + axisV[2]*overlapPoint.y};
            const double *firstNormal = firstFacet.surfaceNormal;
            const double *secondNormal = secondFacet.surfaceNormal;
            double firstDelta[3], secondDelta[3];
            for (int coordinate = 0; coordinate < 3; ++coordinate)
            {
                firstDelta[coordinate] =
                    firstFacet.vertices[0][coordinate]
                    - rayOrigin[coordinate];
                secondDelta[coordinate] =
                    secondFacet.vertices[0][coordinate]
                    - rayOrigin[coordinate];
            }
            const double firstDepth =
                dot3_dev_d(firstDelta, firstNormal)/firstDenominator;
            const double secondDepth =
                dot3_dev_d(secondDelta, secondNormal)/secondDenominator;
            if (fabs(firstDepth - secondDepth) <= depthTolerance)
                continue;
            const int before = firstDepth < secondDepth ? first : second;
            const int after = firstDepth < secondDepth ? second : first;
            const unsigned int afterMask = 1u << (after & 31);
            unsigned int *edgeWord = precedes
                                   + before*wordCount
                                   + (after >> 5);
            atomicOr(edgeWord, afterMask);
            atomicAdd(indegree + after, 1);
        }
        __syncwarp(groupMask);
    }
    __syncthreads();

    if (lane == 0)
    {
        unsigned int emitted[GPU_TRACE_LARGE_SORT_WORDS] = {};
        unsigned int available[GPU_TRACE_LARGE_SORT_WORDS] = {};
        for (int candidate = 0; candidate < count; ++candidate)
            if (indegree[candidate] == 0)
                available[candidate >> 5] |= 1u << (candidate & 31);
        int output = 0;
        while (output < count)
        {
            int selected = count;
            for (int word = 0; word < wordCount; ++word)
            {
                if (available[word] != 0)
                {
                    selected = (word << 5) + __ffs(available[word]) - 1;
                    break;
                }
            }
            if (selected == count)
            {
                for (int candidate = 0; candidate < count; ++candidate)
                {
                    if (!(emitted[candidate >> 5]
                          & (1u << (candidate & 31))))
                        candidateIds[output++] = static_cast<uint32_t>(
                            ordered[candidate].facetId);
                }
                break;
            }
            const unsigned int selectedMask = 1u << (selected & 31);
            emitted[selected >> 5] |= selectedMask;
            available[selected >> 5] &= ~selectedMask;
            candidateIds[output++] =
                static_cast<uint32_t>(ordered[selected].facetId);
            const unsigned int *row =
                precedes + selected*wordCount;
            for (int word = 0; word < wordCount; ++word)
            {
                unsigned int successors = row[word];
                while (successors != 0)
                {
                    const int bit = __ffs(successors) - 1;
                    const int after = (word << 5) + bit;
                    successors &= successors - 1;
                    if (--indegree[after] == 0)
                        available[word] |= 1u << bit;
                }
            }
        }
        state = 2;
    }
    __syncthreads();
}

template <bool PrefilterFirst, bool TriangleFacets>
__global__ void trace_prepare_kernel(const GpuTraceBeamRecord *beams,
                                     const double (*beamVertices)[3],
                                     const GpuTraceFacetRecord *facets,
                                     const GpuTraceItemRecord *items,
                                     const int *itemIndices,
                                     uint32_t *candidateIds,
                                     double geometryScale,
                                     int itemCount)
{
    const int launchIndex = (int)blockIdx.x;
    if (launchIndex >= itemCount)
        return;
    const int itemIndex = itemIndices == nullptr
        ? launchIndex : itemIndices[launchIndex];
    const GpuTraceItemRecord item = items[itemIndex];
    const GpuTraceBeamRecord &beam = beams[itemIndex];
    __shared__ int itemSorted;
    __shared__ GpuTraceFacetDepth ordered[GPU_TRACE_SORT_CANDIDATES];
    __shared__ double sortAxes[9];
    __shared__ unsigned int precedes[GPU_TRACE_SORT_CANDIDATES];
    __shared__ double planeDenominators[GPU_TRACE_SORT_CANDIDATES];
    __shared__ float projectedBounds[GPU_TRACE_SORT_CANDIDATES][4];
    constexpr int projectedEntries = TriangleFacets
        ? GPU_TRACE_SORT_CANDIDATES*GPU_TRACE_TRIANGLE_VERTICES : 1;
    __shared__ GpuTracePoint2 cachedProjected[projectedEntries];
    __shared__ int keptCount;

    if (PrefilterFirst)
    {
        if (beam.nVertices <= 0
            || beam.nVertices > GPU_TRACE_COMPACT_VERTICES)
        {
            if (threadIdx.x == 0)
                itemSorted = 0;
            return;
        }

        for (int candidate = (int)threadIdx.x;
             candidate < item.count; candidate += (int)blockDim.x)
        {
            const int resultIndex = item.offset + candidate;
            const int facetId = static_cast<int>(candidateIds[resultIndex]);
            const GpuTraceFacetRecord &facet = facets[facetId];
            const double *normal = beam.location == 0
                                 ? facet.normalIn : facet.normalOut;
            const double dp0 = dot3_dev_d(beam.dir, normal);
            if (fabs(dp0) < 0x1p-40)
            {
                candidateIds[resultIndex] = 0;
                continue;
            }

            const double ax = fabs(normal[0]);
            const double ay = fabs(normal[1]);
            const double az = fabs(normal[2]);
            const int drop = (ax > ay && ax > az) ? 0 : ((ay > az) ? 1 : 2);

            double bMinU = 1e300, bMaxU = -1e300;
            double bMinV = 1e300, bMaxV = -1e300;
            for (int vertex = 0; vertex < beam.nVertices; ++vertex)
            {
                const double *source =
                    beamVertices[beam.vertexOffset + vertex];
                double p[3] = {source[0], source[1], source[2]};
                const double t =
                    (dot3_dev_d(p, normal) + normal[3]) / dp0;
                p[0] -= beam.dir[0] * t;
                p[1] -= beam.dir[1] * t;
                p[2] -= beam.dir[2] * t;

                const double u = drop == 0 ? p[1] : p[0];
                const double v = drop == 2 ? p[1] : p[2];
                add_bounds_d(u, v, bMinU, bMaxU, bMinV, bMaxV);
            }

            const double *facetBounds = facet.bounds[drop];
            const double margin = facet.margin;
            const bool overlap = !(bMaxU < facetBounds[0] - margin
                                || facetBounds[1] < bMinU - margin
                                || bMaxV < facetBounds[2] - margin
                                || facetBounds[3] < bMinV - margin);
            candidateIds[resultIndex] = overlap
                ? static_cast<uint32_t>(facetId) + 1u : 0u;
        }
        __syncthreads();

        if (threadIdx.x == 0)
        {
            int kept = 0;
            for (int candidate = 0; candidate < item.count; ++candidate)
            {
                const uint32_t encoded =
                    candidateIds[item.offset + candidate];
                if (encoded != 0)
                    candidateIds[item.offset + kept++] = encoded - 1u;
            }
            keptCount = kept;
        }
        __syncthreads();

        trace_sort_facets_warp_dev<TriangleFacets>(
            beam, facets, candidateIds + item.offset,
            keptCount, geometryScale, ordered, sortAxes, precedes,
            planeDenominators, projectedBounds, cachedProjected, itemSorted);
        __syncthreads();

        if (threadIdx.x == 0 && itemSorted && item.count > 0)
        {
            uint32_t metadata = GPU_TRACE_SORTED_MARKER
                              | (static_cast<uint32_t>(keptCount)
                                 << GPU_TRACE_COUNT_SHIFT);
            if (keptCount > 0)
                metadata |= candidateIds[item.offset]
                          & GPU_TRACE_FACET_ID_MASK;
            candidateIds[item.offset] = metadata;
        }
        return;
    }

    trace_sort_facets_warp_dev<TriangleFacets>(
        beam, facets, candidateIds + item.offset,
        item.count, geometryScale, ordered, sortAxes, precedes,
        planeDenominators, projectedBounds, cachedProjected, itemSorted);

    for (int candidate = (int)threadIdx.x;
        candidate < item.count; candidate += (int)blockDim.x)
    {
        const int resultIndex = item.offset + candidate;
        if (!itemSorted)
            continue;

        const int facetId = static_cast<int>(candidateIds[resultIndex]);
        const GpuTraceFacetRecord &facet = facets[facetId];
        const double *normal = beam.location == 0
                             ? facet.normalIn : facet.normalOut;
        const double dp0 = dot3_dev_d(beam.dir, normal);
        if (fabs(dp0) < 0x1p-40)
        {
            candidateIds[resultIndex] = 0;
            continue;
        }

        const double ax = fabs(normal[0]);
        const double ay = fabs(normal[1]);
        const double az = fabs(normal[2]);
        const int drop = (ax > ay && ax > az) ? 0 : ((ay > az) ? 1 : 2);

        double bMinU = 1e300, bMaxU = -1e300;
        double bMinV = 1e300, bMaxV = -1e300;

        for (int vertex = 0; vertex < beam.nVertices; ++vertex)
        {
            const double *source =
                beamVertices[beam.vertexOffset + vertex];
            double p[3] = {source[0], source[1], source[2]};
            const double t = (dot3_dev_d(p, normal) + normal[3]) / dp0;
            p[0] -= beam.dir[0] * t;
            p[1] -= beam.dir[1] * t;
            p[2] -= beam.dir[2] * t;

            const double u = drop == 0 ? p[1] : p[0];
            const double v = drop == 2 ? p[1] : p[2];
            add_bounds_d(u, v, bMinU, bMaxU, bMinV, bMaxV);
        }

        const double *facetBounds = facet.bounds[drop];
        const double margin = facet.margin;
        const bool overlap = !(bMaxU < facetBounds[0] - margin
                               || facetBounds[1] < bMinU - margin
                               || bMaxV < facetBounds[2] - margin
                               || facetBounds[3] < bMinV - margin);
        candidateIds[resultIndex] = overlap
            ? static_cast<uint32_t>(facetId) + 1u : 0u;
    }
    __syncthreads();

    if (threadIdx.x == 0 && itemSorted && item.count > 0)
    {
        int kept = 0;
        for (int candidate = 0; candidate < item.count; ++candidate)
        {
            const uint32_t encoded = candidateIds[item.offset + candidate];
            if (encoded != 0)
                candidateIds[item.offset + kept++] = encoded - 1u;
        }
        uint32_t metadata = GPU_TRACE_SORTED_MARKER
                          | (static_cast<uint32_t>(kept)
                             << GPU_TRACE_COUNT_SHIFT);
        if (kept > 0)
            metadata |= candidateIds[item.offset] & GPU_TRACE_FACET_ID_MASK;
        candidateIds[item.offset] = metadata;
    }
}

__device__ __forceinline__ void trace_add_skip_vertex_bounds(
    const float vertex[3],
    const float direction[3],
    const float normal[4],
    float denominator,
    int drop,
    float &minU,
    float &maxU,
    float &minV,
    float &maxV)
{
    const float x0 = vertex[0];
    const float y0 = vertex[1];
    const float z0 = vertex[2];
    const float t = (x0*normal[0] + y0*normal[1]
                   + z0*normal[2] + normal[3]) / denominator;
    float u;
    float v;
    if (drop == 0)
    {
        u = y0 - direction[1]*t;
        v = z0 - direction[2]*t;
    }
    else if (drop == 1)
    {
        u = x0 - direction[0]*t;
        v = z0 - direction[2]*t;
    }
    else
    {
        u = x0 - direction[0]*t;
        v = y0 - direction[1]*t;
    }
    minU = fminf(minU, u);
    maxU = fmaxf(maxU, u);
    minV = fminf(minV, v);
    maxV = fmaxf(maxV, v);
}

template <bool TriangleFacets>
__global__ void trace_prepare_large_kernel(
    const GpuTraceBeamRecord *beams,
    const double (*beamVertices)[3],
    const GpuTraceFacetRecord *facets,
    const GpuTraceItemRecord *items,
    const int *itemIndices,
    uint32_t *candidateIds,
    double geometryScale,
    int itemCount,
    int markSkips)
{
    const int launchIndex = (int)blockIdx.x;
    if (launchIndex >= itemCount)
        return;
    const int itemIndex = itemIndices == nullptr
        ? launchIndex : itemIndices[launchIndex];
    const GpuTraceItemRecord item = items[itemIndex];
    if (item.count <= GPU_TRACE_SORT_CANDIDATES)
        return;
    const GpuTraceBeamRecord &beam = beams[itemIndex];
    __shared__ int itemSorted;
    __shared__ GpuTraceFacetDepth
        ordered[GPU_TRACE_LARGE_SORT_CANDIDATES];
    __shared__ double sortAxes[9];
    __shared__ unsigned int
        precedes[GPU_TRACE_LARGE_SORT_CANDIDATES
                 * GPU_TRACE_LARGE_SORT_WORDS];
    __shared__ int indegree[GPU_TRACE_LARGE_SORT_CANDIDATES];
    __shared__ int intervalEnds[GPU_TRACE_LARGE_SORT_CANDIDATES];
    __shared__ double
        planeDenominators[GPU_TRACE_LARGE_SORT_CANDIDATES];
    constexpr int cachedVertices = TriangleFacets
        ? GPU_TRACE_TRIANGLE_VERTICES : GPU_TRACE_COMPACT_VERTICES;
    constexpr int projectedEntries = 32*cachedVertices;
    __shared__ GpuTracePoint2 cachedProjected[projectedEntries];
    __shared__ float skipDirection[3];
    __shared__ float skipVertices[GPU_TRACE_COMPACT_VERTICES][3];

    trace_sort_facets_block_dev<TriangleFacets>(
        beam, facets, candidateIds + item.offset,
        item.count, geometryScale, ordered, sortAxes,
        precedes, indegree, intervalEnds, planeDenominators,
        cachedProjected, itemSorted);
    if (!itemSorted)
        return;

    if (markSkips && beam.nVertices > 0
        && beam.nVertices <= GPU_TRACE_COMPACT_VERTICES)
    {
        for (int coordinate = static_cast<int>(threadIdx.x);
             coordinate < 3; coordinate += static_cast<int>(blockDim.x))
            skipDirection[coordinate] =
                static_cast<float>(beam.dir[coordinate]);
        for (int element = static_cast<int>(threadIdx.x);
             element < beam.nVertices*3;
             element += static_cast<int>(blockDim.x))
        {
            skipVertices[element/3][element%3] = static_cast<float>(
                beamVertices[beam.vertexOffset + element/3][element%3]);
        }
        __syncthreads();
        for (int candidate = static_cast<int>(threadIdx.x);
             candidate < item.count;
             candidate += static_cast<int>(blockDim.x))
        {
            const int resultIndex = item.offset + candidate;
            const int facetId = static_cast<int>(
                candidateIds[resultIndex] & GPU_TRACE_FACET_ID_MASK);
            const GpuTraceFacetRecord &facet = facets[facetId];
            const double *normalD = beam.location == 0
                                  ? facet.normalIn : facet.normalOut;
            const float normal[4] = {
                static_cast<float>(normalD[0]),
                static_cast<float>(normalD[1]),
                static_cast<float>(normalD[2]),
                static_cast<float>(normalD[3])};
            const float denominator = skipDirection[0]*normal[0]
                                    + skipDirection[1]*normal[1]
                                    + skipDirection[2]*normal[2];
            const float absDenominator = fabsf(denominator);
            if (absDenominator <= 256.0f*FLT_EPSILON)
                continue;
            const float ax = fabsf(normal[0]);
            const float ay = fabsf(normal[1]);
            const float az = fabsf(normal[2]);
            const int drop = (ax > ay && ax > az)
                ? 0 : ((ay > az) ? 1 : 2);
            float minU = FLT_MAX, maxU = -FLT_MAX;
            float minV = FLT_MAX, maxV = -FLT_MAX;
            if (beam.nVertices <= 4)
            {
                #pragma unroll
                for (int vertex = 0; vertex < 4; ++vertex)
                {
                    if (vertex < beam.nVertices)
                        trace_add_skip_vertex_bounds(
                            skipVertices[vertex], skipDirection, normal,
                            denominator, drop,
                            minU, maxU, minV, maxV);
                }
            }
            else
            {
                for (int vertex = 0; vertex < beam.nVertices; ++vertex)
                {
                    trace_add_skip_vertex_bounds(
                        skipVertices[vertex], skipDirection, normal,
                        denominator, drop,
                        minU, maxU, minV, maxV);
                }
            }
            const double *bounds = facet.bounds[drop];
            const double roundingGuard = 64.0*FLT_EPSILON
                * fmax(1.0, fabs(geometryScale))
                * (1.0 + 1.0/static_cast<double>(absDenominator));
            const double margin = facet.margin + roundingGuard;
            const bool overlap = !(
                   static_cast<double>(maxU) < bounds[0] - margin
                || bounds[1] < static_cast<double>(minU) - margin
                || static_cast<double>(maxV) < bounds[2] - margin
                || bounds[3] < static_cast<double>(minV) - margin);
            if (!overlap)
                candidateIds[resultIndex] |= GPU_TRACE_SKIP_MARKER;
        }
        __syncthreads();
    }

    if (threadIdx.x == 0 && itemSorted && item.count > 0)
    {
        uint32_t metadata = GPU_TRACE_SORTED_MARKER
                          | (static_cast<uint32_t>(item.count)
                             << GPU_TRACE_COUNT_SHIFT);
        const uint32_t first = candidateIds[item.offset];
        metadata |= first & GPU_TRACE_FACET_ID_MASK;
        if ((first & GPU_TRACE_SKIP_MARKER) != 0)
            metadata |= GPU_TRACE_FIRST_SKIP_MARKER;
        candidateIds[item.offset] = metadata;
    }
}

static int gpu_trace_threshold()
{
    const char *value = std::getenv("MBS_GPU_TRACE_MIN_CANDIDATES");
    if (!value || !*value)
        return 1024;
    char *end = nullptr;
    long parsed = std::strtol(value, &end, 10);
    if (!end || *end != '\0' || parsed < 1)
        return 1024;
    return (int)parsed;
}

static int gpu_trace_large_threads()
{
    const char *value = std::getenv("MBS_GPU_TRACE_LARGE_THREADS");
    if (!value || !*value)
        return 64;
    char *end = nullptr;
    const long parsed = std::strtol(value, &end, 10);
    if (!end || *end != '\0'
        || (parsed != 64 && parsed != 128 && parsed != 256))
        return 64;
    return static_cast<int>(parsed);
}

bool GpuTracePrefilterBeamFacets(const Beam &beam,
                                 const Facet *facets,
                                 const IntArray &facetIds,
                                 std::vector<unsigned char> &mayIntersect)
{
    IntArray mutableFacetIds = facetIds;
    bool sortedOnGpu = false;
    GpuTraceBeamFacets item;
    item.beam = &beam;
    item.facetIds = &mutableFacetIds;
    item.mayIntersect = &mayIntersect;
    item.sortedOnGpu = &sortedOnGpu;
    std::vector<GpuTraceBeamFacets> items(1, item);
    double minCoord[3] = {1e300, 1e300, 1e300};
    double maxCoord[3] = {-1e300, -1e300, -1e300};
    const auto include = [&minCoord, &maxCoord](const Point3f &point) {
        for (int coordinate = 0; coordinate < 3; ++coordinate)
        {
            minCoord[coordinate] = std::min(
                minCoord[coordinate],
                static_cast<double>(point.coordinates[coordinate]));
            maxCoord[coordinate] = std::max(
                maxCoord[coordinate],
                static_cast<double>(point.coordinates[coordinate]));
        }
    };
    for (int vertex = 0; vertex < beam.nVertices; ++vertex)
        include(beam.arr[vertex]);
    for (size_t index = 0; index < facetIds.size; ++index)
    {
        const Facet &facet = facets[facetIds.arr[index]];
        for (int vertex = 0; vertex < facet.nVertices; ++vertex)
            include(facet.arr[vertex]);
    }
    const double dx = maxCoord[0] - minCoord[0];
    const double dy = maxCoord[1] - minCoord[1];
    const double dz = maxCoord[2] - minCoord[2];
    const double geometryScale = std::sqrt(dx*dx + dy*dy + dz*dz);
    const bool prepared = GpuTracePrepareBeamFacetBatch(
        facets, geometryScale, items);
    if (!prepared || !sortedOnGpu)
        return prepared;

    mayIntersect.assign(facetIds.size, 0);
    for (size_t original = 0; original < facetIds.size; ++original)
    {
        for (size_t sorted = 0; sorted < mutableFacetIds.size; ++sorted)
        {
            if (facetIds.arr[original] == mutableFacetIds.arr[sorted])
            {
                mayIntersect[original] = 1;
                break;
            }
        }
    }
    return true;
}

bool GpuTracePrepareBeamFacetBatch(const Facet *facets,
                                   double geometryScale,
                                   const std::vector<GpuTraceBeamFacets> &items)
{
    g_traceLastFailureWasOutOfMemory = false;
    size_t total = 0;
    size_t totalVertices = 0;
    for (size_t itemIdx = 0; itemIdx < items.size(); ++itemIdx)
    {
        GpuTraceBeamFacets item = items[itemIdx];
        item.mayIntersect->clear();
        if (item.sortedOnGpu != nullptr)
            *item.sortedOnGpu = false;
        total += item.facetIds->size;
        totalVertices += static_cast<size_t>(std::max(
            0, std::min(item.beam->nVertices,
                        GPU_TRACE_COMPACT_VERTICES)));
    }

    if (items.empty() || total < (size_t)gpu_trace_threshold())
        return false;
    if (!ensure_trace_stream())
        return false;
    const cudaStream_t stream = trace_stream();
    const GpuTraceInputLayout layout =
        trace_input_layout(items.size(), totalVertices, total);
    if (!ensure_trace_input_capacity(layout.bytes)
        || !ensure_trace_host_input_capacity(layout.bytes))
        return false;

    unsigned char *hostInput = g_traceWorkspace.hostInput;
    GpuTraceBeamRecord *hostBeams =
        reinterpret_cast<GpuTraceBeamRecord *>(hostInput + layout.beams);
    GpuTraceItemRecord *hostItems =
        reinterpret_cast<GpuTraceItemRecord *>(hostInput + layout.items);
    int *hostItemIndices =
        reinterpret_cast<int *>(hostInput + layout.itemIndices);
    double (*hostVertices)[3] = reinterpret_cast<double (*)[3]>(
        hostInput + layout.vertices);
    uint32_t *hostCandidateIds =
        reinterpret_cast<uint32_t *>(hostInput + layout.candidates);

    size_t out = 0;
    size_t vertexOut = 0;
    size_t smallItemCount = 0;
    size_t largeItemCount = 0;
    int maxFacetId = -1;
    for (size_t itemIdx = 0; itemIdx < items.size(); ++itemIdx)
    {
        const Beam &beam = *items[itemIdx].beam;
        GpuTraceBeamRecord &record = hostBeams[itemIdx];
        record.vertexOffset = static_cast<int>(vertexOut);
        record.nVertices = beam.nVertices;
        record.location = beam.location == Location::In ? 0 : 1;
        record.dir[0] = beam.direction.coordinates[0];
        record.dir[1] = beam.direction.coordinates[1];
        record.dir[2] = beam.direction.coordinates[2];
        const int copiedVertices = std::min(
            beam.nVertices, GPU_TRACE_COMPACT_VERTICES);
        for (int v = 0; v < copiedVertices; ++v)
            for (int k = 0; k < 3; ++k)
                hostVertices[vertexOut + static_cast<size_t>(v)][k] =
                    beam.arr[v].coordinates[k];
        vertexOut += static_cast<size_t>(std::max(0, copiedVertices));

        GpuTraceBeamFacets item = items[itemIdx];
        hostItems[itemIdx].offset = (int)out;
        hostItems[itemIdx].count = (int)item.facetIds->size;
        if (item.facetIds->size <= GPU_TRACE_SORT_CANDIDATES)
            hostItemIndices[smallItemCount++] = static_cast<int>(itemIdx);
        else
            hostItemIndices[items.size() - 1 - largeItemCount++] =
                static_cast<int>(itemIdx);
        for (size_t facetIdx = 0; facetIdx < item.facetIds->size; ++facetIdx)
        {
            const int facetId = item.facetIds->arr[facetIdx];
            maxFacetId = std::max(maxFacetId, facetId);
            hostCandidateIds[out++] = static_cast<uint32_t>(facetId);
        }
    }

    if (out == 0)
        return false;
    if (maxFacetId < 0
        || static_cast<uint32_t>(maxFacetId) > GPU_TRACE_FACET_ID_MASK)
        return false;
    if (!upload_trace_facets_if_needed(facets, maxFacetId, geometryScale))
        return false;
    const bool triangleFacets =
        g_traceWorkspace.copiedFacetsAreTriangles;

    const bool traceTiming = gpu_trace_timing_enabled()
                          && ensure_trace_timing_events();
    const auto h2dStarted = std::chrono::high_resolution_clock::now();
    if (!trace_cuda_ok(cudaMemcpyAsync(g_traceWorkspace.input, hostInput,
                                      layout.bytes, cudaMemcpyHostToDevice,
                                      stream),
                       "copy packed trace input"))
        return false;
    const auto h2dFinished = std::chrono::high_resolution_clock::now();
    GpuTraceBeamRecord *deviceBeams =
        reinterpret_cast<GpuTraceBeamRecord *>(
            g_traceWorkspace.input + layout.beams);
    GpuTraceItemRecord *deviceItems =
        reinterpret_cast<GpuTraceItemRecord *>(
            g_traceWorkspace.input + layout.items);
    const int *deviceItemIndices = reinterpret_cast<const int *>(
        g_traceWorkspace.input + layout.itemIndices);
    const double (*deviceVertices)[3] =
        reinterpret_cast<const double (*)[3]>(
            g_traceWorkspace.input + layout.vertices);
    uint32_t *deviceCandidateIds = reinterpret_cast<uint32_t *>(
        g_traceWorkspace.input + layout.candidates);
    const bool fullSort = gpu_trace_full_sort();
    const bool largeSkipFlags = fullSort && gpu_trace_large_skip_flags();
    if (traceTiming)
    {
        cudaEventRecord(g_traceWorkspace.timingStart, stream);
    }
    if (fullSort)
    {
        if (smallItemCount > 0 && triangleFacets)
        {
            trace_prepare_kernel<false, true><<<
                (int)smallItemCount, 32, 0, stream>>>(
                deviceBeams, deviceVertices, g_traceWorkspace.facets,
                deviceItems, deviceItemIndices, deviceCandidateIds,
                geometryScale, (int)smallItemCount);
        }
        else if (smallItemCount > 0)
        {
            trace_prepare_kernel<false, false><<<
                (int)smallItemCount, 32, 0, stream>>>(
                deviceBeams, deviceVertices, g_traceWorkspace.facets,
                deviceItems, deviceItemIndices, deviceCandidateIds,
                geometryScale, (int)smallItemCount);
        }
        if (!trace_cuda_ok(cudaGetLastError(), "trace_prepare_kernel"))
            return false;
        if (traceTiming)
            cudaEventRecord(g_traceWorkspace.timingSmallDone, stream);

        const int *largeItemIndices = deviceItemIndices
            + items.size() - largeItemCount;
        if (largeItemCount > 0 && triangleFacets)
        {
            trace_prepare_large_kernel<true><<<
                (int)largeItemCount, gpu_trace_large_threads(), 0, stream>>>(
                deviceBeams, deviceVertices, g_traceWorkspace.facets,
                deviceItems, largeItemIndices, deviceCandidateIds,
                geometryScale, (int)largeItemCount,
                largeSkipFlags ? 1 : 0);
        }
        else if (largeItemCount > 0)
        {
            trace_prepare_large_kernel<false><<<
                (int)largeItemCount, gpu_trace_large_threads(), 0, stream>>>(
                deviceBeams, deviceVertices, g_traceWorkspace.facets,
                deviceItems, largeItemIndices, deviceCandidateIds,
                geometryScale, (int)largeItemCount,
                largeSkipFlags ? 1 : 0);
        }
        if (!trace_cuda_ok(cudaGetLastError(),
                           "trace_prepare_large_kernel"))
            return false;
        if (traceTiming)
            cudaEventRecord(g_traceWorkspace.timingLargeDone, stream);
    }
    else if (gpu_trace_prefilter_first())
    {
        if (triangleFacets)
            trace_prepare_kernel<true, true><<<
                (int)items.size(), 32, 0, stream>>>(
                deviceBeams, deviceVertices, g_traceWorkspace.facets,
                deviceItems, nullptr, deviceCandidateIds,
                geometryScale, (int)items.size());
        else
            trace_prepare_kernel<true, false><<<
                (int)items.size(), 32, 0, stream>>>(
                deviceBeams, deviceVertices, g_traceWorkspace.facets,
                deviceItems, nullptr, deviceCandidateIds,
                geometryScale, (int)items.size());
    }
    else
    {
        if (triangleFacets)
            trace_prepare_kernel<false, true><<<
                (int)items.size(), 32, 0, stream>>>(
                deviceBeams, deviceVertices, g_traceWorkspace.facets,
                deviceItems, nullptr, deviceCandidateIds,
                geometryScale, (int)items.size());
        else
            trace_prepare_kernel<false, false><<<
                (int)items.size(), 32, 0, stream>>>(
                deviceBeams, deviceVertices, g_traceWorkspace.facets,
                deviceItems, nullptr, deviceCandidateIds,
                geometryScale, (int)items.size());
    }
    if (!fullSort
        && !trace_cuda_ok(cudaGetLastError(), "trace_prepare_kernel"))
        return false;
    if (traceTiming && !fullSort)
    {
        cudaEventRecord(g_traceWorkspace.timingSmallDone, stream);
        cudaEventRecord(g_traceWorkspace.timingLargeDone, stream);
    }

    const auto d2hStarted = std::chrono::high_resolution_clock::now();
    if (!trace_cuda_ok(cudaMemcpyAsync(hostCandidateIds,
                                      deviceCandidateIds,
                                      out * sizeof(uint32_t),
                                      cudaMemcpyDeviceToHost,
                                      stream),
                       "copy sorted trace candidate ids"))
        return false;
    if (!trace_cuda_ok(cudaStreamSynchronize(stream),
                       "wait for trace preparation"))
        return false;
    const auto d2hFinished = std::chrono::high_resolution_clock::now();
    if (traceTiming)
    {
        float smallMs = 0.0f;
        float largeMs = 0.0f;
        cudaEventElapsedTime(&smallMs, g_traceWorkspace.timingStart,
                             g_traceWorkspace.timingSmallDone);
        cudaEventElapsedTime(&largeMs, g_traceWorkspace.timingSmallDone,
                             g_traceWorkspace.timingLargeDone);
        const double h2dMs = std::chrono::duration<double, std::milli>(
            h2dFinished - h2dStarted).count();
        const double d2hMs = std::chrono::duration<double, std::milli>(
            d2hFinished - d2hStarted).count();
        std::fprintf(stderr,
                     "GPU trace timing: items=%zu candidates=%zu bytes=%zu "
                     "h2d_ms=%.6f small_ms=%.6f large_ms=%.6f "
                     "d2h_wait_ms=%.6f\n",
                     items.size(), out, layout.bytes, h2dMs,
                     static_cast<double>(smallMs),
                     static_cast<double>(largeMs), d2hMs);
    }

    for (size_t itemIdx = 0; itemIdx < items.size(); ++itemIdx)
    {
        GpuTraceBeamFacets item = items[itemIdx];
        const size_t offset = (size_t)hostItems[itemIdx].offset;
        const size_t originalCount = item.facetIds->size;
        const uint32_t metadata = originalCount == 0 ? 0
            : hostCandidateIds[offset];
        const bool sorted = originalCount == 0
            || (hostCandidateIds[offset]
                & GPU_TRACE_SORTED_MARKER) != 0;
        if (item.sortedOnGpu != nullptr)
            *item.sortedOnGpu = sorted;
        if (!sorted)
            continue;

        const size_t kept = originalCount == 0 ? 0
            : static_cast<size_t>((metadata & GPU_TRACE_COUNT_MASK)
                                  >> GPU_TRACE_COUNT_SHIFT);
        const size_t maximumKept = fullSort
            ? GPU_TRACE_LARGE_SORT_CANDIDATES
            : GPU_TRACE_SORT_CANDIDATES;
        if (kept > originalCount || kept > maximumKept)
            return false;
        const bool hasSkipFlags = largeSkipFlags
                               && originalCount > GPU_TRACE_SORT_CANDIDATES;
        if (hasSkipFlags)
            item.mayIntersect->assign(kept, 1);
        for (size_t facetIdx = 0; facetIdx < kept; ++facetIdx)
        {
            const uint32_t encoded = facetIdx == 0
                ? metadata : hostCandidateIds[offset + facetIdx];
            const uint32_t facetId = encoded & GPU_TRACE_FACET_ID_MASK;
            item.facetIds->arr[facetIdx] = static_cast<int>(facetId);
            if (hasSkipFlags)
            {
                const bool skipped = facetIdx == 0
                    ? (encoded & GPU_TRACE_FIRST_SKIP_MARKER) != 0
                    : (encoded & GPU_TRACE_SKIP_MARKER) != 0;
                if (skipped)
                    (*item.mayIntersect)[facetIdx] = 0;
            }
        }
        item.facetIds->size = kept;
        if (!hasSkipFlags)
            item.mayIntersect->clear();
    }
    return true;
}

bool GpuTraceLastFailureWasOutOfMemory()
{
    return g_traceLastFailureWasOutOfMemory;
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

bool GpuTracePrepareBeamFacetBatch(const Facet */*facets*/,
                                   double /*geometryScale*/,
                                   const std::vector<GpuTraceBeamFacets> &items)
{
    for (size_t i = 0; i < items.size(); ++i)
    {
        items[i].mayIntersect->assign(items[i].facetIds->size, 1);
        if (items[i].sortedOnGpu != nullptr)
            *items[i].sortedOnGpu = false;
    }
    return false;
}

bool GpuTraceLastFailureWasOutOfMemory()
{
    return false;
}

void GpuTraceInvalidateFacetCache()
{
}

#endif
