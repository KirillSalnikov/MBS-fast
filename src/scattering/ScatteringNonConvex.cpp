#include "ScatteringNonConvex.h"

#include "macro.h"
#include <tgmath.h>
#include <assert.h>
#include <iostream>
#include <cstdlib>
#include <algorithm>
#include <cmath>
#include <chrono>
#include <cstring>
#include <cstdint>
#include <utility>
#include <vector>
#include <iomanip>
#include <map>
#ifdef _OPENMP
#include <omp.h>
#endif

#include "BigIntegerLibrary.hh"
#ifdef USE_CUDA
#include "GpuTraceSupport.h"
#endif

//#ifdef _DEBUG // DEB
#include "Tracer.h"
//#endif

#ifdef USE_CUDA
namespace {
std::atomic<int> &GpuTraceBatchBeamLimitState()
{
    static std::atomic<int> limit([]() {
        const char *value = std::getenv("MBS_GPU_TRACE_BATCH_BEAMS");
        if (value && *value)
        {
            char *end = nullptr;
            const long parsed = std::strtol(value, &end, 10);
            if (end && *end == '\0' && parsed >= 1)
                return static_cast<int>(std::min<long>(parsed, 16384));
        }

        int workers = 1;
#ifdef _OPENMP
        const char *openmp = std::getenv("MBS_GPU_TRACE_OPENMP");
        if (openmp && openmp[0] == '1' && openmp[1] == '\0')
            workers = std::max(1, omp_get_max_threads());
#endif
        // Each CUDA grid reserves per-thread stack backing. Keep the total
        // number of simultaneously submitted blocks bounded when tracing
        // orientations from several OpenMP workers. Complex non-convex
        // particles can exhaust CUDA kernel-stack backing well before global
        // VRAM looks full. A 512-beam process-wide budget remained stable on
        // a 12 GiB consumer GPU with eight concurrent tracing workers.
        return std::max(64, (512 + workers - 1) / workers);
    }());
    return limit;
}

int GpuTraceBatchBeamLimit()
{
    return GpuTraceBatchBeamLimitState().load(std::memory_order_relaxed);
}

int GpuTraceReduceBatchBeamLimit(size_t failedBatchSize)
{
    const int reduced = static_cast<int>(
        std::max<size_t>(32, failedBatchSize / 2));
    std::atomic<int> &limit = GpuTraceBatchBeamLimitState();
    int current = limit.load(std::memory_order_relaxed);
    while (reduced < current
           && !limit.compare_exchange_weak(current, reduced,
                                           std::memory_order_relaxed))
    {
    }
    return limit.load(std::memory_order_relaxed);
}

bool GpuTracePrepareBeamFacetBatchAdaptive(
    const Facet *facets,
    double geometryScale,
    const std::vector<GpuTraceBeamFacets> &items)
{
    if (GpuTracePrepareBeamFacetBatch(facets, geometryScale, items))
        return true;
    if (!GpuTraceLastFailureWasOutOfMemory() || items.size() <= 32)
        return false;

    const int reducedLimit = GpuTraceReduceBatchBeamLimit(items.size());
    static std::atomic<int> lastWarnedLimit(0);
    const int previousWarning = lastWarnedLimit.exchange(
        reducedLimit, std::memory_order_relaxed);
    if (previousWarning != reducedLimit)
    {
        std::cerr
            << "WARNING: CUDA trace batch exhausted device memory; "
            << "retrying in smaller pieces and reducing subsequent batches to "
            << reducedLimit << " beams."
            << std::endl;
    }

    const size_t middle = items.size() / 2;
    const std::vector<GpuTraceBeamFacets> first(
        items.begin(), items.begin() + middle);
    const std::vector<GpuTraceBeamFacets> second(
        items.begin() + middle, items.end());
    const bool firstPrepared = GpuTracePrepareBeamFacetBatchAdaptive(
        facets, geometryScale, first);
    const bool secondPrepared = GpuTracePrepareBeamFacetBatchAdaptive(
        facets, geometryScale, second);
    return firstPrepared && secondPrepared;
}

bool GpuTraceAllowOpenMP()
{
    const char *value = std::getenv("MBS_GPU_TRACE_OPENMP");
    return !(value && value[0] == '0' && value[1] == '\0');
}

bool GpuTraceVerifySort()
{
    const char *value = std::getenv("MBS_GPU_TRACE_VERIFY_SORT");
    return value && value[0] == '1' && value[1] == '\0';
}

bool GpuTraceFullSortEnabled()
{
    const char *value = std::getenv("MBS_GPU_TRACE_FULL_SORT");
    return !(value && value[0] == '0' && value[1] == '\0');
}

bool GpuTraceSortCacheEnabled()
{
    const char *value = std::getenv("MBS_GPU_TRACE_SORT_CACHE");
    return GpuTraceFullSortEnabled()
        && !(value && value[0] == '0' && value[1] == '\0');
}

bool GpuTraceRuntimeEnabled()
{
#ifdef _OPENMP
    if (omp_get_max_threads() > 1 && !GpuTraceAllowOpenMP())
    {
        static std::atomic<bool> warned(false);
        if (!warned.exchange(true, std::memory_order_relaxed))
        {
            std::cerr
                << "WARNING: --gpu-trace-prefilter is disabled with "
                << "--threads > 1. Use --threads 1, or set "
                << "MBS_GPU_TRACE_OPENMP=1 with "
                << "--allow-experimental-environment after validating "
                << "GPU memory use."
                << std::endl;
        }
        return false;
    }
#endif
    return true;
}
}
#endif

namespace {
bool PgLikePartShadow()
{
    const char *v = std::getenv("MBS_AGGREGATE_PART_SHADOW");
    return v && *v;
}

double HullCross(const Point2f &o, const Point2f &a, const Point2f &b)
{
    return (a.x - o.x)*(b.y - o.y) - (a.y - o.y)*(b.x - o.x);
}

double CpuTraceProjectedPrefilterMargin()
{
    static const double margin = []() {
        const char *value = std::getenv("MBS_TRACE_CPU_PREFILTER_MARGIN");
        if (!value || !*value)
            return -1.0;
        char *end = nullptr;
        double parsed = std::strtod(value, &end);
        if (!end || *end != '\0' || parsed < 0.0)
            return -1.0;
        return parsed;
    }();
    return margin;
}

bool FacetReachesForwardHalfSpace(const Facet &facet,
                                  const Point3f &origin,
                                  const Point3f &forward,
                                  double tolerance)
{
    // A facet may cross the source plane even when its arithmetic center is
    // behind it.  Reject it only when every vertex is safely behind the plane;
    // otherwise later direction and projected-AABB tests decide intersection.
    for (int vertex = 0; vertex < facet.nVertices; ++vertex)
    {
        if (DotProduct(forward, facet.arr[vertex] - origin) >= -tolerance)
            return true;
    }
    return false;
}

double ProjectedCross(const Point2f &a, const Point2f &b, const Point2f &c)
{
    return (b.x - a.x)*(c.y - a.y) - (b.y - a.y)*(c.x - a.x);
}

struct ProjectedPolygon
{
    Point2f vertices[2*MAX_VERTEX_NUM];
    size_t size = 0;

    void Clear()
    {
        size = 0;
    }

    void Push(const Point2f &point)
    {
        if (size >= 2*MAX_VERTEX_NUM)
            throw std::runtime_error(
                "projected facet intersection exceeds geometry limit");
        vertices[size++] = point;
    }

    void Assign(const ProjectedPolygon &other)
    {
        size = other.size;
        for (size_t index = 0; index < size; ++index)
            vertices[index] = other.vertices[index];
    }
};

double ProjectedSignedArea(const ProjectedPolygon &polygon)
{
    double twiceArea = 0.0;
    for (size_t i = 0; i < polygon.size; ++i)
    {
        const Point2f &a = polygon.vertices[i];
        const Point2f &b = polygon.vertices[(i + 1) % polygon.size];
        twiceArea += a.x*b.y - a.y*b.x;
    }
    return 0.5*twiceArea;
}

void IntersectProjectedFacets(const ProjectedPolygon &subject,
                              const ProjectedPolygon &clip,
                              double areaTolerance,
                              ProjectedPolygon &overlap)
{
    ProjectedPolygon first;
    ProjectedPolygon second;
    first.Assign(subject);
    ProjectedPolygon *input = &first;
    ProjectedPolygon *result = &second;
    const double orientation = ProjectedSignedArea(clip) >= 0.0 ? 1.0 : -1.0;
    for (size_t edge = 0; edge < clip.size && input->size != 0; ++edge)
    {
        const Point2f &a = clip.vertices[edge];
        const Point2f &b = clip.vertices[(edge + 1) % clip.size];
        result->Clear();
        Point2f previous = input->vertices[input->size - 1];
        double previousSide = orientation*ProjectedCross(a, b, previous);
        bool previousInside = previousSide >= -areaTolerance;
        for (size_t index = 0; index < input->size; ++index)
        {
            const Point2f &current = input->vertices[index];
            const double currentSide = orientation*ProjectedCross(a, b, current);
            const bool currentInside = currentSide >= -areaTolerance;
            if (currentInside != previousInside)
            {
                const double denominator = previousSide - currentSide;
                if (std::fabs(denominator) > DBL_MIN)
                {
                    const double t = previousSide/denominator;
                    result->Push(Point2f(
                        previous.x + t*(current.x - previous.x),
                        previous.y + t*(current.y - previous.y)));
                }
            }
            if (currentInside)
                result->Push(current);
            previous = current;
            previousSide = currentSide;
            previousInside = currentInside;
        }
        std::swap(input, result);
    }
    overlap.Assign(*input);
}

bool ProjectedOverlapPoint(const ProjectedPolygon &first,
                           const ProjectedPolygon &second,
                           double areaTolerance,
                           Point2f &point)
{
    ProjectedPolygon overlap;
    IntersectProjectedFacets(first, second, areaTolerance, overlap);
    if (overlap.size < 3
        || std::fabs(ProjectedSignedArea(overlap)) <= areaTolerance)
        return false;

    point = Point2f();
    for (size_t index = 0; index < overlap.size; ++index)
    {
        const Point2f &vertex = overlap.vertices[index];
        point.x += vertex.x;
        point.y += vertex.y;
    }
    point.x /= overlap.size;
    point.y /= overlap.size;
    return true;
}

}

#ifdef _DEBUG
using namespace std;
ofstream trackMapFile("tracks_deb.dat", ios::out);
#endif

using namespace std;

ScatteringNonConvex::ScatteringNonConvex(Particle *particle, Light *incidentLight,
                                         bool isOpticalPath, int nActs)
    : Scattering(particle, incidentLight, isOpticalPath, nActs)
{
}

bool ScatteringNonConvex::ScatterLight(double /*beta*/, double /*gamma*/,
                                       std::vector<Beam> &scaterredBeams)
{
    const auto traceStarted = std::chrono::high_resolution_clock::now();
    // m_particle->Rotate(beta, gamma, 0);
    scaterredBeams.reserve(scaterredBeams.size()
                           + std::max(8192, 4 * m_particle->nFacets));
    m_traceSortCalls = 0;
    m_traceSortCandidates = 0;
    m_traceSortOverlapCalls = 0;
    if (!m_visibilityCacheBuilt)
        BuildFacetVisibilityCache();
    BuildFacetNormalCache();
    if (m_traceCpuProjectedPrefilter)
        BuildFacetProjectionCache();
#ifdef USE_CUDA
    if (m_gpuTracePrefilter)
        GpuTraceInvalidateFacetCache();
#endif
    const auto initialStarted = std::chrono::high_resolution_clock::now();
    if (!SplitLightToBeams())
        return false;
    const auto splitStarted = std::chrono::high_resolution_clock::now();
    const bool ok = SplitBeams(scaterredBeams);
    if (m_tracePrefilterStats)
    {
        const auto finished = std::chrono::high_resolution_clock::now();
        std::cout << "Trace phase stats: setup="
                  << std::chrono::duration<double>(initialStarted
                                                    - traceStarted).count()
                  << " initial="
                  << std::chrono::duration<double>(splitStarted
                                                    - initialStarted).count()
                  << " recursive="
                  << std::chrono::duration<double>(finished
                                                    - splitStarted).count()
                  << " scatter_total="
                  << std::chrono::duration<double>(finished
                                                    - traceStarted).count()
                  << std::endl;
    }
    return ok;
}

bool ScatteringNonConvex::PushBeamsToTree(int facetId,
                                          const PolygonArray &polygons,
                                          Beam &inBeam, Beam &outBeam)
{
    auto id = RecomputeTrackId(0, facetId);
    inBeam.SetTracingParams(facetId, 0, Location::In);
    outBeam.SetTracingParams(facetId, 0, Location::Out);
    inBeam.id = id;
    outBeam.id = id;

    for (unsigned j = 0; j < polygons.size; ++j)
    {
        Beam in = inBeam;
        Beam out = outBeam;
        const Polygon &pol = polygons.arr[j];
        in = pol;
        out = pol;
#ifdef _DEBUG // DEB
        in.pols.push_back(pol);
        out.pols.push_back(pol);
#endif
        Point3f p = pol.Center();
        in.opticalPath = 0;
        out.opticalPath = 0;
        double path = splitting.ComputeIncidentOpticalPath(m_incidentDir, p);
#ifdef _DEBUG // DEB
        in.ops.push_back(path);
        out.ops.push_back(path);
#endif
        in.AddOpticalPath(path);
        out.AddOpticalPath(path);
        UpdateTraceReference(in);
        UpdateTraceReference(out);
        if (const char *initialDebugPath = std::getenv("MBS_INITIAL_BEAMR_DEBUG"))
        {
            if (*initialDebugPath)
            {
                std::ofstream dbg(initialDebugPath, std::ios::app);
                dbg << std::scientific << std::setprecision(15)
                    << facetId << " " << j << " " << out.nVertices
                    << " " << out.Area()
                    << " " << out.direction.cx << " " << out.direction.cy << " " << out.direction.cz
                    << " " << m_facets[facetId].ex_normal.cx << " " << m_facets[facetId].ex_normal.cy << " " << m_facets[facetId].ex_normal.cz
                    << " " << out.polarizationBasis.cx << " " << out.polarizationBasis.cy << " " << out.polarizationBasis.cz
                    << " " << real(out.J.m11) << " " << imag(out.J.m11)
                    << " " << real(out.J.m12) << " " << imag(out.J.m12)
                    << " " << real(out.J.m21) << " " << imag(out.J.m21)
                    << " " << real(out.J.m22) << " " << imag(out.J.m22)
                    << " " << out.opticalPath;
                for (int vi = 0; vi < out.nVertices; ++vi)
                    dbg << " " << out.arr[vi].cx << " " << out.arr[vi].cy << " " << out.arr[vi].cz;
                dbg << std::endl;
            }
        }
#ifdef _DEBUG // DEB
        in.dirs.push_back(in.direction);
        out.dirs.push_back(out.direction);
#endif
        if (!Scattering::PushBeamToTree(in, facetId, 0, Location::In)
            || !Scattering::PushBeamToTree(out, facetId, 0, Location::Out))
            return false;
        AddProjectedIncidentEnergy(facetId, out);
    }

    return true;
}

bool ScatteringNonConvex::SplitByFacet(const IntArray &facetIDs, int facetIndex)
{
    m_polygonResultBuffer.Clear();
    IntersectWithFacet(facetIDs, facetIndex, m_polygonResultBuffer);

    if (m_polygonResultBuffer.size != 0)
    {
        int id = facetIDs.arr[facetIndex];
        Beam inBeam, outBeam;

#ifdef _DEBUG // DEB
        auto id0 = RecomputeTrackId(0, id);
        if (id0 == 45)
            int fff = 0;
#endif
        SetIncidentBeamOpticalParams(id, inBeam, outBeam);
        if (!PushBeamsToTree(id, m_polygonResultBuffer, inBeam, outBeam))
            return false;
    }

    return true;
}

bool ScatteringNonConvex::SplitLightToBeams()
{
    m_incidentEnergy = 0;
    m_treeSize = 0;
    m_beamTree.clear();
    ResetTraceReference();

    IntArray facetIDs;
    SelectVisibleFacetsForLight(facetIDs);

    for (size_t i = 0; i < facetIDs.size; ++i)
    {
        if (!SplitByFacet(facetIDs, (int)i))
            return false;
    }
    return true;
//#ifdef _DEBUG // DEB
//	ofstream logfile("logscat1.txt", ios::out);
//	for (int i = 0; i < m_treeSize; ++i)
//	{
//		logfile << m_beamTree[i].id << " ";
//	}
//#endif
}

void ScatteringNonConvex::SelectVisibleFacetsForLight(IntArray &facetIDs)
{
    FindVisibleFacetsForLight(facetIDs);
    SortFacets_faster(m_incidentDir, facetIDs);
}

void ScatteringNonConvex::IntersectWithFacet(const IntArray &facetIds, int prevFacetNum,
                                             PolygonArray &resFacets)
{
    int id = facetIds.arr[prevFacetNum];

    if (prevFacetNum == 0 || m_facets[id].isVisibleOut)
    {
        if (PgLikePartShadow())
            CutPolygonByAggregateParts(m_facets[id], id, m_facets[id].ex_normal, m_incidentDir, resFacets);
        else
            resFacets.Push(m_facets[id]);
    }
    else // facet is probably shadowed by others
    {
        const Facet &facet = m_facets[id];
        const Point3f &normal = facet.ex_normal;

        if (PgLikePartShadow())
            CutPolygonByAggregateParts(facet, id, normal, m_incidentDir, resFacets);
        else
            CutPolygonByFacets(facet, facetIds, prevFacetNum, normal, normal,
                               m_incidentDir, resFacets);
    }
}

void ScatteringNonConvex::SelectVisibleFacets(const Beam &beam, IntArray &facetIDs)
{
    FindVisibleFacets(beam, facetIDs);
    SortVisibleFacets(beam, facetIDs);
}

void ScatteringNonConvex::SortVisibleFacets(const Beam &beam,
                                             IntArray &facetIDs)
{
    Point3f dir = beam.direction;
    dir.d_param = m_facets[beam.lastFacetId].in_normal.d_param;
    SortFacets_faster(dir, facetIDs);
}


void ScatteringNonConvex::CutPolygonByAggregateParts(const Polygon &pol, int facetId,
                                                     const Vector3f &polNormal,
                                                     const Vector3f &dir,
                                                     PolygonArray &pols)
{
    pols.Push(pol);

    if (!m_particle->isAggregated || m_particle->nFacetsInPart <= 0)
        return;

    int targetBegin = 0, targetEnd = 0;
    m_particle->GetParticalFacetIdRangeByFacetId(facetId, targetBegin, targetEnd);
    const int facetsPerPart = m_particle->nFacetsInPart;
    const int nParts = m_particle->nFacets / facetsPerPart;

    Point3f origin = pol.arr[0];
    Point3f axisU = pol.arr[1] - pol.arr[0];
    Normalize(axisU);
    Point3f axisV = CrossProduct(polNormal, axisU);
    Normalize(axisV);
    const double denom = DotProduct(dir, polNormal);
    if (std::fabs(denom) < EPS_PROJECTION)
        return;

    for (int part = 0; part < nParts; ++part)
    {
        const int begin = part * facetsPerPart;
        const int end = begin + facetsPerPart;
        if (begin == targetBegin)
            continue;

        bool mayBlock = false;
        std::vector<Point2f> projected;
        projected.reserve(facetsPerPart * 6);

        double targetMinX = 1.0e100, targetMinY = 1.0e100;
        double targetMaxX = -1.0e100, targetMaxY = -1.0e100;
        for (int vi = 0; vi < pol.nVertices; ++vi)
        {
            Point3f rel = pol.arr[vi] - origin;
            const double x = DotProduct(rel, axisU);
            const double y = DotProduct(rel, axisV);
            targetMinX = std::min(targetMinX, x); targetMaxX = std::max(targetMaxX, x);
            targetMinY = std::min(targetMinY, y); targetMaxY = std::max(targetMaxY, y);
        }

        for (int fid = begin; fid < end; ++fid)
        {
            const Facet &f = m_facets[fid];
            if (DotProduct(f.ex_normal, dir) >= -EPS_PROJECTION)
                continue;

            bool faceInFront = false;
            double faceMinX = 1.0e100, faceMinY = 1.0e100;
            double faceMaxX = -1.0e100, faceMaxY = -1.0e100;
            std::vector<Point2f> faceProjected;
            faceProjected.reserve(f.nVertices);

            for (int vi = 0; vi < f.nVertices; ++vi)
            {
                const Point3f &src = f.arr[vi];
                const double t = (DotProduct(src, polNormal) + polNormal.d_param) / denom;
                const char *depthMode = std::getenv("MBS_AGGREGATE_PART_SHADOW_DEPTH");
                const double depthTolerance = geometry_length_tolerance(
                    m_geometryScale);
                if (depthMode && depthMode[0] == '1')
                {
                    if (t > depthTolerance)
                        faceInFront = true;
                }
                else if (depthMode && depthMode[0] == '2')
                {
                    if (t < -depthTolerance)
                        faceInFront = true;
                }
                else if (DotProduct(src - pol.arr[0], polNormal)
                         > geometry_length_tolerance(
                               m_geometryScale))
                {
                    faceInFront = true;
                }
                Point3f pp = src - dir * t;
                Point3f rel = pp - origin;
                const double x = DotProduct(rel, axisU);
                const double y = DotProduct(rel, axisV);
                faceMinX = std::min(faceMinX, x); faceMaxX = std::max(faceMaxX, x);
                faceMinY = std::min(faceMinY, y); faceMaxY = std::max(faceMaxY, y);
                faceProjected.push_back(Point2f(x, y));
            }

            const double bboxTolerance = geometry_length_tolerance(
                m_geometryScale);
            const bool bboxOverlap = faceMaxX + bboxTolerance >= targetMinX
                                  && faceMinX - bboxTolerance <= targetMaxX
                                  && faceMaxY + bboxTolerance >= targetMinY
                                  && faceMinY - bboxTolerance <= targetMaxY;
            if (!faceInFront || !bboxOverlap)
                continue;

            mayBlock = true;
            projected.insert(projected.end(), faceProjected.begin(), faceProjected.end());
        }

        if (!mayBlock || projected.size() < 3)
            continue;

        std::sort(projected.begin(), projected.end(),
            [](const Point2f &a, const Point2f &b) {
                return a.x < b.x || (a.x == b.x && a.y < b.y);
            });
        std::vector<Point2f> unique;
        unique.reserve(projected.size());
        for (const Point2f &p2 : projected)
        {
            const double tolerance = geometry_length_tolerance(
                m_geometryScale);
            if (unique.empty()
                || std::fabs(static_cast<double>(unique.back().x) - p2.x)
                       > tolerance
                || std::fabs(static_cast<double>(unique.back().y) - p2.y)
                       > tolerance)
                unique.push_back(p2);
        }
        if (unique.size() < 3)
            continue;

        std::vector<Point2f> hull;
        for (const Point2f &p2 : unique)
        {
            while (hull.size() >= 2
                   && HullCross(hull[hull.size()-2], hull[hull.size()-1], p2)
                      <= 0.0)
                hull.pop_back();
            hull.push_back(p2);
        }
        const size_t lower = hull.size();
        for (int i = (int)unique.size() - 2; i >= 0; --i)
        {
            const Point2f &p2 = unique[i];
            while (hull.size() > lower
                   && HullCross(hull[hull.size()-2], hull[hull.size()-1], p2)
                      <= 0.0)
                hull.pop_back();
            hull.push_back(p2);
        }
        if (hull.size() > 1)
            hull.pop_back();
        if (hull.size() < 3)
            continue;
        if (hull.size() > MAX_VERTEX_NUM)
            throw std::runtime_error(
                "aggregate shadow hull exceeds the 64-vertex geometry limit");

        Polygon clip;
        for (const Point2f &p2 : hull)
            clip.AddVertex(origin + axisU * p2.x + axisV * p2.y);
        if (DotProduct(clip.Normal(), polNormal) < 0)
            Polygon::InverseVertexOrder(clip);

        m_polygonBuffer.Clear();
        while (pols.size != 0)
        {
            const Polygon &subj = pols.Pop();
            Difference(subj, polNormal, clip, polNormal, dir, m_polygonBuffer);
        }
        if (m_polygonBuffer.size == 0)
            break;
        for (unsigned i = 0; i < m_polygonBuffer.size; ++i)
            pols.Push(m_polygonBuffer.arr[i]);
    }
}

void ScatteringNonConvex::CutPolygonByFacets(const Polygon &pol,
                                             const IntArray &facetIds, size_t size,
                                             const Vector3f &polNormal,
                                             const Vector3f &clipNormal,
                                             const Vector3f &dir,
                                             PolygonArray &pols)
{
    pols.Push(pol);

    // cut facet projections out of polygon one by one
    for (unsigned i = 0; i < size; ++i)
    {
        int id = facetIds.arr[i];
        m_polygonBuffer.Clear();

        while (pols.size != 0)
        {
            const Polygon &subj = pols.Pop();
            const Polygon &clip = m_facets[id];

            /// REF: объединить 2 первых аргумента и 2 вторых
            Difference(subj, polNormal, clip, clipNormal, dir, m_polygonBuffer);
        }

        if (m_polygonBuffer.size != 0)
        {
            for (unsigned i = 0; i < m_polygonBuffer.size; ++i)
            {
                pols.Push(m_polygonBuffer.arr[i]);
            }
        }
        else // beam has layed on the facet totally
        {
            break;
        }
    }
}

void ScatteringNonConvex::CutExternalBeam(const Beam &beam,
                                          std::vector<Beam> &scaterredBeams,
                                          const IntArray *readyFacetIds)
{
    const Point3f &n1 = m_facets[beam.lastFacetId].ex_normal;
    const Point3f &n2 = m_facets[beam.lastFacetId].in_normal;
    const bool debugCut = []() {
        const char *value = std::getenv("MBS_CUT_DEBUG");
        return value && value[0] == '1' && value[1] == 0;
    }();
    const double sourceArea = debugCut ? beam.Area() : 0.0;

    IntArray localFacetIds;
    const IntArray *facetIds = readyFacetIds;
    if (facetIds == nullptr)
    {
        SelectVisibleFacets(beam, localFacetIds);
        facetIds = &localFacetIds;
    }

#ifdef _DEBUG // DEB
    if (beam.id == 66557)
        int fff = 0;
#endif
    m_polygonResultBuffer.Clear();
    CutPolygonByFacets(beam, *facetIds, facetIds->size, n1, n2,
                       -beam.direction, m_polygonResultBuffer);

    Beam tmp = beam;
    double path = splitting.ComputeOutgoingOpticalPath(tmp); // добираем оптический путь
    tmp.opticalPath += path;
#ifdef _DEBUG // DEB
    tmp.ops.push_back(path);
#endif
    double cutArea = 0.0;
    for (unsigned i = 0; i < m_polygonResultBuffer.size; ++i)
    {
        tmp.SetPolygon(m_polygonResultBuffer.arr[i]);
        if (debugCut)
            cutArea += tmp.Area();
        scaterredBeams.push_back(tmp);
    }
    if (debugCut)
    {
        std::cout << "MBS_CUT_DEBUG lastFacet=" << beam.lastFacetId
                  << " nActs=" << beam.nActs
                  << " loc=" << (beam.location == Location::Out ? "Out" : "In")
                  << " sourceArea=" << sourceArea
                  << " cutArea=" << cutArea
                  << " candidates=" << facetIds->size
                  << " pieces=" << m_polygonResultBuffer.size
                  << std::endl;
    }
}

void ScatteringNonConvex::SortFacets_faster(const Point3f &beamDir,
                                             IntArray &facetIDs)
{
    if (m_tracePrefilterStats)
    {
        ++m_traceSortCalls;
        m_traceSortCandidates += facetIDs.size;
    }
    if (facetIDs.size < 2)
        return;

    struct FacetDepth
    {
        int facetId;
        double depth;
        double maximumDepth;
    };
    FacetDepth ordered[MAX_FACET_NUM];
    for (size_t i = 0; i < facetIDs.size; ++i)
    {
        const int id = facetIDs.arr[i];
        const Polygon &facet = m_facets[id];
        double key = facet.arr[0].cx * beamDir.cx
                   + facet.arr[0].cy * beamDir.cy
                   + facet.arr[0].cz * beamDir.cz;
        double maximumKey = key;

        for (int vertex = 1; vertex < facet.nVertices; ++vertex)
        {
            const Point3f &p = facet.arr[vertex];
            double projection = p.cx * beamDir.cx
                              + p.cy * beamDir.cy
                              + p.cz * beamDir.cz;
            if (projection < key)
                key = projection;
            if (projection > maximumKey)
                maximumKey = projection;
        }
        ordered[i] = {id, key, maximumKey};
    }

    // Candidate lists are short in normal tracing. Insertion sort avoids a
    // heap allocation and has lower overhead than a generic stable sort here.
    // Exact double ordering remains scale invariant.
    for (size_t i = 1; i < facetIDs.size; ++i)
    {
        const FacetDepth value = ordered[i];
        size_t position = i;
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

    // Floating-point rotations perturb theoretically coplanar facets by a few ulps.
    // Give every such group a deterministic order because subtraction of
    // zero-area clipping remnants is otherwise order-dependent. Particle-file
    // facets are canonically indexed by the loader, so this id is independent
    // of their serialization order.
    const double depthTolerance = geometry_length_tolerance(
        m_geometryScale);
    for (size_t begin = 0; begin < facetIDs.size;)
    {
        size_t end = begin + 1;
        while (end < facetIDs.size
               && ordered[end].depth - ordered[begin].depth <= depthTolerance)
            ++end;
        std::sort(ordered + begin, ordered + end,
            [](const FacetDepth &a, const FacetDepth &b) {
                return a.facetId < b.facetId;
            });
        begin = end;
    }

    bool hasOverlappingDepthIntervals = false;
    for (size_t first = 0; first < facetIDs.size && !hasOverlappingDepthIntervals;
         ++first)
    {
        if (first + 1 < facetIDs.size
            && ordered[first + 1].depth
               <= ordered[first].maximumDepth + depthTolerance)
            hasOverlappingDepthIntervals = true;
    }
    if (!hasOverlappingDepthIntervals)
    {
        for (size_t i = 0; i < facetIDs.size; ++i)
            facetIDs.arr[i] = ordered[i].facetId;
        return;
    }
    if (m_tracePrefilterStats)
        ++m_traceSortOverlapCalls;

    Point3f direction(beamDir.cx, beamDir.cy, beamDir.cz);
    const double directionLength = Length(direction);
    if (!(directionLength > DBL_MIN))
    {
        for (size_t i = 0; i < facetIDs.size; ++i)
            facetIDs.arr[i] = ordered[i].facetId;
        return;
    }
    direction = direction/directionLength;

    Point3f reference;
    const double ax = std::fabs(direction.cx);
    const double ay = std::fabs(direction.cy);
    const double az = std::fabs(direction.cz);
    if (ax <= ay && ax <= az)
        reference = Point3f(1, 0, 0);
    else if (ay <= az)
        reference = Point3f(0, 1, 0);
    else
        reference = Point3f(0, 0, 1);
    Point3f axisU = CrossProduct(direction, reference);
    Normalize(axisU);
    Point3f axisV = CrossProduct(direction, axisU);
    Normalize(axisV);

    const size_t count = facetIDs.size;
    std::vector<ProjectedPolygon> projected(count);
    for (size_t index = 0; index < count; ++index)
    {
        const Polygon &facet = m_facets[ordered[index].facetId];
        for (int vertex = 0; vertex < facet.nVertices; ++vertex)
            projected[index].Push(Point2f(
                DotProduct(facet.arr[vertex], axisU),
                DotProduct(facet.arr[vertex], axisV)));
    }

    std::vector<unsigned char> precedes(count*count, 0);
    std::vector<int> indegree(count, 0);
    const double areaTolerance = geometry_area_tolerance(
        m_geometryScale);
    for (size_t first = 0; first < count; ++first)
    {
        for (size_t second = first + 1; second < count; ++second)
        {
            if (ordered[second].depth
                > ordered[first].maximumDepth + depthTolerance)
                break;
            Point2f overlapPoint;
            if (!ProjectedOverlapPoint(projected[first], projected[second],
                                       areaTolerance, overlapPoint))
                continue;

            const Point3f rayOrigin = axisU*overlapPoint.x
                                    + axisV*overlapPoint.y;
            const Polygon &firstFacet = m_facets[ordered[first].facetId];
            const Polygon &secondFacet = m_facets[ordered[second].facetId];
            const Point3f &firstNormal =
                m_facets[ordered[first].facetId].ex_normal;
            const Point3f &secondNormal =
                m_facets[ordered[second].facetId].ex_normal;
            const double firstDenominator = DotProduct(direction, firstNormal);
            const double secondDenominator = DotProduct(direction, secondNormal);
            if (std::fabs(firstDenominator) <= EPS_PROJECTION
                || std::fabs(secondDenominator) <= EPS_PROJECTION)
                continue;
            const double firstDepth = DotProduct(
                firstFacet.arr[0] - rayOrigin, firstNormal)/firstDenominator;
            const double secondDepth = DotProduct(
                secondFacet.arr[0] - rayOrigin, secondNormal)/secondDenominator;
            if (std::fabs(firstDepth - secondDepth) <= depthTolerance)
                continue;

            const size_t before = firstDepth < secondDepth ? first : second;
            const size_t after = firstDepth < secondDepth ? second : first;
            const size_t relation = before*count + after;
            if (!precedes[relation])
            {
                precedes[relation] = 1;
                ++indegree[after];
            }
        }
    }

    std::vector<unsigned char> emitted(count, 0);
    size_t output = 0;
    while (output < count)
    {
        size_t selected = count;
        for (size_t candidate = 0; candidate < count; ++candidate)
        {
            if (!emitted[candidate] && indegree[candidate] == 0)
            {
                selected = candidate;
                break;
            }
        }
        if (selected == count)
        {
            // A cycle can only arise from a numerically unresolved contact.
            // Preserve the deterministic scale-invariant fallback order.
            for (size_t candidate = 0; candidate < count; ++candidate)
                if (!emitted[candidate])
                    facetIDs.arr[output++] = ordered[candidate].facetId;
            break;
        }
        emitted[selected] = 1;
        facetIDs.arr[output++] = ordered[selected].facetId;
        for (size_t after = 0; after < count; ++after)
            if (precedes[selected*count + after])
                --indegree[after];
    }
}

bool ScatteringNonConvex::CutBeamByFacet(const Polygon &intersection,
                                         int facetId, Beam &beam,
                                         Polygon &reachedIntersection,
                                         PolygonArray &result)
{
    reachedIntersection.Clear();
    const Location &loc = beam.location;
    const Facet &beamFacet = m_facets[beam.lastFacetId];
    const double sourceArea = beam.Area();
    Polygon reachedAperture;
    // Build the reached region once on the source plane.  Difference returns
    // both its complement and the exact overlap, so child creation and parent
    // subtraction cannot use different projected polygons.
    Difference(beam, beamFacet.normal[loc], intersection,
               beamFacet.normal[loc], -beam.direction, result,
               &reachedAperture);

    if (reachedAperture.nVertices < MIN_VERTEX_NUM)
    {
        result.Clear();
        return false;
    }

    double remainingArea = 0.0;
    for (unsigned i = 0; i < result.size; ++i)
        remainingArea += result.arr[i].Area();
    const double overlapArea = reachedAperture.Area();
    const double areaImbalance = std::fabs(
        sourceArea - remainingArea - overlapArea);
    const double relativeScale = std::max(
        sourceArea, std::max(remainingArea, overlapArea));
    const double areaTolerance = std::max(
        1024.0*geometry_area_tolerance(m_geometryScale),
        1.0e-9*relativeScale);
    if (areaImbalance > areaTolerance)
    {
        throw std::runtime_error(
            "beam partition does not conserve projected area: imbalance="
            + std::to_string(areaImbalance)
            + ", source=" + std::to_string(sourceArea)
            + ", remaining=" + std::to_string(remainingArea)
            + ", overlap=" + std::to_string(overlapArea)
            + ", tolerance=" + std::to_string(areaTolerance));
    }

    const Facet &targetFacet = m_facets[facetId];
    Point3f targetProjection[MAX_VERTEX_NUM];
    if (!ProjectToFacetPlane(reachedAperture, beam.direction,
                             targetFacet.normal[loc], targetProjection))
    {
        throw std::runtime_error(
            "reached aperture cannot be projected to the target facet");
    }
    SetOutputPolygon(targetProjection, reachedAperture.nVertices,
                     reachedIntersection);
    if (reachedIntersection.nVertices < MIN_VERTEX_NUM)
    {
        // The source partition is still physical, but a nearly tangent target
        // can shrink its child below the representable geometry threshold.
        m_traceCutoffStatistics->smallFragmentSimplifications.fetch_add(
            1, std::memory_order_relaxed);
    }

    // Preserve the historical internal-visibility behavior: the reached
    // child is constrained to the source aperture, but this parent continues
    // to the remaining candidate facets without subtraction.
    if (loc == Location::In && beamFacet.isVisibleIn)
    {
        result.Clear();
        return false;
    }

    return true;
}

bool ScatteringNonConvex::IsOutgoingBeam(Beam &incidentBeam)
{
    return (incidentBeam.location == Location::Out
            && incidentBeam.nVertices != 0); // OPT: replace each other
}

bool ScatteringNonConvex::SplitBeams(std::vector<Beam> &scaterredBeams)
{
    bool ok = true;
    int processedBeams = 0;
    size_t traceCandidateTests = 0;
    size_t traceGpuPrefilterSkipped = 0;
    size_t traceCpuPrefilterSkipped = 0;
    size_t traceIntersectCalls = 0;
    size_t traceGpuSortCalls = 0;
    size_t traceCpuSortFallbacks = 0;
    size_t traceBatchCalls = 0;
    size_t traceBatchItems = 0;
    size_t traceGpuCandidateItems = 0;
    size_t traceGpuCandidates = 0;
    size_t traceGpuSingleCandidateItems = 0;
    size_t traceGpuComplexCandidateItems = 0;
    size_t traceGpuCandidatesLe24 = 0;
    size_t traceGpuCandidates25To32 = 0;
    size_t traceGpuCandidates33To64 = 0;
    size_t traceGpuCandidates65To128 = 0;
    size_t traceGpuCandidates129To256 = 0;
    size_t traceGpuCandidatesGt256 = 0;
    size_t traceGpuBeamVertices3 = 0;
    size_t traceGpuBeamVertices4 = 0;
    size_t traceGpuBeamVertices5To8 = 0;
    size_t traceGpuBeamVertices9To16 = 0;
    size_t traceGpuBeamVerticesGt16 = 0;
    size_t traceGpuFallbackCandidateLimit = 0;
    size_t traceGpuFallbackBeamLimit = 0;
    size_t traceGpuFallbackFacetLimit = 0;
    size_t traceGpuFallbackUnexplained = 0;
#ifdef USE_CUDA
    size_t traceGpuVerifiedSorts = 0;
    size_t traceGpuSortMismatches = 0;
    size_t traceGpuSkipFlagsVerified = 0;
    size_t traceGpuSkipFalseNegatives = 0;
    size_t traceGpuPrefilterGpuOnly = 0;
    size_t traceGpuPrefilterCpuOnly = 0;
    size_t traceSortCacheGlobalHits = 0;
    size_t traceSortCacheBatchHits = 0;
    size_t traceSortCacheCandidatesSaved = 0;
#endif
    double traceSelectTime = 0.0;
    size_t traceSelectTimingCalls = 0;
    size_t traceSelectTimingSamples = 0;
    double traceAabbTime = 0.0;
    size_t traceAabbTimingCalls = 0;
    size_t traceAabbTimingSamples = 0;
    double traceIntersectTime = 0.0;
    size_t traceIntersectTimingCalls = 0;
    size_t traceIntersectTimingSamples = 0;
    double traceSplitTime = 0.0;
    size_t traceSplitTimingCalls = 0;
    size_t traceSplitTimingSamples = 0;
    double traceOutgoingTime = 0.0;
    size_t traceOutgoingTimingCalls = 0;
    size_t traceOutgoingTimingSamples = 0;
    double traceBatchCollectTime = 0.0;
    double traceVisibleTime = 0.0;
    size_t traceVisibleTimingCalls = 0;
    size_t traceVisibleTimingSamples = 0;
    double traceGpuPrepareTime = 0.0;
    double traceFallbackSortTime = 0.0;
    double tracePreparedProcessTime = 0.0;
    const auto traceNow = []() {
        return std::chrono::high_resolution_clock::now();
    };
    const auto traceSeconds = [](const std::chrono::high_resolution_clock::time_point &a,
                                 const std::chrono::high_resolution_clock::time_point &b) {
        return std::chrono::duration<double>(b - a).count();
    };
    const auto traceTimingSample = [](size_t &calls, size_t &samples) {
        const bool sample = (calls++ & 63u) == 0;
        if (sample)
            ++samples;
        return sample;
    };
    auto processBeam = [this, &scaterredBeams, &ok,
                        &traceCandidateTests,
                        &traceGpuPrefilterSkipped,
                        &traceCpuPrefilterSkipped,
                        &traceIntersectCalls,
                        &traceSelectTime,
                        &traceSelectTimingCalls,
                        &traceSelectTimingSamples,
                        &traceAabbTime,
                        &traceAabbTimingCalls,
                        &traceAabbTimingSamples,
                        &traceIntersectTime,
                        &traceIntersectTimingCalls,
                        &traceIntersectTimingSamples,
                        &traceSplitTime,
                        &traceSplitTimingCalls,
                        &traceSplitTimingSamples,
                        &traceOutgoingTime,
                        &traceOutgoingTimingCalls,
                        &traceOutgoingTimingSamples,
                        &traceNow,
                        &traceSeconds,
                        &traceTimingSample](Beam &beam,
                                      const IntArray *readyFacetIds,
                                      const std::vector<unsigned char> *mayIntersect,
                                      bool gpuPrefiltered) -> bool
    {
        if (!IsTerminalAct(beam)) // REF, OPT: перенести проверку во все места, где пучок закидывается в дерево, чтобы пучки заранее не закидывались в него
        {
            IntArray localFacetIds;
            const IntArray *facetIds = readyFacetIds;
            if (facetIds == nullptr)
            {
                if (m_tracePrefilterStats
                    && traceTimingSample(traceSelectTimingCalls,
                                         traceSelectTimingSamples))
                {
                    auto t0 = traceNow();
                    SelectVisibleFacets(beam, localFacetIds);
                    traceSelectTime += traceSeconds(t0, traceNow());
                }
                else
                {
                    SelectVisibleFacets(beam, localFacetIds);
                }
                facetIds = &localFacetIds;
            }

            bool isDivided = false;

            for (unsigned i = 0; (i < facetIds->size) && !isDivided; ++i)// OPT: move this loop to SplitBeamByFacet
            {
                if (m_tracePrefilterStats)
                    ++traceCandidateTests;
                if (mayIntersect != nullptr && !mayIntersect->empty() && !(*mayIntersect)[i])
                {
                    if (m_tracePrefilterStats)
                        ++traceGpuPrefilterSkipped;
                    continue;
                }
                int facetId = facetIds->arr[i];
                if (m_traceCpuProjectedPrefilter && !gpuPrefiltered)
                {
                    bool mayIntersectProjected;
                    if (m_tracePrefilterStats)
                    {
                        const bool sampleTiming =
                            (traceAabbTimingCalls++ & 63u) == 0;
                        if (sampleTiming)
                        {
                            auto t0 = traceNow();
                            mayIntersectProjected =
                                MayBeamIntersectFacetProjected(beam, facetId);
                            traceAabbTime += traceSeconds(t0, traceNow());
                            ++traceAabbTimingSamples;
                        }
                        else
                        {
                            mayIntersectProjected =
                                MayBeamIntersectFacetProjected(beam, facetId);
                        }
                    }
                    else
                    {
                        mayIntersectProjected = MayBeamIntersectFacetProjected(beam, facetId);
                    }
                    if (!mayIntersectProjected)
                    {
                        if (m_tracePrefilterStats)
                            ++traceCpuPrefilterSkipped;
                        continue;
                    }
                }

                Polygon intersection;
                if (m_tracePrefilterStats)
                    ++traceIntersectCalls;
                if (m_tracePrefilterStats
                    && traceTimingSample(traceIntersectTimingCalls,
                                         traceIntersectTimingSamples))
                {
                    auto tIntersect0 = traceNow();
                    Intersect(facetId, beam, intersection);
                    traceIntersectTime += traceSeconds(tIntersect0, traceNow());
                }
                else
                {
                    Intersect(facetId, beam, intersection);
                }

                if (intersection.nVertices >= MIN_VERTEX_NUM)
                {
                    if (m_tracePrefilterStats
                        && traceTimingSample(traceSplitTimingCalls,
                                             traceSplitTimingSamples))
                    {
                        auto tSplit0 = traceNow();
                        isDivided = SplitBeamByFacet(intersection, facetId, beam, ok);
                        traceSplitTime += traceSeconds(tSplit0, traceNow());
                    }
                    else
                    {
                        isDivided = SplitBeamByFacet(intersection, facetId, beam, ok);
                    }

                    if (!ok)
                    {
                        return false;
                    }
                    if (!isDivided && beam.location == Location::In && IsTracePruned(beam))
                    {
                        beam.nVertices = 0;
                        break;
                    }
                }
            }

            if (IsOutgoingBeam(beam))
            {	// посылаем обрезанный всеми гранями внешний пучок на сферу
                double path;
                if (m_tracePrefilterStats
                    && traceTimingSample(traceOutgoingTimingCalls,
                                         traceOutgoingTimingSamples))
                {
                    auto tOut0 = traceNow();
                    path = splitting.ComputeOutgoingOpticalPath(beam); // добираем оптический путь
                    traceOutgoingTime += traceSeconds(tOut0, traceNow());
                }
                else
                {
                    path = splitting.ComputeOutgoingOpticalPath(beam); // добираем оптический путь
                }
                beam.opticalPath += path;
#ifdef _DEBUG // DEB
    beam.ops.push_back(path);
    if (beam.lastFacetId == 6)
        int fff = 0;
#endif
                scaterredBeams.push_back(beam);
            }
        }
        else if (beam.location == Location::Out)
        {
            if (m_tracePrefilterStats
                && traceTimingSample(traceOutgoingTimingCalls,
                                     traceOutgoingTimingSamples))
            {
                auto tOut0 = traceNow();
                CutExternalBeam(beam, scaterredBeams, readyFacetIds);
                traceOutgoingTime += traceSeconds(tOut0, traceNow());
            }
            else
            {
                CutExternalBeam(beam, scaterredBeams, readyFacetIds);
            }
        }
        return true;
    };
#ifdef USE_CUDA
    struct TraceSortKey
    {
        int lastFacetId;
        int location;
        uint64_t direction[3];
        uint64_t candidateHash;
        size_t candidateCount;

        bool operator<(const TraceSortKey &other) const
        {
            if (lastFacetId != other.lastFacetId)
                return lastFacetId < other.lastFacetId;
            if (location != other.location)
                return location < other.location;
            for (int coordinate = 0; coordinate < 3; ++coordinate)
            {
                if (direction[coordinate] != other.direction[coordinate])
                    return direction[coordinate]
                         < other.direction[coordinate];
            }
            if (candidateHash != other.candidateHash)
                return candidateHash < other.candidateHash;
            return candidateCount < other.candidateCount;
        }

    };
    struct TraceSortCacheEntry
    {
        std::vector<int> original;
        std::vector<int> sorted;
    };
    struct TraceSortRepresentatives
    {
        size_t first = static_cast<size_t>(-1);
        std::vector<size_t> collisions;
    };
    struct TraceBatchItem
    {
        Beam beam;
        IntArray facetIds;
        std::vector<unsigned char> mayIntersect;
        bool hasFacetIds = false;
        bool sortedOnGpu = false;
        bool sortedFromCache = false;
        bool prefilteredOnGpu = false;
        size_t originalFacetCount = 0;
        size_t sortRepresentative = static_cast<size_t>(-1);
    };
    const bool traceSortCacheEnabled = GpuTraceSortCacheEnabled();
    const auto makeTraceSortKey = [](const Beam &beam,
                                     const IntArray &facetIds) {
        TraceSortKey key;
        key.lastFacetId = beam.lastFacetId;
        key.location = beam.location == Location::In ? 0 : 1;
        for (int coordinate = 0; coordinate < 3; ++coordinate)
        {
            std::memcpy(&key.direction[coordinate],
                        &beam.direction.coordinates[coordinate],
                        sizeof(uint64_t));
        }
        uint64_t hash = 1469598103934665603ULL;
        for (size_t index = 0; index < facetIds.size; ++index)
        {
            hash ^= static_cast<uint32_t>(facetIds.arr[index]);
            hash *= 1099511628211ULL;
        }
        key.candidateHash = hash;
        key.candidateCount = facetIds.size;
        return key;
    };
    const auto equalsCandidateList = [](const IntArray &facetIds,
                                        const std::vector<int> &values) {
        if (facetIds.size != values.size())
            return false;
        for (size_t index = 0; index < facetIds.size; ++index)
            if (facetIds.arr[index] != values[index])
                return false;
        return true;
    };
    const auto assignCandidateList = [](IntArray &facetIds,
                                        const std::vector<int> &values) {
        facetIds.size = values.size();
        for (size_t index = 0; index < values.size(); ++index)
            facetIds.arr[index] = values[index];
    };
    std::vector<TraceBatchItem> traceBatch;
    std::vector<GpuTraceBeamFacets> traceGpuItems;
    std::vector<size_t> traceGpuBatchIndices;
    std::vector<IntArray> traceGpuOriginalFacetIds;
    std::map<TraceSortKey, std::vector<TraceSortCacheEntry>> traceSortCache;
    std::map<TraceSortKey, TraceSortRepresentatives>
        traceBatchSortRepresentatives;
#endif
//#ifdef _DEBUG // DEB
//    ofstream logfile("logscat.txt", ios::out);
//    int count = 0;
//#endif
    while (m_treeSize != 0)
    {
#ifdef USE_CUDA
        if (m_gpuTracePrefilter && GpuTraceRuntimeEnabled())
        {
            std::chrono::high_resolution_clock::time_point batchCollectStarted;
            if (m_tracePrefilterStats)
                batchCollectStarted = traceNow();
            const size_t batchLimit = (size_t)GpuTraceBatchBeamLimit();
            if (traceBatch.empty())
            {
                traceBatch.resize(batchLimit);
                traceGpuItems.reserve(batchLimit);
            }
            size_t batchSize = 0;
            while (m_treeSize != 0 && batchSize < batchLimit)
            {
                if (m_traceMaxBeams > 0 && ++processedBeams > m_traceMaxBeams)
                {
                    m_traceCutoffStatistics->configuredBeamLimitHits.fetch_add(
                        1, std::memory_order_relaxed);
                    throw TraceLimitExceeded(
                        "tracing processed more than --trace-max-beams="
                        + std::to_string(m_traceMaxBeams)
                        + " branches; raise the limit or use 0 to disable it");
                }
                TraceBatchItem &item = traceBatch[batchSize++];
                item.facetIds.size = 0;
                item.mayIntersect.clear();
                item.hasFacetIds = false;
                item.sortedOnGpu = false;
                item.sortedFromCache = false;
                item.prefilteredOnGpu = false;
                item.originalFacetCount = 0;
                item.sortRepresentative = static_cast<size_t>(-1);
                item.beam = std::move(m_beamTree.back());
                m_beamTree.pop_back();
                m_treeSize = (int)m_beamTree.size();
                if (!IsTerminalAct(item.beam)
                    || item.beam.location == Location::Out)
                {
                    if (m_tracePrefilterStats
                        && traceTimingSample(traceVisibleTimingCalls,
                                             traceVisibleTimingSamples))
                    {
                        const auto visibleStarted = traceNow();
                        FindVisibleFacets(item.beam, item.facetIds);
                        traceVisibleTime += traceSeconds(
                            visibleStarted, traceNow());
                    }
                    else
                    {
                        FindVisibleFacets(item.beam, item.facetIds);
                    }
                    item.hasFacetIds = true;
                    item.originalFacetCount = item.facetIds.size;
                }
            }
            if (m_tracePrefilterStats)
            {
                ++traceBatchCalls;
                traceBatchItems += batchSize;
                traceBatchCollectTime += traceSeconds(
                    batchCollectStarted, traceNow());
            }

            traceGpuItems.clear();
            traceGpuBatchIndices.clear();
            traceGpuOriginalFacetIds.clear();
            traceBatchSortRepresentatives.clear();
            for (size_t i = 0; i < batchSize; ++i)
            {
                TraceBatchItem &batchItem = traceBatch[i];
                if (batchItem.hasFacetIds && batchItem.facetIds.size != 0)
                {
                    if (m_tracePrefilterStats)
                    {
                        ++traceGpuCandidateItems;
                        traceGpuCandidates += batchItem.facetIds.size;
                        if (batchItem.facetIds.size == 1)
                            ++traceGpuSingleCandidateItems;
                        else
                            ++traceGpuComplexCandidateItems;
                        if (batchItem.facetIds.size <= 24)
                            ++traceGpuCandidatesLe24;
                        else if (batchItem.facetIds.size <= 32)
                            ++traceGpuCandidates25To32;
                        else if (batchItem.facetIds.size <= 64)
                            ++traceGpuCandidates33To64;
                        else if (batchItem.facetIds.size <= 128)
                            ++traceGpuCandidates65To128;
                        else if (batchItem.facetIds.size <= 256)
                            ++traceGpuCandidates129To256;
                        else
                            ++traceGpuCandidatesGt256;
                        if (batchItem.beam.nVertices == 3)
                            ++traceGpuBeamVertices3;
                        else if (batchItem.beam.nVertices == 4)
                            ++traceGpuBeamVertices4;
                        else if (batchItem.beam.nVertices <= 8)
                            ++traceGpuBeamVertices5To8;
                        else if (batchItem.beam.nVertices <= 16)
                            ++traceGpuBeamVertices9To16;
                        else
                            ++traceGpuBeamVerticesGt16;
                    }
                    bool sortHandled = false;
                    if (traceSortCacheEnabled && batchItem.facetIds.size > 24)
                    {
                        const TraceSortKey key = makeTraceSortKey(
                            batchItem.beam, batchItem.facetIds);
                        const auto cached = traceSortCache.find(key);
                        if (cached != traceSortCache.end())
                        {
                            for (const TraceSortCacheEntry &entry
                                 : cached->second)
                            {
                                if (!equalsCandidateList(batchItem.facetIds,
                                                         entry.original))
                                    continue;
                                assignCandidateList(batchItem.facetIds,
                                                    entry.sorted);
                                batchItem.sortedOnGpu = true;
                                batchItem.sortedFromCache = true;
                                ++traceSortCacheGlobalHits;
                                traceSortCacheCandidatesSaved +=
                                    batchItem.originalFacetCount;
                                sortHandled = true;
                                break;
                            }
                        }
                        if (!sortHandled)
                        {
                            TraceSortRepresentatives &representatives =
                                traceBatchSortRepresentatives[key];
                            if (representatives.first
                                != static_cast<size_t>(-1))
                            {
                                const size_t representative =
                                    representatives.first;
                                if (batchItem.facetIds.size
                                    == traceBatch[representative].facetIds.size)
                                {
                                    bool same = true;
                                    for (size_t index = 0;
                                         index < batchItem.facetIds.size;
                                         ++index)
                                    {
                                        if (batchItem.facetIds.arr[index]
                                            != traceBatch[representative]
                                                   .facetIds.arr[index])
                                        {
                                            same = false;
                                            break;
                                        }
                                    }
                                    if (same)
                                    {
                                        batchItem.sortRepresentative =
                                            representative;
                                        sortHandled = true;
                                    }
                                }
                            }
                            for (size_t collisionIndex = 0;
                                 !sortHandled
                                 && collisionIndex
                                    < representatives.collisions.size();
                                 ++collisionIndex)
                            {
                                const size_t representative =
                                    representatives.collisions[collisionIndex];
                                if (batchItem.facetIds.size
                                    != traceBatch[representative].facetIds.size)
                                    continue;
                                bool same = true;
                                for (size_t index = 0;
                                     index < batchItem.facetIds.size; ++index)
                                {
                                    if (batchItem.facetIds.arr[index]
                                        != traceBatch[representative]
                                               .facetIds.arr[index])
                                    {
                                        same = false;
                                        break;
                                    }
                                }
                                if (same)
                                {
                                    batchItem.sortRepresentative =
                                        representative;
                                    sortHandled = true;
                                }
                            }
                            if (!sortHandled)
                            {
                                if (representatives.first
                                    == static_cast<size_t>(-1))
                                    representatives.first = i;
                                else
                                    representatives.collisions.push_back(i);
                            }
                        }
                    }
                    if (sortHandled)
                        continue;

                    GpuTraceBeamFacets item;
                    item.beam = &batchItem.beam;
                    item.facetIds = &batchItem.facetIds;
                    item.mayIntersect = &batchItem.mayIntersect;
                    item.sortedOnGpu = &batchItem.sortedOnGpu;
                    traceGpuItems.push_back(item);
                    traceGpuBatchIndices.push_back(i);
                    if (traceSortCacheEnabled || GpuTraceVerifySort())
                        traceGpuOriginalFacetIds.push_back(
                            batchItem.facetIds);
                }
            }

            if (!traceGpuItems.empty())
            {
                if (m_tracePrefilterStats)
                {
                    const auto gpuPrepareStarted = traceNow();
                    GpuTracePrepareBeamFacetBatchAdaptive(
                        m_facets, m_geometryScale, traceGpuItems);
                    traceGpuPrepareTime += traceSeconds(
                        gpuPrepareStarted, traceNow());
                }
                else
                {
                    GpuTracePrepareBeamFacetBatchAdaptive(
                        m_facets, m_geometryScale, traceGpuItems);
                }
            }

            if (traceSortCacheEnabled)
            {
                for (size_t gpuIndex = 0;
                     gpuIndex < traceGpuItems.size(); ++gpuIndex)
                {
                    TraceBatchItem &batchItem =
                        traceBatch[traceGpuBatchIndices[gpuIndex]];
                    if (!batchItem.sortedOnGpu
                        || batchItem.originalFacetCount <= 24)
                        continue;
                    const IntArray &original =
                        traceGpuOriginalFacetIds[gpuIndex];
                    TraceSortCacheEntry entry;
                    entry.original.assign(original.arr,
                                          original.arr + original.size);
                    entry.sorted.assign(batchItem.facetIds.arr,
                                        batchItem.facetIds.arr
                                            + batchItem.facetIds.size);
                    traceSortCache[makeTraceSortKey(batchItem.beam, original)]
                        .push_back(std::move(entry));
                }
                for (size_t i = 0; i < batchSize; ++i)
                {
                    TraceBatchItem &batchItem = traceBatch[i];
                    if (batchItem.sortRepresentative
                        == static_cast<size_t>(-1))
                        continue;
                    const TraceBatchItem &representative =
                        traceBatch[batchItem.sortRepresentative];
                    if (!representative.sortedOnGpu)
                        continue;
                    batchItem.facetIds = representative.facetIds;
                    batchItem.sortedOnGpu = true;
                    batchItem.sortedFromCache = true;
                    ++traceSortCacheBatchHits;
                    traceSortCacheCandidatesSaved +=
                        batchItem.originalFacetCount;
                }
            }

            if (GpuTraceVerifySort())
            {
                for (size_t itemIndex = 0;
                     itemIndex < traceGpuItems.size(); ++itemIndex)
                {
                    GpuTraceBeamFacets &gpuItem = traceGpuItems[itemIndex];
                    if (gpuItem.sortedOnGpu == nullptr
                        || !*gpuItem.sortedOnGpu)
                        continue;
                    ++traceGpuVerifiedSorts;

                    if (gpuItem.mayIntersect != nullptr
                        && gpuItem.mayIntersect->size()
                           == gpuItem.facetIds->size)
                    {
                        for (size_t candidate = 0;
                             candidate < gpuItem.facetIds->size; ++candidate)
                        {
                            if ((*gpuItem.mayIntersect)[candidate] != 0)
                                continue;
                            ++traceGpuSkipFlagsVerified;
                            const int facetId =
                                gpuItem.facetIds->arr[candidate];
                            if (MayBeamIntersectFacetProjected(
                                    *gpuItem.beam, facetId))
                            {
                                ++traceGpuSkipFalseNegatives;
                                if (traceGpuSkipFalseNegatives <= 8)
                                {
                                    std::cout
                                        << "Trace GPU skip false negative: facet="
                                        << facetId << " candidates="
                                        << gpuItem.facetIds->size << std::endl;
                                }
                            }
                        }
                    }

                    IntArray expected = traceGpuOriginalFacetIds[itemIndex];
                    SortVisibleFacets(*gpuItem.beam, expected);
                    bool retained[MAX_FACET_NUM] = {};
                    for (size_t index = 0; index < gpuItem.facetIds->size;
                         ++index)
                    {
                        const int facetId = gpuItem.facetIds->arr[index];
                        if (facetId >= 0 && facetId < MAX_FACET_NUM)
                            retained[facetId] = true;
                    }
                    const bool gpuPrefiltered =
                        !GpuTraceFullSortEnabled()
                        || traceGpuOriginalFacetIds[itemIndex].size <= 24;
                    for (size_t index = 0;
                         gpuPrefiltered
                         && index < traceGpuOriginalFacetIds[itemIndex].size;
                         ++index)
                    {
                        const int facetId =
                            traceGpuOriginalFacetIds[itemIndex].arr[index];
                        const bool gpuKeeps = facetId >= 0
                            && facetId < MAX_FACET_NUM
                            && retained[facetId];
                        const bool cpuKeeps =
                            MayBeamIntersectFacetProjected(
                                *gpuItem.beam, facetId);
                        if (gpuKeeps && !cpuKeeps)
                        {
                            ++traceGpuPrefilterGpuOnly;
                            if (traceGpuPrefilterGpuOnly <= 4)
                            {
                                std::cout
                                    << "Trace GPU prefilter GPU-only: facet="
                                    << facetId << " original="
                                    << traceGpuOriginalFacetIds[itemIndex].size
                                    << std::endl;
                            }
                        }
                        else if (!gpuKeeps && cpuKeeps)
                        {
                            ++traceGpuPrefilterCpuOnly;
                            if (traceGpuPrefilterCpuOnly <= 4)
                            {
                                std::cout
                                    << "Trace GPU prefilter CPU-only: facet="
                                    << facetId << " original="
                                    << traceGpuOriginalFacetIds[itemIndex].size
                                    << std::endl;
                            }
                        }
                    }
                    IntArray filteredExpected;
                    for (size_t index = 0; index < expected.size; ++index)
                    {
                        const int facetId = expected.arr[index];
                        if (facetId >= 0 && facetId < MAX_FACET_NUM
                            && retained[facetId])
                            filteredExpected.Add(facetId);
                    }

                    bool same = filteredExpected.size
                             == gpuItem.facetIds->size;
                    size_t firstMismatch = 0;
                    while (same && firstMismatch < filteredExpected.size)
                    {
                        if (filteredExpected.arr[firstMismatch]
                            != gpuItem.facetIds->arr[firstMismatch])
                            same = false;
                        else
                            ++firstMismatch;
                    }
                    if (!same)
                    {
                        ++traceGpuSortMismatches;
                        if (traceGpuSortMismatches <= 8)
                        {
                            std::cout
                                << "Trace GPU sort mismatch: original="
                                << traceGpuOriginalFacetIds[itemIndex].size
                                << " kept=" << gpuItem.facetIds->size
                                << " first=" << firstMismatch
                                << " cpu="
                                << (firstMismatch < filteredExpected.size
                                    ? filteredExpected.arr[firstMismatch] : -1)
                                << " gpu="
                                << (firstMismatch < gpuItem.facetIds->size
                                    ? gpuItem.facetIds->arr[firstMismatch] : -1)
                                << std::endl;
                        }
                    }
                }
            }

            std::chrono::high_resolution_clock::time_point preparedProcessStarted;
            if (m_tracePrefilterStats)
                preparedProcessStarted = traceNow();
            for (size_t i = 0; i < batchSize; ++i)
            {
                TraceBatchItem &batchItem = traceBatch[i];
                batchItem.prefilteredOnGpu = batchItem.sortedOnGpu
                    && (!GpuTraceFullSortEnabled()
                        || batchItem.originalFacetCount <= 24);
                if (batchItem.hasFacetIds && batchItem.facetIds.size != 0)
                {
                    if (batchItem.sortedOnGpu)
                    {
                        if (!batchItem.sortedFromCache)
                            ++traceGpuSortCalls;
                    }
                    else
                    {
                        ++traceCpuSortFallbacks;
                        if (m_tracePrefilterStats)
                        {
                            if (batchItem.facetIds.size > 24)
                            {
                                ++traceGpuFallbackCandidateLimit;
                            }
                            else if (batchItem.beam.nVertices <= 0
                                     || batchItem.beam.nVertices > 8)
                            {
                                ++traceGpuFallbackBeamLimit;
                            }
                            else
                            {
                                bool oversizedFacet = false;
                                for (size_t facetIndex = 0;
                                     facetIndex < batchItem.facetIds.size;
                                     ++facetIndex)
                                {
                                    const int facetId =
                                        batchItem.facetIds.arr[facetIndex];
                                    if (m_facets[facetId].nVertices < 3
                                        || m_facets[facetId].nVertices > 8)
                                    {
                                        oversizedFacet = true;
                                        break;
                                    }
                                }
                                if (oversizedFacet)
                                    ++traceGpuFallbackFacetLimit;
                                else
                                    ++traceGpuFallbackUnexplained;
                            }
                            const auto fallbackStarted = traceNow();
                            SortVisibleFacets(batchItem.beam, batchItem.facetIds);
                            traceFallbackSortTime += traceSeconds(
                                fallbackStarted, traceNow());
                        }
                        else
                        {
                            SortVisibleFacets(batchItem.beam, batchItem.facetIds);
                        }
                    }
                }
                if (!processBeam(batchItem.beam,
                                 batchItem.hasFacetIds ? &batchItem.facetIds : nullptr,
                                 batchItem.mayIntersect.empty() ? nullptr : &batchItem.mayIntersect,
                                 batchItem.prefilteredOnGpu))
                    return false;
            }
            if (m_tracePrefilterStats)
                tracePreparedProcessTime += traceSeconds(
                    preparedProcessStarted, traceNow());
            continue;
        }
#endif
        if (m_traceMaxBeams > 0 && ++processedBeams > m_traceMaxBeams)
        {
            m_traceCutoffStatistics->configuredBeamLimitHits.fetch_add(
                1, std::memory_order_relaxed);
            throw TraceLimitExceeded(
                "tracing processed more than --trace-max-beams="
                + std::to_string(m_traceMaxBeams)
                + " branches; raise the limit or use 0 to disable it");
        }
        Beam beam = m_beamTree.back();
        m_beamTree.pop_back();
        m_treeSize = (int)m_beamTree.size();
#ifdef _DEBUG // DEB
//        logfile << count << " " << m_treeSize << " " << beam.id << std::endl;
//        logfile.flush();
//    ++count;
        if (beam.id == 171)
            int ffgf = 0;
#endif
        if (!processBeam(beam, nullptr, nullptr, false))
            return false;
    }

    if (m_tracePrefilterStats)
    {
        const auto traceEstimatedTime = [](double sampledTime,
                                           size_t calls,
                                           size_t samples) {
            return samples == 0 ? 0.0 : sampledTime
                * static_cast<double>(calls)
                / static_cast<double>(samples);
        };
        const double traceSelectEstimatedTime = traceEstimatedTime(
            traceSelectTime, traceSelectTimingCalls, traceSelectTimingSamples);
        const double traceAabbEstimatedTime = traceEstimatedTime(
            traceAabbTime, traceAabbTimingCalls, traceAabbTimingSamples);
        const double traceIntersectEstimatedTime = traceEstimatedTime(
            traceIntersectTime,
            traceIntersectTimingCalls,
            traceIntersectTimingSamples);
        const double traceSplitEstimatedTime = traceEstimatedTime(
            traceSplitTime, traceSplitTimingCalls, traceSplitTimingSamples);
        const double traceOutgoingEstimatedTime = traceEstimatedTime(
            traceOutgoingTime,
            traceOutgoingTimingCalls,
            traceOutgoingTimingSamples);
        const double traceVisibleEstimatedTime = traceEstimatedTime(
            traceVisibleTime, traceVisibleTimingCalls, traceVisibleTimingSamples);
        std::cout << "Trace prefilter stats: candidates=" << traceCandidateTests
                  << " gpu_skip=" << traceGpuPrefilterSkipped
                  << " cpu_aabb_skip=" << traceCpuPrefilterSkipped
                  << " intersect_calls=" << traceIntersectCalls
                  << std::endl;
        std::cout << "Trace timing stats: select_sort=" << traceSelectEstimatedTime
                  << " aabb=" << traceAabbEstimatedTime
                  << " intersect=" << traceIntersectEstimatedTime
                  << " split=" << traceSplitEstimatedTime
                  << " outgoing=" << traceOutgoingEstimatedTime
                  << std::endl;
        std::cout << "Trace timing samples: select_sort="
                  << traceSelectTimingSamples << "/" << traceSelectTimingCalls
                  << " aabb=" << traceAabbTimingSamples
                  << "/" << traceAabbTimingCalls
                  << " intersect=" << traceIntersectTimingSamples
                  << "/" << traceIntersectTimingCalls
                  << " split=" << traceSplitTimingSamples
                  << "/" << traceSplitTimingCalls
                  << " outgoing=" << traceOutgoingTimingSamples
                  << "/" << traceOutgoingTimingCalls
                  << " visible=" << traceVisibleTimingSamples
                  << "/" << traceVisibleTimingCalls
                  << std::endl;
        std::cout << "Trace sort stats: calls=" << m_traceSortCalls
                  << " candidates=" << m_traceSortCandidates
                  << " overlap_calls=" << m_traceSortOverlapCalls
                  << " gpu_calls=" << traceGpuSortCalls
                  << " cpu_fallbacks=" << traceCpuSortFallbacks
                  << std::endl;
        std::cout << "Trace batch stats: batches=" << traceBatchCalls
                  << " items=" << traceBatchItems
                  << " collect=" << traceBatchCollectTime
                  << " visible=" << traceVisibleEstimatedTime
                  << " gpu_prepare=" << traceGpuPrepareTime
                  << " fallback_sort=" << traceFallbackSortTime
                  << " prepared_process=" << tracePreparedProcessTime
                  << std::endl;
        std::cout << "Trace GPU item stats: nonempty="
                  << traceGpuCandidateItems
                  << " candidates=" << traceGpuCandidates
                  << " single=" << traceGpuSingleCandidateItems
                  << " complex=" << traceGpuComplexCandidateItems
                  << std::endl;
#ifdef USE_CUDA
        if (traceSortCacheEnabled)
        {
            std::cout << "Trace sort cache stats: entries="
                      << traceSortCache.size()
                      << " global_hits=" << traceSortCacheGlobalHits
                      << " batch_hits=" << traceSortCacheBatchHits
                      << " candidates_saved="
                      << traceSortCacheCandidatesSaved
                      << std::endl;
        }
#endif
        std::cout << "Trace GPU eligibility stats: candidates<=24="
                  << traceGpuCandidatesLe24
                  << " candidates25-32=" << traceGpuCandidates25To32
                  << " candidates33-64=" << traceGpuCandidates33To64
                  << " candidates65-128=" << traceGpuCandidates65To128
                  << " candidates129-256=" << traceGpuCandidates129To256
                  << " candidates>256=" << traceGpuCandidatesGt256
                  << " beam_vertices=3=" << traceGpuBeamVertices3
                  << " beam_vertices=4=" << traceGpuBeamVertices4
                  << " beam_vertices5-8=" << traceGpuBeamVertices5To8
                  << " beam_vertices9-16=" << traceGpuBeamVertices9To16
                  << " beam_vertices>16=" << traceGpuBeamVerticesGt16
                  << std::endl;
        std::cout << "Trace GPU fallback reasons: candidate_limit="
                  << traceGpuFallbackCandidateLimit
                  << " beam_limit=" << traceGpuFallbackBeamLimit
                  << " facet_limit=" << traceGpuFallbackFacetLimit
                  << " unexplained=" << traceGpuFallbackUnexplained
                  << std::endl;
#ifdef USE_CUDA
        if (GpuTraceVerifySort())
        {
            std::cout << "Trace GPU sort verification: compared="
                      << traceGpuVerifiedSorts
                      << " mismatches=" << traceGpuSortMismatches
                      << " prefilter_gpu_only="
                      << traceGpuPrefilterGpuOnly
                      << " prefilter_cpu_only="
                      << traceGpuPrefilterCpuOnly
                      << " skip_flags_verified="
                      << traceGpuSkipFlagsVerified
                      << " skip_false_negatives="
                      << traceGpuSkipFalseNegatives
                      << std::endl;
        }
#endif
    }

    return true;
}

bool ScatteringNonConvex::SetOpticalBeamParams(const Facet &facet, const Beam &incidentBeam,
                                               Beam &inBeam, Beam &outBeam)
{
    const Point3f &dir = incidentBeam.direction;
    const Point3f &normal = facet.ex_normal;

    bool hasOutBeam = true;
    splitting.ComputeCosA(dir, normal);
    splitting.ComputeSplittingParams(incidentBeam.direction, normal);

    if (splitting.IsNormalIncidence()) // normal incidence
    {
        splitting.ComputeNormalBeamParams(incidentBeam, inBeam, outBeam);
    }
    else // regular incidence
    {
        // The splitting routines need the incident optical state, not its
        // aperture polygon.  Copying the complete Beam here duplicated the
        // polygon for every facet hit in the hottest tracing path.
        Beam incBeam;
        incBeam.SetLight(incidentBeam);
        incBeam.J = incidentBeam.J;
        incBeam.front = incidentBeam.front;
        incBeam.location = incidentBeam.location;
        incBeam.opticalPath = incidentBeam.opticalPath;
#ifdef _DEBUG // DEB
        incBeam.ops = incidentBeam.ops;
#endif

        if (incidentBeam.location == Location::In)
        {
            ComputePolarisationParams(incBeam.direction, normal, incBeam);
//			incBeam.direction = -incBeam.direction;

            hasOutBeam = !splitting.IsCompleteReflection();

            if (hasOutBeam)
            {
                splitting.ComputeRegularBeamsParams(normal, incBeam,
                                                      inBeam, outBeam);
            }
            else // complete internal reflection incidence
            {
                splitting.ComputeCRBeamParams(normal, incBeam, inBeam);
            }
        }
        else // beam is external
        {
            inBeam.J = incidentBeam.J;
            splitting.ComputeCosA(dir, -normal);

            const Point3f &facetNormal = facet.in_normal;
            ComputePolarisationParams(-incBeam.direction, facetNormal, incBeam);
            splitting.ComputeRegularBeamParamsExternal(facetNormal, incBeam,
                                                         inBeam, outBeam);
//			if (m_isOpticalPath)
            {
                double path = splitting.ComputeSegmentOpticalPath(incidentBeam,
                                                                    inBeam.Center());
#ifdef _DEBUG // DEB
    inBeam.ops = incidentBeam.ops;
    outBeam.ops = incidentBeam.ops;
    inBeam.ops.push_back(path);
    outBeam.ops.push_back(path);
#endif
                path += incidentBeam.opticalPath;
                inBeam.AddOpticalPath(path);
                outBeam.AddOpticalPath(path);
            }
        }
    }

    return hasOutBeam;
}

void ScatteringNonConvex::FindVisibleFacetsForLight(IntArray &facetIDs)
{
    for (int i = 0; i < m_particle->nFacets; ++i)
    {
        double cosA = DotProduct(m_incidentDir, m_facets[i].in_normal);

        if (cosA > EPS_PROJECTION) // beam incidents to this facet
        {
            facetIDs.Add(i);
        }
    }
}

bool ScatteringNonConvex::IsVisibleFacet(int facetID, const Beam &beam)
{
    const Point3f &beamNormal = -m_facets[beam.lastFacetId].normal[!beam.location];
    const Point3f &beamCenter = m_facets[beam.lastFacetId].center;
    const double tolerance = geometry_length_tolerance(
        m_geometryScale);
    return FacetReachesForwardHalfSpace(m_facets[facetID], beamCenter,
                                        beamNormal, tolerance);
}

bool ScatteringNonConvex::MayBeamIntersectFacetProjected(const Beam &beam, int facetId) const
{
    if (beam.nVertices <= 0)
        return false;

    const Point3f &normal = (beam.location == Location::In)
                          ? m_facets[facetId].in_normal
                          : m_facets[facetId].ex_normal;
    const double dp = DotProduct(beam.direction, normal);
    if (std::fabs(dp) < EPS_PROJECTION)
        return false;
    const double invDp = 1.0 / dp;

    const int locInt = beam.location == Location::In ? 0 : 1;
    const int drop = m_facetProjectionDrop[locInt][facetId];

    double bMinU = DBL_MAX, bMaxU = -DBL_MAX;
    double bMinV = DBL_MAX, bMaxV = -DBL_MAX;

    if (drop == 0)
    {
        for (int i = 0; i < beam.nVertices; ++i)
        {
            const Point3f &p = beam.arr[i];
            const double t = (DotProduct(p, normal) + normal.d_param) * invDp;
            const double u = p.cy - beam.direction.cy * t;
            const double v = p.cz - beam.direction.cz * t;
            if (u < bMinU) bMinU = u;
            if (u > bMaxU) bMaxU = u;
            if (v < bMinV) bMinV = v;
            if (v > bMaxV) bMaxV = v;
        }
    }
    else if (drop == 1)
    {
        for (int i = 0; i < beam.nVertices; ++i)
        {
            const Point3f &p = beam.arr[i];
            const double t = (DotProduct(p, normal) + normal.d_param) * invDp;
            const double u = p.cx - beam.direction.cx * t;
            const double v = p.cz - beam.direction.cz * t;
            if (u < bMinU) bMinU = u;
            if (u > bMaxU) bMaxU = u;
            if (v < bMinV) bMinV = v;
            if (v > bMaxV) bMaxV = v;
        }
    }
    else
    {
        for (int i = 0; i < beam.nVertices; ++i)
        {
            const Point3f &p = beam.arr[i];
            const double t = (DotProduct(p, normal) + normal.d_param) * invDp;
            const double u = p.cx - beam.direction.cx * t;
            const double v = p.cy - beam.direction.cy * t;
            if (u < bMinU) bMinU = u;
            if (u > bMaxU) bMaxU = u;
            if (v < bMinV) bMinV = v;
            if (v > bMaxV) bMaxV = v;
        }
    }

    const double *bounds = m_facetProjectionBounds[drop][facetId];
    const double environmentMargin = CpuTraceProjectedPrefilterMargin();
    const double margin = m_traceCpuProjectedPrefilterMargin >= 0.0
        ? m_traceCpuProjectedPrefilterMargin
        : (environmentMargin >= 0.0
            ? environmentMargin
            : geometry_length_tolerance(m_geometryScale));
    return !(bMaxU < bounds[0] - margin || bounds[1] < bMinU - margin ||
             bMaxV < bounds[2] - margin || bounds[3] < bMinV - margin);
}

void ScatteringNonConvex::BuildFacetVisibilityCache()
{
    m_visibilityCacheOwner = nullptr;
    for (int locInt = 0; locInt < 2; ++locInt)
    {
        Location loc = locInt == 0 ? Location::In : Location::Out;
        for (int sourceFacet = 0; sourceFacet < m_particle->nFacets; ++sourceFacet)
        {
            size_t &count = m_visibleFacetCacheSize[locInt][sourceFacet];
            count = 0;

            int begin = 0;
            int end = m_particle->nFacets;
            if (m_particle->isAggregated && loc == Location::In)
                m_particle->GetParticalFacetIdRangeByFacetId(sourceFacet, begin, end);

            const Point3f &beamNormal = -m_facets[sourceFacet].normal[!loc];
            const Point3f &beamCenter = m_facets[sourceFacet].center;
            const double tolerance = geometry_length_tolerance(
                m_geometryScale);
            for (int targetFacet = begin; targetFacet < end; ++targetFacet)
            {
                if (targetFacet == sourceFacet)
                    continue;

                if (FacetReachesForwardHalfSpace(m_facets[targetFacet],
                                                 beamCenter, beamNormal,
                                                 tolerance))
                    m_visibleFacetCache[locInt][sourceFacet][count++] = targetFacet;
            }
        }
    }
    m_visibilityCacheBuilt = true;
}

void ScatteringNonConvex::BuildFacetNormalCache()
{
    for (int locInt = 0; locInt < 2; ++locInt)
    {
        for (int facetId = 0; facetId < m_particle->nFacets; ++facetId)
        {
            const Point3f &normal = (locInt == 0)
                                  ? m_facets[facetId].in_normal
                                  : m_facets[facetId].ex_normal;
            m_facetNormalComponents[locInt][facetId][0] = normal.cx;
            m_facetNormalComponents[locInt][facetId][1] = normal.cy;
            m_facetNormalComponents[locInt][facetId][2] = normal.cz;
        }
    }
}

void ScatteringNonConvex::BuildFacetProjectionCache()
{
    for (int locInt = 0; locInt < 2; ++locInt)
    {
        for (int facetId = 0; facetId < m_particle->nFacets; ++facetId)
        {
            const Point3f &normal = (locInt == 0)
                                  ? m_facets[facetId].in_normal
                                  : m_facets[facetId].ex_normal;
            const double ax = std::fabs(normal.cx);
            const double ay = std::fabs(normal.cy);
            const double az = std::fabs(normal.cz);
            m_facetProjectionDrop[locInt][facetId] =
                (ax > ay && ax > az) ? 0 : ((ay > az) ? 1 : 2);
        }
    }
    for (int facetId = 0; facetId < m_particle->nFacets; ++facetId)
    {
        const Facet &facet = m_facets[facetId];
        double minX = DBL_MAX, maxX = -DBL_MAX;
        double minY = DBL_MAX, maxY = -DBL_MAX;
        double minZ = DBL_MAX, maxZ = -DBL_MAX;
        for (int i = 0; i < facet.nVertices; ++i)
        {
            const Point3f &p = facet.arr[i];
            if (p.cx < minX) minX = p.cx;
            if (p.cx > maxX) maxX = p.cx;
            if (p.cy < minY) minY = p.cy;
            if (p.cy > maxY) maxY = p.cy;
            if (p.cz < minZ) minZ = p.cz;
            if (p.cz > maxZ) maxZ = p.cz;
        }
        m_facetProjectionBounds[0][facetId][0] = minY;
        m_facetProjectionBounds[0][facetId][1] = maxY;
        m_facetProjectionBounds[0][facetId][2] = minZ;
        m_facetProjectionBounds[0][facetId][3] = maxZ;
        m_facetProjectionBounds[1][facetId][0] = minX;
        m_facetProjectionBounds[1][facetId][1] = maxX;
        m_facetProjectionBounds[1][facetId][2] = minZ;
        m_facetProjectionBounds[1][facetId][3] = maxZ;
        m_facetProjectionBounds[2][facetId][0] = minX;
        m_facetProjectionBounds[2][facetId][1] = maxX;
        m_facetProjectionBounds[2][facetId][2] = minY;
        m_facetProjectionBounds[2][facetId][3] = maxY;
    }
}

void ScatteringNonConvex::PrepareForParallelTrace()
{
    if (!m_visibilityCacheBuilt)
        BuildFacetVisibilityCache();
}

void ScatteringNonConvex::CopyVisibilityCacheFrom(const ScatteringNonConvex &source)
{
    m_visibilityCacheBuilt = source.m_visibilityCacheBuilt;
    m_visibilityCacheOwner = nullptr;
    if (!m_visibilityCacheBuilt)
        return;
    const ScatteringNonConvex &cache = source.m_visibilityCacheOwner
                                     ? *source.m_visibilityCacheOwner
                                     : source;
    std::memcpy(m_visibleFacetCache, cache.m_visibleFacetCache,
                sizeof(m_visibleFacetCache));
    std::memcpy(m_visibleFacetCacheSize, cache.m_visibleFacetCacheSize,
                sizeof(m_visibleFacetCacheSize));
}

void ScatteringNonConvex::FindVisibleFacets(const Beam &beam, IntArray &facetIds)
{
    const int locInt = beam.location == Location::In ? 0 : 1;
    const ScatteringNonConvex &cache =
        m_visibilityCacheOwner ? *m_visibilityCacheOwner : *this;
    const size_t count = cache.m_visibleFacetCacheSize[locInt][beam.lastFacetId];
    const double (*facetNormals)[3] =
        m_facetNormalComponents[beam.location == Location::In ? 1 : 0];
    const double bx = beam.direction.cx;
    const double by = beam.direction.cy;
    const double bz = beam.direction.cz;
    for (size_t idx = 0; idx < count; ++idx)
    {
        int i = cache.m_visibleFacetCache[locInt][beam.lastFacetId][idx];
        const double *facetNormal = facetNormals[i];
        const double cosFB = bx*facetNormal[0]
                           + by*facetNormal[1]
                           + bz*facetNormal[2];

        if (cosFB >= -EPS_PROJECTION) // conservative candidate filter
        {
            facetIds.Add(i);
        }
    }
}

/* TODO: Разобраться с параметром 'n' (кол-во вн. столкновений)
 при заданных траекториях, возможно он не нужен т.к. заранее известен путь */
bool ScatteringNonConvex::PushBeamPartsToTree(const Beam &beam,
                                              const PolygonArray &parts)
{
    Beam tmp = beam; // OPT: try to replace 'tmp' to 'beam'

    for (size_t i = 0; i < parts.size; ++i)
    {
        tmp = parts.arr[i];

        if (!Scattering::PushBeamToTree(tmp, tmp.lastFacetId,
                                        tmp.nActs, tmp.location))
            return false;
    }

    return true;
}

template<class T>
bool ScatteringNonConvex::PushBeamToTree(Beam &beam, const Beam &oldBeam,
                                         const T &newId, int facetId,
                                         Location loc)
{
    beam.id = newId;
    beam.locations = oldBeam.locations;
#ifdef _DEBUG // DEB
    beam.dirs = oldBeam.dirs;
    if (beam.id == 66557)
        int fff = 0;
#endif
    return Scattering::PushBeamToTree(beam, facetId, oldBeam.nActs+1, loc);
}

bool ScatteringNonConvex::SplitBeamByFacet(const Polygon &intersection,
                                           int facetId, Beam &beam, bool &ok)
{
    auto newId = RecomputeTrackId(beam.id, facetId);
    Facet &facet = m_facets[facetId];

    m_polygonBuffer.Clear();
    Polygon reachedIntersection;
    const bool cutParent = CutBeamByFacet(intersection, facetId, beam,
                                          reachedIntersection,
                                          m_polygonBuffer);
    if (reachedIntersection.nVertices >= MIN_VERTEX_NUM)
    {
        Beam inBeam, outBeam;
        inBeam.SetPolygon(reachedIntersection);
        outBeam.SetPolygon(reachedIntersection);
#ifdef _DEBUG // DEB
//if (beam.lastFacetId==0 && facetId==6)
//if (beam.trackId.toLong()==9633 && facetID==0)
//    int f =0;
        inBeam.pols = beam.pols;
        outBeam.pols = beam.pols;
        inBeam.pols.push_back(reachedIntersection);
        outBeam.pols.push_back(reachedIntersection);

        if (newId == 3496)
            int fff = 0;
#endif

        const bool hasOutBeam = SetOpticalBeamParams(facet, beam, inBeam,
                                                     outBeam);
        ok = PushBeamToTree(inBeam, beam, newId, facetId, Location::In);
        if (!ok)
            return false;

        if (hasOutBeam)
        {
            ok = PushBeamToTree(outBeam, beam, newId, facetId, Location::Out);
            if (!ok)
                return false;
        }
    }

    if (!cutParent)
        return false;

    if (m_polygonBuffer.size == 0)
    {
        beam.nVertices = 0;
        return false;
    }

    bool isDivided = m_polygonBuffer.size > CLIP_RESULT_SINGLE;

    if (isDivided)
    {	// beam had divided by facet
        if (m_polygonBuffer.size == 2)
        {
            double a0 = m_polygonBuffer.arr[0].Area();
            double a1 = m_polygonBuffer.arr[1].Area();

            const double areaTolerance = geometry_area_tolerance(
                m_geometryScale);
            if (a0 <= areaTolerance || a1 <= areaTolerance)
            {
                beam = a0 >= a1 ? m_polygonBuffer.arr[0]
                                : m_polygonBuffer.arr[1];
                isDivided = false;
                m_traceCutoffStatistics->smallFragmentSimplifications.fetch_add(
                    1, std::memory_order_relaxed);
            }
            else
            {
                double r = a0/a1;

                if (std::isfinite(restriction) && r >= restriction)
                {
                    beam = m_polygonBuffer.arr[0];
                    isDivided = false;
                    m_traceCutoffStatistics->smallFragmentSimplifications.fetch_add(
                        1, std::memory_order_relaxed);
                }
                else if (std::isfinite(restriction) && r < 1.0/restriction)
                {
                    beam = m_polygonBuffer.arr[1];
                    isDivided = false;
                    m_traceCutoffStatistics->smallFragmentSimplifications.fetch_add(
                        1, std::memory_order_relaxed);
                }
                else
                {
                    ok = PushBeamPartsToTree(beam, m_polygonBuffer);

                    if (!ok)
                    {
                        return false;
                    }

                    beam.nVertices = 0;
                }
            }

//            double a0 = sqrt(m_polygonBuffer.arr[0].Area());
//            double a1 = sqrt(m_polygonBuffer.arr[1].Area());

//            if (a0 < m_wave && a1 > m_wave)
//            {
//                beam = m_polygonBuffer.arr[1];
//                isDivided = false;
//            }

//            if (a0 > m_wave && a1 < m_wave)
//            {
//                beam = m_polygonBuffer.arr[0];
//                isDivided = false;
//            }

//            if (a0 < m_wave && a1 < m_wave)
//            {
//                ok = PushBeamPartsToTree(beam, m_polygonBuffer);

//                if (!ok)
//                {
//                    return false;
//                }

//                beam.nVertices = 0;
//            }
        }
        else
        {
            ok = PushBeamPartsToTree(beam, m_polygonBuffer);

            if (!ok)
            {
                return false;
            }

            beam.nVertices = 0;
        }
    }
    else if (m_polygonBuffer.size == CLIP_RESULT_SINGLE)
    {
        beam = m_polygonBuffer.arr[0];
//#ifdef _DEBUG // DEB
//    if (beam.id == 477)
//        int fff = 0;
//#endif
    }

    return isDivided;
}

bool ScatteringNonConvex::ScatterLight(double beta, double gamma,
                                       const std::vector<std::vector<int>> &tracks,
                                       std::vector<Beam> &scaterredBeams)
{
//	m_particle->Rotate(beta, gamma, 0);

//	for (const std::vector<int> &track : tracks)
//	{
//		int facetID = track.at(0);

//		bool isIncident;
//		TraceFirstBeamFixedFacet(facetID, isIncident);

//		if (!isIncident)
//		{
//			continue;
//		}

//		for (size_t i = 1; i < track.size(); ++i)
//		{
//			int facetID = track.at(i);

//			std::vector<Beam> buffer; // для прошедших пучков (не дублированных)

//			while (m_treeSize != 0)
//			{
//				Beam beam = m_beamTree[--m_treeSize];

//				IntArray facetIDs;
//				SelectVisibleFacets(beam, facetIDs);
//				int index = FindFacetId(facetID, facetIDs);

//				if (index != -1)
//				{
//					bool isDivided;
//					SplitBeamByFacet(beam, facetID, isDivided);

//					Polygon intersected;
//					Intersect(facetID, beam, intersected);

//					if (intersected.size >= MIN_VERTEX_NUM)
//					{
//						Beam inBeam, outBeam;
//						inBeam.SetPolygon(intersected);
//						outBeam.SetPolygon(intersected);

//						Facet &facet = m_facets[facetID];
//						bool hasOutBeam = SetOpticalBeamParams(facet, beam,
//															   inBeam, outBeam);
//						PushBeamsToBuffer(facetID, beam, hasOutBeam,
//										  inBeam, outBeam, buffer);
//					}
//				}
//			}

//			if (buffer.empty())
//			{
//				break;
//			}

//			for (const Beam &b : buffer)
//			{	// добавляем прошедшие пучки
//				assert(m_treeSize < MAX_BEAM_REFL_NUM);
//				m_beamTree[m_treeSize++] = b;
//			}
//		}

//		while (m_treeSize != 0)
//		{
//			scaterredBeams.push_back(m_beamTree[--m_treeSize]);
//		}
//	}
//}

//void ScatteringNonConvex::TraceFirstBeamFixedFacet(int facetID, bool &isIncident)
//{
//	isIncident = false;

//	IntArray facetIDs;
//	SelectVisibleFacetsForLight(facetIDs);

//	int index = FindFacetId(facetID, facetIDs);

//	for (int i = 0; (facetIDs.arr[i] < index) || (i < facetIDs.size); ++i)
//	{
//		int id = facetIDs.arr[i];
//		CutBeamByFacet();
//	}
//	if (index != -1)
//	{
//		SplitByFacet(facetIDs, index);
//		isIncident = true;
//	}
    (void)beta;
    (void)gamma;
    (void)tracks;
    (void)scaterredBeams;
    return false;
}
