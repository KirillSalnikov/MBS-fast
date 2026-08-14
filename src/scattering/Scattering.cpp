#include "Scattering.h"

#include <float.h>
#include <assert.h>
#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <stdexcept>

#include "geometry_lib.h"

using namespace std;
using ::complex;

static bool DisableTraceTrackIds()
{
    static const bool disabled = []() {
        const char *value = std::getenv("MBS_DISABLE_TRACK_IDS");
        return value != nullptr && value[0] == '1' && value[1] == '\0';
    }();
    return disabled;
}

static bool ForceTraceTrackIds()
{
    static const bool forced = []() {
        const char *value = std::getenv("MBS_FORCE_TRACK_IDS");
        return value != nullptr && value[0] == '1' && value[1] == '\0';
    }();
    return forced;
}

Scattering::Scattering(Particle *particle, Light *incidentLight, bool isOpticalPath,
                       int nActs)
    : m_particle(particle),
      m_geometryScale(particle->MaximalDimention()),
      splitting(isOpticalPath),
      m_incidentLight(incidentLight),
      m_nActs(nActs),
      m_incidentEnergy(0.0)
{
    m_traceCutoffStatistics = std::make_shared<TraceCutoffStatistics>();
    m_facets = m_particle->facets;

    m_incidentDir = m_incidentLight->direction;
    m_incidentDir.d_param = m_incidentLight->direction.d_param;

    m_polarBasis = m_incidentLight->polarizationBasis;

    ::complex ri = m_particle->GetRefractiveIndex();
    splitting.ComputeRiParams(ri);
}

Scattering::~Scattering()
{
}

bool Scattering::EnsureBeamTree()
{
    if (m_beamTree.capacity() == 0)
    {
        const char *value = std::getenv("MBS_TRACE_TREE_RESERVE");
        long reserve = 65536;
        if (value && *value)
        {
            char *end = nullptr;
            long parsed = std::strtol(value, &end, 10);
            if (end && *end == '\0' && parsed > 0)
                reserve = parsed;
        }
        m_beamTree.reserve((size_t)reserve);
    }
    return true;
}

void Scattering::CopyRuntimeOptionsFrom(const Scattering &source)
{
    m_wave = source.m_wave;
    restriction = source.restriction;
    m_traceCutoffJRel = source.m_traceCutoffJRel;
    m_traceCutoffAreaRel = source.m_traceCutoffAreaRel;
    m_traceCutoffImportanceRel = source.m_traceCutoffImportanceRel;
    m_traceMaxBeams = source.m_traceMaxBeams;
    m_cutoffProfileName = source.m_cutoffProfileName;
    m_traceCutoffStatistics = source.m_traceCutoffStatistics;
    m_gpuTracePrefilter = source.m_gpuTracePrefilter;
    m_traceCpuProjectedPrefilter = source.m_traceCpuProjectedPrefilter;
    m_traceCpuProjectedPrefilterMargin = source.m_traceCpuProjectedPrefilterMargin;
    m_tracePrefilterStats = source.m_tracePrefilterStats;
    m_trackIdsRequired = source.m_trackIdsRequired;
}

IdType Scattering::Scattering::RecomputeTrackId(const IdType &oldId, int facetId)
{
    if (DisableTraceTrackIds() && !m_trackIdsRequired)
        return 0;
    if (!m_trackIdsRequired && !ForceTraceTrackIds()
        && std::fabs(imag(m_particle->GetRefractiveIndex())) <= DBL_EPSILON)
        return 0;

    return (oldId + (facetId + 1)) * (m_particle->nFacets + 1);
}

bool Scattering::PushBeamToTree(Beam &beam, int facetId, int level, Location location)
{
    beam.SetTracingParams(facetId, level, location);
#ifdef _DEBUG // DEB
    beam.dirs.push_back(beam.direction);
#endif
    if (location == Location::In && IsTracePruned(beam))
    {
        return true;
    }
    if (m_treeSize >= MAX_BEAM_REFL_NUM)
    {
        m_traceCutoffStatistics->hardBeamLimitHits.fetch_add(
            1, std::memory_order_relaxed);
        throw std::runtime_error(
            "beam tree reached the compiled 65536-branch hard limit; "
            "the orientation result would be incomplete");
    }
    if (m_traceMaxBeams > 0 && m_treeSize >= m_traceMaxBeams)
    {
        m_traceCutoffStatistics->configuredBeamLimitHits.fetch_add(
            1, std::memory_order_relaxed);
        throw std::runtime_error(
            "beam tree reached --trace-max-beams="
            + std::to_string(m_traceMaxBeams)
            + "; raise the limit or use 0 to disable it");
    }
    if (!EnsureBeamTree())
        return false;

    m_beamTree.push_back(beam);
    m_treeSize = (int)m_beamTree.size();
    return true;
}

void Scattering::SetIncidentBeamOpticalParams(unsigned facetId,
                                              Beam &inBeam, Beam &outBeam)
{
    const Point3f &normal = m_facets[facetId].in_normal;
    splitting.ComputeCosA(m_incidentDir, normal);

    if (!splitting.IsNormalIncidence()) // regular incidence
    {
        Beam fakeIncidentBeam;
        fakeIncidentBeam.SetLight(*m_incidentLight);
        const Point3f &facetNormal = m_facets[facetId].in_normal;
        ComputePolarisationParams(-fakeIncidentBeam.direction, facetNormal,
                                  fakeIncidentBeam);
        splitting.ComputeRegularBeamParamsExternal(facetNormal,
                                                     fakeIncidentBeam,
                                                     inBeam, outBeam);
    }
    else // normal incidence
    {
        splitting.ComputeNormalBeamParamsExternal(*m_incidentLight,
                                                    inBeam, outBeam);
    }
}

void Scattering::ComputePolarisationParams(const Vector3f &dir,
                                           const Vector3f &facetNormal, Beam &beam)
{
    Point3f newBasis = CrossProduct(facetNormal, dir);
    Normalize(newBasis);
    beam.RotateJMatrix(newBasis);
    beam.polarizationBasis = newBasis;
}

void Scattering::SplitLightToBeams(int facetId, Beam &inBeam, Beam &outBeam)
{
    SetPolygonByFacet(facetId, inBeam); // REF: try to exchange this to inBeam = m_facets[facetId]
    SetPolygonByFacet(facetId, outBeam);
    SetIncidentBeamOpticalParams(facetId, inBeam, outBeam);

//	if (m_isOpticalPath)
    {
        Point3f p = inBeam.Center();
        double path = splitting.ComputeIncidentOpticalPath(m_incidentDir, p);
        inBeam.opticalPath = 0;
        outBeam.opticalPath = 0;
        inBeam.AddOpticalPath(path);
        outBeam.AddOpticalPath(path);
        outBeam.opticalPath += splitting.ComputeOutgoingOpticalPath(outBeam);
    }
}

void Scattering::AddProjectedIncidentEnergy(int facetId, const Polygon &lightedPolygon)
{
    const Point3f &normal = m_facets[facetId].in_normal;
    double cosA = DotProduct(m_incidentDir, normal);
    m_incidentEnergy += lightedPolygon.Area() * cosA;
}

// TODO: пофиксить
bool Scattering::ScatterLight(double /*beta*/, double /*gamma*/, const std::vector<std::vector<int>> &/*tracks*/,
                              std::vector<Beam> &/*outBeams*/)
{
//	m_particle->Rotate(beta, gamma, 0);

//	for (unsigned int i = 0; i < tracks.size(); ++i)
//	{
//		int facetId = tracks.at(i).at(0);
//		const Point3f &extNormal = m_facets[facetId].ex_normal;

//		double cosIN = DotProduct(m_incidentDir, extNormal);

//		if (cosIN < EPS_COS_90) /// beam is not incident to this facet
//		{
//			continue;
//		}

//		std::vector<Beam> outBuff;
//		Beam incidentBeam;

//		// first incident beam
//		{
//			Beam outBeam;
//			SplitLightToBeams(facetId, incidentBeam, outBeam);
//			outBuff.push_back(outBeam);
//		}

//		unsigned int size = tracks.at(i).size();

//		try // internal beams
//		{
//			for (unsigned int j = 1; j < size; ++j)
//			{
//				facetId = tracks.at(i).at(j);

//				Beam inBeam;
//				SplitSecondaryBeams(incidentBeam, facetId, inBeam, outBuff);

//				incidentBeam = inBeam;
//			}
//		}
//		catch (const std::exception &)
//		{
//			continue;
//		}

//		outBeams.push_back(outBuff.back());
    //	}
    return false;
}

void Scattering::OrderVertices2f(std::vector<Point2f> &vertices,
                                 Polygon &orderedPolygon)
{
    if (vertices.empty())
        return;

    std::sort(vertices.begin(), vertices.end(),
              [](const Point2f &a, const Point2f &b)
              {
                  return a.x < b.x || (a.x == b.x && a.y < b.y);
              });

    std::vector<Point2f> hull;
    hull.reserve(vertices.size());

    auto cross = [](const Point2f &o, const Point2f &a, const Point2f &b)
    {
        const double ax = static_cast<double>(a.x) - o.x;
        const double ay = static_cast<double>(a.y) - o.y;
        const double bx = static_cast<double>(b.x) - o.x;
        const double by = static_cast<double>(b.y) - o.y;
        return ax*by - ay*bx;
    };

    for (const Point2f &p : vertices)
    {
        while (hull.size() >= 2
               && cross(hull[hull.size()-2], hull[hull.size()-1], p) <= 0)
        {
            hull.pop_back();
        }
        hull.push_back(p);
    }

    const size_t lowerSize = hull.size();
    for (int i = (int)vertices.size() - 2; i >= 0; --i)
    {
        const Point2f &p = vertices[i];
        while (hull.size() > lowerSize
               && cross(hull[hull.size()-2], hull[hull.size()-1], p) <= 0)
        {
            hull.pop_back();
        }
        hull.push_back(p);
    }

    if (hull.size() > 1)
        hull.pop_back();

    double signedArea = 0.0;
    for (size_t i = 0; i < hull.size(); ++i)
    {
        const Point2f &a = hull[i];
        const Point2f &b = hull[(i + 1) % hull.size()];
        signedArea += a.x*b.y - b.x*a.y;
    }

    // Keep the shadow aperture normal along -z, matching Natalia's N = -k.
    if (signedArea > 0.0)
        std::reverse(hull.begin(), hull.end());

    for (const Point2f &p : hull)
        orderedPolygon.AddVertex(Point3f(p.x, p.y, 0));
}

void Scattering::ProjectParticleToXY(std::vector<Point2f> &projected)
{
    Point3f n(0, 0, 1, 0);
    n.d_param = m_geometryScale;

    for (int i = 0; i < m_particle->nFacets; i++)
    {
        auto &f = m_particle->facets[i];

        if (DotProduct(f.in_normal, -n) < EPS_GRAZING_INCIDENCE)
        {
            for (int j = 0; j < f.nVertices; j++)
            {
//				auto p = ProjectPointToPlane(f.arr[j], -m_incidentDir, n);
//				projected.push_back(Point2f(-p.coordinates[0], -p.coordinates[1]));
                double tmp = (n.d_param - DotProduct(n, f.arr[j]));
                auto p = f.arr[j] + n*tmp;
                projected.push_back(Point2f(p.coordinates[0], p.coordinates[1]));
            }
        }
    }
}

void Scattering::RemoveDublicatedVertices2f(const std::vector<Point2f> &projected,
                                          std::vector<Point2f> &cleared)
{
    for (size_t i = 0; i < projected.size(); ++i)
    {
        bool isUnique = true;

        for (size_t j = i + 1; j < projected.size(); ++j)
        {
            const double scale = std::max({
                std::fabs(static_cast<double>(projected[i].x)),
                std::fabs(static_cast<double>(projected[i].y)),
                std::fabs(static_cast<double>(projected[j].x)),
                std::fabs(static_cast<double>(projected[j].y))});
            const double tolerance = geometry_length_tolerance(scale);
            if (std::fabs(static_cast<double>(projected[i].x) - projected[j].x)
                    <= tolerance
                && std::fabs(static_cast<double>(projected[i].y) - projected[j].y)
                    <= tolerance)
            {
                isUnique = false;
            }
        }

        if (isUnique)
        {
            cleared.push_back(projected[i]);
        }
    }
}

void Scattering::FormShadowBeam(std::vector<Beam> &scaterredBeams)
{
    std::vector<Point2f> projected;
    ProjectParticleToXY(projected);

    std::vector<Point2f> projectedCleared;
    RemoveDublicatedVertices2f(projected, projectedCleared);

    Beam shadowBeam;
    OrderVertices2f(projectedCleared, shadowBeam);

    Matrix2x2c jones;
    jones.m11 = -jones.m11;
    jones.m22 = -jones.m22;
    shadowBeam.SetMatrix(jones);

    shadowBeam.direction = m_incidentLight->direction;
    shadowBeam.polarizationBasis = m_incidentLight->polarizationBasis;
    // All physical beams used the same artificial 2*FAR_ZONE_DISTANCE phase.
    // Omitting that common phase avoids reducing a very large argument in
    // sin/cos and leaves every Mueller matrix unchanged.
    shadowBeam.opticalPath = 0.0;
    shadowBeam.lastFacetId = __INT_MAX__;
    scaterredBeams.push_back(shadowBeam);
}

bool Scattering::IsTerminalAct(const Beam &beam)
{
    // Test the relative cutoff before selecting and intersecting candidate
    // facets.  For an outgoing branch the caller still performs the final
    // external clipping, matching the historical termination semantics.
    return beam.nActs >= m_nActs || IsTracePruned(beam);
}

void Scattering::ResetTraceReference()
{
    m_traceRefJNorm = 0;
    m_traceRefArea = 0;
    m_traceRefImportance = 0;
}

void Scattering::UpdateTraceReference(const Beam &beam)
{
    double jn = beam.J.Norm();
    double ar = beam.Area();
    double importance = jn * ar;
    if (jn > m_traceRefJNorm) m_traceRefJNorm = jn;
    if (ar > m_traceRefArea) m_traceRefArea = ar;
    if (importance > m_traceRefImportance) m_traceRefImportance = importance;
}

bool Scattering::IsTracePruned(const Beam &beam) const
{
    if (beam.nActs == 0)
        return false;
    const bool useJ = (m_traceCutoffJRel > 0 && m_traceRefJNorm > 0);
    const bool useArea = (m_traceCutoffAreaRel > 0 && m_traceRefArea > 0);
    const bool useImportance = (m_traceCutoffImportanceRel > 0
                                && m_traceRefImportance > 0);
    if (!useJ && !useArea && !useImportance)
        return false;

    m_traceCutoffStatistics->evaluated.fetch_add(1, std::memory_order_relaxed);

    double jn = 0;
    double ar = 0;
    if (useJ || useImportance)
        jn = beam.J.Norm();
    if (useArea || useImportance)
        ar = beam.Area();

    const bool smallJ = useJ && jn < m_traceCutoffJRel * m_traceRefJNorm;
    const bool smallArea = useArea && ar < m_traceCutoffAreaRel * m_traceRefArea;
    const bool smallImportance = useImportance
        && jn * ar < m_traceCutoffImportanceRel * m_traceRefImportance;
    if (!(smallJ || smallArea || smallImportance))
        return false;

    m_traceCutoffStatistics->rejected.fetch_add(1, std::memory_order_relaxed);
    if (smallJ)
        m_traceCutoffStatistics->rejectedJones.fetch_add(
            1, std::memory_order_relaxed);
    if (smallArea)
        m_traceCutoffStatistics->rejectedArea.fetch_add(
            1, std::memory_order_relaxed);
    if (smallImportance)
        m_traceCutoffStatistics->rejectedImportance.fetch_add(
            1, std::memory_order_relaxed);
    return true;
}

std::string Scattering::TraceCutoffReport() const
{
    return FormatTraceCutoffReport(*m_traceCutoffStatistics,
                                   m_cutoffProfileName,
                                   m_traceCutoffJRel,
                                   m_traceCutoffAreaRel,
                                   m_traceCutoffImportanceRel,
                                   restriction,
                                   m_traceMaxBeams);
}

void Scattering::Difference(const Polygon &subject, const Vector3f &subjNormal,
                         const Polygon &clip, const Vector3f &clipNormal,
                         const Vector3f &clipDir, PolygonArray &difference,
                         Polygon *overlap) const
{
    if (overlap != nullptr)
        overlap->Clear();
    if (subject.nVertices < MIN_VERTEX_NUM)
        return;
    if (clip.nVertices < MIN_VERTEX_NUM)
    {
        difference.Push(subject);
        return;
    }

    Point3f clipProjection[MAX_VERTEX_NUM];
    bool isProjected = ProjectToFacetPlane(clip, clipDir, subjNormal,
                                           clipProjection);

    if (!isProjected)
    {
        difference.Push(subject);
        return;
    }

    // Clipped target patches inherit the winding of the incident aperture.
    // Enforce the Difference contract after projection instead of relying on
    // the target facet or aperture to have a particular relative orientation.
    Polygon projectedClip;
    if (overlap != nullptr)
    {
        SetConvexOutputPolygon(clipProjection, clip.nVertices, clipNormal,
                               projectedClip);
    }
    else
    {
        SetOutputPolygon(clipProjection, clip.nVertices, projectedClip);
        if (projectedClip.nVertices >= MIN_VERTEX_NUM
            && DotProduct(projectedClip.Normal(), clipNormal) < 0.0)
        {
            Polygon::InverseVertexOrder(projectedClip);
        }
    }
    if (projectedClip.nVertices < MIN_VERTEX_NUM)
    {
        difference.Push(subject);
        return;
    }
    const Point3f effectiveClipNormal = projectedClip.Normal();

    const double projectionDenominator = DotProduct(clipDir, subjNormal);
    const double depthTolerance = geometry_depth_tolerance(m_geometryScale);
    bool hasForwardDepth = false;
    for (int i = 0; i < clip.nVertices; ++i)
    {
        const double depth = (DotProduct(clip.arr[i], subjNormal)
                              + subjNormal.d_param) / projectionDenominator;
        if (depth <= depthTolerance)
        {
            hasForwardDepth = true;
            break;
        }
    }
    if (!hasForwardDepth)
    {
        difference.Push(subject);
        return;
    }

    Polygon retained = subject;
    Point3f edgeEnd = projectedClip.arr[projectedClip.nVertices - 1];
    const double normalLength = Length(effectiveClipNormal);
    for (int edgeIndex = 0; edgeIndex < projectedClip.nVertices; ++edgeIndex)
    {
        if (retained.nVertices < MIN_VERTEX_NUM)
            break;

        const Point3f edgeStart = edgeEnd;
        edgeEnd = projectedClip.arr[edgeIndex];
        const Point3f edge = edgeEnd - edgeStart;
        const double edgeLength = Length(edge);

        double insideScalar[MAX_VERTEX_NUM];
        double outsideScalar[MAX_VERTEX_NUM];
        double maximumOffset = edgeLength;
        for (int vertex = 0; vertex < retained.nVertices; ++vertex)
        {
            const Point3f offset = retained.arr[vertex] - edgeStart;
            maximumOffset = std::max(maximumOffset, Length(offset));
            const double side = DotProduct(CrossProduct(edge, offset),
                                           effectiveClipNormal);
            insideScalar[vertex] = side;
            outsideScalar[vertex] = -side;
        }
        const double sideTolerance = GEOMETRY_RELATIVE_EPSILON
            * edgeLength * maximumOffset * normalLength;

        Point3f insideVertices[MAX_VERTEX_NUM];
        Point3f outsideVertices[MAX_VERTEX_NUM];
        const int insideSize = clip_polygon_by_nonnegative_scalar(
            retained.arr, insideScalar, retained.nVertices, sideTolerance,
            insideVertices);
        const int outsideSize = clip_polygon_by_nonnegative_scalar(
            retained.arr, outsideScalar, retained.nVertices, sideTolerance,
            outsideVertices);

        Polygon outside;
        SetOutputPolygon(outsideVertices, outsideSize, outside);
        if (outside.nVertices >= MIN_VERTEX_NUM)
            difference.Push(outside);

        SetOutputPolygon(insideVertices, insideSize, retained);
    }
    if (overlap != nullptr)
        *overlap = retained;
}

bool Scattering::ProjectToFacetPlane(const Polygon &polygon, const Vector3f &dir,
                                  const Point3f &normal, Point3f *projection) const
{
    const double dp = DotProduct(dir, normal);
    const double productScale = Length(dir)*Length(normal);
    if (productScale <= DBL_MIN
        || std::fabs(dp) <= geometry_parallel_tolerance(productScale))
    {
        return false; /// beam is parallel to facet
    }

    for (int i = 0; i < polygon.nVertices; ++i)
    {
        const Point3f &p = polygon.arr[i];
        const double t = (DotProduct(p, normal) + normal.d_param) / dp;
        projection[i] = Point3f(p.cx - t*dir.cx,
                                p.cy - t*dir.cy,
                                p.cz - t*dir.cz);
    }

    return true;
}

/// NOTE: вершины пучка и грани должны быть ориентированы в одном направлении
void Scattering::Intersect(int facetID, const Beam &beam, Polygon &intersection) const
{
    intersection.nVertices = 0;
    Point3f outputPoints[MAX_VERTEX_NUM];
    const Facet &facet = m_particle->facets[facetID];
    // REF: перенести в случай невыпуклых частиц
    const Point3f &normal = facet.in_normal;

    const Point3f &normal1 = (beam.location == Location::In) ? facet.in_normal
                                                            : facet.ex_normal;
    bool isProjected = ProjectToFacetPlane(beam, beam.direction, normal1,
                                           outputPoints);
    if (!isProjected)
    {
        return;
    }
    if (beam.nVertices <= 0)
    {
        intersection.nVertices = 0;
        return;
    }
    if (beam.nVertices > MAX_VERTEX_NUM)
        throw std::runtime_error("beam has an invalid vertex count");

    const Point3f normalToFacet = -normal;
    Point3f *outputPtr = outputPoints;
    int outputSize = beam.nVertices;

    Point3f buffer[MAX_VERTEX_NUM];
    Point3f *bufferPtr = buffer;
    int bufferSize;

    int facetSize = facet.nVertices;
    if (facetSize <= 0)
    {
        intersection.nVertices = 0;
        return;
    }
    if (facetSize > MAX_VERTEX_NUM)
        throw std::runtime_error("facet has an invalid vertex count");

    Point3f p1, p2;
    Point3f startPoint, endPoint;
    bool isInsideE, isInsideS;

    p2 = facet.arr[facetSize-1];

    for (int i = 0; i < facetSize; ++i)
    {
        p1 = p2;
        p2 = facet.arr[i];

        bufferSize = outputSize;
        outputSize = 0;

        Point3f *temp = outputPtr;
        outputPtr = bufferPtr;
        bufferPtr = temp;

        startPoint = bufferPtr[bufferSize-1];
        isInsideS = is_inside_i(startPoint, p1, p2, normalToFacet);

        bool isIntersected;

        for (int j = 0; j < bufferSize; ++j)
        {
            endPoint = bufferPtr[j];
            isInsideE = is_inside_i(endPoint, p1, p2, normalToFacet);

            if (isInsideE)
            {
                if (!isInsideS)
                {
                    const Point3f x = intersect_i(startPoint, endPoint,
                                                  p1, p2, normalToFacet,
                                                  isIntersected);
                    if (isIntersected)
                    {
                        if (outputSize >= MAX_VERTEX_NUM)
                            throw std::runtime_error(
                                "beam/facet intersection exceeds the 64-vertex geometry limit");
                        outputPtr[outputSize++] = x;
                    }
                }

                if (outputSize >= MAX_VERTEX_NUM)
                    throw std::runtime_error(
                        "beam/facet intersection exceeds the 64-vertex geometry limit");
                outputPtr[outputSize++] = endPoint;
            }
            else if (isInsideS)
            {
                const Point3f x = intersect_i(startPoint, endPoint,
                                              p1, p2, normalToFacet,
                                              isIntersected);
                if (isIntersected)
                {
                    if (outputSize >= MAX_VERTEX_NUM)
                        throw std::runtime_error(
                            "beam/facet intersection exceeds the 64-vertex geometry limit");
                    outputPtr[outputSize++] = x;
                }
            }

            startPoint = endPoint;
            isInsideS = isInsideE;
        }

        if (outputSize < MIN_VERTEX_NUM)
        {
            return;
        }
    }

    SetOutputPolygon(outputPtr, outputSize, intersection);
    if (intersection.nVertices < MIN_VERTEX_NUM)
        return;

    // A nonconvex target facet may cross the aperture plane of the current
    // beam.  Keep only the part reached for a nonnegative ray parameter.
    // SplitBeamByFacet subsequently intersects its back-projection with the
    // source aperture and uses that one region for both children and parent
    // subtraction.
    if (beam.lastFacetId >= 0 && beam.lastFacetId < m_particle->nFacets)
    {
        const Point3f &beamPlane =
            m_facets[beam.lastFacetId].normal[beam.location];
        const double denominator = DotProduct(beam.direction, beamPlane);
        const double denominatorScale = Length(beam.direction)*Length(beamPlane);
        if (denominatorScale > DBL_MIN
            && std::fabs(denominator)
                > geometry_parallel_tolerance(denominatorScale))
        {
            const double denominatorSign = denominator >= 0.0 ? 1.0 : -1.0;
            double forwardScalar[MAX_VERTEX_NUM];
            for (int i = 0; i < intersection.nVertices; ++i)
            {
                forwardScalar[i] =
                    (DotProduct(intersection.arr[i], beamPlane)
                     + beamPlane.d_param) * denominatorSign;
            }

            Point3f forwardIntersection[MAX_VERTEX_NUM];
            const double depthTolerance = geometry_depth_tolerance(
                m_geometryScale);
            const int forwardSize = clip_polygon_by_nonnegative_scalar(
                intersection.arr, forwardScalar, intersection.nVertices,
                depthTolerance, forwardIntersection);
            SetOutputPolygon(forwardIntersection, forwardSize, intersection);
        }
    }
}

void Scattering::SetOutputPolygon(const Point3f *outputPoints, int outputSize,
                                  Polygon &polygon) const
{
    polygon.nVertices = 0;
    if (outputSize <= 0)
    {
        polygon.nVertices = 0;
        return;
    }
    if (outputSize > MAX_VERTEX_NUM)
        throw std::runtime_error(
            "clipped polygon exceeds the 64-vertex geometry limit");

    const double geometryScale = m_geometryScale;
    Point3f p0 = outputPoints[outputSize-1];
    for (int i = 0; i < outputSize; ++i)
    {
        if (geometry_points_distinct(outputPoints[i], p0, geometryScale))
        {
            polygon.arr[polygon.nVertices++] = outputPoints[i];
        }

        p0 = outputPoints[i];
    }

    if (polygon.nVertices < MIN_VERTEX_NUM
        || polygon.Area() <= geometry_area_tolerance(
               m_geometryScale))
        polygon.nVertices = 0;
}

void Scattering::SetConvexOutputPolygon(const Point3f *outputPoints,
                                        int outputSize,
                                        const Point3f &normal,
                                        Polygon &polygon) const
{
    polygon.Clear();
    if (outputSize < MIN_VERTEX_NUM)
        return;
    if (outputSize > MAX_VERTEX_NUM)
        throw std::runtime_error(
            "convex polygon exceeds the 64-vertex geometry limit");

    const double ax = std::fabs(normal.cx);
    const double ay = std::fabs(normal.cy);
    const double az = std::fabs(normal.cz);
    const int drop = ax >= ay && ax >= az ? 0 : (ay >= az ? 1 : 2);

    struct ProjectedVertex
    {
        double x;
        double y;
        Point3f point;
    };
    std::vector<ProjectedVertex> vertices;
    vertices.reserve(outputSize);
    for (int i = 0; i < outputSize; ++i)
    {
        const Point3f &point = outputPoints[i];
        if (drop == 0)
            vertices.push_back({point.cy, point.cz, point});
        else if (drop == 1)
            vertices.push_back({point.cx, point.cz, point});
        else
            vertices.push_back({point.cx, point.cy, point});
    }
    std::sort(vertices.begin(), vertices.end(),
              [](const ProjectedVertex &left,
                 const ProjectedVertex &right)
              {
                  return left.x < right.x
                      || (left.x == right.x && left.y < right.y);
              });

    const double lengthTolerance = geometry_length_tolerance(m_geometryScale);
    std::vector<ProjectedVertex> unique;
    unique.reserve(vertices.size());
    for (const ProjectedVertex &vertex : vertices)
    {
        if (unique.empty()
            || std::fabs(vertex.x - unique.back().x) > lengthTolerance
            || std::fabs(vertex.y - unique.back().y) > lengthTolerance)
        {
            unique.push_back(vertex);
        }
    }
    if (unique.size() < MIN_VERTEX_NUM)
        return;

    const auto cross = [](const ProjectedVertex &origin,
                          const ProjectedVertex &first,
                          const ProjectedVertex &second)
    {
        return (first.x - origin.x)*(second.y - origin.y)
             - (first.y - origin.y)*(second.x - origin.x);
    };
    const double crossTolerance = geometry_area_tolerance(m_geometryScale);
    std::vector<ProjectedVertex> hull;
    hull.reserve(unique.size() + 1);
    for (const ProjectedVertex &vertex : unique)
    {
        while (hull.size() >= 2
               && cross(hull[hull.size()-2], hull.back(), vertex)
                      <= crossTolerance)
        {
            hull.pop_back();
        }
        hull.push_back(vertex);
    }
    const size_t lowerSize = hull.size();
    for (int i = static_cast<int>(unique.size()) - 2; i >= 0; --i)
    {
        const ProjectedVertex &vertex = unique[i];
        while (hull.size() > lowerSize
               && cross(hull[hull.size()-2], hull.back(), vertex)
                      <= crossTolerance)
        {
            hull.pop_back();
        }
        hull.push_back(vertex);
    }
    if (hull.size() > 1)
        hull.pop_back();
    if (hull.size() < MIN_VERTEX_NUM)
        return;

    Point3f hullPoints[MAX_VERTEX_NUM];
    for (size_t i = 0; i < hull.size(); ++i)
        hullPoints[i] = hull[i].point;
    SetOutputPolygon(hullPoints, static_cast<int>(hull.size()), polygon);
    if (polygon.nVertices >= MIN_VERTEX_NUM
        && DotProduct(polygon.Normal(), normal) < 0.0)
    {
        Polygon::InverseVertexOrder(polygon);
    }
}

/** NOTE: Result beams are ordered in inverse direction */
void Scattering::SetPolygonByFacet(int facetId, Polygon &polygon) const
{
    const Polygon &facet = m_facets[facetId];
    int size = facet.nVertices;
    polygon.nVertices = size;
    --size;

    for (int i = 0; i <= size; ++i)
    {
        polygon.arr[i] = facet.arr[size-i];
    }
}

double Scattering::GetIncedentEnergy() const
{
    return m_incidentEnergy;
}

double Scattering::ComputeInternalOpticalPath(const Beam &beam,
                                              const Point3f sourcePoint,
                                              const vector<int> &track)
{
    double path = 0;
    Point3f dir = -beam.direction; // back direction
    Location loc = Location::Out;
    Location nextLoc;

    Point3f p1 = sourcePoint;
    Point3f p2;

    for (int i = track.size()-1; i > 0; --i)
    {
        nextLoc = beam.GetLocationByActNumber(i-1);

        Point3f &exNormal = m_facets[track[i]].ex_normal;
        dir = splitting.ChangeBeamDirection(dir, exNormal, loc, nextLoc);

        Point3f &inNormal = m_facets[track[i-1]].in_normal;
        if (!ProjectPointToPlane(p1, dir, inNormal, p2))
        {
            // A valid traced segment cannot be parallel to the previous
            // interface.  Stop the absorption-path reconstruction instead of
            // injecting an infinite point and NaNs into every Mueller element.
            return path;
        }
        double len = Length(p2 - p1);

        // Natalia_PO uses the geometric internal segment length here. The
        // effective-index scaling was tested, but it changes the back cone
        // interference pattern and does not match Natalia's active code path.

//#ifdef _DEBUG // DEB
//        Point3f dddd = inNormal;
//        dddd.d_param = -dddd.d_param;
//        Point3f p22 = ProjectPointToPlane(p1, dir, dddd);
//        double len1 = Length(p1 - p22);
//        len1 *= sqrt(real(splitting.GetRi()));
//        path1 += len1;
//#endif
        path += len;

        p1 = p2;
        loc = nextLoc;
    }

//#ifdef _DEBUG // DEB
//	path *= real(splitting.GetRi());
//    Point3f nFar1 = m_incidentDir;
//    Point3f nFar2 = -beam.direction;
//    double dd1 = splitting.FAR_ZONE_DISTANCE + DotProductD(p2, nFar1);
//    double dd2 = fabs(DotProductD(sourcePoint, nFar2) + splitting.FAR_ZONE_DISTANCE);
//    path += dd1;
//    path += dd2;
//	if (fabs(path - beam.opticalPath) > 1)
//		int ff = 0;
//#endif
    return path;
}
