#include "Scattering.h"

#include <float.h>
#include <assert.h>
#include <algorithm>
#include <cstdlib>

#include "geometry_lib.h"

#define NORM_CEIL	FLT_EPSILON + 1

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

Scattering::Scattering(Particle *particle, Light *incidentLight, bool isOpticalPath,
                       int nActs)
    : m_particle(particle),
      splitting(isOpticalPath),
      m_incidentLight(incidentLight),
      m_nActs(nActs)
{
    m_facets = m_particle->facets;

    m_incidentDir = m_incidentLight->direction;
    m_incidentDir.d_param = m_incidentLight->direction.d_param;

    m_polarBasis = m_incidentLight->polarizationBasis;

    ::complex ri = m_particle->GetRefractiveIndex();
    splitting.ComputeRiParams(ri);
    const double riContrast = abs((ri - ::complex(1.0, 0.0)) / (ri + ::complex(1.0, 0.0)));
    EPS_BEAM_ENERGY = (0.25*0.25*m_particle->Area())/M_PI*riContrast*riContrast*1e-7;
}

Scattering::~Scattering()
{
}

bool Scattering::EnsureBeamTree()
{
    if (m_beamTree.capacity() == 0)
        m_beamTree.reserve(4096);
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
    m_gpuTracePrefilter = source.m_gpuTracePrefilter;
}

IdType Scattering::Scattering::RecomputeTrackId(const IdType &oldId, int facetId)
{
    if (DisableTraceTrackIds())
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
        return false;
    }
    if (m_traceMaxBeams > 0 && m_treeSize >= m_traceMaxBeams)
    {
        return false;
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

void Scattering::ComputeFacetEnergy(int facetId, const Polygon &lightedPolygon)
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
                  if (fabs(a.x - b.x) > 1e-7)
                      return a.x < b.x;
                  return a.y < b.y;
              });

    std::vector<Point2f> hull;
    hull.reserve(vertices.size());

    auto cross = [](const Point2f &o, const Point2f &a, const Point2f &b)
    {
        return (a.x - o.x)*(b.y - o.y) - (a.y - o.y)*(b.x - o.x);
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
        orderedPolygon.AddVertex(Point3f(p.x, p.y, 10000));
}

void Scattering::ProjectParticleToXY(std::vector<Point2f> &projected)
{
    Point3f n(0, 0, 1, 10000);
    n.d_param = m_particle->MaximalDimention();

    for (int i = 0; i < m_particle->nFacets; i++)
    {
        auto &f = m_particle->facets[i];

        if (DotProduct(f.in_normal, -n) < EPS_COS_90)
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
    for (int i = 0; i < projected.size(); ++i)
    {
        bool isUnique = true;

        for (int j = i + 1; j < projected.size(); ++j)
        {
            if (projected[i].IsEqualTo(projected[j], 0.0001))
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
    shadowBeam.opticalPath = 2.0 * splitting.FAR_ZONE_DISTANCE;
    shadowBeam.lastFacetId = __INT_MAX__;
    scaterredBeams.push_back(shadowBeam);
}

bool Scattering::IsTerminalAct(const Beam &beam)
{
    return (beam.nActs >= m_nActs) || (beam.J.Norm() < EPS_BEAM_ENERGY);
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

    double jn = 0;
    double ar = 0;
    if (useJ || useImportance)
        jn = beam.J.Norm();
    if (useArea || useImportance)
        ar = beam.Area();

    return (useJ && jn < m_traceCutoffJRel * m_traceRefJNorm)
        || (useArea && ar < m_traceCutoffAreaRel * m_traceRefArea)
        || (useImportance && jn * ar < m_traceCutoffImportanceRel * m_traceRefImportance);
}

void Scattering::Difference(const Polygon &subject, const Vector3f &subjNormal,
                         const Polygon &clip, const Vector3f &clipNormal,
                         const Vector3f &clipDir, PolygonArray &difference) const
{
    __m128 _clip[MAX_VERTEX_NUM];
    bool isProjected = ProjectToFacetPlane(clip, clipDir, subjNormal, _clip);

    if (!isProjected)
    {
        difference.Push(subject);
        return;
    }

    __m128 _clip_normal = _mm_setr_ps(clipNormal.cx, clipNormal.cy, clipNormal.cz, 0.0);

    int clipSize = clip.nVertices;
    __m128 _diff_pol[MAX_VERTEX_NUM];

    __m128 _subject[MAX_VERTEX_NUM];
    __m128 _buffer[MAX_VERTEX_NUM];

    for (int i = 0; i < subject.nVertices; ++i)
    {
        _subject[i] = _mm_load_ps(subject.arr[i].coordinates);
    }

    __m128 *_subj = _buffer;
    __m128 *_buff = _subject;
    int bufSize = subject.nVertices;

    __m128 _first_p, _second_p;
    bool isInFirst, isInSecond;

    __m128 _p2 = _clip[clipSize-1];

    for (int i = 0; i < clip.nVertices; ++i)
    {
        int difSize = 0;

        __m128 _p1 = _p2;
        _p2 = _clip[i];

        __m128 *_tmp = _buff;
        _buff = _subj;
        _subj = _tmp;

        int subjSize = bufSize;
        bufSize = 0;

        _first_p = _subj[subjSize-1];
        isInFirst = is_inside_i(_first_p, _p1, _p2, _clip_normal);

        bool isIntersected;

        for (int j = 0; j < subjSize; ++j)
        {
            _second_p = _subj[j];
            isInSecond = is_inside_i(_second_p, _p1, _p2, _clip_normal);

            if (isInSecond)
            {
                if (!isInFirst)
                {
                    __m128 x = intersect_i(_first_p, _second_p, _p1, _p2,
                                           _clip_normal, isIntersected);

                    if (isIntersected && is_layOnLine_i(x, _first_p, _second_p))
                    {
                        _diff_pol[difSize++] = x;
                        _buff[bufSize++] = x;
                    }
                }

                _buff[bufSize++] = _second_p;
            }
            else
            {
                if (isInFirst)
                {
                    __m128 x = intersect_i(_first_p, _second_p, _p1, _p2,
                                           _clip_normal, isIntersected);

                    if (isIntersected && is_layOnLine_i(x, _first_p, _second_p))
                    {
                        _diff_pol[difSize++] = x;
                        _buff[bufSize++] = x;
                    }
                }

                _diff_pol[difSize++] = _second_p;
            }

            _first_p = _second_p;
            isInFirst = isInSecond;
        }

        if (difSize >= MIN_VERTEX_NUM)
        {
            Polygon resPolygon;
            SetOutputPolygon(_diff_pol, difSize, resPolygon);

            if (resPolygon.nVertices >= MIN_VERTEX_NUM)
            {
//                if (resPolygon.nVertices > 9)
//                {
//                    std::cout << std::endl << resPolygon.nVertices;
//                }

                difference.Push(resPolygon);
            }
        }
    }
}

bool Scattering::ProjectToFacetPlane(const Polygon &polygon, const Vector3f &dir,
                                  const Point3f &normal, __m128 *_projection) const
{
    __m128 _normal = _mm_load_ps(normal.coordinates);
    __m128 _direction = _mm_load_ps(dir.coordinates);

    __m128 _d_param = _mm_set_ps1(normal.d_param);
    __m128 _dp0 = _mm_dp_ps(_direction, _normal, MASK_FULL);

    __m128 _sign_mask = _mm_set1_ps(-0.f);
    __m128 _abs_dp = _mm_andnot_ps(_sign_mask, _dp0);

    if (_abs_dp[0] < EPS_PROJECTION)
    {
        return false; /// beam is parallel to facet
    }

    for (int i = 0; i < polygon.nVertices; ++i)
    {
        const Point3f &p = polygon.arr[i];
        __m128 _point = _mm_load_ps(p.coordinates);
        __m128 _dp1 = _mm_dp_ps(_point, _normal, MASK_FULL);
        __m128 add = _mm_add_ps(_dp1, _d_param);
        __m128 _t = _mm_div_ps(add, _dp0);
        __m128 mul = _mm_mul_ps(_t, _direction);

        _projection[i] = _mm_sub_ps(_point, mul);
    }

    return true;
}

/// NOTE: вершины пучка и грани должны быть ориентированы в одном направлении
void Scattering::Intersect(int facetID, const Beam &beam, Polygon &intersection) const
{
    __m128 _output_points[MAX_VERTEX_NUM];
    const Facet &facet = m_particle->facets[facetID];
    // REF: перенести в случай невыпуклых частиц
    const Point3f &normal = facet.in_normal;

    const Point3f &normal1 = (beam.location == Location::In) ? facet.in_normal
                                                            : facet.ex_normal;
    bool isProjected = ProjectToFacetPlane(beam, beam.direction, normal1,
                                           _output_points);
    if (!isProjected)
    {
        return;
    }

    __m128 _normal_to_facet = _mm_setr_ps(-normal.cx, -normal.cy, -normal.cz, 0.0);
    __m128 *_output_ptr = _output_points;
    int outputSize = beam.nVertices;

    __m128 _buffer[MAX_VERTEX_NUM];
    __m128 *_buffer_ptr = _buffer;
    int bufferSize;

    int facetSize = facet.nVertices;

    __m128 _p1, _p2; // vertices of facet
    __m128 _s_point, _e_point;	// points of projection
    bool isInsideE, isInsideS;

    Point3f p2 = facet.arr[facetSize-1];
    _p2 = _mm_load_ps(p2.coordinates);

    for (int i = 0; i < facetSize; ++i)
    {
        _p1 = _p2;
        p2 = facet.arr[i];
        _p2 = _mm_load_ps(p2.coordinates);

        bufferSize = outputSize;
        outputSize = 0;

        __m128 *_temp = _output_ptr;
        _output_ptr = _buffer_ptr;
        _buffer_ptr = _temp;

        _s_point = _buffer_ptr[bufferSize-1];
        isInsideS = is_inside_i(_s_point, _p1, _p2, _normal_to_facet);

        bool isIntersected;

        for (int j = 0; j < bufferSize; ++j)
        {
            _e_point = _buffer_ptr[j];
            isInsideE = is_inside_i(_e_point, _p1, _p2, _normal_to_facet);

            if (isInsideE)
            {
                if (!isInsideS)
                {
                    __m128 x = intersect_i(_s_point, _e_point, _p1, _p2,
                                           _normal_to_facet, isIntersected);
                    if (isIntersected)
                    {
                        _output_ptr[outputSize++] = x;
                    }
                }

                _output_ptr[outputSize++] = _e_point;
            }
            else if (isInsideS)
            {
                __m128 x = intersect_i(_s_point, _e_point, _p1, _p2,
                                       _normal_to_facet, isIntersected);
                if (isIntersected)
                {
                    _output_ptr[outputSize++] = x;
                }
            }

            _s_point = _e_point;
            isInsideS = isInsideE;
        }
    }

    SetOutputPolygon(_output_ptr, outputSize, intersection);
}

void Scattering::SetOutputPolygon(__m128 *_output_points, int outputSize,
                                  Polygon &polygon) const
{
    Point3f p;

    __m128 eps = _mm_load_ps1(&EPS_MERGE);
    __m128 sign_mask = _mm_set1_ps(-0.f);

    __m128 p0 = _output_points[outputSize-1];

    for (int i = 0; i < outputSize; ++i)
    {
        __m128 difference = _mm_sub_ps(_output_points[i], p0);
        __m128 abs = _mm_andnot_ps(sign_mask, difference);
        __m128 cmp = _mm_cmplt_ps(eps, abs);

        int res = _mm_movemask_ps(cmp) & 0b111;

        if (res != 0)
        {
            p.cx = _output_points[i][0];
            p.cy = _output_points[i][1];
            p.cz = _output_points[i][2];
            polygon.arr[polygon.nVertices++] = p;
        }

        p0 = _output_points[i];
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
    double path1 = 0;
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
        p2 = ProjectPointToPlane(p1, dir, inNormal);
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
