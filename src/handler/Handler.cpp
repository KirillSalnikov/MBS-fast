#include "Handler.h"
#include "HandlerPO_fast.h"

namespace
{
bool HasRequestedTrackIds(const Tracks *tracks)
{
    if (!tracks)
        return false;
    for (const TrackGroup &group : *tracks)
    {
        if (group.size != 0)
            return true;
    }
    return false;
}
}
#include "Intersection.h"

#include "Mueller.hpp"
#include <iostream>
#include <algorithm>
#ifdef _OPENMP
#include <omp.h>
#endif

using namespace std;
using ::complex;

namespace
{
complex StableComplexSinhc(const complex &z)
{
    if (abs(z) < 1e-4)
    {
        const complex z2 = z*z;
        return 1.0 + z2/6.0 + z2*z2/120.0;
    }

    const double a = real(z);
    const double b = imag(z);
    return complex(std::sinh(a)*std::cos(b),
                   std::cosh(a)*std::sin(b)) / z;
}

complex StableComplexEdgeQuotient(const complex &coefficient,
                                  double first, double second,
                                  double waveIndex)
{
    const double delta = second - first;
    const double middle = 0.5*(first + second);
    const complex halfDeltaPhase(-0.5*waveIndex*imag(coefficient)*delta,
                                  0.5*waveIndex*real(coefficient)*delta);
    const complex middlePhase(-waveIndex*imag(coefficient)*middle,
                               waveIndex*real(coefficient)*middle);
    return exp(middlePhase) * complex(0.0, waveIndex*delta)
         * StableComplexSinhc(halfDeltaPhase);
}

bool SmallComplexPhasePolygonFactor(const double *vx, const double *vy, int nv,
                                    const complex &a, const complex &b,
                                    double waveIndex, complex &factor)
{
    if (nv < MIN_VERTEX_NUM)
        return false;

    const complex qx(-waveIndex*imag(a), waveIndex*real(a));
    const complex qy(-waveIndex*imag(b), waveIndex*real(b));
    double maxPhase = 0.0;
    for (int i = 0; i < nv; ++i)
        maxPhase = std::max(maxPhase, abs(qx*vx[i] + qy*vy[i]));
    if (maxPhase > 1e-3)
        return false;

    double twiceArea = 0.0;
    double firstX6 = 0.0, firstY6 = 0.0;
    double secondX12 = 0.0, secondY12 = 0.0, secondXY24 = 0.0;
    double minX = DBL_MAX, maxX = -DBL_MAX;
    double minY = DBL_MAX, maxY = -DBL_MAX;
    for (int i = 0; i < nv; ++i)
    {
        const int j = (i + 1 < nv) ? i + 1 : 0;
        const double xi = vx[i], yi = vy[i];
        const double xj = vx[j], yj = vy[j];
        const double cross = xi*yj - xj*yi;
        twiceArea += cross;
        firstX6 += (xi + xj)*cross;
        firstY6 += (yi + yj)*cross;
        secondX12 += (xi*xi + xi*xj + xj*xj)*cross;
        secondY12 += (yi*yi + yi*yj + yj*yj)*cross;
        secondXY24 += (2.0*xi*yi + xi*yj + xj*yi + 2.0*xj*yj)*cross;
        minX = std::min(minX, xi); maxX = std::max(maxX, xi);
        minY = std::min(minY, yi); maxY = std::max(maxY, yi);
    }
    const double areaScale = std::max(
        (maxX - minX)*(maxY - minY), DBL_MIN);
    if (std::fabs(twiceArea) <= 0x1p-40*areaScale)
    {
        factor = 1.0;
        return true;
    }

    const double invArea = 2.0/twiceArea;
    const double meanX = firstX6*invArea/6.0;
    const double meanY = firstY6*invArea/6.0;
    const double meanXX = secondX12*invArea/12.0;
    const double meanYY = secondY12*invArea/12.0;
    const double meanXY = secondXY24*invArea/24.0;
    const complex meanExponent = qx*meanX + qy*meanY;
    const complex meanExponent2 = qx*qx*meanXX
                                + 2.0*qx*qy*meanXY + qy*qy*meanYY;
    factor = 1.0 + meanExponent + 0.5*meanExponent2;
    return true;
}
}

Handler::Handler(Particle *particle, Light *incidentLight, int nTheta,
                 double wavelength)
    : m_incidentLight(incidentLight),
      nTheta(nTheta),
      m_sphere(0.0, 0.0, 0, 0),
      m_particle(particle),
      m_geometryScale(particle->MaximalDimention()),
      m_wavelength(wavelength),
      m_hasAbsorption(false),
      m_normIndex(1)
{
//	m_wavelength = 0.532;
    m_waveIndex = M_2PI/m_wavelength;
    m_wi2 = m_waveIndex*m_waveIndex;

    ::complex one(0, -1);
    m_complWave = (one * m_wavelength) / SQR(M_2PI);
    m_invComplWave = -one/m_wavelength;

    // Coefficient threshold for the analytic zero-denominator limit. Physical
    // edge lengths use geometry_length_tolerance() separately below.
    m_eps1 = 1e9*DBL_EPSILON;
    m_eps2 = 1e6*DBL_EPSILON;
    m_eps3 = 1e1;

}

void Handler::HandleBeams(std::vector<Beam> &/*beams*/, double /*sinZenith*/)
{
}

void Handler::SetTracks(Tracks *tracks)
{
//	if (!m_tracks)
//	{
//		std::cerr << "Tracks are not set" << std::endl;
//		throw std::exception();
//	}

    m_tracks = tracks;
    if (m_scattering)
        m_scattering->SetTrackIdsRequired(
            m_hasAbsorption || HasRequestedTrackIds(m_tracks));
}

void Handler::WriteMatricesToFile(string &/*destName*/, double /*nrg*/)
{
}

void Handler::WriteTotalMatricesToFile(const string &/*destName*/)
{

}

void Handler::SetNormIndex(double value)
{
    m_normIndex = value;
}

void Handler::SetSinZenith(double value)
{
    m_sinZenith = value;
}

void Handler::SetAbsorptionAccounting(bool value)
{
    m_hasAbsorption = value;
    if (m_scattering)
        m_scattering->SetTrackIdsRequired(
            m_hasAbsorption || HasRequestedTrackIds(m_tracks));
    m_ri = m_particle->GetRefractiveIndex();
    m_riIm = imag(m_ri);
    m_cAbs = -m_riIm*m_waveIndex;
    m_absMag = -m_waveIndex*m_riIm;
}

void Handler::SetAbsorptionPointCount(int value)
{
    m_absorptionPointCount = (value == -1 || value > 0) ? value : 1;
}

void Handler::SetScatteringSphere(const ScatteringRange &/*grid*/)
{

}

void Handler::SetScattering(Scattering *scattering)
{
    m_scattering = scattering;
    if (m_scattering)
        m_scattering->SetTrackIdsRequired(
            m_hasAbsorption || HasRequestedTrackIds(m_tracks));
}

void Handler::ExtropolateOpticalLenght(Beam &beam, const std::vector<int> &tr)
{
    std::vector<double> lengths;

    for (int i = 0; i < beam.nVertices; ++i)
    {
        double d = m_scattering->ComputeInternalOpticalPath(
                    beam, beam.arr[i], tr);
        lengths.push_back(d);
    }

//	Point3d lens = ComputeLengthIndices(beam, info);

//	for (int i = 3; i < beam.nVertices; ++i)
//	{
//		Point3d newP = ChangeCoordinateSystem(hor, ver, n, beam.arr[i]) - center;
//		double newL = lens.z + lens.x*newP.x + lens.y*newP.y;
//		double err = fabs((lengths[i] - newL)/lengths[i])*100;
//		if (err > 5)
//			m_logFile << Polygon(beam) << "Area: " << beam.Area() << ' '
//					  << i << ", " << "Error: " << err << std::endl;
//	}
}

void Handler::ApplyAbsorption(Beam &beam)
{
    vector<int> tr;
    Tracks::RecoverTrack(beam, m_particle->nFacets, tr);

    double absorption = 0.0;
    int nSamples = 0;

    auto addAbsorptionAt = [&](const Point3f &point)
    {
        double path = m_scattering->ComputeInternalOpticalPath(beam, point, tr);
        absorption += (path > DBL_EPSILON) ? exp(m_cAbs*path) : 1.0;
        ++nSamples;

#ifdef _DEBUG // DEB
        double ddd = fabs(path - beam.opticalPath);
//	m_logFile << ddd << endl;
        if (fabs(path - beam.opticalPath) >= 10e-4)
            int ggg = 0;
#endif
    };

    if (m_absorptionPointCount == 1 || beam.nVertices <= 0)
    {
        addAbsorptionAt(beam.Center());
    }
    else
    {
        int count = (m_absorptionPointCount == -1)
            ? beam.nVertices
            : std::min(m_absorptionPointCount, beam.nVertices);

        for (int i = 0; i < count; ++i)
        {
            int idx = (i * beam.nVertices) / count;
            addAbsorptionAt(beam.arr[idx]);
        }
    }

    if (nSamples > 0)
        beam.J *= absorption / nSamples;
}

Point3d Handler::ChangeCoordinateSystem(const Point3d& hor, const Point3d& ver,
                                        const Point3d& normal,
                                        const Point3d& point) const
{
    // расчёт коор-т в СК наблюдателя
    const Point3d p_pr = point - normal*DotProductD(normal, point);

    return Point3d(DotProductD(p_pr, hor), DotProductD(p_pr, ver), 0);
}

void Handler::ComputeCoordinateSystemAxes(const Point3d& normal,
                                          Point3d &hor, Point3d &ver) const
{
    const double xyNorm = hypot(normal.x, normal.y);
    if (xyNorm <= 64.0*DBL_EPSILON)
    {
        hor = Point3d(0, -normal.z, 0);
        ver = Point3d(1, 0, 0);
    }
    else
    {
        hor = Point3d(normal.y/xyNorm, -normal.x/xyNorm, 0);
        ver = CrossProductD(normal, hor);
    }
}

void Handler::ComputeLengthIndices(const Beam &beam, BeamInfo &info)
{
    auto *lens = info.opticalLengths;

    Point3d p1 = ChangeCoordinateSystem(info.horAxis, info.verAxis, info.normald,
                                        beam.arr[0]) - info.projectedCenter;
    Point3d p2 = ChangeCoordinateSystem(info.horAxis, info.verAxis, info.normald,
                                        beam.arr[1]) - info.projectedCenter;
    Point3d p3 = ChangeCoordinateSystem(info.horAxis, info.verAxis, info.normald,
                                        beam.arr[2]) - info.projectedCenter;

//	double den = p1.x*(p2.y - p3.y) - p2.x*(p1.y + p3.y) + p3.x*(p1.y - p2.y);
    double den = (p1.x*p2.y-p1.x*p3.y-p2.x*p1.y+p2.x*p3.y+p3.x*p1.y-p3.x*p2.y);

    const double coordinateScale = std::max({
        std::fabs(p1.x), std::fabs(p1.y),
        std::fabs(p2.x), std::fabs(p2.y),
        std::fabs(p3.x), std::fabs(p3.y)});
    const double determinantTolerance = geometry_area_tolerance(coordinateScale);

    if (std::isfinite(den) && fabs(den) > determinantTolerance)
    {
#ifdef _DEBUG // DEB
        double aa = lens[0]*(p2.x*p3.y - p3.x*p2.y);
        double bb = lens[1]*(p1.x*p3.y + p3.x*p1.y);
        double cc = lens[2]*(p1.x*p2.y - p2.x*p1.y);
#endif
        info.lenIndices.z = (lens[0]*p2.x*p3.y-lens[0]*p3.x*p2.y-
                lens[1]*p1.x*p3.y+lens[1]*p3.x*p1.y+
                lens[2]*p1.x*p2.y-lens[2]*p2.x*p1.y) / den;

        info.lenIndices.x = (lens[0]*p2.y - lens[0]*p3.y -
                lens[1]*p1.y + lens[1]*p3.y +
                lens[2]*p1.y - lens[2]*p2.y) / den;

        info.lenIndices.y = -(lens[0]*p2.x - lens[0]*p3.x -
                lens[1]*p1.x + lens[1]*p3.x +
                lens[2]*p1.x - lens[2]*p2.x) / den;

#ifdef _DEBUG // DEB
//        if (fabs(info.lenIndices.z) > 700000)
//            int fff = 0;
#endif
    }
    else
    {
        info.isBad = true;
        m_isBadBeam = true;
        ++m_nBadBeams;
    }
}

Tracks *Handler::GetTracks() const
{
    return m_tracks;
}

Point3d Handler::ChangeCoordinateSystem(const Point3d& normal,
                                        const Point3d &point) const
{
    Point3d hor;  // условная горизонталь СК экрана в СК тела
    Point3d ver;  // третья ось (условная вертикаль СК экрана)
    ComputeCoordinateSystemAxes(normal, hor, ver);
    return ChangeCoordinateSystem(hor, ver, normal, point);
}

::complex Handler::DiffractInclineAbs(const BeamInfo &info, const Beam &beam,
                             const Point3d &direction) const
{
    if (beam.nVertices < MIN_VERTEX_NUM || beam.nVertices > MAX_VERTEX_NUM)
        return complex(0, 0);

    const Point3f &dir = beam.direction;
    Point3d k_k0 = -direction + Point3d(dir.cx, dir.cy, dir.cz);

    // Project k_k0 onto aperture plane via dot products (avoid ChangeCoordinateSystem)
    double pt_x = DotProductD(k_k0, info.horAxis)
                - DotProductD(info.normald, k_k0) * DotProductD(info.normald, info.horAxis);
    double pt_y = DotProductD(k_k0, info.verAxis)
                - DotProductD(info.normald, k_k0) * DotProductD(info.normald, info.verAxis);

    const ::complex A(pt_x, info.lenIndices.x*m_riIm);
    const ::complex B(pt_y, info.lenIndices.y*m_riIm);

    double vx[MAX_VERTEX_NUM] = {};
    double vy[MAX_VERTEX_NUM] = {};
    double coordinateScale = 0.0;
    for (int i = 0; i < beam.nVertices; ++i)
    {
        const Point3f &vertex = beam.arr[i];
        const Point3d p(vertex.cx, vertex.cy, vertex.cz);
        const Point3d projected = p - info.normald*DotProductD(info.normald, p);
        vx[i] = DotProductD(projected, info.horAxis) - info.projectedCenter.x;
        vy[i] = DotProductD(projected, info.verAxis) - info.projectedCenter.y;
        coordinateScale = std::max(coordinateScale,
                                   std::max(std::fabs(vx[i]), std::fabs(vy[i])));
    }

    complex smallPhaseFactor;
    if (SmallComplexPhasePolygonFactor(vx, vy, beam.nVertices, A, B,
                                       m_waveIndex, smallPhaseFactor))
    {
        const complex areaLimit =
            (m_legacySign ? m_invComplWave : -m_invComplWave)*info.area;
        double absorp = exp(m_absMag*info.lenIndices.z);
        return areaLimit*smallPhaseFactor*absorp;
    }

    ::complex s(0, 0);
    double p1x = vx[beam.nVertices-1];
    double p1y = vy[beam.nVertices-1];
    const double edgeTolerance = geometry_length_tolerance(coordinateScale);

    if (abs(B) > abs(A))
    {
        for (int i = 0; i < beam.nVertices; ++i)
        {
            const double p2x = vx[i];
            const double p2y = vy[i];

            double dx = p1x - p2x;
            if (fabs(dx) <= edgeTolerance)
            {
                p1x = p2x;
                p1y = p2y;
                continue;
            }

            const double ai = (p1y - p2y)/dx;
            const ::complex Ci = A + ai*B;
            const ::complex tmp = StableComplexEdgeQuotient(
                Ci, p1x, p2x, m_waveIndex);

            const double bi = p1y - ai*p1x;
            const ::complex phaseB(-m_waveIndex*imag(B)*bi,
                                    m_waveIndex*real(B)*bi);
            const ::complex contribution = exp(phaseB)*tmp;

            if (!std::isfinite(real(contribution))
                || !std::isfinite(imag(contribution)))
            {
                p1x = p2x; p1y = p2y;
                continue;
            }

            s += contribution;

            p1x = p2x; p1y = p2y;
        }

        s /= B;
    }
    else
    {
        for (int i = 0; i < beam.nVertices; ++i)
        {
            const double p2x = vx[i];
            const double p2y = vy[i];

            double dy = p1y - p2y;
            if (fabs(dy) <= edgeTolerance)
            {
                p1x = p2x;
                p1y = p2y;
                continue;
            }

            const double ci = (p1x - p2x)/dy;
            const ::complex Ei = A*ci + B;
            const ::complex tmp = StableComplexEdgeQuotient(
                Ei, p1y, p2y, m_waveIndex);

            const double di = p1x - ci*p1y;
            const ::complex phaseA(-m_waveIndex*imag(A)*di,
                                    m_waveIndex*real(A)*di);
            const ::complex contribution = exp(phaseA)*tmp;

            if (!std::isfinite(real(contribution))
                || !std::isfinite(imag(contribution)))
            {
                p1x = p2x; p1y = p2y;
                continue;
            }

            s += contribution;
            p1x = p2x; p1y = p2y;
        }

        s /= -A;
    }

    double absorp = exp(m_absMag*info.lenIndices.z);
    return m_complWave * s * absorp;
}

double Handler::BeamCrossSection(const Beam &beam) const
{
    const double eps = 1e7*DBL_EPSILON;

    Point3f normal = m_particle->facets[beam.lastFacetId].ex_normal; // normal of last facet of beam
    double cosA = DotProduct(normal, beam.direction);
    double e = fabs(cosA);

    if (e < eps)
    {
        return 0;
    }

    double area = beam.Area();
    double len = Length(normal);
    return (e*area)/len;
}

::complex Handler::DiffractIncline(const BeamInfo &info, const Beam &beam,
                                 const Point3d &direction) const
{
    const Point3f &dir = beam.direction;
    Point3d k_k0 = -direction + Point3d(dir.cx, dir.cy, dir.cz);

    Point3d	pt_proj = ChangeCoordinateSystem(info.horAxis, info.verAxis,
                                             info.normald, k_k0);
    const double A = pt_proj.x;
    const double B = pt_proj.y;

    double absA = fabs(A);
    double absB = fabs(B);
#ifdef _DEBUG // DEB
    if (beam.id == 423)
        int fff = 0;
#endif
    // This is the slow fallback for polygons that do not fit BeamEdgeData.
    // Use the same dimensionless small-phase limit as the regular fast path.
    std::vector<double> vx(beam.nVertices), vy(beam.nVertices);
    for (int i = 0; i < beam.nVertices; ++i)
    {
        const Point3d p = ChangeCoordinateSystem(
            info.horAxis, info.verAxis, info.normald, beam.arr[i])
            - info.projectedCenter;
        vx[i] = p.x;
        vy[i] = p.y;
    }
    double smallPhaseReal = 0.0, smallPhaseImag = 0.0;
    if (small_phase_polygon_factor(
            vx.data(), vy.data(), beam.nVertices, A, B, m_waveIndex,
            smallPhaseReal, smallPhaseImag))
    {
        const ::complex areaLimit =
            (m_legacySign ? m_invComplWave : -m_invComplWave) * info.area;
        return apply_small_phase_factor(
            areaLimit, smallPhaseReal, smallPhaseImag);
    }

    ::complex s(0, 0);

    int begin, startIndex, endIndex;

//    if (info.order)
    {
//        begin = 0;
//        startIndex = beam.nVertices-1;
//        endIndex = -1;
    }
//    else
    {
        begin = beam.nVertices-1;
        startIndex = 0;
        endIndex = beam.nVertices;
    }

    Point3d p1 = ChangeCoordinateSystem(info.horAxis, info.verAxis,
                                        info.normald, beam.arr[begin]) - info.projectedCenter;
    Point3d p2;

    const double edgeTolerance = geometry_length_tolerance(
        m_geometryScale);
    if (absB > absA)
    {
        for (int i = startIndex; i != endIndex; i = /*info.order ? i-1 :*/ i+1)
        {
            p2 = ChangeCoordinateSystem(info.horAxis, info.verAxis,
                                        info.normald, beam.arr[i]) - info.projectedCenter;

            if (fabs(p1.x - p2.x) < edgeTolerance)
            {
                p1 = p2;
                continue;
            }

            const double ai = (p1.y - p2.y)/(p1.x - p2.x);
            const double Ci = A+ai*B;

            ::complex tmp;

            if (fabs(Ci) < m_eps1)
            {
                tmp = StableComplexEdgeQuotient(
                    ::complex(Ci, 0.0), p1.x, p2.x, m_waveIndex);
            }
            else
            {
                double kCi = m_waveIndex*Ci;
                tmp = (exp_im(kCi*p2.x) - exp_im(kCi*p1.x))/Ci;
            }

            const double bi = p1.y - ai*p1.x;
            s += exp_im(m_waveIndex*B*bi) * tmp;

            p1 = p2;
        }

        s /= B;
    }
    else
    {
        for (int i = startIndex; i != endIndex; i = /*info.order ? i-1 :*/ i+1)
        {
            p2 = ChangeCoordinateSystem(info.horAxis, info.verAxis,
                                        info.normald, beam.arr[i]) - info.projectedCenter;

            if (fabs(p1.y - p2.y) < edgeTolerance)
            {
                p1 = p2;
                continue;
            }

            const double ci = (p1.x - p2.x)/(p1.y - p2.y);
            const double Ei = A*ci+B;

            ::complex tmp;

            if (fabs(Ei) < m_eps1)
            {
                tmp = StableComplexEdgeQuotient(
                    ::complex(Ei, 0.0), p1.y, p2.y, m_waveIndex);
            }
            else
            {
                double kEi = m_waveIndex*Ei;
                tmp = (exp_im(kEi*p2.y) - exp_im(kEi*p1.y))/Ei;
            }

            const double di = p1.x - ci*p1.y;
            s += exp_im(m_waveIndex*A*di) * tmp;

            p1 = p2;
        }

        s /= -A;
    }

    return m_complWave * s;
}

void Handler::PrecomputeEdgeData(const BeamInfo &info, const Beam &beam,
                                  BeamEdgeData &edgeData) const
{
    edgeData.nVertices = 0;
    edgeData.valid = false;
    if (beam.nVertices < 3 || beam.nVertices >= BeamEdgeData::MAX_EDGES)
        return;

    edgeData.nVertices = beam.nVertices;
    edgeData.valid = true;

    int nv = beam.nVertices;

    // Project all vertices into aperture 2D system ONCE
    int begin = nv - 1;
    for (int i = 0; i < nv; ++i)
    {
        int idx = (i == 0) ? begin : (i - 1);
        Point3f v = beam.arr[idx];
        Point3d p(v.cx, v.cy, v.cz);
        Point3d p_pr = p - info.normald * DotProductD(info.normald, p);
        edgeData.x[i] = DotProductD(p_pr, info.horAxis) - info.projectedCenter.x;
        edgeData.y[i] = DotProductD(p_pr, info.verAxis) - info.projectedCenter.y;
    }

    // Precompute per-edge data (slopes, intercepts). Edge validity is a
    // geometric decision and therefore scales with the particle.
    const double edgeTolerance = geometry_length_tolerance(
        m_geometryScale);
    for (int i = 0; i < nv; ++i)
    {
        int inext = (i + 1) % nv;
        double dxi = edgeData.x[inext] - edgeData.x[i];
        double dyi = edgeData.y[inext] - edgeData.y[i];
        edgeData.dx[i] = dxi;
        edgeData.dy[i] = dyi;

        // For absB > absA branch: ai = (p1y-p2y)/(p1x-p2x) where p1=v[i], p2=v[next]
        // ai = (y[i]-y[next])/(x[i]-x[next]) = -dy/(-dx) = dy/dx
        edgeData.edge_valid_x[i] = (fabs(dxi) >= edgeTolerance);
        if (edgeData.edge_valid_x[i])
        {
            edgeData.slope_yx[i] = dyi / dxi;  // NOT negated: (y[next]-y[i])/(x[next]-x[i])
            // but old code uses (p1y-p2y)/(p1x-p2x) = (-dy)/(-dx) = dy/dx. Same!
            edgeData.intercept_y[i] = edgeData.y[i] - edgeData.slope_yx[i] * edgeData.x[i];
        }

        // For absA >= absB branch: ci = (p1x-p2x)/(p1y-p2y) = dx/dy
        edgeData.edge_valid_y[i] = (fabs(dyi) >= edgeTolerance);
        if (edgeData.edge_valid_y[i])
        {
            edgeData.slope_xy[i] = dxi / dyi;
            edgeData.intercept_x[i] = edgeData.x[i] - edgeData.slope_xy[i] * edgeData.y[i];
        }
    }
}

::complex Handler::DiffractInclineFast(const BeamInfo &info, const BeamEdgeData &ed,
                                      const Point3d &beamDir,
                                      const Point3d &direction) const
{
    // Compute A, B directly via dot products (no ChangeCoordinateSystem)
    Point3d k_k0(-direction.x + beamDir.x,
                 -direction.y + beamDir.y,
                 -direction.z + beamDir.z);
    const double A = DotProductD(k_k0, info.horAxis)
                   - DotProductD(info.normald, k_k0) * DotProductD(info.normald, info.horAxis);
    const double B = DotProductD(k_k0, info.verAxis)
                   - DotProductD(info.normald, k_k0) * DotProductD(info.normald, info.verAxis);

    double absA = fabs(A);
    double absB = fabs(B);

    double smallPhaseReal = 0.0, smallPhaseImag = 0.0;
    if (small_phase_polygon_factor(
            ed.x, ed.y, ed.nVertices, A, B, m_waveIndex,
            smallPhaseReal, smallPhaseImag))
    {
        const ::complex areaLimit =
            (m_legacySign ? m_invComplWave : -m_invComplWave) * info.area;
        return apply_small_phase_factor(
            areaLimit, smallPhaseReal, smallPhaseImag);
    }

    ::complex s(0, 0);
    const int nv = ed.nVertices;

    // Use precomputed 2D vertices - no ChangeCoordinateSystem calls.
    const double edgeTolerance = geometry_length_tolerance(
        m_geometryScale);
    if (absB > absA)
    {
        double p1x = ed.x[0], p1y = ed.y[0];
        for (int i = 1; i <= nv; ++i)
        {
            double p2x = ed.x[i % nv], p2y = ed.y[i % nv];

            double dx = p1x - p2x;
            if (fabs(dx) < edgeTolerance) { p1x = p2x; p1y = p2y; continue; }

            double ai = (p1y - p2y) / dx;
            double Ci = A + ai * B;

            ::complex tmp;
            if (fabs(Ci) < m_eps1)
                tmp = StableComplexEdgeQuotient(
                    ::complex(Ci, 0.0), p1x, p2x, m_waveIndex);
            else
            {
                double kCi = m_waveIndex * Ci;
                tmp = (exp_im(kCi*p2x) - exp_im(kCi*p1x)) / Ci;
            }

            double bi = p1y - ai * p1x;
            s += exp_im(m_waveIndex * B * bi) * tmp;

            p1x = p2x; p1y = p2y;
        }
        s /= B;
    }
    else
    {
        double p1x = ed.x[0], p1y = ed.y[0];
        for (int i = 1; i <= nv; ++i)
        {
            double p2x = ed.x[i % nv], p2y = ed.y[i % nv];

            double dy = p1y - p2y;
            if (fabs(dy) < edgeTolerance) { p1x = p2x; p1y = p2y; continue; }

            double ci = (p1x - p2x) / dy;
            double Ei = A * ci + B;

            ::complex tmp;
            if (fabs(Ei) < m_eps1)
                tmp = StableComplexEdgeQuotient(
                    ::complex(Ei, 0.0), p1y, p2y, m_waveIndex);
            else
            {
                double kEi = m_waveIndex * Ei;
                tmp = (exp_im(kEi*p2y) - exp_im(kEi*p1y)) / Ei;
            }

            double di = p1x - ci * p1y;
            s += exp_im(m_waveIndex * A * di) * tmp;

            p1x = p2x; p1y = p2y;
        }
        s /= -A;
    }

    return m_complWave * s;
}
