#include "Handler.h"
#include "HandlerPO_fast.h"

#include "Mueller.hpp"
#include <iostream>
#include <limits>
#include <iomanip>

using namespace std;

Handler::Handler(Particle *particle, Light *incidentLight, int nTheta,
                 double wavelength)
    : m_incidentLight(incidentLight),
      m_particle(particle),
      m_wavelength(wavelength),
      m_hasAbsorption(false),
      m_normIndex(1),
      m_nBadBeams(0),
      m_sphere(0.0, 0.0, 0, 0),
      nTheta(nTheta)
{
//	m_wavelength = 0.532;
    m_waveIndex = M_2PI/m_wavelength;
    m_wi2 = m_waveIndex*m_waveIndex;

    complex one(0, -1);
    m_complWave = (one * m_wavelength) / SQR(M_2PI);
    m_invComplWave = -one/m_wavelength;

    m_eps1 = 1e9*DBL_EPSILON;
//    m_eps2 = 1e6*DBL_EPSILON;
    m_eps2 = __FLT_EPSILON__*2;
    m_eps3 = 1e1;

    m_logFile.open("log1.txt", std::ios::out);
    m_logFile << std::setprecision(std::numeric_limits<long double>::digits10 + 1);
}

void Handler::HandleBeams(std::vector<Beam> &/*beams*/, double sinZenith)
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
}

void Handler::WriteMatricesToFile(string &/*destName*/, double nrg)
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
    m_ri = m_particle->GetRefractiveIndex();
    m_riIm = imag(m_ri);
    m_cAbs = -m_riIm*m_waveIndex;
    m_absMag = -m_waveIndex*m_riIm;
}

void Handler::SetScatteringSphere(const ScatteringRange &grid)
{

}

void Handler::SetScattering(Scattering *scattering)
{
    m_scattering = scattering;
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

    Vector3f _n = beam.Normal();
    Point3d n = Point3d(_n.cx, _n.cy, _n.cz);

    Point3d hor;
    Point3d ver;
    ComputeCoordinateSystemAxes(n, hor, ver);

    Point3f cntr = beam.Center();
    Point3d center = ChangeCoordinateSystem(hor, ver, n,
                                            Point3d(cntr.cx, cntr.cy, cntr.cz));
    double ls[3];
    ls[0] = lengths[0];
    ls[1] = lengths[1];
    ls[2] = lengths[2];

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

//	double opAbs = CalcOpticalPathAbsorption(beam);
    double path = m_scattering->ComputeInternalOpticalPath(beam, beam.Center(), tr);

#ifdef _DEBUG // DEB
    double ddd = fabs(path - beam.opticalPath);
//	m_logFile << ddd << endl;
    if (fabs(path - beam.opticalPath) >= 10e-4)
        int ggg = 0;
#endif

    if (path > DBL_EPSILON)
    {
        double abs = exp(m_cAbs*path);
        beam.J *= abs;
    }
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
    if (fabs(normal.z) > 1-DBL_EPSILON)
    {
        hor = Point3d(0, -normal.z, 0);
        ver = Point3d(1, 0, 0);
    }
    else
    {
        const double tmp = sqrt(SQR(normal.x) + SQR(normal.y));
        hor = Point3d(normal.y/tmp, -normal.x/tmp, 0);
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

    if (fabs(den) > 1e-3)
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

complex Handler::DiffractInclineAbs(const BeamInfo &info, const Beam &beam,
                             const Point3d &direction) const
{
    const Point3f &dir = beam.direction;
    Point3d k_k0 = -direction + Point3d(dir.cx, dir.cy, dir.cz);

    // Project k_k0 onto aperture plane via dot products (avoid ChangeCoordinateSystem)
    double pt_x = DotProductD(k_k0, info.horAxis)
                - DotProductD(info.normald, k_k0) * DotProductD(info.normald, info.horAxis);
    double pt_y = DotProductD(k_k0, info.verAxis)
                - DotProductD(info.normald, k_k0) * DotProductD(info.normald, info.verAxis);

    const complex A(pt_x, info.lenIndices.x*m_riIm);
    const complex B(pt_y, info.lenIndices.y*m_riIm);

    if (abs(A) < m_eps2 && abs(B) < m_eps2)
    {
        double absorp = exp(m_absMag*info.lenIndices.z);
        return (m_legacySign ? m_invComplWave : -m_invComplWave) * info.area * absorp;
    }

    complex s(0, 0);

    int begin, startIndex, endIndex;

    // Always iterate forward (order-dependent iteration was broken:
    // infinite loop + buffer overflow when info.order==true).
    // Same approach as DiffractIncline where order logic is commented out.
    {
        begin = beam.nVertices-1;
        startIndex = 0;
        endIndex = beam.nVertices;
    }

    // Project first vertex once
    Point3f v0 = beam.arr[begin];
    Point3d pv0(v0.cx, v0.cy, v0.cz);
    Point3d pv0_pr = pv0 - info.normald * DotProductD(info.normald, pv0);
    double p1x = DotProductD(pv0_pr, info.horAxis) - info.projectedCenter.x;
    double p1y = DotProductD(pv0_pr, info.verAxis) - info.projectedCenter.y;

    if (abs(B) > abs(A))
    {
        double reB = real(B), imB = imag(B);

        for (int i = startIndex; i != endIndex; i = i+1)
        {
            // Project vertex inline (avoid ChangeCoordinateSystem overhead)
            Point3f vi = beam.arr[i];
            Point3d pvi(vi.cx, vi.cy, vi.cz);
            Point3d pvi_pr = pvi - info.normald * DotProductD(info.normald, pvi);
            double p2x = DotProductD(pvi_pr, info.horAxis) - info.projectedCenter.x;
            double p2y = DotProductD(pvi_pr, info.verAxis) - info.projectedCenter.y;

            double dx = p1x - p2x;
            if (fabs(dx) < m_eps1) { p1x = p2x; p1y = p2y; continue; }

            const double ai = (p1y - p2y)/dx;

            if (fabs(ai) > m_eps3)
            {
                p1x = p2x; p1y = p2y;
                continue;
            }

            complex Ci = A + ai*B;

            complex tmp;
            double absCi = abs(Ci);

            if (absCi < m_eps1)
            {
                double mul = p2x*p2x - p1x*p1x;
                tmp = complex(-m_wi2*real(Ci)*mul/2.0,
                              m_waveIndex*(p2x - p1x) + m_wi2*imag(Ci)*mul/2.0);
            }
            else if (absCi > m_eps3)
            {
                p1x = p2x; p1y = p2y;
                continue;
            }
            else
            {
                double kReCi = m_waveIndex*real(Ci);
                double kImCi = -m_waveIndex*imag(Ci);
                // Fused exp_im(kRe*p)*exp(kIm*p) = exp(kIm*p)*(cos(kRe*p) + i*sin(kRe*p))
                double s2, c2, s1, c1;
                fast_sincos(kReCi*p2x, s2, c2);
                fast_sincos(kReCi*p1x, s1, c1);
                double e2 = exp(kImCi*p2x);
                double e1 = exp(kImCi*p1x);
                tmp = complex(e2*c2 - e1*c1, e2*s2 - e1*s1)/Ci;
            }

            const double bi = p1y - ai*p1x;
            double kBi = m_waveIndex*bi;
            // Fused: exp_im(kBi*reB) * exp(-kBi*imB)
            double sb, cb;
            fast_sincos(kBi*reB, sb, cb);
            double eb = exp(-kBi*imB);
            complex phase_b(eb*cb, eb*sb);
            complex tmp2 = phase_b * tmp;

            if (isnan(real(tmp2)))
            {
                p1x = p2x; p1y = p2y;
                continue;
            }

            s += tmp2;

            p1x = p2x; p1y = p2y;
        }

        s /= B;
    }
    else
    {
        double reA = real(A), imA = imag(A);

        for (int i = startIndex; i != endIndex; i = i+1)
        {
            // Project vertex inline
            Point3f vi = beam.arr[i];
            Point3d pvi(vi.cx, vi.cy, vi.cz);
            Point3d pvi_pr = pvi - info.normald * DotProductD(info.normald, pvi);
            double p2x = DotProductD(pvi_pr, info.horAxis) - info.projectedCenter.x;
            double p2y = DotProductD(pvi_pr, info.verAxis) - info.projectedCenter.y;

            double dy = p1y - p2y;
            if (fabs(dy) < m_eps1) { p1x = p2x; p1y = p2y; continue; }

            const double ci = (p1x - p2x)/dy;

            if (fabs(ci) > m_eps3)
            {
                p1x = p2x; p1y = p2y;
                continue;
            }

            const complex Ei = A*ci + B;

            complex tmp;
            double absEi = abs(Ei);

            if (absEi < m_eps1)
            {
                double mul = p2y*p2y - p1y*p1y;
                tmp = complex(-m_wi2*real(Ei)*mul/2.0,
                              m_waveIndex*(p2y - p1y) + m_wi2*imag(Ei)*mul/2.0);
            }
            else if (absEi > m_eps3)
            {
                p1x = p2x; p1y = p2y;
                continue;
            }
            else
            {
                double kReEi = m_waveIndex*real(Ei);
                double kImEi = -m_waveIndex*imag(Ei);
                // Fused exp_im(kRe*p)*exp(kIm*p)
                double s2, c2, s1, c1;
                fast_sincos(kReEi*p2y, s2, c2);
                fast_sincos(kReEi*p1y, s1, c1);
                double e2 = exp(kImEi*p2y);
                double e1 = exp(kImEi*p1y);
                tmp = complex(e2*c2 - e1*c1, e2*s2 - e1*s1)/Ei;
            }

            const double di = p1x - ci*p1y;
            double kDi = m_waveIndex*di;
            // Fused: exp_im(kDi*reA) * exp(-kDi*imA)
            double sd, cd;
            fast_sincos(kDi*reA, sd, cd);
            double ed = exp(-kDi*imA);
            complex phase_d(ed*cd, ed*sd);
            complex tmp2 = phase_d * tmp;

            if (isnan(real(tmp2)))
            {
                p1x = p2x; p1y = p2y;
                continue;
            }

            s += tmp2;
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

complex Handler::DiffractIncline(const BeamInfo &info, const Beam &beam,
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
    if (absA < m_eps2 && absB < m_eps2)
    {
        return (m_legacySign ? m_invComplWave : -m_invComplWave) * info.area;
    }

    complex s(0, 0);

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

    if (absB > absA)
    {
        for (int i = startIndex; i != endIndex; i = /*info.order ? i-1 :*/ i+1)
        {
            p2 = ChangeCoordinateSystem(info.horAxis, info.verAxis,
                                        info.normald, beam.arr[i]) - info.projectedCenter;

            if (fabs(p1.x - p2.x) < m_eps1)
            {
                p1 = p2;
                continue;
            }

            const double ai = (p1.y - p2.y)/(p1.x - p2.x);
            const double Ci = A+ai*B;

            complex tmp;

            if (fabs(Ci) < m_eps1)
            {
                tmp = complex(-m_wi2*Ci*(p2.x*p2.x - p1.x*p1.x)/2.0,
                              m_waveIndex*(p2.x - p1.x));
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

            if (fabs(p1.y - p2.y) < m_eps1)
            {
                p1 = p2;
                continue;
            }

            const double ci = (p1.x - p2.x)/(p1.y - p2.y);
            const double Ei = A*ci+B;

            complex tmp;

            if (fabs(Ei) < m_eps1)
            {
                tmp = complex(-m_wi2*Ei*(p2.y*p2.y - p1.y*p1.y)/2.0,
                              m_waveIndex*(p2.y - p1.y));
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
    edgeData.nVertices = beam.nVertices;
    edgeData.valid = (beam.nVertices >= 3 && beam.nVertices < BeamEdgeData::MAX_EDGES);
    if (!edgeData.valid) return;

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

    // Precompute per-edge data (slopes, intercepts)
    for (int i = 0; i < nv; ++i)
    {
        int inext = (i + 1) % nv;
        double dxi = edgeData.x[inext] - edgeData.x[i];
        double dyi = edgeData.y[inext] - edgeData.y[i];
        edgeData.dx[i] = dxi;
        edgeData.dy[i] = dyi;

        // For absB > absA branch: ai = (p1y-p2y)/(p1x-p2x) where p1=v[i], p2=v[next]
        // ai = (y[i]-y[next])/(x[i]-x[next]) = -dy/(-dx) = dy/dx
        edgeData.edge_valid_x[i] = (fabs(dxi) >= m_eps1);
        if (edgeData.edge_valid_x[i])
        {
            edgeData.slope_yx[i] = dyi / dxi;  // NOT negated: (y[next]-y[i])/(x[next]-x[i])
            // but old code uses (p1y-p2y)/(p1x-p2x) = (-dy)/(-dx) = dy/dx. Same!
            edgeData.intercept_y[i] = edgeData.y[i] - edgeData.slope_yx[i] * edgeData.x[i];
        }

        // For absA >= absB branch: ci = (p1x-p2x)/(p1y-p2y) = dx/dy
        edgeData.edge_valid_y[i] = (fabs(dyi) >= m_eps1);
        if (edgeData.edge_valid_y[i])
        {
            edgeData.slope_xy[i] = dxi / dyi;
            edgeData.intercept_x[i] = edgeData.x[i] - edgeData.slope_xy[i] * edgeData.y[i];
        }
    }
}

complex Handler::DiffractInclineFast(const BeamInfo &info, const BeamEdgeData &ed,
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

    if (absA < m_eps2 && absB < m_eps2)
        return (m_legacySign ? m_invComplWave : -m_invComplWave) * info.area;

    complex s(0, 0);
    const int nv = ed.nVertices;

    // Use precomputed 2D vertices - no ChangeCoordinateSystem calls
    if (absB > absA)
    {
        double p1x = ed.x[0], p1y = ed.y[0];
        for (int i = 1; i <= nv; ++i)
        {
            double p2x = ed.x[i % nv], p2y = ed.y[i % nv];

            double dx = p1x - p2x;
            if (fabs(dx) < m_eps1) { p1x = p2x; p1y = p2y; continue; }

            double ai = (p1y - p2y) / dx;
            double Ci = A + ai * B;

            complex tmp;
            if (fabs(Ci) < m_eps1)
                tmp = complex(-m_wi2*Ci*(p2x*p2x - p1x*p1x)*0.5,
                              m_waveIndex*(p2x - p1x));
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
            if (fabs(dy) < m_eps1) { p1x = p2x; p1y = p2y; continue; }

            double ci = (p1x - p2x) / dy;
            double Ei = A * ci + B;

            complex tmp;
            if (fabs(Ei) < m_eps1)
                tmp = complex(-m_wi2*Ei*(p2y*p2y - p1y*p1y)*0.5,
                              m_waveIndex*(p2y - p1y));
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
