#include "Splitting.h"

#include <algorithm>
#include <cmath>
#include <float.h>

Splitting::Splitting(bool isOpticalPath)
{
    m_isOpticalPath = isOpticalPath;
}

void Splitting::ComputeRiParams(const complex &ri)
{
    m_ri = ri;
    double re = real(m_ri);
    double im = imag(m_ri);
    m_cRiRe = re*re - im*im;
    m_cRiRe2 = m_cRiRe * m_cRiRe;
    m_cRiIm = 4*re*re*im*im;
}

void Splitting::ComputeSplittingParams(const Point3f &dir, const Point3f &normal)
{
    // Work with the tangential component directly.  The former formulation
    // used r = dir/cosA-normal and s ~ 1/cosA^2, although the factors cancel
    // in every resulting direction.  Near grazing incidence those temporary
    // values overflowed or lost most significant digits.
    r = dir - normal*cosA;
    reRiEff = ComputeEffectiveReRi();
    s = 1.0/reRiEff - Norm(r);
}

bool Splitting::IsCompleteReflection()
{
    return s <= 0.0;
}

bool Splitting::IsNormalIncidence()
{
    return cosA >= EPS_COS_00;
}

bool Splitting::IsIncident()
{
    return cosA > EPS_GRAZING_INCIDENCE;
}

double Splitting::ComputeSegmentOpticalPath(const Beam &beam, const Point3f &facetPoint) const
{
    // beam.front is the phase-plane offset established when the beam is
    // created.  Clipping may replace its aperture polygon and move Center(),
    // but it must not move that phase plane.  Therefore the signed propagation
    // distance follows directly from the plane equation instead of being
    // re-signed relative to the current aperture centroid.
    double path = SignedPhasePlaneDistance(beam.direction, beam.front,
                                           facetPoint);

    if (beam.location == Location::In)
    {
        path *= sqrt(reRiEff);
    }

    return path;
}

// REF: создать отдельные классы RegularSplitting, CRSplitting, NormalSplitting

void Splitting::ComputeCRBeamParams(const Point3f &normal, const Beam &incidentBeam,
                                    Beam &inBeam)
{
    Point3f reflDir = incidentBeam.direction - normal*(2.0*cosA);
    Normalize(reflDir);
    inBeam.SetLight(reflDir, incidentBeam.polarizationBasis);

    complex cv, ch;
    ComputeCRJonesParams(cv, ch);

    inBeam.J = incidentBeam.J;
    inBeam.MultiplyJonesMatrix(cv, ch);

    if (m_isOpticalPath)
    {
        double path = ComputeSegmentOpticalPath(incidentBeam, inBeam.Center());
#ifdef _DEBUG // DEB
        inBeam.ops = incidentBeam.ops;
        inBeam.ops.push_back(path);
#endif
        path += incidentBeam.opticalPath;
        inBeam.AddOpticalPath(path);
    }
}

void Splitting::ComputeRegularJonesParams(const Point3f &normal, const Beam &incidentBeam, Beam &inBeam, Beam &outBeam)
{
    inBeam.J = incidentBeam.J;
    outBeam.J = incidentBeam.J;

    double cosG = DotProduct(normal, outBeam.direction);

    complex tmp0 = m_ri * cosA;
    complex tmp1 = m_ri * cosG;
    complex Tv0 = tmp1 + cosA;
    complex Th0 = tmp0 + cosG;

    complex tmp = 2.0 * tmp0;

    outBeam.MultiplyJonesMatrix(tmp/Tv0, tmp/Th0);

    complex Tv = cosA - tmp1;
    complex Th = tmp0 - cosG;
    inBeam.MultiplyJonesMatrix(Tv/Tv0, Th/Th0);
}

void Splitting::ComputeInternalRefractiveDirection(const Vector3f &tangential,
                                                   const Vector3f &normal,
                                                   Vector3f &dir)
{
    const double sinA2 = std::max(0.0, Norm(tangential));
    const double delta = m_cRiRe - sinA2;
    const double loss = std::sqrt(std::max(0.0, m_cRiIm));
    const double root = std::hypot(delta, loss);

    // q=(delta+sqrt(delta^2+loss^2))/2.  Use the conjugate form when
    // delta<0 to avoid cancellation near the critical angle.
    double q = 0.0;
    if (delta >= 0.0)
        q = 0.5*(delta + root);
    else if (loss > 0.0)
        q = (0.5*loss)*(loss/(root - delta));

    dir = tangential - normal*std::sqrt(std::max(0.0, q));
    Normalize(dir);
}

void Splitting::ComputeCRJonesParams(complex &cv, complex &ch)
{
    const double bf = reRiEff*(1.0 - cosA*cosA) - 1.0;
    double im = (bf > 0) ? sqrt(bf) : 0;

    const complex sq(0, im);
    complex tmp0 = m_ri * cosA;
    complex tmp1 = m_ri * sq;

    cv = (cosA - tmp1)/(tmp1 + cosA);
    ch = (tmp0 - sq)/(tmp0 + sq);
}

void Splitting::ComputeCosA(const Point3f &normal, const Point3f &incidentDir)
{
    cosA = std::max(-1.0, std::min(1.0,
        static_cast<double>(DotProduct(normal, incidentDir))));
}

void Splitting::ComputeRegularBeamsParams(const Point3f &normal,
                                          const Beam &incidentBeam,
                                          Beam &inBeam, Beam &outBeam)
{
    Point3f reflDir = incidentBeam.direction - normal*(2.0*cosA);
    Normalize(reflDir);
    inBeam.SetLight(reflDir, incidentBeam.polarizationBasis);

    Point3f refrDir = r + normal*std::sqrt(std::max(0.0, s));
    Normalize(refrDir);
    outBeam.SetLight(refrDir, incidentBeam.polarizationBasis);

    ComputeRegularJonesParams(normal, incidentBeam, inBeam, outBeam);

    if (m_isOpticalPath)
    {
        double path = ComputeSegmentOpticalPath(incidentBeam, inBeam.Center());
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

void Splitting::ComputeNormalBeamParams(const Beam &incidentBeam,
                                        Beam &inBeam, Beam &outBeam)
{
    const Point3f &dir = incidentBeam.direction;
    inBeam.SetLight(-dir, incidentBeam.polarizationBasis);
    outBeam.SetLight(dir, incidentBeam.polarizationBasis);

    inBeam.J = incidentBeam.J;
    outBeam.J = incidentBeam.J;

    complex temp;

    temp = (2.0 * m_ri)/(1.0 + m_ri); // OPT: вынести целиком
    outBeam.MultiplyJonesMatrix(temp, temp);

    temp = (1.0 - m_ri)/(1.0 + m_ri); // OPT: вынести целиком
    inBeam.MultiplyJonesMatrix(temp, -temp);

    if (m_isOpticalPath)
    {
        double path = ComputeSegmentOpticalPath(incidentBeam, inBeam.Center());
        path += incidentBeam.opticalPath;
#ifdef _DEBUG // DEB
        inBeam.ops = incidentBeam.ops;
        inBeam.ops.push_back(path);
        outBeam.ops = incidentBeam.ops;
        outBeam.ops.push_back(path);
#endif
        inBeam.AddOpticalPath(path);
        outBeam.AddOpticalPath(path);
    }
}

void Splitting::ComputeNormalBeamParamsExternal(const Light &incidentLight,
                                                Beam &inBeam, Beam &outBeam)
{
    inBeam.SetLight(incidentLight);
    outBeam.SetLight(-incidentLight.direction, incidentLight.polarizationBasis);

    inBeam.J.m11 = 2.0/(m_ri + 1.0); // OPT: вынести
    inBeam.J.m22 = inBeam.J.m11;

    outBeam.J.m11 = (m_ri - 1.0)/(m_ri + 1.0);
    outBeam.J.m22 = -outBeam.J.m11;
}

void Splitting::ComputeRegularBeamParamsExternal(const Point3f &facetNormal,
                                                 Beam &incidentBeam,
                                                 Beam &inBeam, Beam &outBeam)
{
    inBeam.polarizationBasis = incidentBeam.polarizationBasis;

    inBeam.J = incidentBeam.J;
    outBeam.J = incidentBeam.J;

    Point3f refrDir;
    const Point3f &dir = incidentBeam.direction;

    const Point3f tangential = dir - facetNormal*cosA;

    refrDir = dir - facetNormal*(2.0*cosA);
    Normalize(refrDir);

    ComputeInternalRefractiveDirection(tangential, -facetNormal, inBeam.direction);

    outBeam.SetLight(refrDir, incidentBeam.polarizationBasis);

    double cosB = DotProduct(facetNormal, inBeam.direction);

    complex tmp0 = m_ri * cosA;
    complex tmp1 = m_ri * cosB;

    complex Tv0 = tmp0 + cosB;
    complex Th0 = tmp1 + cosA;

    complex Tv = tmp0 - cosB;
    complex Th = cosA - tmp1;
    outBeam.MultiplyJonesMatrix(Tv/Tv0, Th/Th0);

    double cos2A = 2.0*cosA;
    inBeam.MultiplyJonesMatrix(cos2A/Tv0, cos2A/Th0);
}

double Splitting::ComputeIncidentOpticalPath(const Point3f &direction,
                                             const Point3f &facetPoint)
{
#ifdef _DEBUG  /* DEB */
//    double dp = DotProduct(direction, facetPoint);
#endif
    return DotProductD(direction, facetPoint);
}

double Splitting::ComputeOutgoingOpticalPath(const Beam &beam)
{
    return beam.front;
}

Point3f Splitting::ChangeBeamDirection(const Vector3f &oldDir,
                                       const Vector3f &normal,
                                       Location oldLoc, Location loc)
{
    Point3f newDir;

    if (oldLoc == Location::Out) // refraction
    {
        ComputeCosA(oldDir, -normal);
        const Vector3f interfaceNormal = -normal;
        const Vector3f tangential = oldDir
            - interfaceNormal*cosA;

        if (oldLoc == loc)
        {
            newDir = oldDir
                - interfaceNormal*(2.0*cosA);
            Normalize(newDir);
        }
        else
        {
            ComputeInternalRefractiveDirection(tangential, normal, newDir);
        }
    }
    else // reflection
    {
        ComputeCosA(oldDir, normal);
        ComputeSplittingParams(oldDir, normal);

        if (oldLoc == loc)
        {
            newDir = oldDir - normal*(2.0*cosA);
        }
        else
        {
            newDir = r + normal*std::sqrt(std::max(0.0, s));
        }

        Normalize(newDir);
    }

    return newDir;
}

complex Splitting::GetRi() const
{
    return m_ri;
}

double Splitting::ComputeEffectiveReRi() const
{
    const double absCos = std::max(std::fabs(cosA), DBL_MIN);
    const double lossTerm = std::sqrt(std::max(0.0, m_cRiIm))/absCos;
    const double root = std::hypot(m_cRiRe, lossTerm);

    if (m_cRiRe >= 0.0)
        return 0.5*(m_cRiRe + root);
    if (lossTerm == 0.0)
        return 0.0;
    return (0.5*lossTerm)*(lossTerm/(root - m_cRiRe));
}
