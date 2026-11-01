#include "Splitting.h"
#include <algorithm>
#include <cassert>
#include <cmath>
#include <complex>
#include <iostream>

using C = std::complex<double>;
C convert(const complex& value) { return C(real(value), imag(value)); }

int main()
{
    double maximumError = 0.;
    double maximumBoundaryError = 0.;
    double maximumPowerError = 0.;
    size_t cases = 0;
    for (double index : {1.001, 1.3116, 2., 4.})
    for (double loss : {0., 1e-8, .01, .1, 1., 10., 100., -1., -10.})
    for (double angle : {0., 1e-3, 15., 30., 45., 60., 75., 85., 89.9, 89.999999})
    {
        const double theta = angle*std::acos(-1.)/180.;
        const double ci = std::cos(theta), si = std::sin(theta);
        const C m(index, loss);
        // Independent Snell-angle form; production uses the normal wavevector.
        const C ct = std::cos(std::asin(C(si, 0)/m));
        const C rp = (m*ci-ct)/(m*ci+ct);
        const C rs = (ci-m*ct)/(ci+m*ct);
        const C tp = 2.*ci/(m*ci+ct);
        const C ts = 2.*ci/(ci+m*ct);
        Splitting split(false);
        split.ComputeRiParams(complex(index, loss));
        Beam incident, transmitted, reflected;
        incident.SetLight(Point3f(si, 0., ci), Point3f(0., 1., 0.));
        incident.J.Fill(complex(0., 0.));
        incident.J.m11 = incident.J.m22 = complex(1., 0.);
        split.ComputeCosA(Point3f(0., 0., 1.), incident.direction);
        if (split.IsNormalIncidence())
            split.ComputeNormalBeamParamsExternal(incident, transmitted, reflected);
        else
            split.ComputeRegularBeamParamsExternal(Point3f(0., 0., 1.),
                                                   incident, transmitted, reflected);
        const C actualRp = convert(reflected.J.m11);
        const C actualRs = convert(reflected.J.m22);
        const C actualTp = convert(transmitted.J.m11);
        const C actualTs = convert(transmitted.J.m22);
        const double error = std::max({std::abs(actualRp-rp),
            std::abs(actualRs-rs), std::abs(actualTp-tp), std::abs(actualTs-ts)});
        maximumError = std::max(maximumError, error);
        assert(error < 5e-12);

        // Tangential E and H must both be continuous at the entry interface.
        // In particular, changing only reflection breaks ts = 1 + rs.
        const double boundaryError = std::max({
            std::abs(actualTs-(1.+actualRs)),
            std::abs(m*ct*actualTs-ci*(1.-actualRs)),
            std::abs(ct*actualTp-ci*(1.-actualRp)),
            std::abs(m*actualTp-(1.+actualRp))});
        maximumBoundaryError = std::max(maximumBoundaryError, boundaryError);
        assert(boundaryError < 5e-12);

        // Power entering the absorbing half-space uses the complex-angle
        // Poynting flux, not the cosine of the real phase-ray direction.
        const double powerS = std::norm(actualRs)
            + std::norm(actualTs)*std::real(m*ct)/ci;
        const double powerP = std::norm(actualRp)
            + std::norm(actualTp)*std::real(m*std::conj(ct))/ci;
        const double powerError = std::max(std::abs(powerS-1.), std::abs(powerP-1.));
        maximumPowerError = std::max(maximumPowerError, powerError);
        assert(powerError < 5e-12);
        assert(std::abs(reflected.direction.cx-si) < 2e-14);
        assert(std::abs(reflected.direction.cz+ci) < 2e-14);
        assert(std::norm(convert(reflected.J.m11)) <= 1.+1e-12);
        assert(std::norm(convert(reflected.J.m22)) <= 1.+1e-12);
        // Verify row multiplication for incident polarization mixing.
        if (!split.IsNormalIncidence()) {
            incident.J.m11=complex(.2,.3); incident.J.m12=complex(.4,-.1);
            incident.J.m21=complex(-.7,.5); incident.J.m22=complex(.9,.2);
            split.ComputeRegularBeamParamsExternal(Point3f(0.,0.,1.),
                                                   incident, transmitted, reflected);
            assert(std::abs(convert(reflected.J.m11)-rp*convert(incident.J.m11)) < 5e-12);
            assert(std::abs(convert(reflected.J.m12)-rp*convert(incident.J.m12)) < 5e-12);
            assert(std::abs(convert(reflected.J.m21)-rs*convert(incident.J.m21)) < 5e-12);
            assert(std::abs(convert(reflected.J.m22)-rs*convert(incident.J.m22)) < 5e-12);
            assert(std::abs(convert(transmitted.J.m11)-tp*convert(incident.J.m11)) < 5e-12);
            assert(std::abs(convert(transmitted.J.m12)-tp*convert(incident.J.m12)) < 5e-12);
            assert(std::abs(convert(transmitted.J.m21)-ts*convert(incident.J.m21)) < 5e-12);
            assert(std::abs(convert(transmitted.J.m22)-ts*convert(incident.J.m22)) < 5e-12);
        }
        ++cases;
    }
    // Nearly vacuum at grazing incidence: subtracting sin(theta)^2 from
    // m^2 after rounding sin(theta)^2 to one creates spurious reflection.
    for (double loss : {1e-30, -1e-30})
    for (double angle : {85., 89.9, 89.999999})
    {
        const double theta = angle*std::acos(-1.)/180.;
        Splitting split(false);
        split.ComputeRiParams(complex(1., loss));
        Beam incident, transmitted, reflected;
        incident.SetLight(Point3f(std::sin(theta), 0., std::cos(theta)),
                          Point3f(0., 1., 0.));
        incident.J.Fill(complex(0., 0.));
        incident.J.m11 = incident.J.m22 = complex(1., 0.);
        split.ComputeCosA(Point3f(0., 0., 1.), incident.direction);
        split.ComputeRegularBeamParamsExternal(Point3f(0., 0., 1.),
                                               incident, transmitted, reflected);
        assert(std::abs(convert(reflected.J.m11)) < 1e-12);
        assert(std::abs(convert(reflected.J.m22)) < 1e-12);
        assert(std::abs(convert(transmitted.J.m11)-1.) < 1e-12);
        assert(std::abs(convert(transmitted.J.m22)-1.) < 1e-12);
    }
    std::cout << cases << " complex Fresnel interface cases and 6 grazing vacuum limits passed;"
              << " maximum amplitude error " << maximumError
              << ", boundary error " << maximumBoundaryError
              << ", power error " << maximumPowerError << '\n';
}
