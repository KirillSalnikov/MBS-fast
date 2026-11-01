#include "Splitting.h"
#include <cassert>
#include <cmath>
#include <complex>
#include <iostream>
using std::isfinite;
#include "../src/cuda/GpuConvexQueueMath.cuh"

using C = std::complex<double>;
C cv(const complex& z) { return C(real(z), imag(z)); }

void checkMatrix(const Matrix2x2c& actual, const Matrix2x2c& expected)
{
    assert(std::abs(cv(actual.m11)-cv(expected.m11)) < 1e-14);
    assert(std::abs(cv(actual.m12)-cv(expected.m12)) < 1e-14);
    assert(std::abs(cv(actual.m21)-cv(expected.m21)) < 1e-14);
    assert(std::abs(cv(actual.m22)-cv(expected.m22)) < 1e-14);
}

int main()
{
    size_t cases = 0;
    for (double c : {1., .8, .001, 1e-6, 1e-8, 1e-10, 1e-12})
    for (double rotation : {0., .6})
    {
        const double z = std::sqrt(1.-rotation*rotation);
        const Point3f normal(rotation, 0., z), tangent(z, 0., -rotation);
        Beam incident, inside, outside;
        incident.SetLight(tangent*std::sqrt(1.-c*c)+normal*c, Point3f(0.,1.,0.));
        incident.J.m11=complex(.2,.3); incident.J.m12=complex(.4,-.1);
        incident.J.m21=complex(-.7,.5); incident.J.m22=complex(.9,.2);
        Matrix2x2c zero; zero.Fill(complex(0.,0.));
        Splitting split(false);
        split.ComputeRiParams(complex(1.,0.));
        split.ComputeCosA(normal, incident.direction);
        assert(split.IsIncident());
        split.ComputeRegularBeamParamsExternal(normal, incident, inside, outside);
        checkMatrix(inside.J, incident.J);
        checkMatrix(outside.J, zero);
        assert(Norm(inside.direction-incident.direction) == 0.);

        split.ComputeSplittingParams(incident.direction, normal);
        assert(!split.IsCompleteReflection());
        split.ComputeRegularBeamsParams(normal, incident, inside, outside);
        checkMatrix(inside.J, zero);
        checkMatrix(outside.J, incident.J);
        assert(Norm(outside.direction-incident.direction) == 0.);
        assert(Norm(split.ChangeBeamDirection(incident.direction, normal,
            Location::In, Location::Out)-incident.direction) == 0.);
        assert(Norm(split.ChangeBeamDirection(incident.direction, -normal,
            Location::Out, Location::In)-incident.direction) == 0.);

        // Exercise the very same host/device split used by resident CUDA tracing.
        namespace G = GpuConvexQueue;
        G::Input input{}; input.ri=1.;
        G::Beam source{}, reflected{}, transmitted{};
        source.direction={incident.direction.cx,incident.direction.cy,incident.direction.cz};
        source.basis={0.,1.,0.};
        source.j[0]={.2,.3}; source.j[1]={.4,-.1};
        source.j[2]={-.7,.5}; source.j[3]={.9,.2};
        reflected.count=3;
        G::Vec n={normal.cx,normal.cy,normal.cz};
        const bool outgoing=G::split(input,source,n,G::dot(n,source.direction),reflected,transmitted);
        assert(outgoing);
        assert(G::finiteBeam(reflected) && G::finiteBeam(transmitted));
        assert(G::norm(G::sub(transmitted.direction,source.direction)) == 0.);
        for (int j=0;j<4;++j) {
            assert(std::abs(C(reflected.j[j].re,reflected.j[j].im)) == 0.);
            assert(std::abs(C(transmitted.j[j].re,transmitted.j[j].im)
                          -C(source.j[j].re,source.j[j].im)) < 1e-14);
        }
        ++cases;
    }
    // The phase direction must use the same stable normal wavevector as the
    // coefficients, also for vanishing absorption and nearly matched indices.
    using L = std::complex<long double>;
    for (double n : {1., 1.+1e-12, 1.+1e-8})
    for (double k : {0., 1e-30, -1e-30, .1})
    for (double c : {.001, 1e-6, 1e-8, 1e-10, 1e-12})
    {
        Splitting split(false); split.ComputeRiParams(complex(n,k));
        Beam incident, transmitted, reflected;
        incident.SetLight(Point3f(std::sqrt(1.-c*c),0.,c),Point3f(0.,1.,0.));
        incident.J.Fill(complex(0.,0.));
        incident.J.m11=incident.J.m22=complex(1.,0.);
        split.ComputeCosA(Point3f(0.,0.,1.),incident.direction);
        split.ComputeRegularBeamParamsExternal(Point3f(0.,0.,1.),incident,transmitted,reflected);
        const L m(n,k),q=std::sqrt((m-1.L)*(m+1.L)+(long double)c*c);
        const double expected=static_cast<double>(q.real()/std::sqrt(1.L-(long double)c*c+q.real()*q.real()));
        assert(std::abs(transmitted.direction.cz-expected) < 2e-12*expected);
        const long double ci=c;
        const L rp=(m*m*ci-q)/(m*m*ci+q),rs=(ci-q)/(ci+q);
        assert(std::abs(cv(reflected.J.m11)-C(rp.real(),rp.imag())) < 2e-12);
        assert(std::abs(cv(reflected.J.m22)-C(rs.real(),rs.imag())) < 2e-12);
        ++cases;
    }
    std::cout << cases << " matched-index CPU/device splits and near-vacuum phase directions passed\n";
}
