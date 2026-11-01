#include "Beam.h"
#include "MullerMatrix.h"
#include <array>
#include <cmath>
#include <complex>
#include <cstdlib>
#include <iostream>
using std::isfinite;
#include "../src/cuda/GpuConvexQueueMath.cuh"

using Z = std::complex<double>;
using Field = std::array<Z, 6>;

static Z asComplex(const complex &z) { return Z(real(z), imag(z)); }

// Reconstruct both incident columns as Cartesian electric fields. This tests
// physical basis invariance independently of the Jones rotation formula.
static Field electricField(const Beam &b)
{
    const Point3f p = CrossProduct(b.polarizationBasis, b.direction);
    const double pp[] = {p.cx, p.cy, p.cz};
    const double ss[] = {b.polarizationBasis.cx, b.polarizationBasis.cy,
                         b.polarizationBasis.cz};
    const complex j[] = {b.J.m11, b.J.m12, b.J.m21, b.J.m22};
    Field e;
    for (int a = 0; a < 3; ++a)
        for (int c = 0; c < 2; ++c)
            e[2*a+c] = pp[a]*asComplex(j[c]) + ss[a]*asComplex(j[2+c]);
    return e;
}

static void check(double error, const char *what)
{
    if (!std::isfinite(error) || error > 4e-15) {
        std::cerr << what << ": " << error << '\n';
        std::exit(1);
    }
}

static std::array<double, 4> stokes(Z p, Z s)
{
    const Z cross = p*std::conj(s);
    // The repository's convention is U=-2 Re(Ep Es*), V=+2 Im(Ep Es*).
    return {{std::norm(p)+std::norm(s), std::norm(p)-std::norm(s),
             -2*cross.real(), 2*cross.imag()}};
}

static void checkMueller(const Matrix2x2c &j)
{
    const double h = std::sqrt(.5);
    const Z input[][2] = {{1., 0.}, {0., 1.}, {h, h}, {h, -h},
                         {h, Z(0, h)}, {h, Z(0, -h)}};
    const MuellerMatrix m(j);
    // Six independent physical states exercise every Mueller column,
    // including the signs of phase-sensitive linear/circular conversions.
    for (const auto &state : input) {
        const auto in = stokes(state[0], state[1]);
        const auto out = stokes(asComplex(j.m11)*state[0]+asComplex(j.m12)*state[1],
                                asComplex(j.m21)*state[0]+asComplex(j.m22)*state[1]);
        for (int row = 0; row < 4; ++row) {
            double predicted = 0;
            for (int col = 0; col < 4; ++col)
                predicted += m(row, col)*in[col];
            check(std::abs(predicted-out[row]), "Jones to Stokes/Mueller");
        }
    }
}

int main()
{
    const double pi = std::acos(-1.);
    int count = 0;
    for (double tilt : {0., .3, 1.2, pi}) {
        for (double azimuth : {0., .7}) {
            const Point3f direction(std::sin(tilt)*std::cos(azimuth),
                                    std::sin(tilt)*std::sin(azimuth), std::cos(tilt));
            const Point3f basis(-std::sin(azimuth), std::cos(azimuth), 0);
            for (double delta : {0., 1e-12, -1e-10, 1e-8, -1e-7, 3e-7,
                                 .1, -.7, pi/2, pi-1e-8, -pi+1e-7, pi}) {
                Beam b;
                b.SetLight(direction, basis);
                b.J.m11 = complex(.2, .3); b.J.m12 = complex(.4, -.1);
                b.J.m21 = complex(-.7, .5); b.J.m22 = complex(.9, .2);
                const Field before = electricField(b);
                const Point3f p = CrossProduct(basis, direction);
                const Point3f next = basis*std::cos(delta) + p*std::sin(delta);

                namespace G = GpuConvexQueue;
                G::Beam gpu{};
                gpu.direction = {direction.cx, direction.cy, direction.cz};
                gpu.basis = {basis.cx, basis.cy, basis.cz};
                const complex j[] = {b.J.m11, b.J.m12, b.J.m21, b.J.m22};
                for (int c = 0; c < 4; ++c) gpu.j[c] = {real(j[c]), imag(j[c])};
                G::rotateBasis(gpu, {next.cx, next.cy, next.cz});

                b.RotateJMatrix(next);
                b.polarizationBasis = next;
                const Field after = electricField(b);
                for (int c = 0; c < 6; ++c)
                    check(std::abs(after[c]-before[c]), "CPU Cartesian field invariance");
                checkMueller(b.J);

                b.polarizationBasis = Point3f(gpu.basis.x, gpu.basis.y, gpu.basis.z);
                b.J.m11 = complex(gpu.j[0].re, gpu.j[0].im);
                b.J.m12 = complex(gpu.j[1].re, gpu.j[1].im);
                b.J.m21 = complex(gpu.j[2].re, gpu.j[2].im);
                b.J.m22 = complex(gpu.j[3].re, gpu.j[3].im);
                const Field gpuField = electricField(b);
                for (int c = 0; c < 6; ++c)
                    check(std::abs(gpuField[c]-before[c]), "CUDA math Cartesian field invariance");
                checkMueller(b.J);
                ++count;
            }
        }
    }
    std::cout << count << " CPU/CUDA basis rotations and six-state Mueller checks passed\n";
}
