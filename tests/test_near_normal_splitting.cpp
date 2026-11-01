#include "Splitting.h"
#include <cassert>
#include <cmath>
#include <complex>
#include <iostream>
using std::isfinite;
#include "../src/cuda/GpuConvexQueueMath.cuh"

static void near(Point3f a,Point3f b,double tolerance=3e-15) {
    assert(std::sqrt(Norm(a-b)) < tolerance);
}
int main() {
    size_t count=0;
    for(double rotation:{0.,.6}) {
        const double z=std::sqrt(1-rotation*rotation);
        const Point3f normal(rotation,0,z),tangent(z,0,-rotation),basis(0,1,0);
        for(double angle:{1e-2,1e-4,1e-5,6e-6,5e-6,1e-6,1e-8,1e-10,4e-12}) {
            Point3f direction=normal*std::cos(angle)+tangent*std::sin(angle);
            for(double loss:{0.,.1}) {
                const double n=1.3116;
                Splitting split(false);split.ComputeRiParams(complex(n,loss));
                split.ComputeCosA(normal,direction);
                // Even when cos(angle) rounds to 1, a nonzero tangential
                // wavevector must survive the normal-incidence dispatch.
                assert(!split.IsNormalIncidence());
                Beam inc,inside,outside;inc.SetLight(direction,basis);
                split.ComputeRegularBeamParamsExternal(normal,inc,inside,outside);
                const std::complex<double> ri(n,loss),q=std::sqrt(ri*ri-std::sin(angle)*std::sin(angle));
                Point3f expected=tangent*std::sin(angle)+normal*q.real();Normalize(expected);
                near(inside.direction,expected);
                near(outside.direction,direction-normal*(2*std::cos(angle)));
                ++count;
            }
            Splitting split(false);split.ComputeRiParams(complex(1.3116,0));
            split.ComputeCosA(normal,direction);split.ComputeSplittingParams(direction,normal);
            assert(!split.IsNormalIncidence()&&!split.IsCompleteReflection());
            Beam inc,inside,outside;inc.SetLight(direction,basis);
            split.ComputeRegularBeamsParams(normal,inc,inside,outside);
            Point3f expected=tangent*(1.3116*std::sin(angle))+normal*std::sqrt(1-1.3116*1.3116*std::sin(angle)*std::sin(angle));
            near(outside.direction,expected);
            near(inside.direction,direction-normal*(2*std::cos(angle)));
            namespace G=GpuConvexQueue;
            G::Input input{};input.ri=1.3116;G::Beam source{},reflected{},transmitted{};
            source.direction={direction.cx,direction.cy,direction.cz};source.basis={0,1,0};source.j[0]=source.j[3]={1,0};
            source.count=reflected.count=3;reflected.vertices[0]={0,0,0};reflected.vertices[1]={1,0,0};reflected.vertices[2]={0,1,0};
            G::Vec normalG={normal.cx,normal.cy,normal.cz};
            assert(G::split(input,source,normalG,G::dot(normalG,source.direction),reflected,transmitted));
            near(Point3f(transmitted.direction.x,transmitted.direction.y,transmitted.direction.z),expected);
            near(Point3f(reflected.direction.x,reflected.direction.y,reflected.direction.z),inside.direction);
            assert(G::finiteBeam(reflected)&&G::finiteBeam(transmitted));
            ++count;
        }
        for(double angle:{0.,0x1p-43}) {
            Splitting split(false);split.ComputeRiParams(complex(1.3116,0));
            split.ComputeCosA(normal,normal*std::cos(angle)+tangent*std::sin(angle));
            assert(split.IsNormalIncidence());
        }
    }
    std::cout<<count<<" near-normal Snell directions and shared CPU/CUDA splits passed\n";
}
