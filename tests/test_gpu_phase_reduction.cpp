#include "../src/cuda/GpuPhaseReduction.cuh"
#include <algorithm>
#include <cstdio>
#include <limits>
#include <random>

// Independent libm reference, including phase wrap boundaries and large
// arguments where casting before range reduction would lose accuracy.
int main()
{
    double worst = 0, worstFloat = 0;
    auto check = [&](double x) {
        double r;
        if (!MbsGpuPhase::Reduce(x, r)) return false;
        const long double exactSin = sinl((long double)x);
        const long double exactCos = cosl((long double)x);
        worst = std::max(worst, (double)hypotl(sinl(r)-exactSin, cosl(r)-exactCos));
        const float f = (float)r;
        worstFloat = std::max(worstFloat,
            (double)hypotl((long double)sinf(f)-exactSin, (long double)cosf(f)-exactCos));
        return std::abs(r) < 3.141593;
    };
    std::mt19937_64 gen(372139);
    std::uniform_real_distribution<double> draw(-1048576.,1048576.);
    for (int i=0; i<1000000; ++i) if (!check(draw(gen))) return 1;
    const long double pi = acosl(-1.L);
    for (int i=-333772; i<=333772; i+=7) {
        const double x = (double)(i*pi);
        for (double v : {std::nextafter(x,-INFINITY),x,std::nextafter(x,INFINITY)})
            if (std::abs(v)<=1048576 && !check(v)) return 2;
    }
    for (double x : {-1048576.,-0.,0.,1048576.}) if (!check(x)) return 3;
    double r;
    for (double x : {std::nextafter(1048576.,INFINITY),
                    std::nextafter(-1048576.,-INFINITY),
                    (double)INFINITY,(double)-INFINITY,std::numeric_limits<double>::quiet_NaN()})
        if (MbsGpuPhase::Reduce(x,r)) return 4;
    std::printf("phase reduction complex error: double=%.9g float=%.9g\n",worst,worstFloat);
    return worst<2e-15 && worstFloat<2e-7 ? 0 : 5;
}
