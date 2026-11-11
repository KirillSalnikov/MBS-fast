#include "cuda/GpuAnalyticFacet.h"
#include <cassert>
#include <limits>
#include <cstdio>
int main()
{
    using AnalyticGpuWorkspace::CanRetain;
    const std::array<size_t,6>large={{1000,1000,100,100,50,50}};
    const std::array<size_t,6>small={{100,100,10,10,5,5}};
    assert(CanRetain(2300*sizeof(double),large,small));
    assert(!CanRetain(2299*sizeof(double),large,small));
    // Retained alone fits but growing the sample accumulators would exceed it.
    const std::array<size_t,6>grow={{100,100,10,10,1000,1000}};
    assert(!CanRetain(3000*sizeof(double),large,grow));
    assert(CanRetain(4200*sizeof(double),large,grow));
    const std::array<size_t,6>huge={{std::numeric_limits<size_t>::max(),0,0,0,0,0}};
    assert(!CanRetain(std::numeric_limits<size_t>::max(),huge,small));
    assert(!CanRetain(0,large,small));
    const std::array<size_t,6>empty={{0,0,0,0,0,0}};
    assert(CanRetain(0,empty,empty));
    using AnalyticGpuWorkspace::ReuseTile;
    const std::array<size_t,6>retained={{2048*2080,2048*2080,2080,2080,0,0}};
    const size_t budget=128*1024*1024;
    assert(ReuseTile(budget,retained,2947,2080,0)==(2048*2080)/2947);
    assert(ReuseTile(budget,retained,8192,2080,0)==520);
    assert(ReuseTile(budget,retained,8193,2080,0)==2080); // Would exceed the 4x launch bound.
    assert(ReuseTile(budget,retained,1024,2080,0)==2080);
    assert(ReuseTile(1024,retained,2947,2080,0)==2080); // Retained memory no longer fits.
    assert(ReuseTile(budget,empty,2947,2080,0)==2080);
    assert(ReuseTile(budget,retained,2947,2080,budget/sizeof(double))==2080);
    std::puts("PASS analytic GPU scratch retention bounds");
}
