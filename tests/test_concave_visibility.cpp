#include "particle/Bullet.h"
#include "particle/ConcaveHexagonal.h"
#include "particle/Hexagonal.h"
#include "particle/HexagonalAggregate.h"

#include <cassert>
#include <cmath>

int main()
{
    const complex refractiveIndex(1.3116, 0.0);

    Hexagonal column(refractiveIndex, 10.0, 7.0);
    assert(std::fabs(column.MaximumEdgeLength() - 7.0) < 1e-12);
    assert(column.MaximalDimention() > column.MaximumEdgeLength());
    column.Resize(24.0);
    const double expectedEdge = 7.0 * 24.0 / std::sqrt(149.0);
    assert(std::fabs(column.MaximumEdgeLength() - expectedEdge) < 1e-5);

    ConcaveHexagonal concave(refractiveIndex, 35.0, 50.0, 55.0079798014);

    for (int i = 0; i < 6; ++i)
    {
        assert(concave.facets[i].isVisibleIn);
        assert(!concave.facets[i].isVisibleOut);
    }

    for (int i = 6; i < 12; ++i)
    {
        assert(!concave.facets[i].isVisibleIn);
        assert(concave.facets[i].isVisibleOut);
    }

    for (int i = 12; i < 18; ++i)
    {
        assert(concave.facets[i].isVisibleIn);
        assert(!concave.facets[i].isVisibleOut);
    }

    Bullet bullet(refractiveIndex, 35.0, 50.0, 20.0);
    for (int i = 0; i < bullet.nFacets; ++i)
    {
        assert(!bullet.facets[i].isVisibleIn);
        assert(!bullet.facets[i].isVisibleOut);
    }

    HexagonalAggregate aggregate(refractiveIndex, 35.0, 50.0, 2);
    const int hiddenOut[] = {1, 2, 3, 7, 10, 11, 12, 15};
    for (int i = 0; i < aggregate.nFacets; ++i)
    {
        bool expectedOut = true;
        for (unsigned j = 0; j < sizeof(hiddenOut) / sizeof(hiddenOut[0]); ++j)
            if (hiddenOut[j] == i)
                expectedOut = false;
        assert(!aggregate.facets[i].isVisibleIn);
        assert(aggregate.facets[i].isVisibleOut == expectedOut);
    }

    return 0;
}
