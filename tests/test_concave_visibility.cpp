#include "particle/ConcaveHexagonal.h"

#include <cassert>

int main()
{
    const complex refractiveIndex(1.3116, 0.0);
    ConcaveHexagonal particle(refractiveIndex, 35.0, 50.0, 55.0079798014);

    for (int i = 0; i < 6; ++i)
    {
        assert(particle.facets[i].isVisibleIn);
        assert(!particle.facets[i].isVisibleOut);
    }

    for (int i = 6; i < 12; ++i)
    {
        assert(!particle.facets[i].isVisibleIn);
        assert(particle.facets[i].isVisibleOut);
    }

    for (int i = 12; i < 18; ++i)
    {
        assert(particle.facets[i].isVisibleIn);
        assert(!particle.facets[i].isVisibleOut);
    }

    return 0;
}
