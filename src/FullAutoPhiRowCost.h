#pragma once
#include <algorithm>
#include <cmath>

// The shared trace runs out to the largest orientation prefix. Diffraction
// work accumulates over the prefixes actually requested on individual rows.
namespace FullAutoPhiRowCost {
inline double Project(double wall, double trace, double diffraction,
                      double peakWork, double meanWork) {
    if (!(wall > 0 && trace > 0 && diffraction > 0) ||
        !std::isfinite(trace) || !std::isfinite(diffraction) ||
        !(peakWork >= meanWork && meanWork > 0))
        return wall * peakWork;
    double diffractionWall = std::min(wall, diffraction);
    double activeTrace = std::min(wall, trace);
    return std::max(activeTrace * peakWork,
                    (wall - diffractionWall) * peakWork + diffractionWall * meanWork);
}
}
