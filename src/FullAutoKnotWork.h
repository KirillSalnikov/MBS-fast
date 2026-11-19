#pragma once
#include <algorithm>
#include <cmath>
#include <cstddef>

namespace FullAutoKnotWork {
// Scheduling proxy for a single locally budgeted residual level. A cold angular
// split is charged for tracing two new production guards plus the whole stencil
// at the measured interpolation prefix. All subsequent checks remain direct.
struct Forecast {
    double growth = 0, split = 0;
    bool promote = false;
};
inline Forecast Compare(double tolerance, double mean, double noise,
                        int observed, int initial, int maximum, std::size_t stencil)
{
    Forecast f;
    if (!(tolerance>0 && mean>0 && noise>0 && mean>=2*noise && mean+noise>tolerance)
        || !std::isfinite(mean+noise) || observed<=0 || initial<=0 || maximum<observed || !stencil)
        return f;
    const double cap=double(maximum)/observed;
    const double margin=tolerance-mean;
    const double needed=margin>0 ? std::min(cap,std::pow(noise/margin,2)) : cap;
    f.growth=observed*std::max(0.,needed-1.);
    f.split=2.*initial+double(stencil)*observed;
    f.promote=f.growth>f.split;
    return f;
}
}
