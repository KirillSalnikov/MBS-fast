#pragma once
#include <cmath>
#include <stdexcept>
#include <vector>

namespace ResidentFacetControl {
inline bool Enabled(const char* value) {
    return value && value[0] == '1' && value[1] == '\0';
}
inline bool Eligible(bool enabled, bool facet, bool backscatter, bool samples,
                     bool gaussian, bool hybrid) {
    return enabled && facet && !backscatter && !samples && !gaussian && !hybrid;
}
// Commit one packet's sampled controls only after every value and destination
// has been checked. Regression coefficients and known means are applied later.
inline void Accumulate(const std::vector<double>& reflection,
                       const std::vector<double>& shadow,
                       std::vector<double>& totalReflection,
                       std::vector<double>& totalShadow,
                       std::vector<double>& prefixReflection,
                       std::vector<double>& prefixShadow, bool collectPrefix) {
    const size_t cells = totalReflection.size();
    if (reflection.size() != cells || shadow.size() != cells || totalShadow.size() != cells
        || (collectPrefix && (prefixReflection.size() != cells || prefixShadow.size() != cells)))
        throw std::runtime_error("analytic facet packet has inconsistent angular dimensions");
    for (size_t i = 0; i < cells; ++i) {
        if (!std::isfinite(reflection[i]) || !std::isfinite(shadow[i])
            || !std::isfinite(totalReflection[i] + reflection[i])
            || !std::isfinite(totalShadow[i] + shadow[i])
            || (collectPrefix && (!std::isfinite(prefixReflection[i] + reflection[i])
                               || !std::isfinite(prefixShadow[i] + shadow[i]))))
            throw std::runtime_error("analytic facet packet contains nonfinite controls");
    }
    for (size_t i = 0; i < cells; ++i) {
        totalReflection[i] += reflection[i];
        totalShadow[i] += shadow[i];
        if (collectPrefix) {
            prefixReflection[i] += reflection[i];
            prefixShadow[i] += shadow[i];
        }
    }
}
}
