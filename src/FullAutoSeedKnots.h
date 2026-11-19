#pragma once
#include <algorithm>
#include <cmath>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace FullAutoSeedKnots {
// Input is an independently chosen subset of the requested physical angles.
inline std::vector<size_t> Parse(const std::string& text, const std::vector<double>& grid) {
    if (grid.size() < 2) throw std::runtime_error("theta proposal needs a grid");
    std::istringstream input(text);
    std::vector<size_t> result;
    for (;;) {
        input >> std::ws;
        if (input.eof()) break;
        double angle;
        if (!(input >> angle) || !std::isfinite(angle)) throw std::runtime_error("invalid theta proposal angle");
        auto pos = std::lower_bound(grid.begin(),grid.end(),angle);
        if (pos==grid.end() || (pos!=grid.begin() && angle-pos[-1] < *pos-angle)) --pos;
        const auto id = static_cast<size_t>(pos-grid.begin());
        if (std::abs(*pos-angle)>1e-10 || (!result.empty() && result.back()>=id))
            throw std::runtime_error("theta proposal must list distinct requested angles in increasing order");
        result.push_back(id);
    }
    if (result.size()<2 || result.front()!=0 || result.back()+1!=grid.size())
        throw std::runtime_error("theta proposal must include both boundaries");
    return result;
}
}
