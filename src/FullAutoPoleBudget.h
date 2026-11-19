#pragma once
#include <algorithm>
#include <cstddef>

namespace FullAutoPoleBudget {
// Bound the extra initial angular work at the exact poles to 1/16 of the
// current working grid; all levels above the base keep their initial budgets.
inline int Count(int initial, int limit, int phi, std::size_t rows, std::size_t poles) {
    if (initial <= 0 || initial > limit || !poles || poles > 2 || poles > rows)
        return initial;
    int multiplier = 1;
    const int azimuthCap = std::min(phi, 8);
    while (multiplier * 2 <= azimuthCap && initial <= limit / (multiplier * 2)
        && std::size_t(multiplier * 2 - 1) <= rows / (16 * poles))
        multiplier *= 2;
    return initial * multiplier;
}
}
