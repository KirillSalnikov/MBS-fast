#include "../src/FullAutoPhiRegret.h"
#include <cassert>

int main() {
    using namespace FullAutoPhiRegret;
    // The two estimates favor different extremes. The middle rule costs
    // at most 20 percent more than the best rule under either estimate.
    std::map<int, Cost> costs{{8, {10, 18}}, {16, {12, 12}}, {32, {15, 10}}};
    assert(Select(costs, 8).phi == 16);
    for (auto& entry : costs) {
        entry.second.row *= 1000;
        entry.second.paired *= .25;
    }
    assert(Select(costs, 8).phi == 16);
    costs = {{8, {2, 3}}, {16, {2, 3}}};
    assert(Select(costs, 16).phi == 16);
    costs[8].paired = 0;
    const Choice fallback = Select(costs, 16);
    assert(fallback.phi == 16 && !std::isfinite(fallback.relative));
}
