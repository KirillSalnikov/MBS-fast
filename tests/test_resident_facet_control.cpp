#include "ResidentFacetControl.h"
#include <cassert>
#include <cstdio>
#include <limits>

int main() {
    using namespace ResidentFacetControl;
    for (const char* flag : {"", "0", "true", "01", "1x"}) assert(!Enabled(flag));
    assert(!Enabled(nullptr) && Enabled("1"));
    assert(Eligible(true, true, false, false, false, false));
    assert(!Eligible(false, true, false, false, false, false));
    assert(!Eligible(true, false, false, false, false, false));
    assert(!Eligible(true, true, true, false, false, false));
    assert(!Eligible(true, true, false, true, false, false));
    assert(!Eligible(true, true, false, false, true, false));
    assert(!Eligible(true, true, false, false, false, true));
    std::vector<double> totalR(2), totalS(2), prefixR(2), prefixS(2);
    const std::vector<double> firstR{.25, .5}, firstS{.75, 1};
    const std::vector<double> secondR{2, 3}, secondS{4, 5};
    Accumulate(firstR, firstS, totalR, totalS, prefixR, prefixS, true);
    Accumulate(secondR, secondS, totalR, totalS, prefixR, prefixS, false);
    assert(totalR == std::vector<double>({2.25, 3.5}));
    assert(totalS == std::vector<double>({4.75, 6}));
    assert(prefixR == firstR && prefixS == firstS);
    const auto savedR = totalR, savedS = totalS, savedPrefixR = prefixR, savedPrefixS = prefixS;
    auto bad = secondR; bad.back() = std::numeric_limits<double>::quiet_NaN();
    bool rejected = false;
    try { Accumulate(bad, secondS, totalR, totalS, prefixR, prefixS, true); }
    catch (const std::runtime_error&) { rejected = true; }
    assert(rejected && totalR == savedR && totalS == savedS
        && prefixR == savedPrefixR && prefixS == savedPrefixS);
    rejected = false;
    try { Accumulate({1}, secondS, totalR, totalS, prefixR, prefixS, true); }
    catch (const std::runtime_error&) { rejected = true; }
    assert(rejected && totalR == savedR && prefixR == savedPrefixR);
    // Overflow in the LAST cell must not mutate the first cell or the prefix.
    totalR.back() = std::numeric_limits<double>::max();
    const auto beforeOverflow = totalR;
    rejected = false;
    try { Accumulate({1, std::numeric_limits<double>::max()}, secondS,
                     totalR, totalS, prefixR, prefixS, true); }
    catch (const std::runtime_error&) { rejected = true; }
    assert(rejected && totalR == beforeOverflow && prefixR == savedPrefixR);
    std::puts("PASS ordinary-control eligibility, single prefix subset, dimension/nonfinite/overflow atomic rejection");
}
