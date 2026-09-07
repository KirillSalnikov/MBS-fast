#include "AdaptivePhi.h"
#include <cassert>
#include <cstdlib>

int main()
{
    const char *names[] = {"MBS_PHI_ADAPT_MIN", "MBS_PHI_ADAPT_MAX",
        "MBS_PHI_ADAPT_M11_TOL", "MBS_PHI_ADAPT_POL_TOL"};
    for (const char *name : names) unsetenv(name);
    AdaptivePhi::Config defaults;
    assert(defaults.first == 75 && defaults.maximum == 9600);
    for (const char *invalid : {"0", "15", "75.5", "nan", "inf", "text", "1000000"})
    {
        setenv(names[0], invalid, 1);
        bool rejected = false;
        try { AdaptivePhi::Config cfg; } catch (const std::runtime_error &) { rejected = true; }
        assert(rejected);
        AdaptivePhi::Config inactive(false); // inactive modes ignore these variables
        assert(inactive.first == 75);
    }
    unsetenv(names[0]);
    for (const char *invalid : {"150", "500", "65537", "-1"})
    {
        setenv(names[1], invalid, 1);
        bool rejected = false;
        try { AdaptivePhi::Config cfg; } catch (const std::runtime_error &) { rejected = true; }
        assert(rejected);
    }
    setenv(names[1], "300", 1);
    AdaptivePhi::Config minimum;
    assert(minimum.maximum == 4*minimum.first);
    setenv(names[2], "nan", 1);
    bool rejected = false;
    try { AdaptivePhi::Config cfg; } catch (const std::runtime_error &) { rejected = true; }
    assert(rejected);
}
