#include "FullAutoBinaryProfile.h"
#include <cassert>
int main() {
    const std::string base="\nbackend-build: CUDA sm_86 fp32-diffraction-storage trace-geometry-fp64";
    const std::string tail=" precise-math\ncompiler: GCC 7.5.0\nopenmp: enabled\nmpi: disabled\n";
    const auto old=base+" phase-build-fp64 phase-trig-reduced-fp32 phase-fallback-fp64"+tail;
    const auto ds=base+" phase-build-double-single phase-trig-double-single phase-trig-fallback-fp64 moments-double-single reduction-double-single"+tail;
    assert(FullAutoDoubleSingleProfileBridge(old,ds));
    assert(!FullAutoDoubleSingleProfileBridge(ds,old));
    assert(!FullAutoDoubleSingleProfileBridge(old,old));
    for (const auto &pair : {std::make_pair("sm_86","sm_120"),std::make_pair("GCC 7.5.0","GCC 13.1"),std::make_pair("fp32-diffraction-storage","fp64-diffraction-storage"),std::make_pair("precise-math","fast-math"),std::make_pair("mpi: disabled","mpi: enabled")}) {
        auto changed=ds; changed.replace(changed.find(pair.first),std::string(pair.first).size(),pair.second);
        assert(!FullAutoDoubleSingleProfileBridge(old,changed));
    }
}
