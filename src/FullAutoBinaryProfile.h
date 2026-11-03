#pragma once
#include <string>

// An explicit numerical receipt may bridge only these two phase profiles.
// Architecture, storage, compiler, OpenMP, MPI and all other version fields
// must still match. This predicate alone never authorizes checkpoint reuse.
inline bool FullAutoDoubleSingleProfileBridge(std::string oldProfile,
                                             std::string newProfile) {
    const std::string reduced = " phase-build-fp64 phase-trig-reduced-fp32 phase-fallback-fp64";
    const std::string ds = " phase-build-double-single phase-trig-double-single phase-trig-fallback-fp64 moments-double-single reduction-double-single";
    const auto a = oldProfile.find(reduced), b = newProfile.find(ds);
    if (a == std::string::npos || b == std::string::npos) return false;
    oldProfile.replace(a, reduced.size(), " phase-certified");
    newProfile.replace(b, ds.size(), " phase-certified");
    return oldProfile == newProfile;
}
