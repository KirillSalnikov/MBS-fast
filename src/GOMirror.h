#pragma once

#include "ArgPP.h"

// Exact reflection/rotation domain for the built-in regular hexagonal column
// and the column with identical coaxial pyramidal cavities at both ends.
// Do not infer mirror symmetry from arbitrary particle-file headers.
inline bool SupportsGOMirror(const ArgPP& args)
{
    if (!args.IsCatched("p")) return false;
    const int type = args.GetIntValue("p", 0);
    if (type != 1 && type != 10) return false;
    // The existing absorption quadrature is not reflection-paired. Until it
    // has its own audit, do not promise identical numerical mirror partners.
    if (args.IsCatched("ri") && args.GetDoubleValue("ri",1)!=0) return false;
    const int beta = args.IsCatched("sym") ? args.GetIntValue("sym", 0) : 2;
    const int gamma = args.IsCatched("sym") ? args.GetIntValue("sym", 1) : 6;
    return (beta == 1 || beta == 2)
        && (gamma == 1 || gamma == 2 || gamma == 3 || gamma == 6);
}
