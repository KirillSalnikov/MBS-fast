#pragma once

#include <cmath>

#ifdef __CUDACC__
#define MBS_PHASE_INLINE __host__ __device__ inline
#else
#define MBS_PHASE_INLINE inline
#endif

namespace MbsGpuPhase {
// Preserve the double-precision phase until after range reduction. Casting a
// large unreduced phase to float loses physically relevant interference bits.
// Two explicitly fused subtractions retain the low part of 2*pi. Outside the
// tested range, the caller must use the original FP64 trigonometry.
MBS_PHASE_INLINE bool Reduce(double phase, double& reduced)
{
    if (!(::fabs(phase) <= 1048576.0))
        return false; // Also declines NaN and infinities.
    const double turns = ::nearbyint(phase * 0.15915494309189533577);
    reduced = ::fma(-turns, 6.28318530717958623200, phase);
    reduced = ::fma(-turns, 2.44929359829470641435e-16, reduced);
    return true;
}
}

#undef MBS_PHASE_INLINE
