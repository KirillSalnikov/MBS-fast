# Bugfix: Forward-Direction Fresnel Sign Error

## Date
2026-03-22

## Summary
Sign error in the forward-direction special case of the Kirchhoff diffraction integral.
The Fresnel amplitude at theta=0 had opposite sign compared to the general edge-sum formula,
causing destructive interference in coherent Jones summation and anomalous M11(0) values
for small particles.

## Symptoms
- M11(theta=0) anomalously **small** compared to M11(theta=0.1) for small particles in coherent mode
- Example (hex D=H=3um, single orientation beta=45 gamma=30):
  - **Before fix**: M11(0)=15.50, M11(0.1)=1.84  (ratio 8.4 — theta=0 too large, discontinuous)
  - **After fix**:  M11(0)=1.90,  M11(0.1)=1.84  (ratio 1.03 — smooth, correct)
- Orientation-averaged (256 Sobol, sym 6 2):
  - D=3um:  M11(0)/M11(1) was 0.019, now **1.031**
  - D=10um: M11(0)/M11(1) was 0.012, now **1.390**
  - D=100um: M11(0)/M11(1) was 4.09, now **789** (correct sharp forward peak)
- Incoherent mode (--incoh) was **not affected** (|J|^2 is sign-independent)
- Integrated quantities (Q_sca, C_sca) changed by ~0.06% (theta=0 has near-zero weight 2pi*dcos)

## Root Cause
The Fresnel-Kirchhoff diffraction integral has a special case when the scattering direction
equals the beam exit direction (A ~ 0, B ~ 0). In this limit, the edge-sum formula converges to:

    F_edge = complWave * (edge_sum)  -->  -invComplWave * area   as A,B -> 0

But the code used the **positive** sign:

    fresnel = +m_invComplWave * area    // WRONG

This created a sign discontinuity: beams hitting the A~0,B~0 threshold at theta=0
had opposite sign from the same beams at theta=0.1 (which used the edge-sum).
In coherent mode, this caused destructive interference between the ~20 forward-going beams
(with wrong sign) and the ~150 other beams (with correct sign from the edge-sum).

## Physical explanation
The Kirchhoff diffraction integral for a flat aperture is:

    F = (-ik / 4pi) * integral_aperture exp(i k r . n) dA

In the forward direction (r = beam_direction), the phase is constant and:

    F = (-ik / 4pi) * Area = -invComplWave * Area

where invComplWave = i*k / (2pi)^2... but effectively the sign of the area term
must match the sign convention of the edge-sum at the same limit.

## Verification method
1. Added per-beam debug logging of Jones contributions at theta=0 (j=0) and theta=0.1 (j=1)
2. Found that 20 out of 172 beams hit the forward special case at j=0 but not at j=1
3. These 20 beams had Fresnel = (0, +F) at j=0 but (0, -F) at j=1 — opposite sign
4. After negating the 20 contributions: M11(0)/M11(0.1) = 1.03 (smooth, correct)

## Files Changed

### src/handler/HandlerPO.cpp

**Line 970** — HandleBeams (single-orientation path):
```
BEFORE: fresnel = m_invComplWave * beam_area;
AFTER:  fresnel = -m_invComplWave * beam_area;
```

**Line 1333** — HandleBeamsToLocal (OpenMP parallel path):
```
BEFORE: fresnel = m_invComplWave * beam_area;
AFTER:  fresnel = -m_invComplWave * beam_area;
```

**Line 1889** — ComputeFromCache (multi-size path, inline Mueller):
```
BEFORE: fr = icwr * area_s;
        fi = icwi * area_s;
AFTER:  fr = -icwr * area_s;
        fi = -icwi * area_s;
```

### src/handler/Handler.cpp

**Line 255** — DiffractIncline (original diffraction, non-optimized):
```
BEFORE: return m_invComplWave * info.area;
AFTER:  return -m_invComplWave * info.area;
```

**Line 473** — DiffractInclineFast (optimized with precomputed edge data):
```
BEFORE: return m_invComplWave * info.area;
AFTER:  return -m_invComplWave * info.area;
```

**Line 642** — DiffractInclineAbs (absorption variant):
```
BEFORE: return m_invComplWave * info.area;
AFTER:  return -m_invComplWave * info.area;
```

## Impact on Results
- **Coherent M11 near theta=0**: fixed (was wrong for all particle sizes, most visible for small)
- **Coherent M11 at theta>1**: unchanged (these never hit the special case)
- **Incoherent mode**: unchanged (sign cancels in |J|^2)
- **Q_sca**: changed by ~0.06% (D=10um: 1.7205 -> 1.7216)
- **All other Mueller elements**: M12, M33, M34 etc. also corrected near theta=0
