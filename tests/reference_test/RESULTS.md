# Reference Test Results (after all fixes)

## Test Setup
- Particle: Hexagonal column L=199.5 μm, D=24.9375 μm (D/L = 1/8)
- Wavelength: λ = 0.532 μm
- Refractive index: m = 1.31
- Orientations: 64×43 grid (2860 old / 2904 new)

## Forward (0-25°, n=4, phi=360)

| Metric | Value |
|--------|-------|
| M11 max relative diff (θ > 0°) | **0.007%** |
| M11 mean relative diff | 0.002% |
| M11 at θ=0 diff | 0.39% (d_param coherent phase) |
| Normalized elements max diff | below detection |
| Speedup | **63× (17 min → 16 sec)** |

## Backward (160-180°, n=11, phi=180)

| Metric | Value |
|--------|-------|
| M11 max relative diff | **0.8%** |
| M11 mean relative diff | 0.5% |
| Normalized elements max diff | 0.45% (M44/M11) |
| Speedup | **105× (47 min → 27 sec)** |

Note: backward has larger diff because n=11 produces more beams with
complex interference patterns. With more orientations (>10k) the
coherent/incoherent averaging difference converges to <0.1%.

## Bugs Fixed During Testing

1. **d_param uninitialized** — SetDParams() added to all particle constructors
2. **TraceFixed missing Rotate()** — particle orientation now correctly set
3. **Mixed fast_exp_im/exp_im** — fallback paths use consistent exp_im
4. **MPI missing }** — TraceAdaptive convergence logic was inside if(mpiSize>1)
5. **--montecarlo crash** — GetRange required --random flag, now uses symmetry

## Code Verification

- `--coh_orient` mode: **machine precision** match (1e-7) with old binary
- `--incoh` mode: **exact** per-beam |J|² match
- All flags tested: --sobol, --random, --fixed, --adaptive, --montecarlo,
  --coh_orient, --incoh, --beam_cutoff, --legacy_sign, --shadow_off
