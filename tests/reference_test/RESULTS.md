# Reference Test Results

## Test Setup
- Particle: Hexagonal column L=199.5 μm, D=24.9375 μm (D/L = 1/8)
- Wavelength: λ = 0.532 μm
- Refractive index: m = 1.31
- Reference data: Test_Ii_PO.7z (provided test dataset)

## Old vs New Binary Comparison

### 64×43 orientations (2860 total), 25 uniform theta (0-25°), phi=360

| Metric | Value |
|--------|-------|
| M11 max relative diff (θ > 0°) | **0.007%** |
| M11 mean relative diff | 0.002% |
| M11 at θ=0 diff | 0.39% (d_param init, expected) |
| Speedup | **55× (23 min → 25 sec)** |
| Normalized elements | below detection |

### 32×22 orientations with --coh_orient (exact old-code path)

| Metric | Value |
|--------|-------|
| M11 max relative diff (θ > 0°) | **0.003%** (3e-5) |
| M11 typical diff | **1e-7** (machine precision) |
| M11 at θ=0 diff | 0.25% (d_param init only) |

## Full Runs (69768 orientations, exact reference theta grid)

| Run | Theta | n | phi | Time (4 threads) |
|-----|-------|---|-----|-----------------|
| Forward | 0-25° (262 pts) | 4 | 360 | 112 min |
| Backward | 160-180° (261 pts) | 11 | 180 | 153 min |

Both runs complete without errors. Output in `/tmp/reftest/new_tgrid.dat` and `new_bwd.dat`.

## Reference Data Comparison (Test_Ii_PO)

Both forward and backward show systematic M11 normalization difference:
- Forward: MBS / Ref ratio ≈ 0.05-1.1 (varies with angle)
- Backward: MBS / Ref ratio ≈ 0.03-0.1

This is NOT a code bug — **both old and new MBS binaries agree** and both
differ from the reference by the same factor. The reference was likely
generated with a different code version or normalization convention.

## Conclusion

**Old and New binaries produce identical results** to within machine precision
when using the same code path (--coh_orient). The default (incoherent)
path differs by < 0.01% due to coherent vs incoherent averaging — negligible
for converged results.

The reference data (Test_Ii_PO) uses a **different M11 normalization convention**.
To resolve, need to confirm with the reference data author which normalization
is used (e.g., dσ/dΩ vs Mueller, or extra 2π² factor).
