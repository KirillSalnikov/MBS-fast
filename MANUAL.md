# MBS-raw: Physical Optics for Ice Crystals — Manual

## Overview

MBS-raw computes light scattering by non-spherical particles using the Physical Optics (PO) method. It traces geometric optics rays through the particle (reflections, refractions), then applies Kirchhoff diffraction to each output beam to obtain the far-field Mueller matrix.

**Optimized version**: 36× faster than original, AVX-512 SIMD, OpenMP parallel, Sobol quasi-random orientations.

---

## Quick Start

```bash
# Simple run: hex column, Sobol orientations, auto theta grid
mbs_po --po --sobol 1024 --auto_tgrid \
    -p 1 10 10 -w 0.532 --ri 1.31 0 -n 12 \
    --grid 0 180 48 180 --close

# Fully adaptive (auto everything)
mbs_po --po --adaptive 0.01 --auto_tgrid \
    -p 1 10 10 -w 0.532 --ri 1.31 0 -n 12 \
    --grid 0 180 48 180 --close
```

---

## Particle Shape (`-p`)

```
-p TYPE HEIGHT DIAMETER
```

| Type | Shape | Description |
|------|-------|-------------|
| 1 | Hexagonal column | H = height, D = diameter |
| 2 | Hexagonal plate | H = height, D = diameter |
| 3 | Bullet rosette | H = height, D = diameter |

**Symmetry**: automatically detected from particle type.
- Hex column/plate: β_sym = π/2 (mirror top↔bottom), γ_sym = π/3 (6-fold rotation)
- Effective 12× symmetry reduction for orientation averaging.

**Example**: 10 μm diameter column with aspect ratio D/L = 0.7:
```
-p 1 14.286 10    # H = D/0.7 = 14.286 μm
```

---

## Physical Parameters

### Wavelength (`-w`)
```
-w WAVELENGTH_UM
```
Wavelength in micrometers. Example: `-w 0.532` (green laser, 532 nm).

### Refractive Index (`--ri`)
```
--ri N_REAL N_IMAG
```
Complex refractive index m = n_re + i·n_im.
- Ice at 532 nm: `--ri 1.31 0`
- Ice at 1064 nm: `--ri 1.30 1.9e-8`
- Absorbing ice (10 μm IR): `--ri 1.197 0.248`

When n_im > 0, absorption is computed via `DiffractInclineAbs` (optimized with `fast_sincos`).

---

## Method Selection

### Physical Optics (`--po`)
```
--po
```
Required flag for PO method. Computes Kirchhoff diffraction on each beam.

### Geometrical Optics (`--go`)
No `--po` flag → GO mode. Rays mapped to far-field bins without diffraction. Faster but less accurate.

---

## Internal Reflections (`-n`)

```
-n MAX_REFLECTIONS
```

Maximum number of internal reflections allowed per ray.

| Value | Description | Use case |
|-------|-------------|----------|
| 1 | Entry refraction only | Testing, fast estimate |
| 5 | Standard | Forward scattering, C_sca (< 1% error for m=1.31) |
| 6 | Default | General purpose |
| 12 | High accuracy | Backscattering M₁₁(180°) (< 5% error) |
| 20 | Maximum | Backscattering convergence (< 3% error) |
| 25 | Overkill | Reference calculations |

**Convergence at x=1000 (hex column D/L=0.7, m=1.31)**:
- C_sca converges at n=6 (error 0.4%)
- M₁₁(180°) backscattering converges slowly: n=6 → 11% error, n=12 → 5%, n=20 → 3%
- Time scales linearly: n=6 → 135s, n=20 → 459s (single thread)
- With OpenMP + AVX-512: n=20 costs only ~10s extra

**Recommendation**: use `-n 12` for general work, `-n 20` if backscattering matters.

---

## Orientation Averaging

### Sobol Quasi-Random (`--sobol`) ⭐ RECOMMENDED

```
--sobol N
```

Generate N Sobol quasi-random orientations. N should be a power of 2 (32, 64, 128, ..., 8192, ...).

**Automatically uses particle symmetry**: β ∈ [0, β_sym], γ ∈ [0, γ_sym].
For hex prism: β ∈ [0, π/2], γ ∈ [0, π/3] → 12× effective orientations.

Sobol convergence rate: O((log N)²/N), much faster than random Monte Carlo O(1/√N).

**How many orientations?**

| x | Forward peak, C_sca (< 2%) | Backscattering M₁₁(180°) (< 5%) |
|---|---|---|
| 10–50 | 128 | 512 |
| 50–200 | 256 | 1024 |
| 200–1000 | 512 | 4096 |
| > 1000 | 1024 | 8192+ |

### Adaptive (`--adaptive`) ⭐ EASIEST

```
--adaptive EPS
```

Automatically determine number of orientations. Starts with 256, doubles until C_sca converges to relative accuracy EPS.

Example: `--adaptive 0.01` → 1% accuracy on C_sca.

Note: this optimizes for C_sca convergence. Backscattering may need more orientations.

### Grid (`--random`)

```
--random N_BETA N_GAMMA
```

Uniform grid in β × γ over the symmetry-reduced domain. Traditional method.
Total orientations = (N_BETA+1) × N_GAMMA.

Example: `--random 20 60` → 1260 orientations.

### From File (`--orientfile`)

```
--orientfile FILENAME
```

Read (β, γ) pairs in radians from text file (one pair per line). Comments with `#`.

---

## Scattering Grid

### Uniform (`--grid`)

```
--grid THETA_MIN THETA_MAX N_PHI N_THETA
```

- THETA_MIN, THETA_MAX: scattering angle range in degrees (usually 0 180)
- N_PHI: number of azimuthal bins (≥ 48 recommended; 6 is insufficient!)
- N_THETA: number of zenith bins

**Important**: N_PHI = 6 (old default) causes 50–200% error in M₁₁ at side/back angles! Use N_PHI ≥ 48.

### Non-Uniform Theta (`--tgrid`)

```
--tgrid FILENAME
```

Read custom theta grid from file (one angle in degrees per line). Use `generate_theta_grid.py` to create optimal grids.

For x=100: 1800 uniform bins → 156 non-uniform bins = **11× speedup**.

### Auto Theta Grid (`--auto_tgrid`) ⭐ RECOMMENDED

```
--auto_tgrid
```

Automatically generate optimal non-uniform theta grid based on size parameter x = πD/λ:
- Fine zone (0 to 5× peak width): step = 0.1 × (180°/x)
- Transition zone: logarithmic spacing to 10°
- Coarse zone (10° to 180°): step 2°

---

## Coherence

### Coherent Summation (default)

All beams from one orientation are summed coherently (Jones amplitudes), then converted to Mueller. This is the correct PO behavior.

### Incoherent (`--incoh`)

```
--incoh
```

Each beam → Mueller independently, then sum. Faster, but loses interference effects. Use for testing.

---

## Multi-Size Computation (`--sizefile`)

```
--sizefile FILENAME
```

Compute multiple size parameters from one beam tracing pass. File contains one x value per line.

Uses `BeamCache`: beam topology is invariant across sizes (same shape, same RI), only vertex positions and phases scale with D.

Output: `M_x10.dat`, `M_x20.dat`, etc.

---

## Performance Options

### Beam Importance Cutoff (`--beam_cutoff`)

```
--beam_cutoff EPS
```

Skip diffraction for beams with |J|² × area < EPS × C_geo (geometric cross-section).

Testing shows MBS-raw's built-in area cutoff in polygon clipping already handles this, so additional cutoff has negligible effect. Default threshold: 1e-12.

### Shadow Beam (`--shadow_off`)

```
--shadow_off
```

Disable shadow beam (Babinet external diffraction). For debugging.

---

## Output

### M.dat

Tab-separated columns:
```
theta  2pi*dcos  M11  M12  M13  M14  M21  M22  M23  M24  M31  M32  M33  M34  M41  M42  M43  M44
```

- theta: scattering angle (degrees)
- 2pi*dcos: solid angle weight for integration. C_sca = Σ M₁₁ × (2π·dcos)
- M_{ij}: orientation-averaged Mueller matrix elements (phi-averaged)

### Diagnostics (stderr)

- Q_sca: scattering efficiency. Warning if Q_sca > 2.5 (PO overestimate)
- Phase timings (with OpenMP): Phase 1 (tracing), Phase 2 (diffraction)
- Beam count, direction evaluations

---

## Build

```bash
# Optimized build (AVX-512, LTO, OpenMP)
make

# Or manually:
g++ -O3 -march=native -std=gnu++11 -funroll-loops -flto \
    -mavx512f -mavx512dq -fopenmp \
    -Isrc -Isrc/math -Isrc/handler -Isrc/common -Isrc/geometry \
    -Isrc/geometry/intrinsic -Isrc/geometry/sse -Isrc/particle \
    -Isrc/scattering -Isrc/tracer -Isrc/splitting -Isrc/bigint \
    -o bin/mbs_po $(find src -name '*.cpp') \
    $(find src/bigint -name '*.cc' 2>/dev/null) -lm -lgomp
```

**Requirements**: GCC ≥ 11, AVX-512 CPU (AMD Zen 4+, Intel Skylake-X+)

For non-AVX-512 CPU: remove `-mavx512f -mavx512dq`, polynomial sincos still works via AVX2.

---

## Performance Summary

**Single thread** (x=20, 256 orient, 48 phi, n=5):

| Version | Time | Speedup |
|---------|------|---------|
| Original | ~320s | 1× |
| Optimized (this version) | 8.9s | **36×** |
| GOAD (Rust, reference) | 11s | — |

**Key optimizations**: PrecomputeEdgeData, vertex-cached sincos, polynomial sincos (~14 digits, AVX-512 FMA), inline hot loop, no heap allocation, branchless edge loop, LTO.

**OpenMP** (6 physical + 6 SMT cores):

| Threads | Speedup |
|---------|---------|
| 6 | 3.65× |
| 12 (SMT) | 4.89× |

---

## Typical Workflows

### Quick test
```bash
mbs_po --po --sobol 256 -p 1 10 10 -w 0.532 --ri 1.31 0 -n 6 \
    --grid 0 180 48 90 --close
```

### Production (forward + side scattering)
```bash
mbs_po --po --sobol 1024 --auto_tgrid \
    -p 1 14.3 10 -w 0.532 --ri 1.31 0 -n 12 \
    --grid 0 180 48 180 --close
```

### Production (with backscattering)
```bash
mbs_po --po --sobol 4096 --auto_tgrid \
    -p 1 14.3 10 -w 0.532 --ri 1.31 0 -n 20 \
    --grid 0 180 48 180 --close
```

### Adaptive (automatic)
```bash
mbs_po --po --adaptive 0.01 --auto_tgrid \
    -p 1 14.3 10 -w 0.532 --ri 1.31 0 -n 12 \
    --grid 0 180 48 180 --close
```

### Size scan
```bash
echo -e "10\n20\n50\n100\n200" > sizes.txt
mbs_po --po --sobol 1024 --auto_tgrid \
    -p 1 14.3 10 -w 0.532 --ri 1.31 0 -n 12 \
    --grid 0 180 48 180 --sizefile sizes.txt --close
```

---

## Known Limitations

1. **Q_sca > 2 at large x**: PO does not enforce the optical theorem. Q_sca grows with x instead of approaching 2. Use IGOM correction or renormalize.

2. **Polarization elements M₃₃, M₃₄, M₄₄**: MBS-raw's RotateJones and GOAD's Karczewski matrix have the same norm but different structure. M₁₁ agrees between codes, but M₃₃/M₃₄/M₄₄ differ.

3. **Small particles (x < 20)**: PO is not valid. Use ADDA or T-matrix methods.

4. **Backscattering convergence**: M₁₁(180°) requires 10–100× more orientations than forward scattering.

5. **Absorption**: DiffractInclineAbs is optimized but not as fast as the non-absorbing path.

---

## Verification

**Backward compatibility** (v2 vs original):
- M₁₁ with same parameters (6 phi bins): **max rdiff = 0.04%** (effectively identical)
- M₁₂, M₂₂: up to 18% difference from BAC-CAB optimization in RotateJones (affects off-diagonal only)
- M₁₁ is physics-invariant under RotateJones changes (same norm ‖R‖ = ‖R_old‖)

**Effect of phi bins**:
- 6 phi bins (old default): M₁₁ overestimated by 37-59% at side/back angles
- 48 phi bins (recommended): correct azimuthal averaging
- This is NOT a code change — it's a parameter choice

**Q_sca normalization**:
- `--sobol` and `--random` modes: Q_sca uses actual projected area (m_incomingEnergy) → **Q ≈ 2.0**
- Python scripts using π(D/2)²: Q inflated by factor A_hex/A_sphere ≈ 1.57 for D/L=0.7
- Always use Q_sca from MBS stderr output (correct normalization)
