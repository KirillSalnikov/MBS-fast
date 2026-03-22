# MBS-raw: Physical Optics for Ice Crystals — Manual

## Overview

MBS-raw computes light scattering by non-spherical particles using the Physical Optics (PO) method. It traces geometric optics rays through the particle (reflections, refractions), then applies Kirchhoff diffraction to each output beam to obtain the far-field Mueller matrix.

**Optimized version**: ~70x faster than original, AVX-512/AVX2 SIMD, OpenMP parallel, Sobol quasi-random orientations. Builds for Intel, AMD Zen 2 (EPYC 7H12), and AMD Zen 4 (Ryzen 7000 / EPYC Genoa).

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

## CLI Reference

All flags are grouped by function. Required flags have no default; optional flags are marked with their defaults.

### Particle Definition

#### `-p TYPE HEIGHT DIAMETER [EXTRA]`

Define particle by type and dimensions (in micrometers).

| Type | Shape | Parameters |
|------|-------|------------|
| 1 | Hexagonal column | H = height, D = diameter |
| 2 | Bullet | H = height, D = diameter (cap angle 62 deg, computed automatically) |
| 3 | Bullet rosette | H, D, [optional cap height; default = D*sqrt(3)*tan(62 deg)/4] |
| 4 | Droxtal | arg3 = sup parameter (32.35 deg, 71.81 deg hardcoded) |
| 10 | Concave hexagonal | H, D, concavity depth |
| 12 | Hexagonal aggregate | H, D, number of elements |
| 999 | CertainAggregate | sup parameter (reads geometry from built-in definition) |

**Symmetry**: automatically detected from particle type.
- Hex column/plate: beta_sym = pi/2 (mirror top<->bottom), gamma_sym = pi/3 (6-fold rotation)
- Effective 12x symmetry reduction for orientation averaging.

**Example**: 10 um diameter column with aspect ratio D/L = 0.7:
```
-p 1 14.286 10    # H = D/0.7 = 14.286 um
```

#### `--pf FILENAME`

Load particle geometry from a file instead of using `-p`. The file defines vertices and facets. The `-p` argument is still parsed for the filename (legacy behavior: the filename is read from the `-p` value).

**Example**:
```bash
mbs_po --pf -p myparticle.dat --ri 1.31 0 -n 12 ...
```

#### `--rs SIZE`

Resize a particle loaded with `--pf`. Scales the particle so that its maximal dimension equals SIZE (in micrometers). Only valid when `--pf` is also specified.

**Example**:
```bash
mbs_po --pf --rs 50.0 -p myparticle.dat --ri 1.31 0 -n 12 ...
```

#### `--forced_convex`

Force the convex particle tracing algorithm, even if the particle would normally be detected as concave. Use when you know the particle is convex but auto-detection fails.

#### `--forced_nonconvex`

Force the non-convex (concave) particle tracing algorithm. Use for concave particles that are incorrectly detected as convex.

---

### Physical Parameters

#### `-w WAVELENGTH_UM`

Wavelength in micrometers. Required for PO diffraction and absorption calculations.

Examples:
- `-w 0.532` (green laser, 532 nm)
- `-w 1.064` (Nd:YAG, 1064 nm)

#### `--ri N_REAL N_IMAG`

Complex refractive index m = n_re + i*n_im.

- Ice at 532 nm: `--ri 1.31 0`
- Ice at 1064 nm: `--ri 1.30 1.9e-8`
- Absorbing ice (10 um IR): `--ri 1.197 0.248`

**Important**: When `--abs` is NOT specified, the imaginary part is set to zero regardless of the value provided. You must use `--abs` to enable absorption.

#### `--abs`

Enable absorption accounting. Requires `-w` to be specified. When enabled, the imaginary part of `--ri` is used; without `--abs`, the imaginary part is forced to zero.

Uses `DiffractInclineAbs` (optimized with `fast_sincos`) for the diffraction integrals.

**Example**:
```bash
mbs_po --po --abs -w 10.0 --ri 1.197 0.248 ...
```

#### `-n MAX_REFLECTIONS`

Maximum number of internal reflections allowed per ray.

| Value | Description | Use case |
|-------|-------------|----------|
| 1 | Entry refraction only | Testing, fast estimate |
| 5 | Standard | Forward scattering, C_sca (< 1% error for m=1.31) |
| 6 | Default | General purpose |
| 12 | High accuracy | Backscattering M11(180 deg) (< 5% error) |
| 20 | Maximum | Backscattering convergence (< 3% error) |
| 25 | Overkill | Reference calculations |

**Convergence at x=1000 (hex column D/L=0.7, m=1.31)**:
- C_sca converges at n=6 (error 0.4%)
- M11(180 deg) backscattering converges slowly: n=6 -> 11% error, n=12 -> 5%, n=20 -> 3%
- Time scales linearly: n=6 -> 135s, n=20 -> 459s (single thread)
- With OpenMP + AVX-512: n=20 costs only ~10s extra

**Recommendation**: use `-n 12` for general work, `-n 20` if backscattering matters.

---

### Method

#### `--po`

Required flag for Physical Optics method. Computes Kirchhoff diffraction on each beam.

#### `--go`

Geometrical Optics mode (also activated by omitting `--po`). Rays are mapped to far-field bins without diffraction. Faster but less accurate.

#### `--all`

Compute all trajectory groups. When `--tr` is specified, by default only the listed trajectory groups are computed; `--all` adds computation of all remaining trajectories as well.

Required for most standard runs when `--tr` is used. Without `--tr`, `--all` creates an empty catch-all group for total scattering.

#### `--incoh`

Incoherent beam summation. Each beam is converted to Mueller independently, then summed. Loses interference effects. Default is coherent summation (Jones amplitudes summed first). Use for testing.

#### `--karczewski`

Use Karczewski polarization matrix instead of RotateJones for beam polarization rotation. Experimental flag. Affects M33, M34, M44 elements but not M11. Both methods give the same |R| norm, so M11 is identical.

#### `--jones`

Output Jones matrices to file. Writes raw Jones matrix data for each scattering direction. Available in fixed-orientation PO mode.

---

### Orientations

One of the following orientation modes must be specified.

#### `--sobol N` (RECOMMENDED)

Generate N Sobol quasi-random orientations. N should be a power of 2 (32, 64, 128, ..., 8192, ...).

**Automatically uses particle symmetry**: beta in [0, beta_sym], gamma in [0, gamma_sym].
For hex prism: beta in [0, pi/2], gamma in [0, pi/3] -> 12x effective orientations.

Sobol convergence rate: O((log N)^2/N), much faster than random Monte Carlo O(1/sqrt(N)).

**How many orientations?**

| x | Forward peak, C_sca (< 2%) | Backscattering M11(180 deg) (< 5%) |
|---|---|---|
| 10-50 | 128 | 512 |
| 50-200 | 256 | 1024 |
| 200-1000 | 512 | 4096 |
| > 1000 | 1024 | 8192+ |

#### `--adaptive EPS` (EASIEST)

Automatically determine number of orientations. Starts with 256, doubles until C_sca converges to relative accuracy EPS.

Example: `--adaptive 0.01` -> 1% accuracy on C_sca.

Note: this optimizes for C_sca convergence. Backscattering may need more orientations.

#### `--random N_BETA N_GAMMA`

Uniform grid in beta x gamma over the symmetry-reduced domain. Traditional method.
Total orientations = (N_BETA+1) x N_GAMMA.

Example: `--random 20 60` -> 1260 orientations.

#### `--montecarlo N`

Monte Carlo random orientations. N is the total number of orientations, sampled randomly within the beta/gamma range (uses random seed). Less efficient than Sobol for the same N.

**Example**:
```bash
mbs_po --po --montecarlo 1000 --all \
    -p 1 10 10 -w 0.532 --ri 1.31 0 -n 12 \
    --grid 0 180 48 180 --close
```

#### `--fixed BETA GAMMA`

Fixed orientation (single orientation, no averaging). BETA and GAMMA are in degrees.

**Example**:
```bash
mbs_po --po --fixed 45 30 \
    -p 1 10 10 -w 0.532 --ri 1.31 0 -n 12 \
    --grid 0 180 48 180 --close
```

#### `--orientfile FILENAME`

Read (beta, gamma) pairs in radians from text file (one pair per line). Comments with `#`.

#### `--sym B G`

Override particle symmetry for Sobol/adaptive orientation generation.
B = beta symmetry divisor, G = gamma symmetry divisor.
- `--sym 6 2`: hex prism (beta in [0, pi/6], gamma in [0, pi]) — 12x reduction
- `--sym 2 6`: same total reduction, different split
- `--sym 1 1`: no symmetry (full sphere)

**Example**: `--sobol 1024 --sym 6 2` generates 1024 Sobol points in the symmetry-reduced domain.

#### `-b MIN MAX`

Override beta range (in degrees). By default, beta runs from 0 to beta_sym (particle symmetry). Use with `--random` or `--montecarlo`.

**Example**: `-b 0 90` restricts beta to [0 deg, 90 deg].

#### `-g MIN MAX`

Override gamma range (in degrees). By default, gamma runs from 0 to gamma_sym (particle symmetry). Use with `--random` or `--montecarlo`.

**Example**: `-g 0 60` restricts gamma to [0 deg, 60 deg].

---

### Scattering Grid

#### `--grid THETA_MIN THETA_MAX N_PHI N_THETA`

Uniform scattering angle grid.

- THETA_MIN, THETA_MAX: scattering angle range in degrees (usually 0 180)
- N_PHI: number of azimuthal bins (>= 48 recommended; 6 is insufficient!)
- N_THETA: number of zenith bins

**Important**: N_PHI = 6 (old default) causes 50-200% error in M11 at side/back angles! Use N_PHI >= 48.

Also accepts 3-parameter form: `--grid RADIUS N_PHI N_THETA` (backscattering cone of given radius around 180 deg).

#### `--tgrid FILENAME`

Read custom theta grid from file (one angle in degrees per line). Use `generate_theta_grid.py` to create optimal grids.

For x=100: 1800 uniform bins -> 156 non-uniform bins = **11x speedup**.

#### `--auto_tgrid` (RECOMMENDED)

Automatically generate optimal non-uniform theta grid based on size parameter x = pi*D/lambda:
- Fine zone (0 to 5x peak width): step = 0.1 x (180 deg/x)
- Transition zone: logarithmic spacing to 10 deg
- Coarse zone (10 deg to 180 deg): step 2 deg

#### `--point`

Compute only the exact backscatter point (theta = 180 deg). Requires `--po`. Currently commented out in the source code (non-functional); retained for future use.

#### `--filter ANGLE`

Scattering angle filter in degrees. Sets a backscattering cone aperture: only scattering within ANGLE degrees of the backscatter direction (180 deg) is accumulated. Used with `--all` or `--random` PO modes.

**Example**: `--filter 5` collects only scattering within 5 deg of exact backscatter.

---

### Multi-Size Computation (`--sizefile`)

```
--sizefile FILENAME
```

Compute multiple size parameters from one beam tracing pass. Uses `BeamCache`: beam topology is invariant across sizes (same shape, same RI), only vertex positions and phases scale with D.

**File format**: plain text, one size parameter (x = pi*D/lambda) value per line.

**Example file** (`sizes.txt`):
```
10
20
50
100
200
500
```

**Output**: `M_x10.dat`, `M_x20.dat`, etc.

**Performance note**: Uses `ComputeFromCache`, which recomputes diffraction integrals for each cached beam at each size. This adds overhead compared to running each size independently. For small particles (small x), separate runs may be faster.

**Recommendation**: Use `--sizefile` primarily for large particles (x > 500) with high reflection count (`-n` > 12), where the ray tracing phase dominates and caching provides the most benefit.

---

### Performance Options

#### `--beam_cutoff EPS`

Skip diffraction for beams with |J|^2 x area < EPS x C_geo (geometric cross-section).

Testing shows MBS-raw's built-in area cutoff in polygon clipping already handles this, so additional cutoff has negligible effect. Default threshold: 1e-12.

#### `--shadow_off`

Disable shadow beam (Babinet external diffraction). For debugging only.

#### `-r RATIO`

Restriction ratio for small beams during intersection. Beams smaller than 1/RATIO of the reference size are discarded. Default: 100.

Lowering this value keeps more small beams (slower but potentially more accurate for complex particles). Raising it discards more small beams (faster).

---

### Output

#### `-o DIRNAME`

Output folder name. Default: `M`.

Supports variable substitution with `%KEY` syntax, where KEY is any CLI argument name and an optional digit prefix selects the argument index. For example, `-o results_%0p_%1p` would substitute the first and second values of `-p` into the folder name.

**Example**: `-o results_hex` writes output to `results_hex/M.dat`.

#### `--close`

Close the program after calculation completes. Without this flag, the program waits for a keypress before exiting (interactive mode). **Required for batch/scripted runs and pipelines.**

#### `--log SECONDS`

Progress output interval in seconds. The program prints progress information (orientation count, elapsed time) to stderr at the specified interval.

**Example**: `--log 10` prints progress every 10 seconds.

#### `--tr FILENAME`

Load trajectory groups from file. A trajectory group specifies which internal reflection sequences to track separately. The file format defines facet index sequences for the particle.

When `--tr` is used without `--all`, only the specified trajectory groups are computed. Add `--all` to also compute all remaining (unmatched) trajectories.

#### `--gr`

Output per-group scattering results. When trajectory groups are defined (via `--tr`), writes separate output for each group. Legacy/diagnostic flag.

#### `--shadow`

Legacy flag. Appears in argument parser but has no documented effect in current code.

---

## Output Files

### M.dat

Tab-separated columns:
```
theta  2pi*dcos  M11  M12  M13  M14  M21  M22  M23  M24  M31  M32  M33  M34  M41  M42  M43  M44
```

- theta: scattering angle (degrees)
- 2pi*dcos: solid angle weight for integration. C_sca = sum of M11 x (2pi*dcos)
- M_{ij}: orientation-averaged Mueller matrix elements (phi-averaged)

### Diagnostics (stderr)

- Q_sca: scattering efficiency. Warning if Q_sca > 2.5 (PO overestimate)
- Phase timings (with OpenMP): Phase 1 (tracing), Phase 2 (diffraction)
- Beam count, direction evaluations

---

## Build

### Intel / AVX-512 (Skylake-X, Ice Lake, Sapphire Rapids, ...)

```bash
bash build.sh            # -> bin/mbs_po
```

Uses `-O3 -march=native -mavx512f -mavx512dq -fopenmp`. Requires GCC >= 9.

### AMD EPYC 7H12 (Zen 2, SP3)

```bash
bash build_epyc.sh       # -> bin/mbs_po_epyc
# Or with clang/AOCC (often 5-10% faster on Zen):
bash build_epyc_clang.sh # -> bin/mbs_po_epyc_clang
```

Uses `-O3 -march=znver2 -mtune=znver2 -fopenmp`. No AVX-512 — all SIMD via AVX2 + FMA3.
`fast_sincos_8x` falls back to 2 x `fast_sincos_4x`; rsqrt uses SSE `rsqrt_ss` + Newton-Raphson.

**NUMA optimization** (important for EPYC, 4+ NUMA domains per socket):

```bash
OMP_PROC_BIND=close OMP_PLACES=cores OMP_NUM_THREADS=64 \
    bin/mbs_po_epyc --po --sobol 4096 ...
```

Without `OMP_PROC_BIND=close` threads migrate between chiplets, losing L3 cache locality.
Use `numactl --interleave=all` if memory bandwidth is the bottleneck (rare for PO).

### AMD Zen 4 (Ryzen 7000, EPYC 9004 Genoa)

```bash
bash build_zen4.sh       # -> bin/mbs_po_zen4
```

Uses `-O3 -march=znver4 -mtune=znver4 -mavx512f -mavx512dq -mavx512vl -fopenmp`.
Zen 4 has AVX-512 with 256-bit execution units (ops decoded as 2x256-bit). Still benefits from wider register file and mask operations.

### Generic AVX2 (any x86-64 CPU with AVX2+FMA)

```bash
# Remove -mavx512f -mavx512dq from build.sh, or:
g++ -O3 -march=haswell -std=gnu++11 -funroll-loops -fopenmp \
    $(find src -not -path '*/bigint/*' -name '*.cpp') \
    $(find src/bigint -name '*.cc') -Isrc -Isrc/math ... -lm -lgomp
```

The code auto-detects AVX-512 at compile time (`#ifdef __AVX512F__`).

---

## Performance Summary

**Single thread** (hex D=H=10, 128 Sobol, 48 phi, 181 theta, n=6):

| Version | Phase 2 time | vs Original |
|---------|-------------|-------------|
| Original (before all optimizations) | ~320s | 1x |
| Optimized + batched sincos (current) | 4.6s | **~70x** |
| EPYC build (AVX2 only) | 4.7s | ~68x |

**Key optimizations**: PrecomputeEdgeData, vertex-cached sincos, batched sincos pre-computation (all thetas at once via AVX-512/AVX2), polynomial sincos (~14 digits), inline hot loop, no heap allocation, theta-coefficients (2 FMA per direction).

**OpenMP scaling** (hex D=H=10, 128 Sobol, current code):

| Threads | Phase 2 | Speedup |
|---------|---------|---------|
| 1 | 4.6s | 1.0x |
| 4 | 1.7s | 2.7x |
| 12 (SMT) | ~0.9s | ~5x |
| 64 (EPYC) | ~0.2s* | ~23x* |

*Estimated; actual EPYC scaling depends on orientation count (need >= 512 for 64 cores).

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

### Fixed orientation (single crystal)
```bash
mbs_po --po --fixed 45 30 \
    -p 1 14.3 10 -w 0.532 --ri 1.31 0 -n 12 \
    --grid 0 180 48 180 --close
```

### Absorbing particle
```bash
mbs_po --po --sobol 1024 --auto_tgrid --abs \
    -p 1 14.3 10 -w 10.0 --ri 1.197 0.248 -n 12 \
    --grid 0 180 48 180 --close
```

### Size scan
```bash
echo -e "10\n20\n50\n100\n200" > sizes.txt
mbs_po --po --sobol 1024 --auto_tgrid \
    -p 1 14.3 10 -w 0.532 --ri 1.31 0 -n 12 \
    --grid 0 180 48 180 --sizefile sizes.txt --close
```

### Custom particle from file
```bash
mbs_po --po --pf --rs 100.0 -p custom_crystal.dat \
    --ri 1.31 0 -n 12 --sobol 1024 --auto_tgrid \
    --grid 0 180 48 180 -w 0.532 --close
```

### Trajectory-resolved scattering
```bash
mbs_po --po --all --tr tracks.dat --gr \
    --sobol 1024 --auto_tgrid \
    -p 1 14.3 10 -w 0.532 --ri 1.31 0 -n 12 \
    --grid 0 180 48 180 --close
```

---

## Complete Flag Summary

| Flag | Arguments | Category | Description |
|------|-----------|----------|-------------|
| `-p` | TYPE H D [extra] | Particle | Particle type and dimensions |
| `--pf` | (none) | Particle | Load particle from file (filename via `-p`) |
| `--rs` | SIZE | Particle | Resize particle to maximal dimension SIZE |
| `--forced_convex` | (none) | Particle | Force convex tracing algorithm |
| `--forced_nonconvex` | (none) | Particle | Force non-convex tracing algorithm |
| `-w` | WAVELENGTH | Physical | Wavelength in micrometers |
| `--ri` | N_RE N_IM | Physical | Complex refractive index |
| `--abs` | (none) | Physical | Enable absorption (requires `-w`) |
| `-n` | N | Physical | Max internal reflections |
| `--po` | (none) | Method | Physical Optics mode |
| `--go` | (none) | Method | Geometrical Optics mode |
| `--all` | (none) | Method | Compute all trajectory groups |
| `--incoh` | (none) | Method | Incoherent beam summation |
| `--karczewski` | (none) | Method | Karczewski polarization matrix (experimental) |
| `--jones` | (none) | Method | Output Jones matrices |
| `--sobol` | N | Orientations | Sobol quasi-random (N = power of 2) |
| `--adaptive` | EPS | Orientations | Adaptive convergence to accuracy EPS |
| `--random` | N_BETA N_GAMMA | Orientations | Uniform grid |
| `--montecarlo` | N | Orientations | Monte Carlo random orientations |
| `--fixed` | BETA GAMMA | Orientations | Single orientation (degrees) |
| `--orientfile` | FILENAME | Orientations | Orientations from file (radians) |
| `--sym` | B G | Orientations | Symmetry override: beta/B, gamma/(2pi/G) |
| `-b` | MIN MAX | Orientations | Beta range override (degrees) |
| `-g` | MIN MAX | Orientations | Gamma range override (degrees) |
| `--grid` | T1 T2 N_PHI N_TH | Scattering | Uniform scattering grid |
| `--tgrid` | FILENAME | Scattering | Non-uniform theta grid from file |
| `--auto_tgrid` | (none) | Scattering | Auto-generate optimal theta grid |
| `--point` | (none) | Scattering | Backscatter point only (disabled) |
| `--filter` | ANGLE | Scattering | Backscatter cone aperture (degrees) |
| `--sizefile` | FILENAME | Multi-size | Multiple x values from file |
| `--beam_cutoff` | EPS | Performance | Beam importance cutoff |
| `--shadow_off` | (none) | Performance | Disable shadow beam |
| `-r` | RATIO | Performance | Small beam restriction ratio (default 100) |
| `-o` | DIRNAME | Output | Output folder (default "M") |
| `--close` | (none) | Output | Exit after calculation (for scripts) |
| `--log` | SECONDS | Output | Progress output interval |
| `--tr` | FILENAME | Output | Load trajectory groups from file |
| `--gr` | (none) | Output | Output per-group results |
| `--shadow` | (none) | Output | Legacy flag (no documented effect) |

---

## Bugfixes

### Forward-Direction Fresnel Sign (2026-03-22)

Fixed sign error in the forward-direction special case of the diffraction integral.
When A,B approach 0 (forward scattering), the edge-sum converges to `-invComplWave * area`,
but the shortcut formula used `+invComplWave * area`.

**Impact**: M11(theta=0) was anomalously small in coherent mode for small particles.
For D=3um hex: M11(0)/M11(0.1) was 0.019, now 1.031 (correct).
Incoherent mode and integrated quantities (Q_sca) were not significantly affected.

See `BUGFIX_forward_direction_sign.md` for full details.

---

## Known Limitations

1. **Q_sca > 2 at large x**: PO does not enforce the optical theorem. Q_sca grows with x instead of approaching 2. Use IGOM correction or renormalize.

2. **Polarization elements M33, M34, M44**: MBS-raw's RotateJones and GOAD's Karczewski matrix have the same norm but different structure. M11 agrees between codes, but M33/M34/M44 differ.

3. **Small particles (x < 20)**: PO is not valid. Use ADDA or T-matrix methods.

4. **Backscattering convergence**: M11(180 deg) requires 10-100x more orientations than forward scattering.

5. **LTO**: `-flto` causes ~20% regression on current code due to GCC inlining heuristics interacting poorly with the hand-inlined hot loop. Do not use.

5. **Absorption**: DiffractInclineAbs is optimized but not as fast as the non-absorbing path.

---

## Verification

**Backward compatibility** (v2 vs original):
- M11 with same parameters (6 phi bins): **max rdiff = 0.04%** (effectively identical)
- M12, M22: up to 18% difference from BAC-CAB optimization in RotateJones (affects off-diagonal only)
- M11 is physics-invariant under RotateJones changes (same norm ||R|| = ||R_old||)

**Effect of phi bins**:
- 6 phi bins (old default): M11 overestimated by 37-59% at side/back angles
- 48 phi bins (recommended): correct azimuthal averaging
- This is NOT a code change -- it's a parameter choice

**Q_sca normalization**:
- `--sobol` and `--random` modes: Q_sca uses actual projected area (m_incomingEnergy) -> **Q ~ 2.0**
- Python scripts using pi*(D/2)^2: Q inflated by factor A_hex/A_sphere ~ 1.57 for D/L=0.7
- Always use Q_sca from MBS stderr output (correct normalization)
