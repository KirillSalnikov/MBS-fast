# MBS-fast Manual

MBS-fast computes light scattering by non-spherical particles with geometrical ray tracing plus Physical Optics (PO) diffraction. The code is optimized for CPU OpenMP/MPI runs and CUDA GPU diffraction runs. This manual describes the current command line interface, build targets, grids, multi-GPU scans, and the main diffraction/Jones/Mueller formulas used by the implementation.

## Contents

| Section | What it covers |
|---|---|
| [Build](#build) | CPU, GPU, precision variants, EPYC Zen targets, AVX2/AVX-512 |
| [Quick runs](#quick-runs) | Minimal CPU/GPU examples |
| [Physical model](#physical-model) | Ray tracing, Jones dyads, diffraction integral, Mueller conversion |
| [Particles and sizes](#particles-and-sizes) | `-p`, `--pf`, `--rs`, `--k_eq`, refractive index |
| [Orientation grids](#orientation-grids) | `--oldauto`, `--random`, Sobol, Euler, pole shortcut |
| [Scattering grids](#scattering-grids) | `--grid`, `--tgrid`, `--nphi`, integration weights |
| [GPU and multi-GPU](#gpu-and-multi-gpu) | CUDA backend, FFT backend, size scans, device assignment |
| [All flags](#all-flags) | Every parsed CLI flag |
| [Environment variables](#environment-variables) | Production and diagnostic environment controls |
| [Output](#output) | `.dat` layout and common diagnostics |

## Build

### Requirements

| Component | Required for | Notes |
|---|---|---|
| GCC >= 9 or Clang >= 14 | CPU and host CUDA code | GCC is the normal path on EPYC |
| OpenMP | CPU parallelism | GCC links with `libgomp` |
| MPI | `cpu/` split build | `mpicxx` is used when available |
| CUDA toolkit | `gpu/` split build | `nvcc`, `libcudart`, `libcufft` |
| NVIDIA driver | GPU runs | Must support the requested `sm_XX` binary |

The split build is the recommended build path. CPU objects are kept under `cpu/build/`; CUDA objects are kept under `gpu/build/`. The shared physics code remains in `src/`.

### CPU build

```bash
# MPI + OpenMP CPU binary
make -C cpu -j

# Run one MPI rank with 64 OpenMP threads
OMP_NUM_THREADS=64 cpu/bin/mbs_po_mpi --po -p 1 125.9 78.09 \
    --ri 1.3116 0 -w 0.532 -n 12 --oldauto 2 \
    --grid 0 180 600 180 --threads 64 --close -o out_cpu

# Run four MPI ranks, each with 16 OpenMP threads
mpirun -np 4 cpu/bin/mbs_po_mpi --po -p 1 125.9 78.09 \
    --ri 1.3116 0 -w 0.532 -n 12 --oldauto 2 \
    --grid 0 180 600 180 --threads 16 --close -o out_cpu_mpi
```

Debug help build:

```bash
make -C cpu debug
cpu/bin/mbs_po_mpi_debug --help-debug
```

### GPU build

```bash
# Default/reference GPU binary: FP64 with precise math
PATH=/usr/local/cuda/bin:$PATH make -C gpu -j

# Explicit variants
make -C gpu fp32      -j   # gpu/bin/mbs_po_gpu_float
make -C gpu fp64      -j   # gpu/bin/mbs_po_gpu_double
make -C gpu fp32_fast -j   # gpu/bin/mbs_po_gpu_float_fast
make -C gpu fp64_fast -j   # gpu/bin/mbs_po_gpu_double_fast
```

`gpu/Makefile` also works when `nvcc` is not in `PATH`: it uses
`$CUDA_PATH/bin/nvcc` (default `CUDA_PATH=/usr/local/cuda`). It reads the host
GCC limit from the installed CUDA headers and selects the newest compatible
versioned compiler, such as `g++-13` for CUDA 12.6. Override the detection with
`CUDA_HOST_CXX=/path/to/g++` when needed. On distributions whose default GCC is
newer than the CUDA limit, install a supported version alongside the system
compiler; changing the system-wide GCC alternative is not required.

| Target | Binary | CUDA precision | Fast math | Typical use |
|---|---|---:|---:|---|
| `make -C gpu` or `make -C gpu fp64` | `gpu/bin/mbs_po_gpu_double` | FP64 | no | Numerical reference and production default |
| `make -C gpu fp32` | `gpu/bin/mbs_po_gpu_float` | FP32 | no | Calibrated FP32 throughput |
| `make -C gpu fp32_fast` | `gpu/bin/mbs_po_gpu_float_fast` | FP32 | yes | Maximum throughput after validation |
| `make -C gpu fp64_fast` | `gpu/bin/mbs_po_gpu_double_fast` | FP64 | yes | FP64 throughput after validation |
| `make -C gpu double_debug` | `gpu/bin/mbs_po_gpu_double_debug` | FP64 | yes | Debug help and diagnostics |

FP32 applies to diffraction-beam storage, Jones sums, and Mueller output.
Visibility/topology tracing, optical paths, absorption, cancellation-prone
polygon moments, and critical phase trigonometry remain FP64. `--version`
reports the storage, critical-phase, math, and target-architecture profile.

The GPU split build enables CUDA diffraction by default. `--gpu` is accepted but optional. Use `--cpu` only when you deliberately want the CPU diffraction backend from a GPU-capable binary.

The Makefile detects the GPU compute capability with `nvidia-smi`. Override when needed:

```bash
# Ampere, e.g. RTX 3080 Ti
make -C gpu double_fast -j GPU_ARCH=86

# Ada, e.g. RTX 4070
make -C gpu float_fast -j GPU_ARCH=89
```

### EPYC Zen5, AVX-512, and AVX2

`scripts/detect_arch_flags.sh` is used by the Makefiles through `ARCH_FLAGS`. On a native machine, the safest high-performance option is usually:

```bash
make -C cpu clean
make -C cpu -j ARCH_FLAGS="-march=native -mtune=native"

make -C gpu clean
make -C gpu double_fast -j ARCH_FLAGS="-march=native -mtune=native"
```

For named AMD targets:

| CPU target | GCC/Clang flags | Notes |
|---|---|---|
| EPYC Zen2 | `ARCH_FLAGS="-march=znver2 -mtune=znver2"` | AVX2/FMA path |
| EPYC Zen3 | `ARCH_FLAGS="-march=znver3 -mtune=znver3"` | AVX2/FMA path |
| EPYC Zen4 | `ARCH_FLAGS="-march=znver4 -mtune=znver4"` | AVX-512 capable on GCC versions that support it |
| EPYC Zen5 | `ARCH_FLAGS="-march=znver5 -mtune=znver5"` | Use with a compiler that knows `znver5` |
| Portable AVX2 | `ARCH_FLAGS="-O3 -mavx2 -mfma"` | Runs on AVX2 machines |
| Explicit AVX-512 | `ARCH_FLAGS="-O3 -mavx512f -mavx512dq -mavx512cd -mavx512bw -mavx512vl"` | Use only on machines that support these flags |

There is no GCC/Clang option named `AVX12`. If the intended target is Zen5 with AVX-512, use `-march=znver5` when the compiler supports it. If the compiler is older, either upgrade GCC/Clang or use the explicit AVX-512 flags above after checking `lscpu | grep avx512`.

### Legacy root Makefile

The root `Makefile` still exists for compatibility:

```bash
make                         # bin/mbs_po, CPU
make gpu_double_fast         # delegates to gpu/double_fast
make gpu_float_fast          # delegates to gpu/float_fast
make USE_CUDA=1              # legacy single-tree CUDA build
```

Prefer the split `cpu/` and `gpu/` builds for normal work.

## Quick runs

```bash
# GPU double, oldauto orientation grid, full 0..180 degree output
gpu/bin/mbs_po_gpu_double_fast --po -p 1 125.9 78.09 \
    --ri 1.3116 0 -w 0.532 -n 12 --oldauto 2 \
    --grid 0 180 600 180 --threads 64 --close -o hex_gpu_double

# Same run, force CPU backend from the GPU binary
gpu/bin/mbs_po_gpu_double_fast --po --cpu -p 1 125.9 78.09 \
    --ri 1.3116 0 -w 0.532 -n 12 --oldauto 2 \
    --grid 0 180 600 180 --threads 64 --close -o hex_cpu_from_gpu_bin

# File particle scaled by equivalent size parameter
gpu/bin/mbs_po_gpu_float_fast --po --pf particle.dat --k_eq 58.81 \
    --ri 1.6 0.002 -w 1.064 -n 14 --oldauto 2 --pole \
    --grid 0 180 600 180 --threads 16 --close -o particle_keq

# Two-point forward/backward check
gpu/bin/mbs_po_gpu_double_fast --po -p 1 125.9 78.09 \
    --ri 1.3116 0 -w 0.532 -n 1 --oldauto 2 --pole \
    --grid 0 180 600 1 --threads 64 --close -o check_0_180
```

## Physical model

### Pipeline

| Stage | Main code path | Result |
|---|---|---|
| Particle setup | `src/particle`, `src/main.cpp` | Facets, normals, area, volume, symmetry |
| Geometrical tracing | `src/scattering`, `src/tracer`, `src/Splitting.cpp` | Output beams with direction, polygon, area, optical path, Jones matrix |
| PO diffraction | `HandlerPO::ApplyDiffraction*`, `Handler::DiffractIncline*`, CUDA equivalents | Far-field Jones amplitude per beam and direction |
| Coherent sum | `HandlerPO::AddToMueller`, GPU fused Mueller kernels | Sum Jones matrices over beams for one orientation |
| Mueller conversion | `src/math/Mueller.cpp` | 4x4 Mueller matrix |
| Orientation average | `HandlerPOTotal`, `TracerPOTotal` | Final randomly oriented Mueller matrix |

### Ray Jones matrices

Each traced beam carries a 2x2 complex Jones matrix `beam.J`. Refraction and reflection multiply it by Fresnel coefficients in `Splitting.cpp`. The two polarization channels are the local vertical/horizontal basis of the incident and outgoing ray.

For a beam with local Jones matrix

```text
J = [ Jvv  Jvh ]
    [ Jhv  Jhh ],
```

the code treats `J` as the coherent amplitude transform accumulated along the ray path. Cutoffs based on `|J|^2`, area, or `|J|^2*area` can remove weak beams before expensive diffraction.

### Dyadic rotation into the scattering basis

Before adding a beam to a detector direction, the beam-local Jones matrix must be expressed in the detector polarization basis. `HandlerPO::RotateJones` and `RotateJonesFast` build this 2x2 rotation/projection using dot products between:

| Vector | Meaning |
|---|---|
| `beam.Direction()` | Beam propagation direction after tracing |
| `direction` | Requested far-field scattering direction |
| `vf` | Forward reference direction |
| `info` polarization vectors | Precomputed beam polarization basis |

The effective far-field Jones contribution has the structure used in `ApplyDiffractionFast`:

```text
J_beam(theta, phi) = F_edge(theta, phi) * R_out(theta, phi) * F_n(theta, phi)
```

where:

| Factor | Code | Meaning |
|---|---|---|
| `F_edge` | `DiffractInclineFast`, `DiffractIncline`, `DiffractInclineAbs` | Scalar Kirchhoff edge diffraction integral for the beam polygon |
| `R_out` | `RotateJones` or `KarczewskiJones` | Jones-basis rotation/projection from beam coordinates to scattering coordinates |
| `F_n` | `ComputeFnJones` | Fresnel/Jones correction for the beam normal and direction |

### Karczewski polarization flag

`--karczewski` replaces the default `RotateJones` basis projection with
`HandlerPO::KarczewskiJones`. It changes only the outgoing polarization
rotation/projection factor `R_out`; tracing, beam areas, optical paths,
Fresnel coefficients, diffraction amplitude, and the scattering grid are left
unchanged.

This flag is mainly a polarization-convention and validation tool, not a
universal accuracy switch. The current implementation is marked experimental
in the code because it follows the Karczewski aperture-frame matrix but does
not fully reproduce the whole GOAD coordinate pipeline. Use it when comparing
polarization-sensitive Mueller elements or when a reference calculation uses
the Karczewski convention.

The transverse norms in this branch are evaluated with `hypot`, and every
basis is validated before normalization. A degenerate pole configuration
falls back to the default `RotateJones`; a non-finite experimental result is
never accumulated into the Mueller matrix. This numerical safeguard does not
turn the incomplete GOAD convention implementation into a physical reference.

Expected behavior:

| Quantity | Expected effect |
|---|---|
| `M11` | Should stay unchanged. `M11 = 0.5 * ||J||_F^2`, and the Karczewski/default rotations preserve the same Frobenius norm. |
| `M12`, `M22` | Usually less sensitive than the lower polarization block, but can still shift if the reference basis convention differs. |
| `M33`, `M34`, `M44` | Can change noticeably because the Jones matrix structure in the scattering basis changes. |
| Integral quantities dominated by intensity | Should not be interpreted as improved just because `--karczewski` is enabled; check the Mueller elements directly. |

Practical rule: if only `M33/M34/M44` disagree while `M11` agrees, run the same
case with and without `--karczewski`. If `--karczewski` moves the result toward
the reference, the mismatch is likely a polarization-basis convention issue.
If `M11` changes, the difference is not explained by this flag and points to a
different geometry, diffraction, cutoff, or normalization problem.

### Kirchhoff diffraction integral

For each traced output beam the illuminated aperture is a polygon. The PO far field is evaluated by a boundary/edge form of the Kirchhoff integral. In code this lives in:

| Function | Use |
|---|---|
| `Handler::DiffractIncline` | CPU scalar non-absorbing edge integral |
| `Handler::DiffractInclineFast` | CPU optimized edge integral with precomputed edge data |
| `Handler::DiffractInclineAbs` | Absorbing variant |
| `GpuDiffraction.cu` kernels | CUDA implementation of the same aperture/edge calculation |

Conceptually the scalar diffraction factor is the aperture phase integral

```text
F(q) = integral_A exp(i k q.r) dA,
```

where `A` is the projected beam polygon, `k = 2*pi/lambda`, and `q` is the difference between the outgoing beam direction and the requested far-field direction. The implementation evaluates the polygon integral through its edges to avoid sampling the aperture area directly. For an edge from vertex `a` to `b`, the contribution is a stable complex exponential difference divided by the corresponding phase denominator; special small-denominator branches use sinc-like limits.

The GPU and CPU branches should use the same geometry, phase convention, and theta grid. Differences normally come from precision (`float` vs `double`), fast math, FFT interpolation, atomics/reduction order, and special pole or endpoint handling.

### Coherent and incoherent Mueller accumulation

Default PO mode is coherent per orientation:

```text
J_total(theta, phi) = sum_beams J_beam(theta, phi)
M(theta, phi)       = Mueller(J_total(theta, phi))
```

With `--incoh`, the conversion happens before beam summation:

```text
M(theta, phi) = sum_beams Mueller(J_beam(theta, phi))
```

The Jones-to-Mueller conversion is implemented in `src/math/Mueller.cpp`. For Jones elements `S1..S4`, the code computes the standard 4x4 real Mueller matrix from bilinear products such as `|S1|^2`, `|S2|^2`, `Re(S1 conj(S2))`, and `Im(S1 conj(S2))`.

The coherent sum is limited to ray paths of one particle at one fixed orientation. For an ensemble of randomly and independently oriented particles, the code first forms `Mueller(J_total)` for each orientation and then averages Mueller matrices with the orientation weights. Jones matrices from distinct orientations cannot be summed without particle positions and relative incident phases. The deprecated `--coherent-orientations` option is therefore disabled; its previous implementation was effectively the ordinary incoherent orientation average.

Consequently, this single-particle model does not reproduce coherent backscattering of a particulate medium caused by interference of reciprocal multiple-scattering paths between distinct particles. That problem requires a multi-particle solver with particle positions and inter-particle field propagation.

`--beam-cutoff*` and `--trace-cutoff*` remove amplitudes before coherent summation and can bias cross terms even when the removed intensities are small. Validate backscatter with a repeat using `--cutoff-profile safe` or `--cutoff-profile off`.

### Orientation averaging

For randomly oriented particles, each orientation produces a Mueller matrix on the scattering grid. `HandlerPOTotal` and `TracerPOTotal` average orientations with the beta/gamma weights selected by the orientation mode. At beta poles, `--pole` uses one gamma value because all gamma rotations are equivalent at the exact pole. In current code, `--pole` keeps beta endpoints instead of midpoint beta sampling, so the beta pole is actually present in the grid.

### Scattering-angle weights

Uniform theta grids output `Nth + 1` rows. `2pi*dcos` is the solid-angle ring weight assigned to the theta row:

```text
dOmega(theta_j) = 2*pi * (cos(theta_left) - cos(theta_right))
```

For endpoint rows the interval is half-width and clipped to the requested range. For a full `0..180` grid the sum of `2pi*dcos` over all rows is `4*pi`.

## Particles and sizes

| Flag | Arguments | Description |
|---|---:|---|
| `-p` | `TYPE L D [extra]` | Built-in particle. Exactly one of `-p` or `--pf` is required. |
| `--pf` | `FILE` | Load particle from file. |
| `--rs` | `SIZE` | Resize file particle to `Dmax = SIZE`. Mutually exclusive with `--k_eq`. |
| `--k_eq` | `X` | Resize so `k_eq = 2*pi*r_eq/lambda`. Requires wavelength. |
| `--ri` | `Re Im` | Complex refractive index. Absorption is enabled automatically when `Im != 0`, or explicitly with `--abs`. |
| `-w` | `LAMBDA` | Wavelength in micrometers. |
| `-n` | `N` | Maximum internal reflection/refraction depth. |

File coordinates and all internal tracing geometry are stored in double
precision. The loader translates coordinates relative to the
order-independent bounding-box centre and removes only scale-relative binary
serialization noise. A global translation changes only a common optical phase;
this canonical form preserves the Mueller matrix, small edges, and facet order
for files with large absolute coordinates.

Built-in particle types:

| Type | Shape | Parameters |
|---:|---|---|
| 1 | Hexagonal column/plate | `L D` |
| 2 | Bullet | `L D` |
| 3 | Bullet rosette | `L D [cap]` |
| 4 | Droxtal | `L D extra`, uses the extra shape parameter |
| 10 | Concave hexagonal | `L D concavity` |
| 12 | Hexagonal aggregate | `L D count` |
| 999 | Built-in aggregate | `extra` |

Equivalent-size scaling:

```text
r_eq = (3 V / (4*pi))^(1/3)
k_eq = 2*pi*r_eq/lambda
scale = (k_eq_target * lambda / (2*pi)) / r_eq_original
```

## Orientation grids

| Mode | Arguments | Best use | Notes |
|---|---:|---|---|
| `--oldauto` | `DIV` | Production regular grid | Grid step follows diffraction-limited angular scale; common values are `2`, `4`, `8`. |
| `--random` | `Nb Ng` | Manual beta/gamma regular grid | Uses symmetry-reduced beta/gamma domain unless overridden. |
| `--fixed` | `BETA GAMMA` | Single-orientation debugging | Angles are degrees. |
| `--orientfile` | `FILE` | Reproducible custom orientations | One `beta_deg gamma_deg` pair in degrees per line. |
| `--sobol` | `N` | Quasi-random orientation average | Good convergence for broad scans. |
| `--sobol_seed` | `N S` | Seeded Sobol/Owen sequence | Useful for repeatable convergence checks. |
| `--sobol_ring` | `Nb Ng` | Sobol beta with uniform gamma rings | Hybrid orientation grid. |
| `--so3_quat` | `N` | Isotropic SO(3) average | Hammersley on the symmetry quotient; laboratory alpha is integrated by scattering phi. |
| `--so3_full_quat` | `N` | Independent SO(3) audit | Direct full-group Shoemake quaternions; deliberately ignores particle symmetry. |
| `--hammersley` | `N` | Low-discrepancy orientations | Debug/experimental. |
| `--lattice` | `N` | Rank-1 lattice orientations | Debug/experimental. |
| `--lattice_z` | `N Z` | Rank-1 lattice with explicit generator | Debug/experimental. |
| `--euler_quad` | `Nb Ng` | Gauss in cos(beta), periodic gamma | High-order quadrature. |
| `--euler_adapt` | `Nb NgMax` | Adaptive gamma count per beta ring | Reduces work near beta poles. |
| `--adaptive` | `EPS` | Adaptive Sobol count | Doubles orientations until relative convergence target. |
| `--auto` | `EPS` | Auto theta, phi, orientation count | Convenience mode. |
| `--autofull` | `EPS` | Auto `n`, theta, phi, orientations | Slower, more complete search. |
| `--oldautofull` | `EPS` | Autofull search with oldauto final grid | Useful when final regular grid is required. |

Orientation modifiers:

| Flag | Arguments | Description |
|---|---:|---|
| `--ring_points` | `N` | Points per diffraction ring for oldauto estimates; default `3`. |
| `--mirror_gamma` | none | For a verified reflection-symmetric particle, sample half the gamma range and restore the omitted Mueller contribution with Stokes parity. Supported by `--so3_quat`. |
| `--so3_mirror_audit` | none | With an even `--so3_quat N`, explicitly trace nested reflected pairs for a direct comparison with `N/2 --mirror_gamma`. Validation only. |
| `--sym` | `Sb Sg` | Override symmetry domain: beta range is `pi/Sb`, gamma range is `2*pi/Sg`. |
| `--b` | `B1 B2` | Manual beta range in degrees for `--random`. |
| `--g` | `G1 G2` | Manual gamma range in degrees for `--random`. |
| `--maxorient` | `N` | Maximum orientations for adaptive modes. |
| `--chunk` | `N` | Orientation/gamma chunk size. |
| `--coh_orient` | none | Disabled legacy mode; cross-orientation coherence is undefined without relative phases. |
| `--pole` | none | Use one gamma value at beta poles. Best with endpoint grids such as `--oldauto`. |
| `--owen_avg` | `K` | Average `K` nested Owen final seeds in `--autofull`; default can be controlled by environment. |
| `--owen_seeds` | `S...` | Explicit Owen seeds for `--autofull` final averaging. |

### Symmetry-reduced SO(3) quadrature

For an isotropic ensemble, the Haar measure in the program's Euler convention
is proportional to

```text
dR = sin(beta) d(beta) d(gamma) d(alpha).
```

`--so3-quaternion N` uses the particle's fundamental domain
`0 <= beta <= beta_sym`, `0 <= gamma < gamma_sym`. Point `i=0,...,N-1` is

```text
u_i = (i + 1/2) / N
v_i = radical_inverse_base_2(i)
beta_i  = acos(1 - (1 - cos(beta_sym)) u_i)
gamma_i = gamma_sym v_i
w_i = 1/N,  sum_i w_i = 1.
```

For a paired mirror validation, add `--so3-mirror-audit` to an even full count
and let `K=N/2`. The audit uses the same formulas with `i=0,...,K-1`, `u_i=(i+
1/2)/K`, and `gamma_i=(gamma_sym/2)v_i`, then explicitly traces both
`(beta_i,gamma_i)` and `(beta_i,gamma_sym-gamma_i)` with weight `1/N`. This
layout nests exactly with `--so3-quaternion K --mirror-gamma`. The modifier is
not a production quadrature rule and deliberately leaves the normal
full-domain Hammersley sequence unchanged.

Thus beta is uniform in `cos(beta)`, not in beta itself. `--sym Sb Sg`
overrides the fundamental domain with `beta_sym=pi/Sb` and
`gamma_sym=2*pi/Sg`; use this only for an actual symmetry of the particle.

With `--mirror-gamma`, the sampled range becomes
`0 <= gamma < gamma_sym/2`. The second half is restored before the alpha
quadrature as

```text
M_full(theta,phi) = [M_half(theta,phi)
                    + P M_half(theta,-phi) P] / 2,
P = diag(1,1,-1,-1).
```

This is valid only when the optical particle, including facet visibility and
material assignment, has the asserted reflection plane and the selected PO
tracer is numerically invariant under that reflection. The program cannot
prove this property from a particle file. Validate every new particle class
and parameter regime against one full-gamma run. `N` remains the number of
actually traced beta/gamma orientations. Use an even full count with the audit
modifier and exactly half that count with mirror; for example, compare
`--so3-quaternion 4096 --so3-mirror-audit` with `--so3-quaternion 2048
--mirror-gamma`. These runs share the same base points, so their difference
measures numerical mirror consistency instead of two-sample Hammersley noise.
Also validate the normal full-domain sequence against an independent dense
result. Check backscatter separately because a small global L2 error can hide
a local pole discrepancy.

The remaining laboratory angle alpha is not retraced. A rotation around the
incident beam is equivalent to changing scattering azimuth phi, so each theta
row is accumulated as

```text
M_avg(theta) = (1/Nphi) sum_phi M(theta,phi) L(-phi).
```

Here `L` is the Stokes reference-plane rotation; its Q/U block contains
`cos(2 phi)` and `sin(2 phi)`. Multiplication is on the right because this
rotation changes the incident polarization basis. At theta=0 and theta=pi the
azimuthal basis is singular: the code averages without `L` and then applies the
exact forward/backward pole form. This is the same implementation used by the
regular variable-phi path.

The mode requires `Nphi >= 2`. Recommended production form:

```bash
mbs_po --method po ... --so3-quaternion 4096 \
    --scattering-diffraction-sampling 2 --latitude-phi-grid
```

Use `--so3-full-quaternion N` only to audit the reduced result. Its `N` counts
full three-dimensional rotations, whereas the reduced mode's `N` counts traced
beta/gamma particle orientations; its phi grid performs the alpha quadrature.

## Scattering grids

| Flag | Arguments | Description |
|---|---:|---|
| `--grid` | `T1 T2 Nphi Nth` | Uniform theta grid from `T1` to `T2` degrees. Output has `Nth + 1` theta rows. |
| `--grid` | `R Nphi Nth` | Backscatter cone grid of radius `R` degrees. |
| `--tgrid` | `FILE` | Non-uniform theta grid in degrees, one row per theta value. |
| `--auto_tgrid` | `EPS` | Adaptive theta grid by bisection. |
| `--auto_phi` | none | Choose `Nphi` automatically from size parameter. |
| `--nphi` | `N` | Override phi count. Highest priority for phi. |
| `--diffraction-limit-grid` | `FACTOR` | Derive uniform theta and equatorial-phi steps from the diffraction estimate. |
| `--latitude-phi-grid` | none | Reduce the phi count in proportion to `sin(theta)` away from the equator. |
| `--filter` | `DEG` | Restrict output to a backscattering cone. Legacy/debug use. |
| `--point` | none | Backscatter point mode. Legacy; avoid for optimized regular-grid production. |

Priority:

| Quantity | Priority |
|---|---|
| Theta grid | `--tgrid` > `--grid` > `--auto_tgrid` > `--auto` > default |
| Phi grid | `--nphi` > `--grid` > `--auto_phi` > `--auto` > default |

For full angular integrals use `--grid 0 180 Nphi Nth`. Endpoint rows are physical half-cells; they must not have zero weight.

### Diffraction-limit and latitude-phi grids

For a PO `--diffraction-grid DIV` run, `--diffraction-limit-grid FACTOR`
derives the scattering grid from the longest edge of any current particle
facet, `lmax`, and wavelength `lambda`. This is not the particle diameter or
the maximum distance between arbitrary vertices:

```text
xi = FACTOR * 0.69 * lambda / lmax
Ntheta = ceil(pi / xi)
Nphi,eq = round_up_6(max(12, ceil(2*pi / xi)))
```

`xi` is in radians and the output contains `Ntheta + 1` theta rows from zero
to pi. `FACTOR=1` uses the diffraction-width estimate literally;
`FACTOR=0.5` approximately halves both angular steps, while a factor above
one makes the grid coarser. This is an a-priori resolution estimate, not an
a-posteriori accuracy guarantee. Validate a production factor against a finer
one for every required Mueller element.

With `--latitude-phi-grid`, each theta row uses

```text
Nphi(theta) = min(Nphi,eq,
                  round_up_6(max(12, ceil(2*pi*sin(theta) / xi))))
```

The equatorial spacing is unchanged, while redundant azimuth samples are
removed near the poles. Counts are rounded to a multiple of six, and a pole
shares its neighboring work group for stable basis handling. Results are
still accumulated into the full rectangular `Nphi,eq` output layout, so the
file format does not change.

The variable latitude grid supports `--method po --diffraction-grid DIV` and
the symmetry-reduced `--so3-quaternion N` path, one particle size per process.
For SO(3), prefer `--scattering-diffraction-sampling Q`; the unified
`--diffraction-sampling Q` is itself a primary regular orientation rule.
`--diffraction-limit-grid` conflicts with explicit
`--nphi`; the latitude grid is not supported by a shared serial multi-size
trace. Run each size in its own process so `lmax` and the grid are recomputed.

```bash
CUDA_VISIBLE_DEVICES=0,1,2,3 \
MBS_GPU_GROUPS=1 MBS_GPU_MULTI_MAX=4 \
gpu/bin/mbs_po_gpu_double --po -p 1 100 70 --ri 1.3116 0 -w 0.532 \
    -n 8 --diffraction-grid 2 --diffraction-limit-grid 0.5 \
    --latitude-phi-grid --threads 32 --close -o results/column
```

## GPU and multi-GPU

### CUDA backend

| Flag | Description |
|---|---|
| `--gpu` | Use CUDA diffraction backend. Optional in `gpu/` binaries because they default to GPU. |
| `--cpu` | Force CPU backend from a GPU-capable binary. |
| `--fft` | Use cuFFT angular phi interpolation backend. Requires GPU. |
| `--threads N` | OpenMP host worker threads. Also controls tracing-side parallelism before GPU diffraction. |

For fixed-orientation CPU PO, the direction-grid calculation is parallelized
over azimuth rows when the grid is large enough to amortize OpenMP startup.
Each worker owns disjoint Jones/Mueller cells, while beam contributions to a
cell are still accumulated in their original order. Results are therefore
independent of `--threads` up to ordinary floating-point output precision.

GPU runtime errors are fatal in the split GPU build. For debugging only:

```bash
MBS_GPU_ALLOW_FALLBACK=1 gpu/bin/mbs_po_gpu_float_fast ...
```

### Why the GPU version is faster

The expensive PO part is the repeated evaluation of the same beam aperture integral for many detector directions and many orientations. The CPU traces rays and prepares compact beam records; the GPU then evaluates diffraction and, in the fast paths, Mueller accumulation on a large regular grid.

| Work item | CPU path | GPU path | Parallel granularity |
|---|---|---|---|
| Ray tracing through facets | OpenMP over orientations/gamma blocks | Mostly still CPU-side; `--gpu_trace` is only a candidate prefilter | Orientation/chunk level |
| Beam packing | CPU host code | CPU prepares `GpuBeam`/packed beam buffers and copies to device | Beam records |
| Edge diffraction | Scalar/vector CPU loops | CUDA kernels in `GpuDiffraction.cu` | Beam x theta x phi x orientation |
| Jones coherent sum | CPU arrays of complex 2x2 matrices | Atomic or no-atomic device kernels | Grid cell and orientation |
| Jones to Mueller | CPU after Jones sum | Separate `mueller_batch_kernel` or fused in diffraction kernel | Grid cell |
| Multi-`k_eq` scans | Repeat diffraction per size | Optional fused multi-`k_eq` CUDA pass | Size x beam x grid |
| FFT phi interpolation | CPU direct phi grid unless disabled | cuFFT low-phi pass plus zero-padding reconstruction | Theta/orientation batches |

The main speedup comes from exposing a very large product of independent operations:

```text
Nwork ~= Norient * Nbeams_per_orientation * Ntheta * Nphi
```

Each beam/direction contribution is mostly independent until the coherent Jones sum. That makes the diffraction stage a good CUDA workload: many threads do the same edge-integral arithmetic over different `(orientation, beam, theta, phi)` combinations.

### Atomic and no-atomic GPU accumulation

There are two families of CUDA accumulation paths in `src/cuda/GpuDiffraction.cu`.

| Mode | Main kernels | How accumulation works | When useful |
|---|---|---|---|
| Atomic Jones accumulation | `diffraction_kernel` | Each thread computes one beam/grid contribution and uses `atomicAdd` into the complex Jones buffer. Afterwards `mueller_batch_kernel` converts the summed Jones matrix to Mueller. | General fallback, supports more combinations such as no-shadow output. |
| No-atomic grid kernels | `diffraction_grid_kernel` | Work is reorganized so a thread owns a grid/output slot and loops over beams, avoiding multiple writers to the same Jones cell. | Full-only output when memory/layout allows it. |
| Fused Mueller kernels | `diffraction_grid_mueller_kernel`, `diffraction_grid_mueller_full_kernel`, `*_full8_kernel`, `*_mixed8_kernel` | The kernel computes coherent Jones locally for an output slot and immediately converts/adds Mueller, avoiding large intermediate Jones buffers. | Fast production path for full-only GPU runs. |
| Staged orientation Mueller | `diffraction_grid_mueller_orient_kernel` + `reduce_mueller_orient_kernel` | First writes per-orientation Mueller, then reduces orientations. | Reduces contention and gives a deterministic reduction shape for large orientation batches. |
| Multi-`k_eq` fused | `diffraction_grid_mueller_multik_kernel` | Computes several nearby equivalent sizes from one packed beam batch. | Shared-batch `--multikeq` scans. |

In the atomic path the atomics are on the real and imaginary components of the 2x2 Jones matrix:

```text
J00.re, J00.im, J01.re, J01.im,
J10.re, J10.im, J11.re, J11.im
```

For `_noshadow` output the same set of atomics is also applied to the no-shadow Jones buffer. This is correct for coherent PO because beams must be summed as Jones amplitudes before conversion to Mueller. The cost is contention when many beams hit the same detector cell.

The no-atomic/fused paths avoid that contention by changing ownership: one CUDA thread or thread group owns an output grid element and accumulates the beams for that element in registers/local variables. This is often faster for full Mueller output because it reduces global-memory atomics and may skip the large intermediate Jones buffer. The tradeoff is that these paths are more specialized and may be disabled for diagnostic modes that require extra outputs.

Useful controls:

| Environment variable | Effect |
|---|---|
| `MBS_GPU_NO_ATOMICS=1` | Force no-atomic path when supported. |
| `MBS_GPU_NO_ATOMICS=0` | Force atomic path for comparison/debugging. |
| `MBS_GPU_FUSED_MUELLER=1` | Prefer fused diffraction-to-Mueller kernels when no-atomic mode is active. |
| `MBS_GPU_STAGE_MUELLER=1` | Use staged per-orientation Mueller plus reduction where supported. |
| `MBS_GPU_NO_VERTEX_CACHE=1` | Disable cached/packed vertex path for debugging. |
| `MBS_GPU_TIMING=1` | Print GPU timing breakdown for count/pack/copy/kernels/d2h/add. |
| `MBS_GPU_BLOCK=N` | Override CUDA block size for kernel tuning. |
| `MBS_ORIENTATION_TIMING=1` | Print summed CPU time for rotation, tracing, and prepared-beam construction. |

The diagnostic target below compares FP64 direct quaternion-vector rotation
with a precomputed 3x3 matrix on the GPU:

```bash
make USE_CUDA=1 gpu_quaternion_probe
CUDA_VISIBLE_DEVICES=0 ./bin/gpu_quaternion_rotation_probe 65536 64 50
```

This probe does not calculate scattering and does not replace profiling of a
complete run. Direct quaternions provide uniform full-SO(3) samples in the
audit mode; by themselves they do not accelerate CPU tracing or CUDA
diffraction. The production `--so3-quaternion` speedup comes from symmetry and
from replacing alpha retracing by the scattering-phi integral.

Do not confuse this with full GPU ray tracing. The production GPU backend is primarily a diffraction and Mueller-accumulation accelerator. `--gpu_trace` only prefilters nonconvex tracing candidates; exact intersections remain CPU-side.

### Multi-size and multi-GPU scans

`MBS_GPU_GROUPS=1` enables the dedicated exact scheduler for a variable
latitude-phi grid. Independent theta-row work groups are assigned dynamically
to visible GPUs. Orientation tracing and prepared beams remain shared in host
memory, while each GPU owns one CUDA workspace. Packed beams, weights, and
offsets are uploaded once per orientation chunk per device and reused for the
subsequent theta groups handled by that worker.

Select devices with `CUDA_VISIBLE_DEVICES`. Limit workers with
`MBS_GPU_MULTI=N` or `MBS_GPU_MULTI_MAX=N`; `MBS_GPU_MULTI=0` leaves one GPU.
This scheduler requires direct CUDA diffraction. FFT interpolation and
theta-zone beam filtering retain the existing path. Speedup is not necessarily
linear because CPU tracing, packing, PCIe transfers, final reduction, and
theta-group load imbalance remain in the wall time.

| Flag | Arguments | Description |
|---|---:|---|
| `--multigrid` | `Dmin Dmax N` | Log-spaced scan in maximum particle dimension. |
| `--multikeq` | `Kmin Kmax N` | Log-spaced scan in equivalent size parameter. |
| `--multikeq_list` | `FILE` | Exact list of `k_eq` values, one per line. |
| `--multigrid_parallel` | `N` | Run scan points as child processes; `0` means auto from the GPU list or CPU core count. |
| `--multigrid_threads` | `N` | OpenMP threads per child process. |
| `--gpu_devices` | `LIST` | Comma-separated CUDA device list, for example `0,1,2,3,4`. |
| `--multikeq_shared_batches` | none | Batch nearby `k_eq` values per GPU child and reuse tracing from the largest value in the batch. |
| `--multikeq_batch_ratio` | `R` | Maximum `kmax/kmin` inside one shared batch; default `1.05`. |

The parallel scheduler is process-based. It removes the scan flag from the
parent command, creates one output subdirectory per size or batch under `-o`,
then launches children with either `--rs SIZE`, `--k_eq K`, or a generated
`--multikeq_list FILE`. Each child writes stdout/stderr to
`<output>/<label>.run.log`.

With `--gpu`, each child is assigned one CUDA device by setting
`CUDA_VISIBLE_DEVICES` and `MBS_GPU_SLOT`. If `--gpu_devices` is omitted, the
code uses the CUDA devices visible to the process. `--multigrid_parallel 0`
means "use as many jobs as useful": for GPU runs this normally becomes the
number of listed/visible devices; for CPU runs it is limited by available host
threads. Use explicit `--multigrid_parallel N` when the machine is shared.

Exact multi-GPU scan, one process per active GPU slot and one independently
traced size per child:

```bash
gpu/bin/mbs_po_gpu_float_fast --po --pf particle.dat \
    --multikeq_list keq.txt --ri 1.6 0.002 -w 1.064 -n 12 \
    --oldauto 2 --grid 0 180 600 180 --fft --close \
    --multigrid_parallel 0 --multigrid_threads 16 --gpu_devices 0,1,2,3,4 \
    -o scan_exact
```

Shared-batch `k_eq` scan:

```bash
gpu/bin/mbs_po_gpu_float_fast --po --pf particle.dat \
    --multikeq_list keq.txt --ri 1.6 0.002 -w 1.064 -n 12 \
    --oldauto 2 --grid 0 180 600 180 --fft --close \
    --multigrid_parallel 0 --multigrid_threads 16 --gpu_devices 0,1,2,3,4 \
    --multikeq_shared_batches --multikeq_batch_ratio 1.05 \
    -o scan_shared
```

`--multikeq_shared_batches` is only for `k_eq` scans. Values close enough by
`--multikeq_batch_ratio` are grouped into one child. The child traces the
largest `k_eq` in the batch and reuses that prepared beam set for smaller
values, so it is faster but not identical to independent oldauto grids. Use a
tighter ratio for validation, and use exact mode when every size must have its
own oldauto reference grid.

Experimental fused multi-`k_eq` GPU diffraction:

```bash
MBS_GPU_MULTI_K_FULL=1 gpu/bin/mbs_po_gpu_float_fast ...
```

For large direct `theta x phi` grids, host RAM can dominate even when GPU memory
is mostly free. The oldauto log prints a line like:

```text
Oldauto/random memory: gamma chunk=... (GPU host-RAM guard, MemAvailable=... MB, VmRSS=... MB, grid-transient=... MB, budget=... MB)
```

If `grid-transient` is comparable to the budget, lowering `--chunk` alone cannot
fix memory use; reduce `Nphi/Nth`, use `--fft`, increase the host RAM budget, or
split the calculation into smaller independent grids.

For oldauto/random checkpoints, the no-shadow Mueller grid is stored only when
no-shadow output is requested. Normal full-output runs therefore avoid keeping
and checkpointing a second full `theta x phi x 4 x 4` matrix, which matters for
large direct grids.

Common scheduler controls:

| Variable | Meaning |
|---|---|
| `MBS_PARALLEL_MEM_FRACTION` | Fraction of current `MemAvailable` divided among GPU child processes; default `0.70`. |
| `MBS_HOST_MEM_BUDGET_MB` | Hard host-memory budget seen by a child. The scheduler sets it automatically unless already set. |
| `MBS_HOST_MEM_FRACTION` | Fraction of current available host RAM allowed for prepared-orientation chunks. |
| `MBS_HOST_MEM_RESERVE_MB` | Host RAM reserve kept free by oldauto/random chunking; default `4096`. |
| `MBS_OLDAUTO_GRID_MEM_SAFETY` | Safety multiplier for direct-grid transient arrays; default `1.25`. |
| `MBS_OLDAUTO_BYTES_PER_GAMMA_MB` | Per-gamma host-memory estimate for prepared beams; default is conservative. |
| `MBS_OLDAUTO_GAMMA_CHUNK` | Default oldauto gamma chunk before the host-memory guard clamps it. |
| `MBS_PARALLEL_SHARED_KEQ=1` | Environment equivalent of `--multikeq_shared_batches`. |
| `MBS_PARALLEL_KEQ_BATCH_RATIO=R` | Environment default for `--multikeq_batch_ratio`. |

## All flags

### Method, particle, and physical parameters

| Flag | Arguments | Description |
|---|---:|---|
| `--po` | none | Physical Optics mode. |
| `--go` | none | Geometrical Optics mode. |
| `-p` | `TYPE L D [extra]` | Built-in particle. |
| `--pf` | `FILE` | Particle from file. |
| `--rs` | `SIZE` | Resize file particle by Dmax. |
| `--k_eq` | `X` | Resize by equivalent size parameter. |
| `--ri` | `Re Im` | Refractive index. |
| `-w` | `LAMBDA` | Wavelength in micrometers. |
| `-n` | `N` | Maximum internal reflections/refractions. |
| `--abs` | none | Enable absorption accounting. Also enabled automatically for nonzero `Im(ri)`. |
| `--abs_points` | `N` or `all` | Absorption sampling: center point or all polygon vertices. |

### Orientation flags

| Flag | Arguments | Description |
|---|---:|---|
| `--fixed` | `BETA GAMMA` | Single orientation in degrees. |
| `--random` | `Nb Ng` | Regular beta/gamma grid. |
| `--sobol` | `N` | Sobol orientations. |
| `--so3_quat` | `N` | Symmetry-reduced SO(3); alpha is integrated by scattering phi. |
| `--so3_mirror_audit` | none | Explicit nested reflected pairs for validating `--mirror_gamma`; even `N` only. |
| `--so3_full_quat` | `N` | Direct full-SO(3) quaternion audit. |
| `--sobol_seed` | `N S` | Sobol with explicit Owen seed. |
| `--sobol_ring` | `Nb Ng` | Sobol beta with gamma rings. |
| `--hammersley` | `N` | Hammersley orientation set. |
| `--lattice` | `N` | Rank-1 lattice orientation set. |
| `--lattice_z` | `N Z` | Rank-1 lattice with explicit generator. |
| `--euler_quad` | `Nb Ng` | Gauss beta x periodic gamma quadrature. |
| `--euler_adapt` | `Nb NgMax` | Gauss beta with adaptive gamma count. |
| `--montecarlo` | `N` | Pseudo-random Monte Carlo orientations. |
| `--adaptive` | `EPS` | Adaptive orientation convergence. |
| `--auto` | `EPS` | Auto theta, phi, and orientation count. |
| `--autofull` | `EPS` | Auto `n`, theta, phi, and orientations. |
| `--oldautofull` | `EPS` | Autofull search with oldauto final grid. |
| `--owen_avg` | `K` | Average final Owen seeds. |
| `--owen_seeds` | `S...` | Explicit final Owen seeds. |
| `--oldauto` | `DIV` | Physics-based regular grid. |
| `--ring_points` | `N` | Points per diffraction ring estimate. |
| `--mirror_gamma` | none | Use mirrored half gamma domain. |
| `--orientfile` | `FILE` | Load beta/gamma orientations in degrees from file. |
| `--b` | `B1 B2` | Beta range for `--random`. |
| `--g` | `G1 G2` | Gamma range for `--random`. |
| `--maxorient` | `N` | Maximum adaptive orientation count. |
| `--chunk` | `N` | Chunk size for orientation/gamma processing. |
| `--coh_orient` | none | Disabled legacy mode; do not use. |
| `--pole` | none | Exact beta-pole gamma shortcut. |
| `--sym` | `Sb Sg` | Override particle symmetry. |

### Scattering grid and acceleration flags

| Flag | Arguments | Description |
|---|---:|---|
| `--grid` | `T1 T2 Nphi Nth` | Uniform theta range. |
| `--grid` | `R Nphi Nth` | Backscatter cone. |
| `--tgrid` | `FILE` | Non-uniform theta grid. |
| `--auto_tgrid` | `EPS` | Adaptive theta grid. |
| `--auto_phi` | none | Automatic phi count. |
| `--nphi` | `N` | Override phi count. |
| `--diffraction-limit-grid` | `FACTOR` | Set theta and equatorial-phi steps to `FACTOR*0.69*lambda/lmax`. |
| `--latitude-phi-grid` | none | Use a row-dependent phi count proportional to `sin(theta)`. |
| `--threads` | `N` | OpenMP worker threads. |
| `--gpu` | none | CUDA backend. |
| `--cpu` | none | Force CPU backend in GPU build. |
| `--fft` | none | cuFFT phi interpolation. |
| `--beam_cutoff` | `EPS` | Set relative beam `J` and area cutoffs. |
| `--beam_cutoff_j` | `EPS` | Skip beams by relative `|J|^2`. |
| `--beam_cutoff_area` | `EPS` | Skip beams by relative beam area. |
| `--beam_cutoff_importance` | `EPS` | Skip beams by relative `|J|^2*area`. |
| `--trace_cutoff` | `EPS` | Set trace pruning `J` and area cutoffs. |
| `--trace_cutoff_j` | `EPS` | Prune trace tree by relative `|J|^2`. |
| `--trace_cutoff_area` | `EPS` | Prune trace tree by relative area. |
| `--trace_cutoff_importance` | `EPS` | Prune trace tree by `|J|^2*area`. |
| `--trace_max_beams` | `N` | Abort the calculation after `N` beam nodes in one orientation; `0` disables. |
| `--gpu_trace` | none | Experimental CUDA prefilter for nonconvex tracing candidates. |
| `--trace_prefilter` | none | Enable CPU projected-AABB prefilter. |
| `--no_trace_prefilter` | none | Disable CPU projected-AABB prefilter. |
| `--trace_prefilter_margin` | `M` | AABB prefilter margin. |
| `--trace_prefilter_stats` | none | Print prefilter counters. |
| `-r` | `RATIO` | Beam area restriction ratio; default `100`. |

### Multi-size and multi-GPU flags

| Flag | Arguments | Description |
|---|---:|---|
| `--multigrid` | `Dmin Dmax N` | Log-spaced scan in maximum particle dimension. |
| `--multikeq` | `Kmin Kmax N` | Log-spaced scan in equivalent size parameter. |
| `--multikeq_list` | `FILE` | Exact `k_eq` values, one per line. |
| `--multigrid_parallel` | `N` | Launch scan points as child processes; `0` chooses an automatic job count. |
| `--multigrid_threads` | `N` | OpenMP threads per child process. |
| `--gpu_devices` | `LIST` | CUDA devices assigned round-robin to child processes. |
| `--multikeq_shared_batches` | none | Group nearby `k_eq` values and reuse tracing from the largest value in each batch. |
| `--multikeq_batch_ratio` | `R` | Maximum `kmax/kmin` within one shared batch; default `1.05`. |

### Output, diagnostics, and legacy flags

| Flag | Arguments | Description |
|---|---:|---|
| `-o` | `NAME` | Output path/name. |
| `--close` | none | Exit after computation. |
| `--log` | `SEC` | Progress logging interval. |
| `--checkpoint` | none | Save/resume long `--orientfile`, `--oldauto`, and `--random` runs. |
| `--save_betas` | none | Save per-beta Mueller files. |
| `--full_only` | none | Write only full Mueller output. This is the default. |
| `--noshadow_output` | none | Also write `_noshadow` Mueller output. |
| `--shadow` | none | Legacy flag; currently no effect. |
| `--shadow_off` | none | Disable shadow beam. Diagnostic use. |
| `--incoh` | none | Incoherent per-beam Mueller sum. |
| `--jones` | none | Write Jones matrices where supported. |
| `--karczewski` | none | Use Karczewski polarization matrix. |
| `--legacy_sign` | none | Use old Fresnel sign convention. |
| `--ot_phase_avg` | none | Average optical-theorem extinction over one far-reference phase period. |
| `--ot_phase_shift` | `F` | Diagnostic optical-theorem phase shift in wavelengths. |
| `--ot_ping` | `D` | Legacy far-screen optical-theorem phase rotation. |
| `--filter` | `DEG` | Restrict output to backscatter cone. |
| `--point` | none | Legacy backscatter point mode. |
| `--tr` | `FILE` | Load trajectory file. |
| `--all` | none | Calculate all loaded trajectories. |
| `--gr` | none | Output trajectory groups. |
| `--forced_convex` | none | Force convex processing. |
| `--forced_nonconvex` | none | Force nonconvex processing. |
| `--help`, `-h` | none | Print help. |
| `--help-debug` | none | Print full debug/experimental help. |

## Environment variables

Every active `MBS_*` variable is written to the run summary. Variables that
alter sampling, geometry, cutoffs, or numerical results are rejected by
default. `--allow-experimental-environment` permits them after explicit
acknowledgement and keeps their names and values in the result log. Production
runs should normally unset them. Prepared-orientation modes begin with a
16-orientation pilot, measure actual nested beam/path storage, and adapt the
next chunk to current available RAM and the limits below.

Production-use variables:

| Variable | Meaning |
|---|---|
| `CUDA_VISIBLE_DEVICES` | Standard CUDA device visibility control. |
| `OMP_NUM_THREADS` | OpenMP default thread count when `--threads` is not used. |
| `MBS_GPU_ALLOW_FALLBACK=1` | Allow old GPU-to-CPU fallback after CUDA failure. Debug only. |
| `MBS_GPU_MEM_FRACTION` | Fraction of free GPU memory allowed for internal buffers. |
| `MBS_GPU_NO_ATOMICS` | Select no-atomic (`1`) or atomic (`0`) GPU accumulation path when possible. |
| `MBS_GPU_FUSED_MUELLER` | Select fused diffraction-to-Mueller kernels when no-atomic mode is active. |
| `MBS_GPU_STAGE_MUELLER` | Enable staged per-orientation Mueller reduction when supported. |
| `MBS_GPU_NO_VERTEX_CACHE` | Disable cached/packed vertex path for debugging. |
| `MBS_GPU_TIMING` | Print CUDA timing breakdowns. |
| `MBS_GPU_BLOCK` | Override CUDA kernel block size. |
| `MBS_ORIENTATION_TIMING` | Print CPU-time breakdown for rotation, tracing, and prepared-beam construction. |
| `MBS_HOST_MEM_FRACTION` | Fraction of current available host RAM for measured prepared-orientation chunking. |
| `MBS_HOST_MEM_RESERVE_MB` | Host memory reserve in MB. |
| `MBS_HOST_MEM_BUDGET_MB` | Hard host memory budget, also set for parallel children. |
| `MBS_OLDAUTO_GRID_MEM_SAFETY` | Safety multiplier for oldauto host-RAM estimate of large `theta x phi` grids; default `1.25`. |
| `MBS_OLDAUTO_BYTES_PER_GAMMA_MB` | Per-gamma prepared-beam memory estimate used by oldauto/random host-RAM guard. |
| `MBS_OLDAUTO_GAMMA_CHUNK` | Default oldauto gamma chunk before automatic memory clamping. |
| `MBS_OLDAUTO_BETA_MIDPOINT=0/1` | Use midpoint beta rings in oldauto/random quadrature. |
| `MBS_OLDAUTO_GAMMA_STAGGER=1` | Stagger gamma samples between beta rings to reduce aliasing in narrow events. |
| `MBS_FFT_PHI_FACTOR` | Override FFT reduced phi factor. |
| `MBS_FFT_THETA_FACTOR` | Override FFT theta batching factor. |
| `MBS_FFT_CHECK=1` | Enable FFT diagnostic checks. |
| `MBS_FFT_ADAPTIVE_PHI=1` | Enable adaptive reduced-phi behavior in FFT backend. |
| `MBS_GPU_MULTI=0` | Disable automatic multi-orientation GPU batching. |
| `MBS_GPU_MULTI_MAX=N` | Cap automatic GPU multi batching. |
| `MBS_GPU_GROUPS=1` | Distribute variable latitude-phi theta groups across visible GPUs. |
| `MBS_GPU_MULTI_K_FULL=1` | Experimental fused multi-`k_eq` diffraction. |
| `MBS_GPU_BEAM_STATS=1` | Print GPU beam packing/count diagnostics. |
| `MBS_SHARED_BETA_GROUP=N` | Override shared-batch beta grouping. |
| `MBS_SHARED_ORIENT_CHUNK=N` | Override shared-batch orientation chunk size. |
| `MBS_PARALLEL_MEM_FRACTION` | Memory fraction used by parallel multi-size scheduler. |
| `MBS_PARALLEL_SHARED_KEQ=1` | Environment equivalent of shared `k_eq` batching. |
| `MBS_PARALLEL_KEQ_BATCH_RATIO=R` | Environment default for shared batch ratio. |

Diagnostic variables exist for FFT refinement, autofull heuristics, trace prefilters, and internal debug ranges. They are intentionally not recommended for normal production runs; inspect `grep -R "getenv(\"MBS_" src` when debugging a specific subsystem.

## Output

Main output is a whitespace-separated `.dat` file:

```text
ScAngle 2pi*dcos M11 M12 M13 M14 M21 M22 M23 M24 M31 M32 M33 M34 M41 M42 M43 M44
```

| Column | Meaning |
|---|---|
| `ScAngle` | Scattering angle theta in degrees. |
| `2pi*dcos` | Solid-angle ring weight for this theta row. |
| `Mij` | Mueller matrix elements after beam and orientation averaging. |

Common companion outputs:

| Output | Created by | Meaning |
|---|---|---|
| `_noshadow` file | `--noshadow_output` | Mueller matrix without shadow/external beam contribution. |
| `_betas/` directory | `--save_betas` | Per-beta Mueller diagnostics. |
| checkpoint files | `--checkpoint` | Resume data for long orientation-grid runs. |
| Jones output | `--jones` | Raw Jones matrices in supported PO modes. |

Warnings about optical-theorem/integral mismatch are diagnostics. For non-absorbing runs the physical absorption is fixed to zero; the printed integral characteristic is a consistency check, not a replacement for the output Mueller matrix.
