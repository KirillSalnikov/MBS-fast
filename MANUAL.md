# MBS-raw: Physical Optics for Ice Crystals — Manual

## Overview

MBS-raw computes light scattering by non-spherical particles using the Physical Optics (PO) method. It traces geometric optics rays through the particle (reflections, refractions), then applies Kirchhoff diffraction to each output beam to obtain the far-field Mueller matrix.

**Optimized version**: ~70x faster than original, AVX-512/AVX2 SIMD, OpenMP parallel, Sobol quasi-random orientations. Builds for Intel, AMD Zen 2 (EPYC 7H12), and AMD Zen 4 (Ryzen 7000 / EPYC Genoa).

Default `make` and `bash build.sh` auto-select CPU flags. On AMD EPYC 7H12
they use `-march=znver2 -mtune=znver2`; on other CPUs they use
`-march=native -mtune=native`. For manual override, set `ARCH_FLAGS=...` or
`CXXFLAGS=...`.

---

## Quick Start

```bash
# Full auto: adaptive orientations + auto grids + auto beam cutoff
mbs_po --po --auto 0.05 \
    -p 1 10 10 -w 0.532 --ri 1.31 0 -n 8 --close

# Full auto including n search
mbs_po --po --autofull 0.05 \
    -p 1 10 10 -w 0.532 --ri 1.31 0 --close

# Manual: fixed Sobol orientations with auto grids
mbs_po --po --sobol 1024 --auto_tgrid --auto_phi \
    -p 1 10 10 -w 0.532 --ri 1.31 0 -n 8 \
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

#### `--k_eq X`

Resize a file particle by equivalent-volume radius instead of by maximal
dimension. This is the preferred size control for irregular particles when
comparing to ADDA.

Definitions:

```
V      = particle volume
r_eq  = (3 V / (4 pi))^(1/3)
k_eq  = 2 pi r_eq / lambda
```

If the input particle has equivalent radius `r_eq0`, wavelength `lambda`, and
the requested value is `K`, the code applies the scale factor

```
s = (K lambda / (2 pi)) / r_eq0
Dmax_new = s Dmax_old
V_new    = s^3 V_old
```

`--k_eq` and `--rs` are mutually exclusive. `--k_eq` requires `-w`.

Example:

```bash
mbs_po --po --oldauto 2 --pole \
    --pf shapeA64_mbs.dat --k_eq 20.76 \
    --ri 1.6 0.002 -w 1.064 -n 8 \
    --tgrid scattering_angles --nphi 600 --gpu --fft --close
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

Depth in random/oldauto logs means the internal reflection order of traced
beams. Depth 0 is the external/shadow contribution, depth 1 is the first
refracted/reflected generation, and so on up to `-n`.

#### `--oldauto DIV`

Physics-based regular orientation grid. It computes the beta/gamma grid from
the diffraction-limited angular scale of the current particle and then divides
that grid by `DIV`.

Let

```
L = Dmax
Delta_theta = 0.69 lambda / L    [radians]
Delta_theta_deg = Delta_theta * 180/pi
Delta_orient = Delta_theta_deg / ring_points
```

For particle symmetry ranges `beta_sym` and `gamma_sym` in degrees:

```
N_beta_full  = ceil(beta_sym  / Delta_orient)
N_gamma_full = ceil(gamma_sym / Delta_orient)
N_beta       = ceil(N_beta_full  / DIV)
N_gamma      = ceil(N_gamma_full / DIV)
```

The actual beta count in the log is `N_beta + 1`; gamma count is `N_gamma`.
For example, hexagonal symmetry has `beta_sym=90 deg`, `gamma_sym=60 deg`.

`DIV=2` is the dense mode used for the Greek shape runs. Larger values are
faster and less accurate. `--oldauto` uses `--nphi` if supplied, otherwise
uses `--grid` phi or the default.

#### `--ring_points N`

Number of orientation samples per diffraction ring used by `--oldauto`.
Default is 3.

Formula:

```
Delta_orient = (0.69 lambda / Dmax) * 180/pi / N
```

Increasing `N` increases the orientation count roughly as `N^2`.

#### `--pole`

Fast pole shortcut for regular beta/gamma grids. At `beta=0` and `beta=pi`
all gamma rotations describe the same physical orientation, so the code traces
one gamma and multiplies its weight by the full gamma count.

This keeps the normalization equivalent to the full grid while avoiding
duplicated pole work.

#### `--mirror_gamma`

Use mirror symmetry in the regular `--oldauto` beta/gamma grid. For particles
with a mirror plane inside the rotational symmetry sector, the fundamental
gamma interval is halved:

```
gamma_max_effective = gamma_sym / 2
```

For a hexagonal particle this changes the oldauto gamma domain from `0..60 deg`
to `0..30 deg`. The oldauto formula is then applied to the reduced interval:

```
N_gamma_full = ceil((gamma_sym/2) / Delta_orient)
N_gamma      = ceil(N_gamma_full / DIV)
```

The quadrature normalization uses the reduced fundamental domain. The omitted
mirror half is restored in the Mueller grid by pairing azimuths `phi` and
`-phi`:

```
M_mirror_averaged(phi) = 0.5 * (M(phi) + P M(-phi) P)
P = diag(1, 1, -1, -1)
```

This preserves the even Mueller components and cancels the components that are
odd under the mirror operation. It is valid when the particle really has the
mirror symmetry in the input geometry.

On a hex plate test, `--oldauto 4` changed the grid from `30 x 20` to
`30 x 10` gamma samples and matched the full-gamma `M11(theta)` curve with
maximum relative difference about `4.4e-6`.

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

Use `--checkpoint` with very long orientation-file or oldauto/random runs to save and resume
chunked/beta-ring progress. Checkpointing is disabled by default to avoid extra I/O in
ordinary runs.

#### `--sym B G`

Override particle symmetry for Sobol/adaptive orientation generation.
B = beta symmetry divisor, G = gamma symmetry divisor.
beta_sym = pi/B, gamma_sym = 2*pi/G.

**WARNING**: Usually NOT needed. Particle symmetry is detected automatically from `-p` type (e.g. hex prism: beta_sym=90°, gamma_sym=60°). Use `--sym` only for custom particles loaded with `--pf`.

- `--sym 2 6`: hex prism (beta in [0, pi/2=90°], gamma in [0, 2pi/6=60°]) — 12x reduction
- `--sym 1 1`: no symmetry (full sphere)

**Note**: `--auto` and `--autofull` use particle auto-symmetry. Do NOT add `--sym`.

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

Automatically generate optimal non-uniform theta grid based on size parameter x = pi*D/lambda.
Five zones with adaptive step sizes:

| Zone | Range | Step | Purpose |
|------|-------|------|---------|
| Fine | 0 to 10×peak_width | peak_width/20 | Resolve forward diffraction peak |
| Transition | to 20° (or 40° if x<30) | geometric ×1.3, max 1° | Smooth transition |
| Medium | to 120° | 0.5° (x<100) or 1° | Side scattering |
| Side/back | 120° to 175° | 1° | Backscattering |
| Near-backscatter | 175° to 180° | 0.25° | LDR, depolarization |

**Total points vs size parameter:**

| x | peak width | N_theta | Fine step | N_phi (auto) |
|---|---|---|---|---|
| 10 | 18.0° | 316 | 0.50° | 48 |
| 50 | 3.6° | 389 | 0.18° | 48 |
| 100 | 1.8° | 383 | 0.09° | 48 |
| 200 | 0.9° | 395 | 0.045° | 60 |
| 500 | 0.36° | 404 | 0.018° | 90 |
| 1000 | 0.18° | 408 | 0.009° | 132 |
| 2000 | 0.09° | 412 | 0.0045° | 180 |
| 3000 | 0.06° | 414 | 0.003° | 222 |

See `fig_auto_grid_points.pdf` and `fig_auto_grid_steps.pdf` for visualizations.

#### `--auto_phi`

Automatically select N_phi based on size parameter. Formula: `N_phi = max(48, 6*ceil(sqrt(x)/1.5))`.

Phi convergence **depends on x**: fixed N_phi=48 gives 2% Q_sca error at x=18 but 13% at x=600.
Auto_phi keeps error uniform at ~2-3% across all sizes.

| x | auto N_phi | Q_sca error |
|---|---|---|
| 18 | 48 | 2.0% |
| 59 | 48 | 2.5% |
| 177 | 66 | 2.2% |
| 591 | 120 | 2.9% |

See `fig_auto_phi.pdf` for convergence plots.

**Recommended**: use both `--auto_tgrid --auto_phi` together.

#### `--point`

Compute only the exact backscatter point (theta = 180 deg). Requires `--po`. Currently commented out in the source code (non-functional); retained for future use.

#### `--filter ANGLE`

Scattering angle filter in degrees. Sets a backscattering cone aperture: only scattering within ANGLE degrees of the backscatter direction (180 deg) is accumulated. Used with `--all` or `--random` PO modes.

**Example**: `--filter 5` collects only scattering within 5 deg of exact backscatter.

---

### Multi-Size Computation

#### `--multigrid DMIN DMAX N`

Generate `N` logarithmically spaced maximal dimensions from `DMIN` to `DMAX`
and compute them in one run.

```
D_i = exp(log(DMIN) + i (log(DMAX)-log(DMIN))/(N-1))
```

With `--oldauto`, the particle is first resized to `DMAX`, the orientation
grid is built for this largest size, and the tracing is performed once. The
prepared beam geometry is then scaled for every `D_i`.

#### `--multikeq KMIN KMAX N`

Same as `--multigrid`, but sizes are equivalent-volume size parameters:

```
k_eq,i = exp(log(KMIN) + i (log(KMAX)-log(KMIN))/(N-1))
```

With `--oldauto`, the common trace is done at `KMAX`.

#### `--multikeq_list FILE`

Read exact `k_eq` values from a text file. This is the safest mode when the
target sizes come from an external ADDA table and are not exactly log-spaced.

File format:

```
# first column is k_eq; optional second column is ignored by the parser
2.27 adda
3.03 adda
3.99 adda
5.23 adda
65 extrap
```

The code uses the maximum listed value as the reference size:

```
K_ref = max(k_eq list)
```

Then it traces once at `K_ref` and diffracts every listed size. Output labels
use the exact values, for example:

```
OUT/OUT_keq2p27.dat
OUT/OUT_keq2p27_noshadow.dat
```

Important: only tracing/preparation is shared. Diffraction is still computed
for every size because the phase changes with size:

```
exp(i k optical_path) and exp(i k r)
```

For a size scale `s = K_i / K_ref`, the prepared geometry is scaled as:

```
coordinates       -> s coordinates
areas             -> s^2 areas
optical lengths   -> s optical lengths
incoming energy   -> s^2 incoming energy
```

The diffraction integral must then be recomputed for the scaled beams.

#### `--multigrid_parallel JOBS`

Run `--multigrid` or `--multikeq` as independent child processes instead of
using shared tracing. This is useful when one GPU should process several
independent sizes concurrently, or when fault isolation matters more than
shared tracing.

Each child receives one size via `--rs` or `--k_eq`.

With CUDA builds, `--multigrid_parallel 0` means auto: use the number of
visible GPU devices. Devices are assigned by setting `CUDA_VISIBLE_DEVICES`
for each child. `--gpu_devices 0,1,2,3,4` restricts or orders the device list.

The parent also divides host RAM between children. It exports
`MBS_HOST_MEM_BUDGET_MB` to each child, and the oldauto/random GPU
auto-chunker uses that as an upper host-memory budget. The default total
fraction is controlled by:

```
MBS_PARALLEL_MEM_FRACTION=0.70
```

#### `--multikeq_shared_batches`

Opt-in trace reuse for multi-GPU `k_eq` scans. Instead of launching one child
per size, the parent groups nearby `k_eq` values into per-GPU batch files and
launches each child with its own `--multikeq_list`. Inside that child,
oldauto traces/prepares once at the largest `k_eq` in the batch and then
recomputes diffraction for every size in the batch.

This is faster for dense size scans because tracing is reused, but it is not
bit-identical to independent oldauto for the smaller sizes in a batch: those
sizes use the orientation grid chosen for the batch reference size. This is
usually a denser grid for the smaller sizes, but it is still a different
quadrature. Therefore the exact default remains one child per size.

#### `--multikeq_batch_ratio R`

Maximum allowed ratio inside an opt-in shared k_eq batch:

```
kmax / kmin <= R
```

Default:

```
R = 1.05
```

Use smaller values, e.g. `1.01` or `1.02`, when the result must stay closer to
independent oldauto grids. Use larger values only for exploratory scans where
trace reuse is more important than strict quadrature matching.

#### `--multigrid_threads N`

OpenMP threads per child process in `--multigrid_parallel`. If omitted, GPU
children default to one thread so several child processes can share CPU cores.

---

### Performance Options

#### `--threads N`

Set the OpenMP worker thread count. If omitted, the program uses physical CPU
cores by default, not hyperthreads.

This affects CPU tracing/preparation and CPU diffraction. With `--gpu`, CUDA
does the diffraction, but CPU threads are still used for ray tracing and beam
preparation. In shared multikeq runs this phase can be large because the trace
is done at the largest listed size.

For one GPU job that should also use CPU tracing:

```bash
--threads 12 --gpu --fft
```

For many GPU jobs sharing the same machine:

```bash
--threads 1 --gpu --fft
```

#### `--gpu`

Use the CUDA backend for Phase 2 diffraction. Requires a binary built with
`USE_CUDA=1`.

Build:

```bash
make USE_CUDA=1
```

For Ampere/Ada cards with older drivers, build with an explicit architecture
to avoid PTX JIT errors:

```bash
make USE_CUDA=1 NVCCFLAGS="-O3 -std=c++11 -arch=sm_86 -U_GNU_SOURCE -D_DEFAULT_SOURCE -D_XOPEN_SOURCE=700"
```

CUDA implementation:

- CPU traces rays and prepares `PreparedOrientation` / `PreparedBeam` data.
- Beam polygon vertices, Jones matrices, projected centers, optical lengths,
  and areas are packed into GPU buffers.
- CUDA kernels evaluate diffraction for many scattering directions and beams.
- Results are reduced into local Mueller arrays and copied back to CPU.

Memory roughly scales as:

```
O(batch_orientations * beams_per_orientation * vertices_per_beam)
+ O(N_phi * N_theta * Mueller/Jones accumulators)
```

The code chooses GPU orientation batches dynamically so a 12 GB card can run
large cases without allocating the full orientation set on the GPU at once.

#### `--gpu_trace`

Enable the experimental CUDA prefilter for non-convex tracing candidates.
This is not full GPU ray tracing. It does not replace the exact CPU polygon
intersection. CUDA only projects beam/facet pairs and rejects pairs whose 2D
bounding boxes cannot overlap; all surviving pairs still go through the
existing CPU `Intersect()` path.

The implementation batches many beams from one tracing stack layer, copies
beam records, facet records, and compact `(beam, facet)` index pairs to the
GPU, then copies one byte per candidate back to the CPU. Workspaces are
thread-local so OpenMP tracing can call the prefilter safely.

Environment controls:

```
MBS_TRACE_CPU_PREFILTER=0          # disable conservative projected CPU prefilter
MBS_TRACE_CPU_PREFILTER_MARGIN=8   # projected bounding-box margin, default 8
MBS_TRACE_TREE_RESERVE=16384       # initial non-convex beam-tree capacity
MBS_GPU_TRACE_BATCH_BEAMS=1024      # beams collected per CUDA prefilter batch
MBS_GPU_TRACE_MIN_CANDIDATES=65536  # below this, skip CUDA and use CPU path
MBS_GPU_TRACE_OPENMP=1              # allow CUDA prefilter inside OpenMP tracing
MBS_FORCE_TRACK_IDS=1              # keep BigInteger track IDs even when Im(ri)=0
```

For non-convex particles the conservative projected CPU prefilter is enabled
by default. Set `MBS_TRACE_CPU_PREFILTER=0` to disable it. The prefilter is a
broad-phase test:
the CPU tracer checks whether the projection of a beam onto a candidate facet
plane can overlap the facet bounding box. If the expanded boxes are disjoint,
the expensive exact polygon intersection is skipped. Increasing the margin makes
the test less aggressive. The default margin is intentionally conservative
(`8` in particle length units); reduce it only after checking the Mueller
matrix against a no-prefilter reference for the specific particle.

`MBS_TRACE_TREE_RESERVE` controls the initial `std::vector<Beam>` capacity for
the non-convex beam stack. The default `16384` reduces reallocations and beam
copies in deep non-convex traces. Lower it only if memory is tight; raise it if
logs show much larger per-orientation beam counts.

For non-absorbing particles (`Im(ri)=0`) the tracer skips BigInteger track-id
updates by default because absorption path recovery does not need them. This is
mathematically neutral for Mueller output and speeds up deep non-convex tracing.
Set `MBS_FORCE_TRACK_IDS=1` for debugging or explicit track-output workflows.

Current status: this is a correctness-checked experimental path, not a default
speed path. It is useful as groundwork for a full GPU tracing layer, but on
small and medium Afine30 tests the CPU tracer is still faster because the
remaining exact intersection and CPU/GPU synchronization dominate.
With OpenMP tracing (`--threads > 1`) the CUDA prefilter is disabled unless
`MBS_GPU_TRACE_OPENMP=1`, because many CPU threads launching small CUDA
batches serialize and can be much slower than the CPU path.

For large non-convex particles, Phase 1 can still be the main bottleneck:
the current GPU backend accelerates diffraction, while geometric ray tracing,
beam splitting, internal-reflection expansion, and exact polygon intersections
remain CPU-side. Full CUDA tracing has not been completed yet.

#### `--fft`

Enable the experimental CUDA FFT phi-interpolation backend. Requires `--gpu`.

This is not aperture pFFT/FMM and not a 3D spatial FFT. It accelerates the
azimuthal scattering grid by computing diffraction on a reduced direct
azimuth grid and reconstructing the requested `N_phi` grid using cuFFT-based
angular Fourier interpolation.

The log line:

```
GPU FFT phi interpolation: direct Nphi=60, output Nphi=600 (factor=10)
```

means direct diffraction was evaluated for 60 azimuth samples and interpolated
to 600 output samples.

Environment controls:

```
MBS_FFT_PHI_FACTOR=auto   # default
MBS_FFT_PHI_FACTOR=10     # force output/direct ratio
MBS_FFT_CHECK=1           # optional diagnostics/check path
MBS_GPU_NO_ATOMICS=1      # use orientation-grid CUDA reduction path
MBS_GPU_MEM_FRACTION=0.8  # fraction of free GPU memory available to batches
MBS_GPU_BLOCK=128|256|512 # CUDA block size for diffraction kernels, default 128
MBS_SHARED_PIPELINE=1     # overlap CPU trace and GPU diffraction in shared multikeq
MBS_SHARED_ORIENT_CHUNK=64 # global scheduler chunk size for shared multikeq
```

By default the optimized PO path writes only the full Mueller matrix. The
legacy `_noshadow` matrix is skipped to avoid an extra no-shadow reduction,
host copy, FFT interpolation, and output file. Use `--noshadow_output` when
that diagnostic matrix is needed.

Accuracy depends on smoothness in phi. It is usually appropriate for
orientation-averaged Mueller matrices with dense `N_phi`; validate against
`--gpu` without `--fft` for new particle classes.

`MBS_GPU_NO_ATOMICS=1` switches the CUDA diffraction path from beam-grid
atomic accumulation to orientation-grid accumulation. For many beams per
orientation this can reduce atomic contention and improve GPU occupancy. It
keeps the same mathematical result; small last-digit differences are expected
from a different summation order.

#### `MBS_SHARED_BETA_GROUP=N`

Environment variable for shared `--oldauto --multikeq*` runs. It groups `N`
beta blocks before the diffraction pass when the legacy beta-block scheduler is
used. Current GPU shared multikeq runs use the global orientation scheduler
below, so this variable is mostly a memory/diagnostic knob rather than the main
batch-size control.

Default:

```
auto: min(nThreads/4, 4), bounded by nBeta
```

With `N>1`, the code traces several beta blocks into one prepared-orientation
batch and then diffracts that larger batch for every size. This reduces the
number of CUDA/cuFFT launches and gives the GPU longer work packets.

Memory tradeoff:

```
prepared memory ~ N * N_gamma * average_beams_per_orientation
```

Recommended only for legacy beta-block tests:

```
MBS_SHARED_BETA_GROUP=1   # smallest memory and latency
MBS_SHARED_BETA_GROUP=2   # conservative larger beta block
MBS_SHARED_BETA_GROUP=4   # longer GPU packets, more RAM and longer tails
```

Large `N` can be slower for non-convex file particles: one hard orientation can
hold the entire beta group and leave the GPU idle.

#### `MBS_SHARED_ORIENT_CHUNK=N`

Environment variable for shared `--oldauto --multikeq*` GPU runs. It enables the
global orientation scheduler:

```
globalIndex = betaIndex * N_gamma + gammaIndex
```

The scheduler traces contiguous chunks of `N` global orientations and immediately
passes each completed chunk to CUDA diffraction. With `MBS_SHARED_PIPELINE=1`,
the CPU traces chunk `k+1` while the GPU diffracts chunk `k`.

Default for `--gpu` shared multikeq:

```
MBS_SHARED_ORIENT_CHUNK=64
```

Memory tradeoff:

```
prepared memory ~ chunk * average_beams_per_orientation
```

Recommended values:

```
MBS_SHARED_ORIENT_CHUNK=64   # current default; good Afine30/Greek starting point
MBS_SHARED_ORIENT_CHUNK=128  # fewer CUDA launches, sometimes faster on large GPUs
MBS_SHARED_ORIENT_CHUNK=32   # lower latency, useful when GPU has long idle gaps
```

On Afine30 with `--grid 0 180 600 181 --gpu --fft --oldauto 2`, chunks of
`32..64` improved the observed scheduler rate from about `4.3 orient/s` to
about `5.5 orient/s` on the 4070 SUPER node and reduced host RAM from roughly
`10 GB` to below `2 GB`.

#### `MBS_SHARED_PIPELINE=1`

Environment variable for shared multi-size GPU runs. It starts CPU tracing of
the next global orientation chunk while the current chunk is being diffracted on
the GPU. Use it together with `MBS_SHARED_ORIENT_CHUNK`.

Recommended for file-particle multikeq runs:

```bash
MBS_SHARED_PIPELINE=1 MBS_SHARED_ORIENT_CHUNK=64 \
./bin/mbs_po_float_fast --pf Afine30.dat --multikeq_list keq_list.txt \
    --ri 1.6 0.002 -n 14 --po --oldauto 2 \
    --beam_cutoff_j 0.001 --beam_cutoff_area 0.002 \
    --trace_cutoff_importance 0.0001 --trace_max_beams 20000 \
    --grid 0 180 600 181 -w 1.064 --gpu --fft --close
```

If GPU utilization oscillates with long zero-utilization gaps, reduce
`MBS_SHARED_ORIENT_CHUNK` to `32`. If GPU launch overhead dominates and host RAM
is available, try `128`.

#### `--beam_cutoff EPS`

Compatibility shorthand for beam cutoffs. It sets the separate J and area
beam thresholds to `EPS`.

Modern code exposes separate tests:

```
--beam_cutoff_j EPS
--beam_cutoff_area EPS
--beam_cutoff_importance EPS
```

The quantities are dimensionless and normalized within the current
orientation/chunk:

```
J_rel    = |J|^2 / max(|J|^2)
A_rel    = area / max(area)
I_rel    = |J|^2 area / max(|J|^2 area)
```

`--beam_cutoff_j EPS` skips beams with:

```
J_rel < EPS
```

`--beam_cutoff_area EPS` skips beams with:

```
A_rel < EPS
```

`--beam_cutoff_importance EPS` skips beams with:

```
I_rel < EPS
```

Use the importance cutoff for a size-independent single threshold; it keeps
beams that are either strong, large-area, or both.

EPS = target accuracy (0 to 1). Set automatically by `--auto`/`--adaptive`. Can also be used standalone with `--random` or `--sobol`.

| EPS | Typical beams skipped | Use case |
|-----|----------------------|----------|
| 0.1 (10%) | ~75% | Fast estimate |
| 0.05 (5%) | ~65% | General purpose |
| 0.01 (1%) | ~45% | Accurate |
| 0.001 (0.1%) | ~5% | Reference |

**Example**:
```bash
mbs_po --po --random 10 30 --beam_cutoff 0.05 \
    -p 1 42.857 30 -w 0.532 --ri 1.31 0 -n 14 \
    --grid 0 180 48 180 --close
```

#### `--trace_cutoff EPS`

Compatibility shorthand for tracing cutoffs. It sets:

```
--trace_cutoff_j EPS
--trace_cutoff_area EPS
```

The trace cutoffs prune internal beam trees before diffraction, so they can
reduce both tracing time and later diffraction work.

#### `--trace_cutoff_j EPS`

Stop tracing a beam branch when its relative Jones intensity is small:

```
|J|^2 / |J_initial|max^2 < EPS
```

This is independent of particle size because both numerator and denominator
scale in the same internal amplitude units.

#### `--trace_cutoff_area EPS`

Stop tracing a beam branch when its relative aperture area is small:

```
area / area_initial,max < EPS
```

This removes geometrically tiny internal beams. Use with care for exact
backscatter studies because small apertures can still contribute sharp
features.

#### `--trace_cutoff_importance EPS`

Stop tracing a beam branch when its relative internal importance is small:

```
(|J|^2 * area) / max_initial(|J|^2 * area) < EPS
```

This is usually the safest speed cutoff when particle size changes, because it
combines field strength and aperture size in one dimensionless quantity. It is
applied while building the non-convex beam tree, including newly split beam
parts after facet clipping, so tiny internal leftovers are not put back into
the tracing stack.

#### `--trace_max_beams N`

Abort tracing one orientation after more than `N` beam nodes. `0` disables the
limit. This is a safety valve for pathological non-convex cases where internal
beam splitting explodes.

Recommended heavy-case starting point:

```bash
--beam_cutoff_importance 0.0001 \
--trace_cutoff_j 0.0001 \
--trace_cutoff_area 0.0001 \
--trace_max_beams 20000
```

#### `--chunk N`

Maximum orientation/gamma chunk size for memory-aware streaming modes. Smaller
chunks reduce RAM pressure and make checkpointing/resume safer; larger chunks
reduce overhead. In GPU oldauto/random runs without beam cutoffs the automatic
default is 64 gamma orientations per beta ring; `--chunk` overrides it. Use a
smaller value if RAM pressure appears, or a larger value if the particle is
simple and memory is still far from the limit.

Memory is approximately:

```
O(chunk * average_beams_per_orientation * PreparedBeam_size)
```

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

#### `--checkpoint`

Enable checkpoint save/resume for `--orientfile` and `--oldauto`/`--random`
runs. For `--orientfile`, the checkpoint is written after each completed
orientation chunk. For `--oldauto`/`--random`, it is written after each completed
beta ring and stores the accumulated Mueller matrices plus energy/extinction
accounting. On restart with the same command and `-o`, the existing output folder
is reused and completed beta rings are skipped. The checkpoint is removed after
successful completion. Disabled by default.

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

### Native AVX2/FMA Baseline

```bash
bash build.sh            # -> bin/mbs_po
```

Uses `-O3 -march=native -mavx2 -mfma -fopenmp`. Requires GCC >= 9.
Use the CPU-specific scripts below for AVX-512 or tuned EPYC builds.

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
./build.sh
# Or compile manually:
g++ -O3 -march=haswell -std=gnu++11 -funroll-loops -fopenmp \
    $(find src -not -path '*/bigint/*' -name '*.cpp') \
    $(find src/bigint -name '*.cc') -Isrc -Isrc/math ... -lm -lgomp
```

The code auto-detects AVX-512 at compile time (`#ifdef __AVX512F__`).

### CUDA Build

```bash
make USE_CUDA=1
```

Optional precision selector:

```bash
make USE_CUDA=1 GPU_PRECISION=double   # default
make USE_CUDA=1 GPU_PRECISION=float    # experimental FP32 CUDA path
make USE_CUDA=1 GPU_FAST_MATH=1        # double kernels with nvcc fast math
```

`GPU_PRECISION=float` changes internal CUDA diffraction buffers/kernels to
single precision while keeping the public CPU/output path in double where
applicable. On RTX 3080 Ti / 4070 SUPER this can speed CUDA diffraction by
about 1.4-1.5x on Afine30 smoke tests, while integrated `Q` and `M11` remain
close. Small polarization components are more sensitive: on one Afine30
control run the largest significant FP32 relative error was in `M14/M41` near
theta=0 deg, about 4.5%. `GPU_FAST_MATH=1` with float was faster again but
increased that polarization error to about 8%, so treat it as a diagnostic or
M11-only speed mode unless validated for the target observable.

Convenience CUDA targets build separate binaries without replacing the default
double binary:

```bash
make cuda_float       # bin/mbs_po_float
make cuda_float_fast  # bin/mbs_po_float_fast
make cuda_double_fast # bin/mbs_po_double_fast
make cuda_variants    # all three
```

For RTX 3080 Ti and RTX 4070 SUPER, this explicit architecture build has been
used successfully:

```bash
make USE_CUDA=1 \
  NVCCFLAGS="-O3 -std=c++11 -arch=sm_86 -U_GNU_SOURCE -D_DEFAULT_SOURCE -D_XOPEN_SOURCE=700"
```

Use `--gpu` at runtime to select CUDA. Add `--fft` to use cuFFT phi
interpolation.

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
mbs_po --po --oldauto 2 --pole \
    --pf shapeA64_mbs.dat --multikeq 2.27 160 15 \
    --ri 1.6 0.002 -w 1.064 -n 8 \
    --tgrid scattering_angles --nphi 600 --gpu --fft --close
```

### Exact k_eq list from ADDA table
```bash
cat > keq_list.txt <<EOF
2.27 adda
3.03 adda
3.99 adda
5.23 adda
EOF

mbs_po --po --oldauto 2 --pole \
    --pf shapeA64_mbs.dat --multikeq_list keq_list.txt \
    --ri 1.6 0.002 -w 1.064 -n 8 \
    --tgrid scattering_angles --nphi 600 --gpu --fft \
    --beam_cutoff_importance 0.0001 \
    --trace_cutoff_j 0.0001 --trace_cutoff_area 0.0001 \
    --trace_max_beams 20000 --close
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
| `--k_eq` | K | Particle | Resize file particle so `2*pi*r_eq/lambda = K` |
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
| `--oldauto` | DIV | Orientations | Physics-based regular grid divided by DIV |
| `--ring_points` | N | Orientations | Points per diffraction ring for oldauto grid |
| `--mirror_gamma` | (none) | Orientations | Use mirror symmetry; oldauto gamma domain is halved |
| `--pole` | (none) | Orientations | Use one gamma at beta poles with full weight |
| `--montecarlo` | N | Orientations | Monte Carlo random orientations |
| `--fixed` | BETA GAMMA | Orientations | Single orientation (degrees) |
| `--orientfile` | FILENAME | Orientations | Orientations from file (radians) |
| `--checkpoint` | (none) | Orientations | Save/resume long `--orientfile` and `--oldauto`/`--random` runs |
| `--chunk` | N | Orientations | Max orientation/gamma chunk size |
| `--sym` | B G | Orientations | Symmetry override: beta/B, gamma/(2pi/G) |
| `-b` | MIN MAX | Orientations | Beta range override (degrees) |
| `-g` | MIN MAX | Orientations | Gamma range override (degrees) |
| `--grid` | T1 T2 N_PHI N_TH | Scattering | Uniform scattering grid |
| `--tgrid` | FILENAME | Scattering | Non-uniform theta grid from file |
| `--auto_tgrid` | (none) | Scattering | Auto-generate optimal theta grid |
| `--auto_phi` | (none) | Scattering | Auto-select phi grid |
| `--nphi` | N | Scattering | Override N_phi |
| `--point` | (none) | Scattering | Backscatter point only (disabled) |
| `--filter` | ANGLE | Scattering | Backscatter cone aperture (degrees) |
| `--multigrid` | DMIN DMAX N | Multi-size | Log-spaced Dmax sizes; shared trace in oldauto |
| `--multikeq` | KMIN KMAX N | Multi-size | Log-spaced k_eq sizes; shared trace in oldauto |
| `--multikeq_list` | FILE | Multi-size | Exact k_eq values; shared trace at max |
| `--multigrid_parallel` | JOBS | Multi-size | Run sizes as child processes |
| `--multigrid_threads` | N | Multi-size | Threads per child process |
| `--gpu_devices` | LIST | Multi-size | CUDA devices for parallel GPU children |
| `--multikeq_shared_batches` | (none) | Multi-size | Group nearby k_eq values per GPU child and reuse tracing |
| `--multikeq_batch_ratio` | R | Multi-size | Max `kmax/kmin` inside a shared k_eq batch |
| `--threads` | N | Performance | OpenMP worker threads |
| `--gpu` | (none) | Performance | CUDA diffraction backend |
| `--fft` | (none) | Performance | cuFFT phi interpolation backend; requires `--gpu` |
| `--beam_cutoff` | EPS | Performance | Shorthand for J/area beam cutoffs |
| `--beam_cutoff_j` | EPS | Performance | Beam skip by relative `|J|^2` |
| `--beam_cutoff_area` | EPS | Performance | Beam skip by relative area |
| `--beam_cutoff_importance` | EPS | Performance | Beam skip by relative `|J|^2*area` |
| `--trace_cutoff` | EPS | Performance | Shorthand for trace J/area cutoffs |
| `--trace_cutoff_j` | EPS | Performance | Trace prune by relative `|J|^2` |
| `--trace_cutoff_area` | EPS | Performance | Trace prune by relative area |
| `--trace_cutoff_importance` | EPS | Performance | Trace prune by relative `|J|^2*area` |
| `--trace_max_beams` | N | Performance | Max beam nodes per orientation |
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

1. **Optical-theorem phase convention**: by default `C_ext_OT` is computed from
   the object-centered spherical-wave forward amplitude `S`:
   `C_ext_OT = lambda * Im(f00 + f11) * w_orient`.  This matches the
   Gao/Yang/Kattawar optical-theorem convention after removing the far-screen
   `exp(i k (r-z))` factor.  `--ot_ping D` is a legacy compatibility rotation
   for old files that stored the far-screen amplitude:
   `Im(F) cos(2 k D) - Re(F) sin(2 k D)`.  It is not a physical particle
   parameter.

2. **Q_sca > 2 at large x**: PO does not enforce the optical theorem. Q_sca grows with x instead of approaching 2. Use IGOM correction or renormalize.

3. **Polarization elements M33, M34, M44**: MBS-raw's RotateJones and GOAD's Karczewski matrix have the same norm but different structure. M11 agrees between codes, but M33/M34/M44 differ.

4. **Small particles (x < 20)**: PO is not valid. Use ADDA or T-matrix methods.

5. **Backscattering convergence**: M11(180 deg) requires 10-100x more orientations than forward scattering.

6. **LTO**: `-flto` causes ~20% regression on current code due to GCC inlining heuristics interacting poorly with the hand-inlined hot loop. Do not use.

7. **Absorption**: DiffractInclineAbs is optimized but not as fast as the non-absorbing path.

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
