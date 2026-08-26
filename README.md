# MBS-fast

MBS-fast calculates light scattering by faceted particles using geometrical ray
tracing with either Physical Optics (PO) diffraction or Geometrical Optics
(GO) output. It supports fixed orientations, several deterministic and
quasi-random orientation averages, adaptive convergence, CPU MPI/OpenMP, CUDA
diffraction, and process-parallel size scans across multiple GPUs.

MBS-fast is an optimized and numerically hardened development of the original
[`MBS-raw`](https://github.com/Heart-Under-Blade/MBS-raw) codebase. Consult that
repository for the upstream implementation and project history.

## Method Reference

When using the MBS-1 method, cite D. N. Timofeev, A. V. Konoshonkin, and
N. V. Kustova, "Modified Beam-Splitting 1 (MBS-1) Algorithm for Solving the
Problem of Light Scattering by Nonconvex Atmospheric Ice Particles,"
*Atmospheric and Oceanic Optics* **31**(6), 642-649 (2018):
[publisher page](https://link.springer.com/article/10.1134/S1024856018060179),
[DOI: 10.1134/S1024856018060179](https://doi.org/10.1134/S1024856018060179).

## Documentation

The generated command help is the authoritative list of accepted options,
aliases, argument counts, and short descriptions:

```bash
cpu/bin/mbs_po_mpi --help
cpu/bin/mbs_po_mpi --help-debug   # legacy, diagnostic, and experimental flags
cpu/bin/mbs_po_mpi --version      # revision and compiled feature set
```

Long-form manuals describe the physical model, coordinate systems, ray
tracing, Fresnel coefficients, Jones/Mueller matrices, diffraction, clipping,
orientation averaging, CPU vectorization, and CUDA implementation:

| Language | Editable source | Typeset manual |
|---|---|---|
| English | [`docs/MANUAL.tex`](docs/MANUAL.tex) | [`docs/MANUAL.pdf`](docs/MANUAL.pdf) |
| Russian | [`docs/MANUAL_RU.tex`](docs/MANUAL_RU.tex) | [`docs/MANUAL_RU.pdf`](docs/MANUAL_RU.pdf) |

The Markdown manuals under `docs/` are searchable companion documents. The
typeset TeX/PDF manuals are the release documentation, while generated
`--help` remains authoritative for the current CLI.

## Requirements

CPU build:

- Linux on x86-64;
- GNU Make, Bash, and GNU coreutils (`timeout` is used by the release tests);
- GCC 9 or newer, or Clang 14 or newer, with C++11 and OpenMP support;
- an MPI C++ wrapper and launcher from the same installation. The supplied
  release gate is tested with Open MPI and normally uses `/usr/bin/mpicxx` and
  `/usr/bin/mpirun` when they are available.

CUDA build additionally requires an NVIDIA driver, CUDA toolkit with `nvcc`,
and cuFFT. `gpu/Makefile` calls `scripts/select_cuda_host_compiler.sh` to select
a CUDA-compatible host compiler when one is installed. CUDA and host-compiler
compatibility is determined by the installed toolkit, not by MBS-fast.

Building the manuals requires XeLaTeX, `fontspec`, `unicode-math`, Latin Modern
Math, Computer Modern Unicode, and the standard LaTeX packages imported by
`docs/MANUAL*.tex`. On Ubuntu, install the Russian manual font with
`sudo apt install fonts-cmu`; run `fc-match 'CMU Serif'` to verify discovery.

## Build

```bash
git clone https://github.com/KirillSalnikov/MBS-fast.git
cd MBS-fast

# CPU, MPI + OpenMP
make -C cpu -j

# CUDA variants
make -C gpu -j                 # FP64, precise math + LTO; default/reference
make -C gpu fp32 -j            # FP32 storage, precise math
make -C gpu fp32_fast -j       # FP32 storage + --use_fast_math
make -C gpu fp64_fast -j       # FP64 storage + --use_fast_math

# Rebuild both typeset manuals (two XeLaTeX passes per language)
make docs
```

The supported build entry points are `cpu/Makefile`, `gpu/Makefile`, and their
delegating targets in the root `Makefile`. Architecture-specific root wrappers
are intentionally not provided. Both split builds call
`scripts/detect_arch_flags.sh`; the CUDA build also calls
`scripts/select_cuda_host_compiler.sh`. Override `ARCH_FLAGS`, `GPU_ARCH`, or
`CUDA_HOST_CXX` directly when automatic detection is not suitable.

The resulting binaries are independent:

| Build | Binary |
|---|---|
| CPU MPI/OpenMP | `cpu/bin/mbs_po_mpi` |
| CUDA float | `gpu/bin/mbs_po_gpu_float` |
| CUDA float fast | `gpu/bin/mbs_po_gpu_float_fast` |
| CUDA double | `gpu/bin/mbs_po_gpu_double` |
| CUDA double fast | `gpu/bin/mbs_po_gpu_double_fast` |

The default CUDA target is FP64 without `--use_fast_math`; it is the numerical
reference build. Link-time optimization is enabled by default because it
reduces production runtime without changing CUDA arithmetic. Set `GPU_LTO=0`
to disable it; LTO and non-LTO objects use separate build directories. The
`fp32`, `fp64`, `fp32_fast`, and `fp64_fast` Make targets
select explicit profiles. FP32 stores diffraction-beam geometry and accumulated
optical quantities in FP32. Visibility/topology tracing, optical paths,
absorption, cancellation-prone polygon moments, and critical phase
trigonometry remain FP64. Validate throughput profiles against CPU double and
CUDA FP64 on representative convex, non-convex, fixed-orientation, and
orientation-averaged cases before production use. Record `--version` output.
After building, verify the selected host and CUDA profiles explicitly:

```bash
gpu/bin/mbs_po_gpu_double --version  # must report fp64-diffraction-storage
gpu/bin/mbs_po_gpu_float --version   # must report fp32-diffraction-storage
```

The build now rejects CUDA host code without exactly one precision profile;
`unknown-precision` is not an accepted production binary.

For scientific reference and final archive calculations, use
`mbs_po_gpu_double`. On consumer GPUs with weak FP64 hardware, start with
`mbs_po_gpu_float`: it keeps visibility, topology, optical paths, absorption,
polygon moments, and phase-sensitive operations in FP64 while using FP32 for
the large diffraction storage and accumulation path. Use
`mbs_po_gpu_float_fast` only after an application-specific FP64 comparison;
fast transcendental approximations can change narrow interference features.

A representative RTX 3080 Ti check (concave hexagonal particle, `k_eq=20`,
64 SO(3) samples, `N_phi=720`, `N_theta=360`) measured median wall times of
4.92 s for precise FP64, 3.45 s for precise mixed FP32, and 3.63 s for FP32
fast math. Mixed FP32 was 1.43x faster with global Mueller and M11 weighted-L2
errors near `1e-6` relative to FP64. This is a selection example, not a
universal tolerance guarantee; repeat it for the production particle and grid.
The CUDA diffraction path also enables warp-per-output beam reduction and its
3D-grid specialization by default. Set `MBS_GPU_WARP_BEAMS=0` or
`MBS_GPU_WARP_GRID_3D=0` only for diagnostic A/B runs.

For a double CUDA build without fast math:

```bash
make -C gpu fp64 -j
```

The CPU build defaults to architecture-specific optimization. It distinguishes
EPYC `9004` (Zen 4) from `9005` (Zen 5), uses AVX-512 when the compiler supports
the matching target, and otherwise falls back to `-march=native`. Such a binary
must run on a compatible CPU. Override `ARCH_FLAGS` for a portable build:

```bash
make -C cpu clean
make -C cpu ARCH_FLAGS='-march=x86-64-v3 -mtune=generic' -j
```

An explicit Zen 5 build, with a compiler that supports `znver5`, is:

```bash
make -C cpu clean
make -C cpu \
    ARCH_FLAGS='-march=znver5 -mtune=znver5 -mavx512f -mavx512dq -mavx512vl' -j
```

The GPU Makefile detects compute capability through `nvidia-smi`; if detection
is unavailable it defaults to `sm_86`. Set the architecture explicitly when
building elsewhere, for example `make -C gpu GPU_ARCH=89 -j`. Object
directories are separated by architecture, precision, fast-math, LTO, and MPI
settings. Run `make -C gpu clean` when changing the CUDA toolkit or host
compiler.

## Quick Start

CPU PO calculation:

```bash
cpu/bin/mbs_po_mpi --method po --backend cpu \
    --particle 1 10 10 --refractive-index 1.3116 0 \
    --wavelength-um 0.532 --max-reflections 8 --sobol 1024 \
    --scattering-grid 0 180 48 180 --threads 16 \
    --output results/po_column --close
```

For a fixed-orientation CPU PO calculation, `--threads N` distributes
independent azimuth rows of the scattering grid across OpenMP workers. Beam
contributions within each grid cell retain their original accumulation order,
so changing the thread count does not change the numerical result. Very small
grids remain serial to avoid OpenMP startup overhead.

The CPU PO diffraction loop also batches four regular theta directions in one
AVX2 edge-integral pass. Small-phase polygons and near-singular edge
denominators automatically retain the stable scalar formulas. This path is
enabled by the normal CPU build and has no runtime flag or approximation
tolerance. Temporary phase and sine/cosine buffers are stored vertex-major with
adjacent theta directions contiguous, so the four-lane edge pass uses direct
SIMD loads instead of scalar gathers. On the representative absorbing concave
case used for profiling, this layout reduced the median time from 12.76 s to
11.37 s (1.12x) relative to the preceding vectorized implementation, with
bit-identical Mueller output.

CUDA PO calculation:

```bash
gpu/bin/mbs_po_gpu_double --method po --backend cuda \
    --particle 1 125.9 78.09 --refractive-index 1.3116 0 \
    --wavelength-um 0.532 --max-reflections 12 \
    --orientation-diffraction-sampling 2 \
    --scattering-grid 0 180 600 180 --threads 16 \
    --output results/po_gpu --close
```

CUDA PO builds automatically enable the validated batched GPU tracing profile
for non-convex particles. It performs projected facet ordering and conservative
candidate filtering on the GPU while retaining exact polygon intersections and
beam splitting on the CPU. Use `--no-gpu-trace-prefilter` for a reference run
without this stage; `--gpu-trace-prefilter` explicitly requests it and is
mainly useful in reproducibility scripts. The initial per-worker batch is
derived from the OpenMP worker count; if CUDA reports stack-related OOM, the
process remembers a smaller successful limit for all subsequent orientations.
Each worker uses its own nonblocking CUDA stream. When trace-limit retries are
enabled, the process also remembers the largest successful retry beam limit,
so later orientations do not repeat a known-insufficient first pass; the
configured retry ceiling and all physical cutoffs remain unchanged. Heavy
trace batches are locked per physical GPU across independent processes, while
workers inside one process still overlap. This prevents CUDA stack-memory OOM
when two calculations happen to trace on the same card at once.

GO is a CPU method, including when invoked from a CUDA-capable binary:

```bash
cpu/bin/mbs_po_mpi --method go --backend cpu \
    --particle 1 100 70 --refractive-index 1.3116 0 \
    --wavelength-um 0.532 --max-reflections 12 \
    --fixed-orientation 37 19 --scattering-grid 0 180 360 180 \
    --output results/go_fixed --close
```

Every calculation requires exactly one method, one particle source, and one
primary orientation mode. Canonical hyphenated names are recommended; old
short and underscore names remain compatible aliases.

## Particle Input

Built-in particle syntax is checked before allocation:

| Type | Shape and parameters after `TYPE` |
|---:|---|
| `1` | hexagonal column/plate: `L D` |
| `2` | bullet: `L D` |
| `3` | bullet rosette: `L D [CAP]` |
| `4` | droxtal: `SCALE` |
| `10` | concave hexagonal: `L D CAVITY_DEG` |
| `12` | built-in two-column aggregate: `L D 2` |
| `999` | fixed built-in aggregate: `SCALE` |

The legacy type-4 form `4 UNUSED UNUSED SCALE` is accepted with a warning.
For type 10, the cavity angle must keep both cavities distinct and inside the
requested `L/D` geometry.

Use `--particle-file FILE` for general faceted geometry. Coordinates are
interpreted in the same length unit as the wavelength, normally micrometers.
They are parsed and retained in double precision. The loader translates them
to the order-independent bounding-box centre and snaps only nonphysical text
serialization noise on a scale-relative binary grid. This removes a physically
irrelevant global translation while preserving small edges and deterministic
facet ordering in files with large absolute coordinates.
[`examples/cube.particle`](examples/cube.particle) is a complete runnable
example. The strict text format is:

```text
CONCAVE_0_OR_1
AGGREGATE_0_OR_1 [FACETS_PER_COMPONENT]
BETA_SYMMETRY_DEG GAMMA_SYMMETRY_DEG

x y z
x y z
x y z

x y z
...
```

Blank lines separate facets; `#` starts a comment. A facet needs at least three
finite vertices in one simple convex boundary. Vertices must wind so every
facet normal points outward. The surface must be closed, consistently oriented,
and manifold. The current limits are 256 facets and 64 vertices per facet.
Aggregate facet count must be divisible by `FACETS_PER_COMPONENT`. Invalid
files stop with a line number and a `Fix:` instruction instead of continuing
with partial geometry.

`--resize-dmax-um` applies only to file particles. `--k-eq` scales either
source to

```text
k_eq = 2*pi*r_eq/lambda,   r_eq = (3*V/(4*pi))^(1/3).
```

Neither single-size option can be combined with a size scan.

## Orientation And Scattering Grids

Primary orientation modes include fixed and regular Euler grids, Monte Carlo,
orientation files, diffraction-derived grids, Sobol/Owen variants,
Hammersley, rank-1 lattices, Gauss-Euler quadratures, and adaptive modes. Run
`--help` for their exact arguments. PO implements every listed production
mode; GO is restricted to fixed, Euler-grid, Monte Carlo, diffraction-grid,
Sobol, and Sobol-seed modes. Unsupported combinations fail during preflight.

Uniform scattering grids have two forms:

```text
--scattering-grid T1 T2 NPHI NTH   # T1 < T2; writes NTH+1 theta rows
--scattering-grid R NPHI NTH       # cone of radius R around backscatter
```

`--theta-grid-file` reads distinct angles in `[0,180]`, one per line.
`--phi-points` overrides `NPHI`. In an orientation average, this azimuthal
quadrature is also the equivalent laboratory Euler-alpha average. Particle
axial symmetry reduces body gamma, not generally laboratory alpha.

`--so3-quaternion N` uses that identity explicitly: it traces `N` Hammersley
points only in the particle's symmetry-reduced beta/gamma domain, with uniform
Haar weights in `cos(beta)`, and evaluates the laboratory-alpha integral on the
scattering-azimuth grid. Use at least two phi points; for production pair it
with `--scattering-diffraction-sampling Q` and `--latitude-phi-grid`. The slower
`--so3-full-quaternion N` samples all three Euler degrees of freedom directly
and is retained as an independent audit mode.

`--so3-mirror-audit` provides an opt-in antithetic validation layout for an
even full-gamma `--so3-quaternion N` run: `N/2` Hammersley points cover the
first half of the gamma domain and their `N/2` reflected partners are traced
explicitly. It is directly comparable with `--so3-quaternion N/2
--mirror-gamma`, which traces the same base points and reconstructs the
partners. Without this audit modifier, the production full-domain Hammersley
sequence is unchanged.

For a particle with an actual reflection plane, `--mirror-gamma` also works
with `--so3-quaternion`: it samples `0 <= gamma < gamma_sym/2` and restores the
omitted half as `M(phi) -> P M(-phi) P`, `P=diag(1,1,-1,-1)`. `N` is always the
number of orientations that are really traced. Thus a mirror run with roughly
half the non-mirror `N` keeps the same gamma-point density. Do not use this
modifier for a chiral, asymmetric, or unverified file-loaded particle. The PO
tracer and its facet-visibility data must also be numerically mirror-invariant;
compare one production-resolution case with a full-gamma run before relying on
the reduction for a new particle class.
Use an even full-run count with `--so3-mirror-audit` and compare it with exactly
half that count under `--mirror-gamma`; this separates mirror-model error from
orientation-sampling error. Also compare the normal full-domain sequence with
an independent high-resolution result before production use. Inspect the
requested backscatter angle separately from a global L2 norm.

`--diffraction-sampling Q` uses one convention for both orientation and
scattering grids. Here `lmax` is the longest edge of any particle facet, not
the particle diameter. With `xi = 0.69 * wavelength / lmax`, their target
angular step is `xi/Q`: `Q=1` uses the diffraction scale and `Q=2` halves it. The
specialized `--orientation-diffraction-sampling Q` and
`--scattering-diffraction-sampling Q` override their respective parts. The
orientation-only form can be combined with an explicit `--scattering-grid`.

Add `--latitude-phi-grid` to use `N_phi(theta) proportional to sin(theta)`.
It keeps the equatorial spacing while avoiding redundant azimuth samples near
the poles. Automatic scattering sampling is supported for one size per
process. Legacy `--diffraction-grid`, `--diffraction-limit-grid`, and
`--ring-points` remain available through `--help-debug` for command-line
compatibility. Runs made before the `lmax` edge-length correction can have a
different automatically selected grid.

## Adaptive Convergence

Adaptive modes are PO-only:

| Mode | Parameters selected |
|---|---|
| `--adaptive-orientations EPS` | Sobol `N` only |
| `--auto-theta-grid EPS` | nonuniform theta points only |
| `--adaptive-phi EPS` / `--adaptive-alpha EPS` | `N_phi` / laboratory-alpha only |
| `--adaptive-reflections EPS` | reflection depth `n` only |
| `--adaptive-euler-grid EPS` | `N_beta`, `N_gamma`, then joint refinement |
| `--auto EPS` | `N_phi`, theta, Sobol `N`; fixed `n` |
| `--autofull EPS` | `n`, `N_phi`, theta, Sobol `N` |
| `--diffraction-autofull EPS` | `n`, `N_phi`, theta, `N_beta`, `N_gamma` |

The convergence metric is mode-specific:

- phi and Euler candidates compare relative `M11`, every other element after
  normalization by `M11`, and integrated `M11`;
- reflection-depth search adds outgoing-energy change;
- theta refinement evaluates actual Mueller rows at one-third and two-thirds
  of each interval against interpolation, and also checks integrated `M11`;
- Sobol orientation refinement checks incoming energy and all 16 Mueller
  elements on control-theta rows selected from the actual grid.

The unified modes tune in the order `n` when enabled, `N_phi`, theta, then
orientations, and repeat until the selected settings and pilot count stabilize.
There are no particle-specific halo angles in the theta refinement.
An explicit `--phi-points N` fixes `N_phi` in a unified mode, while
`--scattering-grid` or `--theta-grid-file` fixes the theta grid. Remove the
corresponding explicit option when that dimension must be tuned. Standalone
`--adaptive-phi`/`--auto-phi` still reject `--phi-points` because those commands
have no other purpose than selecting `N_phi`.

Use [`configs/adaptive.example.conf`](configs/adaptive.example.conf) to store
independent minima, maxima, and tolerances for every searched parameter:

```bash
cpu/bin/mbs_po_mpi --method po --backend cpu \
    --particle 1 100 70 --refractive-index 1.3116 0 \
    --wavelength-um 0.532 --autofull 0.02 \
    --adaptive-config configs/adaptive.example.conf \
    --output results/converged --close
```

The file uses strict `key = value` records with `#` or `;` comments. Explicit
CLI limits override file values; per-parameter file tolerances override the
mode-wide `EPS`. `--convergence-passes` applies to orientation, phi,
reflection, and Euler pass streaks; theta uses its interval criterion directly.

Phi, reflection, theta, Euler, and unified searches fail if their requested
accuracy is not reached before a hard limit. Standalone
`--adaptive-orientations` is intentionally exploratory: at its cap it warns
and writes the best available result. Unified modes use strict orientation
convergence and fail at the cap. Adaptive runs write a
`*_convergence.tsv` report with candidates, errors, worst theta/element, and
the final selection.

## MPI, Size Scans, And Multiple GPUs

The CPU binary can run directly as one MPI rank or under `mpirun`. Set
`--threads` per rank and avoid oversubscribing host cores:

```bash
mpirun -np 4 cpu/bin/mbs_po_mpi --method po --backend cpu \
    --particle 1 100 70 --refractive-index 1.3116 0 \
    --wavelength-um 0.532 --sobol 4096 \
    --scattering-grid 0 180 96 360 --threads 16 \
    --output results/mpi --close
```

Serial PO size-scan support is deliberately limited to paths that consume the
scan in their implementation:

| Scan | Supported serial PO orientation modes |
|---|---|
| `--dmax-grid` | Euler-grid, diffraction-grid, Sobol, lattice, lattice-generator, Euler-quadrature, diffraction-autofull |
| `--k-eq-grid` / `--k-eq-list` | Euler-grid, diffraction-grid, diffraction-autofull |

Every `--dmax-grid` value is an absolute physical maximum dimension and is
converted as `x = pi*Dmax/lambda`; it is not a scale factor relative to the
source particle. A k-eq scan instead uses the source shape volume to map each
requested equivalent-size parameter to a geometrically similar particle.

Use process-parallel scanning for other orientation modes. Parallel Dmax scans
require a file particle because each child applies `--resize-dmax-um`:

```bash
gpu/bin/mbs_po_gpu_float_fast --method po --backend cuda \
    --particle-file examples/cube.particle --refractive-index 1.3116 0 \
    --wavelength-um 0.532 --sobol 4096 \
    --scattering-grid 0 180 96 360 \
    --dmax-grid 5 100 20 --scan-jobs 4 --scan-threads 4 \
    --gpu-devices 0,1,2,3 --output results/dmax_scan --close
```

`--scan-jobs 0` selects jobs automatically. Each child receives one size or
size batch and one visible GPU in round-robin order. This is process-level
parallelism across scan sizes. Separately, supported coherent, averaged CUDA
paths automatically split an orientation batch across all visible GPUs in one
process. Select devices with `CUDA_VISIBLE_DEVICES`, disable that second level
with `MBS_GPU_MULTI=0`, or cap it with `MBS_GPU_MULTI_MAX=N`. Automatic
multi-device batching is disabled under a multi-rank MPI launch unless
`MBS_GPU_MULTI` is set explicitly. `--shared-k-eq-batches` can reuse tracing
for nearby k-eq values; control the batch span with `--k-eq-batch-ratio`.

For `--diffraction-autofull`/`--oldauto` runs with a variable azimuth grid,
set `MBS_GPU_GROUPS=1` to distribute independent theta work groups across the
visible GPUs. This keeps tracing and prepared orientations in one process and
creates only one diffraction workspace per GPU, instead of duplicating the
whole calculation. Select devices with `CUDA_VISIBLE_DEVICES`; set
`MBS_GPU_MULTI=N` or `MBS_GPU_MULTI_MAX=N` to cap the group workers. The group
scheduler is used only by direct CUDA diffraction without FFT interpolation or
theta-zone beam filtering; unsupported combinations retain the existing path.
Within each orientation chunk, every GPU packs and uploads the shared beam data
once and reuses it for subsequent theta groups assigned to that worker.

Do not launch `--scan-jobs` under an MPI job with more than one rank: both are
process-level schedulers, and the program rejects that combination. Use
`--threads` for a serial scan or `--scan-threads` for scan children, not both.

## Output And Errors

`--output PATH/PREFIX` creates the directory `PATH/PREFIX` and writes files
inside it using `PREFIX` as the result-name prefix. If omitted, a timestamped
directory under `results/` is created. If the requested directory already
exists, a new sibling with `(N)` appended is used, except when `--checkpoint`
is resuming that directory. Optional output placeholders use
`%[zero-based-index]OPTION_`, for example `%0p_`; malformed or missing-option
placeholders fail during preflight.

Common outputs are:

| File | Produced by | Contents |
|---|---|---|
| `<prefix>.dat` | PO and GO | angle, solid-angle weight, and all 16 Mueller elements |
| `<prefix>_log.txt` | PO | run configuration, timing, energy, and integral diagnostics |
| `<prefix>_back.dat`, `<prefix>_forward.dat` | GO | exact-pole peak values and normalized diagonal elements |
| `<prefix>_noshadow.dat` | PO with `--no-shadow-output` | Mueller matrix without the shadow/external beam |
| `<prefix>_jones.dat` | fixed PO with `--jones-output` | complex Jones matrices |
| `<prefix>_convergence.tsv` | adaptive modes | candidate values, convergence errors, and final selection |
| `FILE` and `FILE.visibility` from `--save-geometry FILE` | explicit request | the effective particle and per-facet visibility metadata needed for an exact round trip |
| `<label>.run.log` and child directories | process-parallel scans | one child log and result directory per size/batch |

`--jones-output` is only valid for fixed-orientation PO. `--no-shadow-output`
writes both full and no-shadow results, while `--no-shadow-beam` calculates
only without the shadow beam; they cannot be combined. Output files are plain
whitespace-separated text; the first row is the column header.

Normal runs do not write diagnostics in the current working directory. To dump
the first handled PO beam set for debugging, set `MBS_BEAM_LOG` to an explicit
file path; failure to open that optional path is reported as a warning.

Every active `MBS_*` environment variable is recorded in the result log.
Overrides that can change sampling, geometry, cutoffs, or numerical results are
rejected unless `--allow-experimental-environment` is passed explicitly. For a
production result, unset such overrides instead of acknowledging them.

Prepared-orientation chunks no longer use a fixed bytes-per-orientation guess.
The calculation starts with 16 orientations, measures nested beam/path
allocations, and grows the next chunk within current `MemAvailable`,
`MBS_HOST_MEM_FRACTION`, `MBS_HOST_MEM_RESERVE_MB`, and an optional
`MBS_HOST_MEM_BUDGET_MB`. Allocation failures include a corrective action.

Exit status is part of the CLI contract:

| Code | Meaning |
|---:|---|
| `0` | successful calculation or help |
| `2` | command/preflight error |
| `1` | runtime, calculation, CUDA, output, or child-process failure |

User-correctable failures contain `ERROR:` and at least one concrete `Fix:`.
CUDA runtime failure does not silently fall back to CPU unless the diagnostic
environment variable `MBS_GPU_ALLOW_FALLBACK=1` is set explicitly.

## Verification

```bash
# Complete release gate; CUDA is run when nvcc and a GPU are available
make test

# Strict project-warning gate
make test_warnings

# Precise FP64 CPU/CUDA numerical agreement gate on an NVIDIA host
make test_cuda

# CLI schema and validation unit tests
tests/run_cli_tests.sh

# Real CPU calculations plus pairwise selector and cross-domain checks
MBS=cpu/bin/mbs_po_mpi tests/run_release_cli_matrix.sh

# Fast adaptive-path tests
MBS=cpu/bin/mbs_po_mpi tests/run_adaptive_tests.sh

# Numerical regression suite
MBS=cpu/bin/mbs_po_mpi SKIP_BUILD=1 tests/run_tests.sh

# Separate AddressSanitizer + UndefinedBehaviorSanitizer release matrix
make test_sanitize
```

The release CLI matrix requires every option shown by `--help-debug` to appear
in at least one success or expected-error case. It exercises every primary PO
orientation mode, every implemented GO orientation mode, all built-in
particles, valid and malformed file geometry and trajectory input, all
pairwise primary-mode conflicts, serial-scan support, adaptive limits, output
failures, integer-overflow guards, and MPI/scan conflicts. Every expected
rejection must have a nonzero status, `ERROR:`, and an actionable `Fix:` line.
GitHub Actions builds the precise FP64 CUDA target. Runtime CUDA validation
still requires a self-hosted runner with the intended GPU and driver; the CPU
gate cannot validate kernels, device-memory limits, cuFFT, or compute
capability.

The matrix is schema-complete and covers pairwise selectors plus selected
cross-domain combinations; it is not an exhaustive enumeration of every
numeric value or arbitrary combination of all flags. New constraints must be
added both to preflight validation and to a success or expected-error case.

## Release Process

1. Run `make test`, `make test_sanitize`, `make docs`, and the hardware FP64
   CPU/CUDA agreement gate on a compatible NVIDIA system.
2. Confirm that the CPU workflow passes. The CUDA workflow compiles the precise
   FP64 target for every CUDA-related change; its runtime job is available by
   manual dispatch on a self-hosted GPU runner.
3. Build every advertised CUDA precision profile on each supported
   architecture and compare representative convex, non-convex,
   fixed-orientation, and orientation-averaged results.
4. Record the exact commit, `--version` output, compiler, MPI/CUDA versions,
   hardware, precision profile, commands, and tolerances in
   `docs/releases/vX.Y.Z.md`.
5. Create and push an annotated `vX.Y.Z` tag. The `Publish release` workflow
   creates the GitHub release from the matching notes file. GitHub supplies the
   source archives; this repository does not currently publish prebuilt binary
   assets.

## Versioning, Citation, And License

Record `--version` output with every numerical data set. Development builds
should be identified by the full Git commit. Published releases use annotated
semantic-version tags. Cite the MBS-1 method article above and include the
repository URL plus exact tag/commit used for the calculation.

The repository includes the GNU General Public License version 3 text; see
[`LICENSE`](LICENSE). Preserve source-level third-party notices when
redistributing modified versions.
Report reproducible defects through the
[GitHub issue tracker](https://github.com/KirillSalnikov/MBS-fast/issues) with
the full command, binary precision, compiler/CUDA versions, hardware, and the
generated log.
