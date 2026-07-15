# MBS-fast

MBS-fast calculates light scattering by faceted particles using geometrical ray
tracing with either Physical Optics (PO) diffraction or Geometrical Optics
(GO) output. It supports fixed orientations, several deterministic and
quasi-random orientation averages, adaptive convergence, CPU MPI/OpenMP, CUDA
diffraction, and process-parallel size scans across multiple GPUs.

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

The older Markdown manuals under `docs/` are retained as historical material;
they are not the release source for CLI behavior.

## Requirements

CPU build:

- Linux on x86-64;
- GNU Make, Bash, and GNU coreutils (`timeout` is used by the release tests);
- GCC 9 or newer, or Clang 14 or newer, with C++11 and OpenMP support;
- an MPI C++ wrapper and launcher from the same installation. The supplied
  release gate is tested with Open MPI and normally uses `/usr/bin/mpicxx` and
  `/usr/bin/mpirun` when they are available.

CUDA build additionally requires an NVIDIA driver, CUDA toolkit with `nvcc`,
and cuFFT. The build script selects a CUDA-compatible host compiler when one is
installed. CUDA and host-compiler compatibility is determined by the installed
toolkit, not by MBS-fast.

Building the manuals requires XeLaTeX, `fontspec`, `unicode-math`, the
Latin Modern/OpenType fonts used by the documents, and the standard LaTeX
packages imported by `docs/MANUAL*.tex`.

## Build

```bash
git clone https://github.com/KirillSalnikov/MBS-fast.git
cd MBS-fast

# CPU, MPI + OpenMP
make -C cpu -j

# CUDA variants
make -C gpu -j                 # float + --use_fast_math, the default
make -C gpu float -j           # float, without --use_fast_math
make -C gpu double_fast -j     # double GPU diffraction + --use_fast_math

# Rebuild both typeset manuals (two XeLaTeX passes per language)
make docs
```

The resulting binaries are independent:

| Build | Binary |
|---|---|
| CPU MPI/OpenMP | `cpu/bin/mbs_po_mpi` |
| CUDA float | `gpu/bin/mbs_po_gpu_float` |
| CUDA float fast | `gpu/bin/mbs_po_gpu_float_fast` |
| CUDA double fast | `gpu/bin/mbs_po_gpu_double_fast` |

The default CUDA target is float with `--use_fast_math`; it is the throughput
build, not the numerical reference build. Validate it against CPU double and a
CUDA double build on representative convex, non-convex, fixed-orientation, and
orientation-averaged cases before using it for a production data set. Record
the output of `--version` with the results.

For a double CUDA build without fast math:

```bash
make -C gpu GPU_PRECISION=double GPU_FAST_MATH=0 \
    TARGET=bin/mbs_po_gpu_double -j
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
directories are separated by architecture, precision, fast-math, and MPI
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

CUDA PO calculation:

```bash
gpu/bin/mbs_po_gpu_double_fast --method po --backend cuda \
    --particle 1 125.9 78.09 --refractive-index 1.3116 0 \
    --wavelength-um 0.532 --max-reflections 12 --diffraction-grid 2 \
    --scattering-grid 0 180 600 180 --threads 16 \
    --output results/po_gpu --close
```

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
finite, non-collinear vertices. The current limits are 256 facets and 63
vertices per facet. Aggregate facet count must be divisible by
`FACETS_PER_COMPONENT`. Invalid files stop with a line number and a `Fix:`
instruction instead of continuing with partial geometry.

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
| `particle_for_check.dat` | every calculation | the effective, resized particle written inside the result directory |
| `<label>.run.log` and child directories | process-parallel scans | one child log and result directory per size/batch |

`--jones-output` is only valid for fixed-orientation PO. `--no-shadow-output`
writes both full and no-shadow results, while `--no-shadow-beam` calculates
only without the shadow beam; they cannot be combined. Output files are plain
whitespace-separated text; the first row is the column header.

Normal runs do not write diagnostics in the current working directory. To dump
the first handled PO beam set for debugging, set `MBS_BEAM_LOG` to an explicit
file path; failure to open that optional path is reported as a warning.

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
# Complete CPU release gate: build, CLI matrix, adaptive tests, regressions
make test

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
CUDA builds still require verification on a host with the intended GPU and
driver; the CPU gate cannot validate CUDA kernels, device memory limits,
cuFFT, or the target compute capability.

The matrix is schema-complete and covers pairwise selectors plus selected
cross-domain combinations; it is not an exhaustive enumeration of every
numeric value or arbitrary combination of all flags. New constraints must be
added both to preflight validation and to a success or expected-error case.

## Release Checklist

A release is ready only after all of the following are true:

1. Start from a clean clone of the intended tag and verify that `git archive`
   contains the generated CLI sources, adaptive example, tests, manuals, and
   example particle. Do not release from a dirty worktree.
2. Run `make test`, `make test_sanitize`, a warning build with no project-source
   warnings, `make docs`, and the two-rank MPI cases with `mpicxx` and `mpirun`
   from the same MPI installation. Filter third-party MPI-header diagnostics
   separately instead of hiding project warnings.
3. Add a CI workflow that runs the portable CPU release gate from a clean clone
   and retains test logs. Keep the hardware CUDA gate as a separately recorded
   release job when no compatible CI runner is available.
4. Review the English and Russian TeX manuals for matching production CLI and
   physics coverage, rebuild both PDFs, and require zero `Overfull` or LaTeX
   errors. Page count alone is not a substitute for content parity.
5. On each supported CUDA architecture, build every advertised precision mode
   and run the CLI matrix plus numerical comparisons for a convex particle, a
   non-convex fixed orientation, and a non-convex orientation average. Check
   CPU double, CUDA double, and float-fast separately.
6. Record the tested compiler, MPI, CUDA, driver, GPU model/compute capability,
   CPU architecture, precision, fast-math state, commands, tolerances, and
   `--version` output in the release notes.
7. Resolve and document third-party source provenance, including the notices in
   `src/math/compl.cpp`, and choose the exact GPL
   form (`GPL-3.0-only` or `GPL-3.0-or-later`), copyright holders, and notices.
   Add project-specific citation metadata before claiming a citable release.
8. Add a changelog, citation metadata, security/contact policy, and third-party
   notices to the tagged source tree.
9. Publish an annotated semantic-version tag, source archive,
   binary artifacts where applicable, and SHA-256 checksums. Verify the
   archives by rebuilding and rerunning the release gate from extraction.

## Versioning, Citation, And License

Record `--version` output with every numerical data set. Development builds
should be identified by the full Git commit. A published
release should use an annotated semantic-version tag and include source and
binary checksums plus compiler, CUDA, GPU architecture, and precision metadata.
Until project-specific citation metadata or a DOI is published, cite the
repository URL and exact tag/commit used for the calculation.

The repository currently includes the GNU General Public License version 3;
see [`LICENSE`](LICENSE). The exact project licensing declaration and
third-party notices must be finalized as part of the release checklist above.
Report reproducible defects through the
[GitHub issue tracker](https://github.com/KirillSalnikov/MBS-fast/issues) with
the full command, binary precision, compiler/CUDA versions, hardware, and the
generated log.
