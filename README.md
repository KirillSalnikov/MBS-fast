# MBS-fast

Light scattering by faceted particles: **physical optics (PO)** and **geometrical
optics (GO)**, with CPU/OpenMP/MPI and CUDA implementations.

For a new orientation-averaged calculation, start with **`--fullauto EPS`**.
Give the particle, refractive index, wavelength and required accuracy. The native
controller selects its numerical work, saves intermediate matrices and resumes
completed work after interruption. Python is not required for the calculation.

**`--fullauto 0.02` means 2%, not `--fullauto 2`.** This is an estimated numerical
accuracy within the selected optical model. It is not the error of PO or GO
against Maxwell's equations, and not a simultaneous 95% guarantee over a curve.

[Fullauto guide](docs/FULLAUTO.md) · [Подробная методичка](docs/FULLAUTO_RU.md) ·
[Русская HTML-методичка](docs/FULLAUTO_ADAPTIVE_RU.html) ·
[Краткие примеры HTML](docs/MINI_FULLAUTO_RU.html) ·
[Complete reference](docs/MANUAL.md) · [Справочник](docs/MANUAL_RU.md)

## Build once

Linux, GNU Make, a C++11 compiler and OpenMP are required. The CPU target also
uses Open MPI. CUDA builds require a compatible NVIDIA toolkit and driver.

```bash
# CPU: MPI-capable executable; run fullauto as one controller process.
make -C cpu -j2
cpu/bin/mbs_po_mpi --version

# NVIDIA GPU: reference FP64 build; set the actual device architecture.
make -C gpu fp64 GPU_ARCH=70 -j2
# 70: V100. 86: RTX 3080 Ti. 89: RTX 4090.
gpu/bin/mbs_po_gpu_double --version
```

For a consumer GPU, a mixed FP32 build is available. Validate its numerical
profile on representative inputs before a production campaign; orientation
statistics alone do not measure arithmetic error:

```bash
make -C gpu fp32 GPU_ARCH=86 -j2
```

See [build profiles](gpu/README.md) for precise, consumer and reduced phase
variants. More compiler flags are not needed to start fullauto.

## First calculation: one command

A hexagonal column with axial length 5 µm and basal diameter 2 µm, wavelength
0.532 µm, real refractive index 1.3116 and a 2% target:

```bash
cpu/bin/mbs_po_mpi --method po --backend cpu \
  --particle 1 5 2 --ri 1.3116 0 --wavelength-um 0.532 \
  --fullauto 0.02 --output results/column_5um
```

`--particle 1 L D` describes a regular hexagonal column. These are its axial
length and basal diameter, not an imposed maximum vertex-to-vertex size.
Use `--particle-file` and `--resize-dmax-um` to scale an imported particle by Dmax:

```bash
cpu/bin/mbs_po_mpi --method po --backend cpu \
  --particle-file examples/cube.particle --resize-dmax-um 5 \
  --ri 1.3116 0 --wavelength-um 0.532 \
  --fullauto 0.02 --output results/file_particle_5um
```

All geometric lengths must use the same units as the wavelength. Imported
symmetry declarations are not trusted automatically: fullauto starts from the
full orientation domain and verifies supported reductions.

## GPU and forward/backward cones

```bash
# Forward cone, one visible GPU.
CUDA_VISIBLE_DEVICES=0 gpu/bin/mbs_po_gpu_double \
  --method po --backend cuda --particle 1 100 50 \
  --ri 1.3116 0 --wavelength-um 0.532 \
  --fullauto 0.02 --fullauto-theta-range 0 25 \
  --output results/column_100um_forward

# Backward cone, two GPUs sharing a controller and saved work.
CUDA_VISIBLE_DEVICES=0,1 gpu/bin/mbs_po_gpu_double \
  --method po --backend cuda --particle 1 100 50 \
  --ri 1.3116 0 --wavelength-um 0.532 \
  --fullauto 0.02 --fullauto-theta-range 170 180 \
  --fullauto-gpus 0,1 --output results/column_100um_backward
```

GPU numbers are relative to `CUDA_VISIBLE_DEVICES`. Each output directory must
have one controller. Do not start two controllers on the same checkpoint.
A cone calculation does not supply full-sphere scattering integrals.

## What fullauto chooses

| Quantity | Default behavior | To fix or limit it |
|---|---|---|
| Orientation samples | Eight Owen-scrambled Sobol replicas; size-dependent start and local refinement | `--fullauto-initial N` fixes the start; `--max-orientations N` limits growth |
| Internal reflection depth | Paired search starting at 6, up to 18 | `--max-reflections N` fixes the physical target |
| Azimuth samples | Selects the cost/variance tradeoff; Haar alpha is sampled | `--phi-points N` fixes phi |
| Scattering angles | Builds support and check angles, interpolates the requested output grid | `--fullauto-theta-range A B` or `--theta-grid-file FILE` |
| Interpolation | Target `EPS/4`; separate matched-sample residual | A 2% run therefore has a 0.5% interpolation allowance |
| Beam pruning, coherent PO | Balanced estimator with a coupled uncut correction | `--cutoff-profile off` disables pruning; an explicit preset changes this automatic policy |
| Analytic controls | Trains and validates variance-reduction components; rejects unsuitable controls | `--fullauto-analytic-controls off` |
| Symmetry and resources | Verifies supported shape symmetries; tunes execution | Explicit flags take precedence; see the guide |

An explicit starting count is **not** a fixed final orientation count. An
explicit reflection depth **is** fixed: the default depth search is disabled.
Final checks may report that this fixed depth is insufficient; they do not
silently replace it. Reaching a resource cap is not convergence.

Sobol is the current production default and the practical starting choice for
our completed campaigns. Its advantage depends on shape, size and angles;
there is no claim that it beats every quadrature on every particle.

## Read the result and resume

| File | Meaning |
|---|---|
| `fullauto_status.json` | Controller state, selected parameters and last published accuracy |
| `intermediate_current.json` | Latest published matrix snapshot, including unfinished calculations |
| `mueller_fullauto.dat` | Full requested-angle matrix; inspect status before treating it as final |
| `mueller_evaluated.dat` | Directly evaluated angles |
| `mueller_fullauto_annotated.csv` | Matrix with direct/interpolated provenance and error estimates |
| `execution_progress.json` | Current worker activity; activity does not constitute a new accuracy estimate |
| `fullauto_history.jsonl` | Parameter choices, refinements and execution events |
| `controller_command.json`, `jobs/`, `adaptive_checkpoint.tsv` | Saved recipe and reusable work; retain the complete output directory |

To resume, run **the same command with the same binary and output directory**.
Completed blocks are reused automatically. An interrupted unfinished block may
be repeated. Rebuilding a binary changes its identity; cross-version recovery
requires an explicitly tested compatibility receipt, not merely the same CLI.
See [recovery](docs/FULLAUTO.md#results-and-recovery) for services and compatible upgrades.

## Geometrical optics

```bash
cpu/bin/mbs_po_mpi --method go --backend cpu \
  --particle 1 100 50 --ri 1.3116 0 --wavelength-um 0.532 \
  --fullauto 0.02 --output results/column_go
```

GO integrates alpha analytically and estimates finite solid-angle bins. Its
forward/backward delta contributions are separate integrated matrices. It does
not use PO's cubic theta interpolation or analytic facet controls. Optional
[GO mirror reduction](docs/GO_MIRROR.md) and
[audited balanced cutoffs](docs/GO_BALANCED.md) have explicit applicability
checks. They are not automatically interchangeable with PO settings.

## Other workflows

- [Fixed orientations, legacy grids, material/particle options](docs/MANUAL.md).
- [CPU extinction only](docs/EXTINCTION_ONLY.md): compute forward optical-theorem
  extinction without an angular Mueller grid. This is not a full albedo solver.
- [Multisize phase cache](docs/DIFFRACTION_OPT_20260907.md): reuse direction-independent
  beam phases in a direct FP64 size scan. Fullauto runs one size per controller;
  do not combine it with legacy size-scan flags.
- Generated `--help` and `--help-debug` list accepted arguments. The fullauto
  guides explain which options are automatic, fixed, unsupported or diagnostic.

## Verification and release

```bash
make test
make test_sanitize
make docs
# On an NVIDIA host, also run the hardware gate:
make test_cuda
```

Rebuilding the manuals requires Pandoc, XeLaTeX, fontconfig and the CMU fonts.
`make docs` regenerates the fullauto HTML/TeX from the maintained Markdown
guides before compiling the English and Russian PDFs.

Fullauto, GO mirror/balanced, extinction-only and multisize integration checks
are described in the [release preparation notes](docs/releases/v2.0.0-preparation.md).
Record the exact commit, binary `--version`, hardware, commands and tolerances.
A successful build or text merge is not a numerical release validation.

The publish workflow requires an annotated version tag and matching release
notes. Development `main` is identified by its commit; a prepared 2.0 codebase
is not a published `v2.0.0` until the release checks and tag are complete.

## Method, provenance and license

MBS-fast develops the original [MBS-raw](https://github.com/Heart-Under-Blade/MBS-raw)
implementation. For MBS-1 cite D. N. Timofeev, A. V. Konoshonkin and N. V. Kustova,
*Modified Beam-Splitting 1 (MBS-1) Algorithm for Solving the Problem of Light
Scattering by Nonconvex Atmospheric Ice Particles*, Atmospheric and Oceanic
Optics 31(6), 642–649 (2018), [DOI](https://doi.org/10.1134/S1024856018060179),
and the exact MBS-fast tag or commit. [GPLv3 license](LICENSE).
