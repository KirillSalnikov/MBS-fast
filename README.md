# MBS-fast: Physical Optics for Ice Crystals

Fast Kirchhoff diffraction code for light scattering by non-spherical ice particles.
OpenMP + MPI parallel. AVX-512/AVX2 SIMD. Adaptive grids. Auto beam cutoff.

## Requirements

- **GCC >= 9** (or Clang >= 14)
- **OpenMP** (libgomp for GCC, libomp for Clang)
- x86-64 CPU with AVX2+FMA (minimum) or AVX-512 (optimal)
- **MPI** (optional, for cluster runs): `sudo apt install libopenmpi-dev`
- **CUDA toolkit** (optional, for `--gpu`/`--fft`): NVIDIA GPU with CUDA-capable driver

## Build

Preferred split builds:

```bash
# CPU cluster binary: MPI between processes, OpenMP inside each process.
make -C cpu -j
mpirun -np 4 cpu/bin/mbs_po_mpi --po --sobol 4096 ...

# CUDA binary. GPU backend is enabled by default in this binary.
PATH=/usr/local/cuda/bin:$PATH \
make -C gpu -j
gpu/bin/mbs_po_gpu_float_fast --po --fft --autofull 0.05 ...
```

The split builds keep object files out of `src`: CPU objects go to
`cpu/build/`, CUDA objects go to `gpu/build/`. The shared physical model and
CLI implementation still live in `src/`, so the CPU and GPU binaries do not
drift into two separate codebases. In the GPU split build `--gpu` is optional;
use `--cpu` only to force the CPU backend from the GPU-capable binary.

The GPU split build auto-detects the CUDA compute capability with
`nvidia-smi` and falls back to `GPU_ARCH=86` when it cannot query a device.
Override explicitly when needed, for example:

```bash
make -C gpu -j GPU_ARCH=86
```

| Script | Target CPU | Binary |
|--------|-----------|--------|
| `make` | Auto CPU flags (EPYC 7H12 -> Zen 2, otherwise native) | `bin/mbs_po` |
| `bash build.sh` | Auto CPU flags (EPYC 7H12 -> Zen 2, otherwise native) | `bin/mbs_po` |
| `bash build_epyc.sh` | AMD EPYC 7H12 (Zen 2) | `bin/mbs_po_epyc` |
| `bash build_zen4.sh` | AMD Zen 4 (Ryzen 7000 / Genoa) | `bin/mbs_po_zen4` |
| `bash build_epyc_clang.sh` | EPYC + Clang/AOCC | `bin/mbs_po_epyc_clang` |
| `bash build_mpi.sh` | MPI + OpenMP (cluster) | `bin/mbs_po_mpi` |

Default `make`/`build.sh` intentionally avoid hard-coded AVX-512. On AMD EPYC
7H12 they use `-march=znver2 -mtune=znver2`; elsewhere they use
`-march=native -mtune=native`. Override with `ARCH_FLAGS=...` or
`CXXFLAGS=...` if you need a portable binary.

CUDA build:

```bash
make USE_CUDA=1

# For RTX 3080 Ti / 4070 SUPER, avoid unsupported PTX JIT on older drivers:
make USE_CUDA=1 NVCCFLAGS="-O3 -std=c++11 -arch=sm_86 -U_GNU_SOURCE -D_DEFAULT_SOURCE -D_XOPEN_SOURCE=700"

# Optional separate CUDA precision/speed binaries:
make cuda_float       # bin/mbs_po_float
make cuda_float_fast  # bin/mbs_po_float_fast
make cuda_double_fast # bin/mbs_po_double_fast
```

CUDA precision builds:

| Command | Binary | CUDA diffraction precision | Use case |
|---------|--------|----------------------------|----------|
| `make USE_CUDA=1` | `bin/mbs_po` | double | default CUDA build, safest |
| `make cuda_float` | `bin/mbs_po_float` | float / FP32 | test FP32 without fast math |
| `make cuda_float_fast` | `bin/mbs_po_float_fast` | float / FP32 + fast math | production fast GPU run on RTX cards |
| `make cuda_double_fast` | `bin/mbs_po_double_fast` | double + fast math | reference/fallback GPU run |

For RTX 3080 Ti / 4070 SUPER build the release binaries explicitly for
Ampere-compatible `sm_86`:

```bash
PATH=/usr/local/cuda/bin:$PATH \
make cuda_float_fast -j1 \
    NVCCFLAGS="-O3 -std=c++11 -arch=sm_86 -U_GNU_SOURCE -D_DEFAULT_SOURCE -D_XOPEN_SOURCE=700"

PATH=/usr/local/cuda/bin:$PATH \
make cuda_double_fast -j1 \
    NVCCFLAGS="-O3 -std=c++11 -arch=sm_86 -U_GNU_SOURCE -D_DEFAULT_SOURCE -D_XOPEN_SOURCE=700"
```

Legacy root CUDA builds still need `--gpu` to use CUDA and `--fft` to use the
cuFFT phi-interpolation backend. The new `gpu/` split build enables CUDA by
default:

```bash
bin/mbs_po_float_fast --po --pf particle.dat --k_eq 50 \
    --ri 1.6 0.002 -w 1.064 -n 14 \
    --oldauto 2 --grid 0 180 600 181 --gpu --fft --close
```

CUDA runtime errors are fatal by default in the GPU split build. This prevents a
bad GPU binary from silently continuing on the CPU and producing misleading
timings. For debugging only, restore the old fallback behavior with:

```bash
MBS_GPU_ALLOW_FALLBACK=1 gpu/bin/mbs_po_gpu_float_fast ...
```

`--fft` first computes diffraction on a reduced uniform `phi` grid and then
reconstructs the requested `Nphi` grid by Fourier zero-padding. For safer runs
use global full-`Nphi` fallback on suspicious theta rows after the FFT pass:

```bash
bin/mbs_po_float_fast --po --pf particle.dat --k_eq 50 \
    --ri 1.6 0.002 -w 1.064 -n 14 \
    --oldauto 2 --grid 0 180 600 181 --gpu --fft \
    --fft_auto_phi 1e-3 --fft_auto_phi_max_rows 16 --close
```

`--fft_auto_phi EPS` keeps the FFT acceleration for the main pass and then
recomputes only suspicious final theta rows directly at the full requested
`Nphi`. This global post-pass is much cheaper than refining inside every GPU
orientation chunk. The older per-chunk detector remains available for debugging
with `MBS_FFT_ADAPTIVE_PHI=1`.

Production recommendation: use `bin/mbs_po_float_fast` for speed. If a point
looks suspicious, rerun that case with `bin/mbs_po_double_fast`. The FP32 CUDA
edge-integral tolerance can be overridden with `MBS_GPU_EPS1=...`; the default
FP32 minimum is `1e-5`.

The default optimized PO output writes only the full Mueller matrix. Add
`--noshadow_output` to also write the legacy `_noshadow` matrix.

## Quick Start

```bash
# Simplest: full auto (finds N_theta, N_phi, N_orient automatically)
bin/mbs_po --po --auto 0.05 \
    -p 1 100 70 -w 0.532 --ri 1.31 0 -n 8 --close

# Full auto including n search (--autofull)
bin/mbs_po --po --autofull 0.05 \
    -p 1 100 70 -w 0.532 --ri 1.31 0 --close

# Manual Sobol with adaptive theta grid (5% tolerance)
bin/mbs_po --po --sobol 1024 --auto_tgrid 0.05 --auto_phi \
    -p 1 100 70 -w 0.532 --ri 1.31 0 -n 8 --close

# File particle scaled by equivalent size parameter k_eq
bin/mbs_po --po --oldauto 2 --pole \
    --pf shapeA64_mbs.dat --k_eq 20.76 \
    --ri 1.6 0.002 -w 1.064 -n 8 \
    --tgrid scattering_angles --nphi 600 --gpu --fft --close

# Multi-size: trace once at max k_eq, then diffract each listed size
MBS_SHARED_BETA_GROUP=4 MBS_GPU_NO_ATOMICS=1 \
bin/mbs_po --po --oldauto 2 --pole \
    --pf shapeA64_mbs.dat --multikeq_list keq_list.txt \
    --ri 1.6 0.002 -w 1.064 -n 8 \
    --tgrid scattering_angles --nphi 600 --gpu --fft --close
```

**Note**: `--grid` and `--sym` are NOT needed with `--auto`/`--autofull`. Theta grid, phi bins, and particle symmetry are set automatically.

## Choosing n (reflections)

| n | Use case |
|---|----------|
| 6 | Fast estimate, forward scattering |
| 8 | General purpose (Q_sca < 0.1% error) |
| 12-14 | Backscattering, depolarization |
| 16-20 | Reference calculations |
| omit | `--autofull` finds n automatically |

## --auto vs --autofull

| Mode | What it automates | Speed |
|------|-------------------|-------|
| `--auto EPS` | N_theta, N_phi, N_orient, beam cutoff | Fast |
| `--autofull EPS` | Same + n (reflections) | 3-4× slower |

Both converge on 5 control points: Q_sca, M₁₁(22°), M₁₁(46°), M₁₁(90°), M₁₁(180°).
All must be within EPS for 2 consecutive iterations.

Auto beam cutoff: skip beams by dimensionless relative thresholds. The modern
fine-grained flags are `--beam_cutoff_j`, `--beam_cutoff_area`,
`--beam_cutoff_importance`, `--trace_cutoff_j`, `--trace_cutoff_area`, and
`--trace_cutoff_importance`.
`--beam_cutoff EPS` remains a shorthand for the J/area beam tests.

`--beam_cutoff EPS` can also be used standalone with `--random` or `--sobol`.

## Multi-core (OpenMP)

```bash
OMP_NUM_THREADS=12 bin/mbs_po --po --auto 0.05 \
    -p 1 100 70 -w 0.532 --ri 1.31 0 -n 8 --close

# EPYC (64 cores, NUMA-aware)
OMP_PROC_BIND=close OMP_PLACES=cores OMP_NUM_THREADS=64 \
    bin/mbs_po_epyc --po --auto 0.05 \
    -p 1 100 70 -w 0.532 --ri 1.31 0 -n 8 --close
```

## Cluster (MPI + OpenMP)

```bash
# Build
sudo apt install libopenmpi-dev
bash build_mpi.sh

# Run on 8 nodes × 64 cores
mpirun -np 8 --hostfile nodes.txt \
    -x OMP_NUM_THREADS=64 -x OMP_PROC_BIND=close \
    bin/mbs_po_mpi --po --sobol 8192 --auto_tgrid --auto_phi \
    -p 1 1000 700 -w 0.532 --ri 1.31 0 -n 16 --close

# SLURM
sbatch cluster_mpi.sh
```

## Output

Each run produces two Mueller matrix files:
- `M.dat` — full Mueller (with Babinet/shadow beam diffraction)
- `M_noshadow.dat` — without shadow beam (refracted beams only, better for halo studies)

Format: `theta 2pi*dcos M11 M12 ... M44` (18 columns, phi-averaged).

## Validation

Old binary vs MBS-fast verified to **machine precision** (1e-7) with `--coh_orient`.
Default incoherent mode: **0.007% max M11 diff** at 2860 orientations.
Speedup: **55×** (23 min → 25 sec on 4 cores).

See `tests/reference_test/RESULTS.md` and comparison plots in `tests/reference_test/`.

## Documentation

- **docs/MANUAL.pdf** — full CLI reference (35+ flags), build guides, output format
- **docs/bugfix_15x_mismatch.pdf** — investigation of 15× mismatch (d_param init, Rotate fix)
- **docs/figures/** — convergence plots, comparisons
- **tests/reference_test/** — reference test scripts, comparison plots, results
- **tests/run_tests.sh** — regression tests

## Key Features

- **Batched sincos**: all vertex phases pre-computed via AVX-512/AVX2 before theta loop
- **Adaptive theta grid**: recursive bisection — adds points where M₁₁ changes fast (2-3× fewer points than uniform)
- **Auto phi**: N_phi = x/5 + 48 (linear in size parameter, validated by --autofull)
- **Incremental adaptive**: reuses previous orientations (Sobol subset property)
- **5 convergence controls**: Q_sca, M₁₁(22°), M₁₁(46°), M₁₁(90°), M₁₁(180°)
- **Auto beam cutoff**: two dimensionless thresholds (|J|²/max < eps AND area/max < eps). Protects forward peak and strong beams, skips 45-80% of negligible beams
- **CUDA diffraction**: `--gpu`; add `--fft` for cuFFT phi interpolation backend
- **Multi-size**: `--multigrid`, `--multikeq`, `--multikeq_list` — trace once at the largest size, recompute diffraction for all sizes
- **GPU batching controls**: `MBS_SHARED_BETA_GROUP=N`, `MBS_GPU_NO_ATOMICS=1`, `MBS_GPU_MEM_FRACTION=...`
- **Per-beta save**: `--save_betas` — write intermediate Mueller per beta (backup/resume)
- **Opt-in checkpoint**: `--checkpoint` — save/resume long `--orientfile` runs; off by default to avoid extra I/O
- **Dual output**: M.dat (with shadow) + M_noshadow.dat (without) at no extra cost
- **Memory-aware chunking**: auto-sizes orientation batches to fit in RAM
- **MPI + OpenMP hybrid**: distributed across nodes, threaded within node
