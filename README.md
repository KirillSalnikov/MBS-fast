# MBS-fast: Physical Optics for Ice Crystals

Fast Kirchhoff diffraction code for light scattering by non-spherical ice particles.
OpenMP + MPI parallel. AVX-512/AVX2 SIMD. Adaptive grids. Auto beam cutoff.

## Requirements

- **GCC >= 9** (or Clang >= 14)
- **OpenMP** (libgomp for GCC, libomp for Clang)
- x86-64 CPU with AVX2+FMA (minimum) or AVX-512 (optimal)
- **MPI** (optional, for cluster runs): `sudo apt install libopenmpi-dev`

## Build

| Script | Target CPU | Binary |
|--------|-----------|--------|
| `bash build.sh` | Intel AVX-512 | `bin/mbs_po` |
| `bash build_epyc.sh` | AMD EPYC 7H12 (Zen 2) | `bin/mbs_po_epyc` |
| `bash build_zen4.sh` | AMD Zen 4 (Ryzen 7000 / Genoa) | `bin/mbs_po_zen4` |
| `bash build_epyc_clang.sh` | EPYC + Clang/AOCC | `bin/mbs_po_epyc_clang` |
| `bash build_mpi.sh` | MPI + OpenMP (cluster) | `bin/mbs_po_mpi` |
| `make` | Intel (Makefile) | `bin/mbs_po` |

## Quick Start

```bash
# Simplest: full auto (finds N_theta, N_phi, N_orient automatically)
bin/mbs_po --po --auto 0.05 \
    -p 1 100 70 -w 0.532 --ri 1.31 0 -n 8 --close

# Full auto including n search (--autofull)
bin/mbs_po --po --autofull 0.05 \
    -p 1 100 70 -w 0.532 --ri 1.31 0 --close

# Manual Sobol orientations
bin/mbs_po --po --sobol 1024 --auto_tgrid --auto_phi \
    -p 1 100 70 -w 0.532 --ri 1.31 0 -n 8 --close
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

Auto beam cutoff: eps³ × total beam energy. Skips 70-90% of weak beams.

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

## Documentation

- **docs/MANUAL.pdf** — full CLI reference (35+ flags), build guides, output format
- **docs/reports/** — bugfix reports (forward-direction Fresnel sign)
- **docs/figures/** — convergence plots, comparisons
- **MANUAL.md** — manual (markdown)
- **tests/run_tests.sh** — regression tests

## Key Features

- **Batched sincos**: all vertex phases pre-computed via AVX-512/AVX2 before theta loop
- **Auto theta grid**: 7 zones (forward peak, halos 22°/46°, side, backscatter)
- **Auto phi**: N_phi = x/5 + 48 (linear in size parameter, validated by --autofull)
- **Incremental adaptive**: reuses previous orientations (Sobol subset property)
- **5 convergence controls**: Q_sca, M₁₁(22°), M₁₁(46°), M₁₁(90°), M₁₁(180°)
- **Auto beam cutoff**: eps³ × total energy, skips negligible beams
- **Dual output**: M.dat (with shadow) + M_noshadow.dat (without) at no extra cost
- **Memory-aware chunking**: auto-sizes orientation batches to fit in RAM
- **MPI + OpenMP hybrid**: distributed across nodes, threaded within node
