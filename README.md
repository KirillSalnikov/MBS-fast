# MBS-fast: Physical Optics for Ice Crystals

Fast Kirchhoff diffraction code for light scattering by non-spherical ice particles.
OpenMP + MPI parallel. AVX-512/AVX2 SIMD. Adaptive grids.

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
# Full auto: adaptive orientations + auto theta/phi grids
bin/mbs_po --po --auto 0.05 \
    -p 1 100 70 -w 0.532 --ri 1.31 0 -n 14 \
    --grid 0 180 48 181 --sym 6 2 --close

# Fixed Sobol orientations
bin/mbs_po --po --sobol 1024 --auto_tgrid --auto_phi \
    -p 1 100 70 -w 0.532 --ri 1.31 0 -n 14 \
    --grid 0 180 48 181 --sym 6 2 --close
```

## Multi-core (OpenMP)

```bash
OMP_NUM_THREADS=12 bin/mbs_po --po --auto 0.05 ...

# EPYC (64 cores, NUMA-aware)
OMP_PROC_BIND=close OMP_PLACES=cores OMP_NUM_THREADS=64 \
    bin/mbs_po_epyc --po --auto 0.05 ...
```

## Cluster (MPI + OpenMP)

```bash
# Build
sudo apt install libopenmpi-dev
bash build_mpi.sh

# Run on 8 nodes × 64 cores = 512 cores
mpirun -np 8 --hostfile nodes.txt \
    -x OMP_NUM_THREADS=64 -x OMP_PROC_BIND=close \
    bin/mbs_po_mpi --po --sobol 8192 --auto_tgrid --auto_phi \
    -p 1 1000 700 -w 0.532 --ri 1.31 0 -n 16 \
    --grid 0 180 48 181 --sym 6 2 --close

# SLURM cluster
sbatch cluster_mpi.sh
```

See `cluster_mpi.sh` for a ready-to-use SLURM script.

**nodes.txt** format:
```
node01 slots=1
node02 slots=1
...
```

Requirements: SSH without password between nodes, shared filesystem (NFS) or binary copied to each node.

## Output

Each run produces two Mueller matrix files:
- `M.dat` — full Mueller (with Babinet/shadow beam)
- `M_noshadow.dat` — without shadow (refracted beams only, better for halo studies)

## Documentation

- **docs/MANUAL.pdf** — full CLI reference (35+ flags), performance, build guides
- **docs/reports/** — bugfix reports, theta investigation
- **docs/figures/** — all plots (size scans, GOAD comparison, halos)
- **MANUAL.md** — manual (markdown)
- **tests/run_tests.sh** — regression tests

## Performance

| Setup | Time | Speedup |
|-------|------|---------|
| Original code, 1 thread | ~320 s | 1× |
| Optimized, 1 thread | 4.6 s | 70× |
| + OpenMP 12 threads | 1.7 s | 190× |
| + MPI 8 nodes × 64 cores | ~0.03 s* | ~10000×* |

*Estimated for 8192 Sobol on 512 cores. Benchmark: hex D=H=10um, n=6.
