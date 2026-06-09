# MBS-fast

MBS-fast computes light scattering by non-spherical ice particles with geometrical ray tracing plus Physical Optics diffraction. It supports CPU OpenMP/MPI builds and CUDA GPU diffraction builds.

The current full manual is [`MANUAL.md`](MANUAL.md). It covers:

| Topic | Where |
|---|---|
| CPU/GPU compilation, EPYC Zen targets, AVX2/AVX-512 | [`MANUAL.md#build`](MANUAL.md#build) |
| CUDA float/double binaries | [`MANUAL.md#gpu-build`](MANUAL.md#gpu-build) |
| GPU speedup, atomics, fused Mueller kernels | [`MANUAL.md#why-the-gpu-version-is-faster`](MANUAL.md#why-the-gpu-version-is-faster) |
| Diffraction, Jones dyads, Mueller conversion | [`MANUAL.md#physical-model`](MANUAL.md#physical-model) |
| Grids and orientation modes | [`MANUAL.md#orientation-grids`](MANUAL.md#orientation-grids), [`MANUAL.md#scattering-grids`](MANUAL.md#scattering-grids) |
| Multi-GPU and multi-size scans | [`MANUAL.md#gpu-and-multi-gpu`](MANUAL.md#gpu-and-multi-gpu) |
| Every parsed CLI flag | [`MANUAL.md#all-flags`](MANUAL.md#all-flags) |

## Quick Build

```bash
# CPU MPI/OpenMP binary
make -C cpu -j

# CUDA GPU binaries
make -C gpu float_fast  -j
make -C gpu double_fast -j
```

## Quick Run

```bash
gpu/bin/mbs_po_gpu_double_fast --po -p 1 125.9 78.09 \
    --ri 1.3116 0 -w 0.532 -n 12 --oldauto 2 \
    --grid 0 180 600 180 --threads 64 --close -o out_gpu
```

Release binaries show a short help:

```bash
gpu/bin/mbs_po_gpu_double_fast --help
```

Full debug/experimental flag help is available through:

```bash
gpu/bin/mbs_po_gpu_double_fast --help-debug
```
