# MBS-fast

MBS-fast computes light scattering by non-spherical ice particles with geometrical ray tracing plus Physical Optics diffraction. It supports CPU OpenMP/MPI builds and CUDA GPU diffraction builds.

The current full manual is available as Markdown, LaTeX, and PDF:

| Language | Markdown | LaTeX | PDF |
|---|---|---|---|
| English | [`docs/MANUAL.md`](docs/MANUAL.md) | [`docs/MANUAL.tex`](docs/MANUAL.tex) | [`docs/MANUAL.pdf`](docs/MANUAL.pdf) |
| Russian | [`docs/MANUAL_RU.md`](docs/MANUAL_RU.md) | [`docs/MANUAL_RU.tex`](docs/MANUAL_RU.tex) | [`docs/MANUAL_RU.pdf`](docs/MANUAL_RU.pdf) |

It covers:

| Topic | Where |
|---|---|
| CPU/GPU compilation, EPYC Zen targets, AVX2/AVX-512 | [`docs/MANUAL.md#build`](docs/MANUAL.md#build) |
| CUDA float/double binaries | [`docs/MANUAL.md#gpu-build`](docs/MANUAL.md#gpu-build) |
| GPU speedup, atomics, fused Mueller kernels | [`docs/MANUAL.md#why-the-gpu-version-is-faster`](docs/MANUAL.md#why-the-gpu-version-is-faster) |
| Diffraction, Jones dyads, Mueller conversion | [`docs/MANUAL.md#physical-model`](docs/MANUAL.md#physical-model) |
| Grids and orientation modes | [`docs/MANUAL.md#orientation-grids`](docs/MANUAL.md#orientation-grids), [`docs/MANUAL.md#scattering-grids`](docs/MANUAL.md#scattering-grids) |
| Multi-GPU and multi-size scans | [`docs/MANUAL.md#gpu-and-multi-gpu`](docs/MANUAL.md#gpu-and-multi-gpu) |
| Every parsed CLI flag | [`docs/MANUAL.md#all-flags`](docs/MANUAL.md#all-flags) |

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
