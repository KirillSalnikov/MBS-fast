# RTX 5080 baseline — 2026-11-03

Measured on **kika**, NVIDIA GeForce RTX 5080 (16 GB, native `sm_120`),
AMD Ryzen 9 9950X3D, **8 OpenMP threads**, one BLAS thread.
Windows/WSL2 Linux 6.6.114.1; driver 621.01; CUDA Toolkit 13.4.3
(nvcc 13.4.92), GCC 14.3.0. GPU clock/power settings were not changed.

One excluded warmup and **five fresh-process repeats** per case. The score is
external wall time including process startup, preparation and result output.
Every warmup/repeat passed completeness, numerical-reference and solver gates.
See the [fixed workload and protocol](../../README.md).

| Profile / case | Median, s | Min–max, s | MAD, s (relative) | Largest scaled matrix difference from reference |
|---|---:|---:|---:|---:|
| gpu-fp32 | 4.533 | 4.526–5.490 | 0.007 (0.16%) | 4.23709e-07% |
| gpu-fp64 | 8.649 | 8.607–9.206 | 0.042 (0.49%) | 3.72573e-14% |
| cpu | 15.470 | 15.407–15.911 | 0.063 (0.40%) | 1.74375e-13% |

The error column is the maximum across all 181 angles, all 16 matrix entries
and all six invocations, normalized by reference M11 at each angle. It is a
repeatability/regression measure, not an independent physical-accuracy estimate.
No slow repeats were discarded. Full samples, resource telemetry, commands,
compiler versions, binary/source/reference hashes and validation are in
[gpu-fp32.json](gpu-fp32.json), [CSV](gpu-fp32.csv); [gpu-fp64.json](gpu-fp64.json), [CSV](gpu-fp64.csv); [cpu.json](cpu.json), [CSV](cpu.csv).

## Native phase medians

These phases may overlap or be subsets of other phases; do not add them to
reconstruct total wall time.

| Profile / case | Native metric | Median, s |
|---|---|---:|
| gpu-fp32 | diffraction_and_controls_s | 4.1926 |
| gpu-fp32 | gpu_diffraction_including_packing_s | 4.1574 |
| gpu-fp32 | gpu_kernels_s | 4.0895 |
| gpu-fp32 | processing_s | 4.3487 |
| gpu-fp32 | trace_prepare_s | 0.1492 |
| gpu-fp64 | diffraction_and_controls_s | 8.1333 |
| gpu-fp64 | gpu_diffraction_including_packing_s | 8.0965 |
| gpu-fp64 | gpu_kernels_s | 7.9569 |
| gpu-fp64 | processing_s | 8.4448 |
| gpu-fp64 | trace_prepare_s | 0.3000 |
| cpu | diffraction_and_controls_s | 15.0823 |
| cpu | gpu_diffraction_including_packing_s | 0.0000 |
| cpu | processing_s | 15.2641 |
| cpu | trace_prepare_s | 0.1594 |

## Reproduce and compare

```bash
python3 scripts/benchmark_hardware.py run --profile gpu-fp32 --output runs/new_machine
python3 scripts/benchmark_hardware.py compare benchmarks/hardware_v1/results/rtx5080_20261103/gpu-fp32.json runs/new_machine/results.json
```

Keep the committed references, workload, profile and eight-thread setting.
A changed numerical-source or toolchain fingerprint is identified explicitly
as a software/toolchain-and-hardware comparison. MBS and BEM run different
models/workloads and these times must not be ranked against each other.

The measured native revision is `781f1e0`.
Result/README publication commits can follow this revision without changing
the measured numerical sources. Each JSON records the complete source hash.

The frozen FP32-storage curve differs from FP64 by at most 0.000118952% on the M11 scale (full-matrix relative L2 1.27094e-07); CPU versus FP64 GPU: 1.42748e-06% and L2 1.78945e-09.
