# RTX 5080 optimization checks — 2026-11-03

Machine: kika, RTX 5080 16 GB (`sm_120`), Ryzen 9 9950X3D,
8 OpenMP threads, one BLAS thread, WSL2, driver 621.01,
CUDA Toolkit 13.4.3 / nvcc 13.4.92, GCC 14.3.0.

## Validated application timings

Same ka=10 column, 2048 seeded orientations, 64 azimuths, all 181 theta
points and all 16 Mueller elements. One excluded warmup, five fresh-process
repeats, no completed-work reuse. All reference/completeness gates passed.

| Precision profile | Median, s | Min–max, s | MAD, s (relative) |
|---|---:|---:|---:|
| gpu-fp32 | 4.607 | 4.603–4.619 | 0.003 (0.07%) |
| gpu-fp32-reduced | 3.576 | 3.568–3.649 | 0.003 (0.09%) |
| gpu-fp64 | 8.608 | 8.604–9.143 | 0.005 (0.05%) |

The explicit reduced-trigonometry profile is **1.288× faster**
than precise FP32 storage on this case (22.4% less wall time).
This is a precision-profile comparison on the same equipment, not a hardware
speedup score. Precise FP32 and FP64 algorithms/defaults were retained.

Raw reports: [FP32](gpu-fp32.json), [range-reduced FP32](gpu-fp32-reduced.json),
[FP64](gpu-fp64.json); matching CSV files are alongside them. The reduced
profile uses the original frozen FP64 matrix as its regression reference.
The measured native revision is `7ee3eb1`; later tooling/result commits do not
change its numerical source. Exact binary, source and toolchain hashes are in JSON.

## Precision controls

[32 full-curve comparisons](shape_validation.json) cover column, plate,
bullet, bullet rosette, droxtal, concave column, two-component hexagonal
aggregate and irregular built-in aggregate, each at k_eq=0.01, 1, 20 and 100.
Each uses the same 64 seeded Haar orientations and 16 azimuths in both
precision modes. Worst difference over every angle and every matrix element,
normalized by M11 at that angle: **0.0000675%** relative to precise FP32.
The absorbing irregular-particle regression also passes for coherent and
incoherent accumulation, including cancellation-sensitive angular points.

Independent long-double/libm and CUDA phase tests verify range reduction and
the FP64 fallback. These are algorithm-equivalence controls; the small-k_eq
cases are not a claim that physical optics is physically converged there.
The orientation/azimuth counts here are fixed comparison inputs, not a
physical-convergence certificate. The reduced profile is opt-in.

```bash
python3 scripts/install_rtx50.py --cuda-version 13.4.3 --with-reduced-phase
python3 scripts/benchmark_hardware.py run --profile gpu-fp32-reduced --output runs/reduced
python3 scripts/validate_cuda_precision.py --output runs/precision_controls
```

## Launch configuration and remaining limits

[Precise-profile exploration](launch_tuning_precise.json) uses one warmup and
three measured repeats per variant, with default runs at both ends. On this
workload, 128/256-thread blocks, a flat thread layout, warp reduction, staged
Mueller accumulation and disabling compact-beam splitting did not improve
the default 64-thread 3D layout. Chunk 256 reduced kernel time by about 1.5%,
without enough evidence across workloads to change the default chunk.
[The portable tuner also ran on the reduced profile](launch_tuning_reduced.json);
128-thread blocks again provided no clear benefit. Every candidate's complete
matrix passed the numerical gate; raw slow repeats were retained.

[Kernel resource data](kernel_resources.json) are from the installed native
binaries. The selected coherent compact kernel uses 102 registers/thread in
precise FP32, 105 in reduced FP32 and 140 in FP64. With 64-thread blocks,
the CUDA driver occupancy API gives resource ceilings of 33.3% for both FP32
profiles and 25% for FP64 on this GPU. Static resource/driver
occupancy estimates are not measured achieved occupancy or bandwidth.
No Nsight hardware-counter trace was collected. Further kernel changes need
such profiling and validation on additional production workloads; this report
does not claim all possible RTX 5080 optimizations are complete.
