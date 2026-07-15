# CUDA GPU build

This directory builds the CUDA binary. The GPU backend is enabled by default in
this build, so `--backend cuda` is not required. Use `--backend cpu` only when
you deliberately want the CPU backend from the GPU-capable binary. The legacy
`--gpu` and `--cpu` switches remain accepted.

```bash
make -C gpu -j
gpu/bin/mbs_po_gpu_float_fast --method po --backend cuda \
    --particle 1 10 10 --refractive-index 1.3116 0 \
    --wavelength-um 10 --max-reflections 2 --sobol 256 \
    --scattering-grid 0 180 48 180 --threads 4 \
    --output results/gpu_example --close
```

Adaptive pilot diffraction uses the same CUDA Mueller kernel as the final
calculation. Parameter candidates are compared on the full Mueller matrix; the
GPU does not use a separate lower-accuracy convergence criterion.
`N_phi` is the laboratory-alpha quadrature implemented through scattering
azimuth. Body-axis symmetry is already applied to gamma and cannot generally
reduce `N_phi` for tilted particles; `--alpha-points` and `--adaptive-alpha`
are descriptive aliases, not extra approximations.

The Makefile auto-detects the NVIDIA compute capability with `nvidia-smi`.
If detection is unavailable it builds for `sm_86`. That fallback is not a
claim that `sm_86` is the best or supported target for every GPU. Query the
deployment GPU and override the architecture explicitly, for example:

```bash
make -C gpu -j GPU_ARCH=86
```

Object directories include architecture, precision, fast-math, and MPI state.
Run `make -C gpu clean` when changing the CUDA toolkit or host compiler.

Other precision variants:

```bash
make -C gpu float       # gpu/bin/mbs_po_gpu_float
make -C gpu float_fast  # gpu/bin/mbs_po_gpu_float_fast
make -C gpu double_fast # gpu/bin/mbs_po_gpu_double_fast
```

Objects are written under `gpu/build/`, so this build does not conflict with
the CPU MPI/OpenMP build.

For supported averaged coherent PO paths, one process automatically partitions
an orientation batch across visible GPUs. Use `CUDA_VISIBLE_DEVICES=0,1`, set
`MBS_GPU_MULTI=0` to disable this behavior, or set `MBS_GPU_MULTI_MAX=N` to cap
the device count. Process-parallel size scans are a separate mechanism driven
by `--scan-jobs` and `--gpu-devices`.

The default target uses float and CUDA fast math. Treat CPU double and CUDA
double results as the numerical reference, and calibrate float-fast on the
same particle, orientation mode, angular grid, and convergence limits before
production use.

CUDA runtime errors stop the process by default. This is intentional: otherwise
an invalid GPU build can silently fall back to CPU and ruin speed measurements.
For debugging only, set `MBS_GPU_ALLOW_FALLBACK=1` to restore the old fallback
behavior.
