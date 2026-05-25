# CUDA GPU build

This directory builds the CUDA binary. The GPU backend is enabled by default in
this build, so `--gpu` is not required. Use `--cpu` only when you deliberately
want to run the CPU backend from the GPU-capable binary.

```bash
make -C gpu -j
gpu/bin/mbs_po_gpu_float_fast --po --fft --autofull 0.05 ...
```

Other precision variants:

```bash
make -C gpu float       # gpu/bin/mbs_po_gpu_float
make -C gpu float_fast  # gpu/bin/mbs_po_gpu_float_fast
make -C gpu double_fast # gpu/bin/mbs_po_gpu_double_fast
```

Objects are written under `gpu/build/`, so this build does not conflict with
the CPU MPI/OpenMP build.
