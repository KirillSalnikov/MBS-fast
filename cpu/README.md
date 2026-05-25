# CPU MPI/OpenMP build

This directory builds the CPU-only MPI + OpenMP binary. CUDA files are not
compiled, and `--gpu` is intentionally unavailable.

```bash
make -C cpu -j
mpirun -np 4 cpu/bin/mbs_po_mpi --po --sobol 4096 ...
```

Objects are written under `cpu/build/`, so this build does not conflict with
the CUDA build.
