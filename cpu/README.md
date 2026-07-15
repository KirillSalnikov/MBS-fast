# CPU MPI/OpenMP build

This directory builds the CPU-only MPI + OpenMP binary. CUDA files are not
compiled, and `--backend cuda` is intentionally unavailable.

```bash
make -C cpu -j
/usr/bin/mpirun --oversubscribe -np 4 cpu/bin/mbs_po_mpi \
    --method po --backend cpu --particle 1 10 10 \
    --refractive-index 1.3116 0 --wavelength-um 10 \
    --max-reflections 1 --sobol 64 \
    --scattering-grid 0 180 12 24 --threads 1 \
    --output results/cpu_mpi_example --close
```

Objects are written under `cpu/build/`, so this build does not conflict with
the CUDA build. Use `mpicxx` and `mpirun` from the same MPI installation; set
`CXX=/path/to/mpicxx` at build time when the system default is not the intended
wrapper.
