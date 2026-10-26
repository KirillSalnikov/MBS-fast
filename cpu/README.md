# CPU MPI/OpenMP build

This directory builds the CPU-only MPI + OpenMP binary. CUDA files are not
compiled, and `--backend cuda` is intentionally unavailable.

For automatic orientation averaging, start one controller:

```bash
make -C cpu -j2
cpu/bin/mbs_po_mpi --method po --backend cpu \
    --particle 1 5 2 --ri 1.3116 0 --wavelength-um 0.532 \
    --fullauto 0.02 --output results/cpu_fullauto
```

`0.02` requests a 2% numerical tolerance. The controller selects the orientation
counts, reflection depth and angular work; repeat the same command to resume.
Do not launch multiple fullauto controllers on one output directory with MPI.
See the [fullauto guide](../docs/FULLAUTO.md) and
[Russian guide](../docs/FULLAUTO_RU.md).

For a directly specified orientation/angular grid, MPI remains available:

```bash
make -C cpu -j2
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
