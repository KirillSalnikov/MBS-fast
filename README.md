# MBS-raw: Physical Optics for Ice Crystals

Fast Kirchhoff diffraction code for light scattering by non-spherical ice particles.
~70x faster than original, AVX-512/AVX2 SIMD, OpenMP parallel.

## Requirements

- **GCC >= 9** (or Clang >= 14)
- **OpenMP** (libgomp for GCC, libomp for Clang)
- x86-64 CPU with AVX2+FMA (minimum) or AVX-512 (optimal)

## Build

### Intel (Skylake-X, Ice Lake, Sapphire Rapids, ...)

```bash
bash build.sh              # -> bin/mbs_po
```

### AMD EPYC 7H12 (Zen 2, SP3)

```bash
bash build_epyc.sh         # -> bin/mbs_po_epyc
```

No AVX-512. SIMD via AVX2+FMA3. For Clang/AOCC:

```bash
sudo apt install libomp-16-dev   # if using clang-16
bash build_epyc_clang.sh         # -> bin/mbs_po_epyc_clang
```

### AMD Zen 4 (Ryzen 7000, EPYC 9004 Genoa)

```bash
bash build_zen4.sh         # -> bin/mbs_po_zen4
```

AVX-512 supported on Zen 4 (256-bit execution units).

### With make (if installed)

```bash
make                       # Intel, uses Makefile
make -f Makefile.epyc      # EPYC 7H12
```

### Generic (any AVX2+FMA CPU)

```bash
g++ -O3 -march=haswell -std=gnu++11 -funroll-loops -fopenmp \
    $(find src -not -path '*/bigint/*' -name '*.cpp') \
    $(find src/bigint -name '*.cc') \
    -Isrc -Isrc/math -Isrc/handler -Isrc/common -Isrc/geometry \
    -Isrc/geometry/intrinsic -Isrc/geometry/sse -Isrc/particle \
    -Isrc/scattering -Isrc/tracer -Isrc/splitting -Isrc/bigint \
    -o bin/mbs_po -lm -lgomp
```

## Quick Start

```bash
# Hex column D=H=10um, 1024 Sobol orientations, auto theta grid
bin/mbs_po --po --sobol 1024 --auto_tgrid \
    -p 1 10 10 -w 0.532 --ri 1.31 0 -n 12 \
    --grid 0 180 48 180 --close

# Adaptive convergence (auto orientations)
bin/mbs_po --po --adaptive 0.01 --auto_tgrid \
    -p 1 10 10 -w 0.532 --ri 1.31 0 -n 12 \
    --grid 0 180 48 180 --close
```

## Multi-core (OpenMP)

```bash
# Intel / desktop
OMP_NUM_THREADS=12 bin/mbs_po --po --sobol 1024 ...

# EPYC 7H12 (64 cores, NUMA-aware)
OMP_PROC_BIND=close OMP_PLACES=cores OMP_NUM_THREADS=64 \
    bin/mbs_po_epyc --po --sobol 4096 ...
```

## Documentation

- **MANUAL.md** / **MANUAL.pdf** — full CLI reference, all flags, performance notes
- **BUGFIX_forward_direction_sign.md** — forward-direction Fresnel sign fix
- **tests/run_tests.sh** — regression tests

## Performance

| CPU | Build | Phase 2 (1 thread) | Speedup |
|-----|-------|--------------------|---------|
| Intel (AVX-512) | `build.sh` | 4.6 s | ~70x |
| EPYC 7H12 (AVX2) | `build_epyc.sh` | 4.7 s | ~68x |
| Original code | — | ~320 s | 1x |

Benchmark: hex D=H=10um, 128 Sobol, 48 phi x 181 theta, n=6.
