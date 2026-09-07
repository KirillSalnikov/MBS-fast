# FP64 shared multi-size diffraction optimization

The serial direct CUDA size scan now selects the fused multi-k implementation
automatically for 2..32 sizes in FP64 builds. Other modes preserve their prior
default selection; `MBS_GPU_MULTI_K_FULL=0` explicitly restores the old per-size
path. Failed fused processing continues through the existing per-size fallback.

The new `prepare_multik_jones_kernel` evaluates the direction-independent
optical-path phase and absorption-adjusted Jones matrix once per beam/size.
The angular kernel reuses these eight FP64 components instead of repeating
the same `sincos`, absorption path traversal, exponentials, and Jones phasing
at every theta/phi cell. Polygon diffraction, angular sampling, coherent beam
summation, Mueller accumulation, and all physics cutoffs remain unchanged.

The cache is a persistent CUDA workspace allocation, refreshed on every chunk
and resized as necessary. Its size is `nSizes * nBeams * 8 * sizeof(GpuReal)`.
If available memory cannot accommodate it within `MBS_GPU_MEM_FRACTION`, the
uncached kernel executes. The fused output allocation also checks available
VRAM before proceeding, because the single-size batch planner does not account
for its larger Mueller output. Repacking the workspace invalidates the
single-size prepared-beam cache before any possible fallback.

`MBS_GPU_MULTI_K_PHASE_CACHE=0` disables only the new cache while retaining
multi-k, allowing an isolated before/after comparison.

Validation command (use an idle GPU):

```sh
CUDA_VISIBLE_DEVICES=0 python3 tests/regression_multik_phase_cache.py \
  --binary gpu/bin/mbs_po_gpu_double --particle /path/to/particle.dat
```

This integration test compares cached, uncached, and legacy paths for three
sizes, multiple orientation chunks, and both zero and nonzero imaginary
refractive index. It validates every Mueller element against M11 with a
1e-9 relative scale tolerance and retains all logs and matrices.

Production-scale probes use the Sahara p0000 (68 facets) and p0070 (148 facets)
meshes, 8x8 orientations, 181 theta points, 600 phi samples, lambda=0.532 um,
m=1.4+0i, n=8, and the current queue cutoffs. For four radii 4, 35, 67,
100 um, diffraction times on a V100 were:

| Shape | Legacy per-size | Multi-k uncached | Multi-k cached |
|---|---:|---:|---:|
| p0000 | 10.05 s | 4.83 s | 4.45 s |
| p0070 | 20.63 s | 9.66 s | 8.79 s |

The new phase cache contributes about 8..10 percent above the pre-existing
multi-k optimization. Total measured acceleration against the old default
is 2.26..2.35x. Maximum pointwise M11 differences against legacy were below
4e-11 percent in these probes. These are implementation equivalence tests,
not a convergence claim for the full orientation/particle ensemble.

The final build was also checked on the full p0000 block of 32 radii (4..35 um):
legacy direct diffraction 81.37 s, uncached multi-k 38.26 s, cached multi-k
34.80 s. Total speedup is 2.338x; phase caching itself gives 1.099x.
Maximum pointwise M11 error is 3.766e-11 percent, and maximum absolute
Mij/M11 difference is 6.312e-13.

The integration regression above passed on p0070, including imaginary index
0.01 and multiple orientation chunks. Maximum all-element error divided by
M11 was 1.152e-12, below the 1e-9 acceptance threshold. The uncached control
also passed. The memory guard is inspected and compiled; an artificial CUDA
OOM was not forced on the shared production server.

First phase-cache build SHA256 (retained as gpu/bin/mbs_po_gpu_double_phasecache_v1):
`ee61760d4967839b465427ef34074620293b92a3dbb666def0355f356fa1ebc5`.
Source and binary are isolated in MBS-fast-diffraction-opt-20260907. The
ongoing campaign still uses its original binary; benchmark runs resumed
the original GPU0 process after finishing.

## Second optimization: streaming polygon phases

The multi-k integral can now carry each edge's endpoint phasor into the next
edge instead of indexing two 32-entry vertex arrays. The first vertex is
retained for closing the polygon. If the list of valid edges skips a vertex,
its phase is explicitly recomputed; no adjacency assumption changes the
mathematics. Small-phase moment evaluation and the near-singular edge quotient
retain their existing implementations, and edges accumulate in their original
order. No FFT, phase approximation, or reduction of the sampling grids is used.

Phase-cached multi-k uses this path by default, with 128 CUDA threads per block
unless MBS_GPU_BLOCK is explicitly set. MBS_GPU_MULTI_K_STREAM=0 restores the
prior vertex-array integral. If phase caching is disabled or cannot allocate
its buffer, the old uncached multi-k path remains available.

| Shape, four radii, 64 orientations | Vertex arrays / block 64 | Streaming / block 64 | Streaming / block 128 |
|---|---:|---:|---:|
| p0000 | 4.39 s | 3.97 s | 3.78 s |
| p0070 | 8.82 s | 7.79 s | 7.43 s |

All output Mueller values in these comparisons were identical at the written
FP64 precision. Additional acceleration is 1.16..1.19x over the phase-cache
version. Streaming with block 256 was slower (5.47 and 10.78 s).

A separate experiment tiled 2 or 4 sizes inside each thread. It was slower
(6.30/8.67 s and 12.66/17.28 s respectively), and was removed from the final
source. The data is retained outside the checkout in diffraction_tile_bench_20260907.

Final 32-size p0000 confirmation (same 64 orientations, 181x600 angles):
vertex arrays 34.59 s, streaming default 29.22 s. The second optimization
gives 1.184x over the paired previous implementation and total 2.785x over
the original 81.37 s per-size direct reference. All 32 output matrices matched
the paired vertex-array results exactly as written; maximum M11 difference
against the original reference remains 3.766e-11 percent.

The expanded regression passed both p0070 and examples/cube.particle, with
imaginary indices 0 and 0.01. It exercises multiple orientation chunks,
three sizes, explicit cache-off, explicit vertex-array, and automatic streaming
paths. Maximum all-element differences divided by M11 against legacy were
1.152e-12 (dust) and 4.716e-13 (cube).

Final streaming-build SHA256:
`6ee02be7158ad095700b24bbd9edb905112e0f900dc29d569d17d5f0a9ae7350`.
Final probe artifacts are in ../diffraction_stream32_bench_20260907/ and
../diffraction_stream_regression_20260907/. The production queue was resumed
on its existing binary after these temporary tests.
