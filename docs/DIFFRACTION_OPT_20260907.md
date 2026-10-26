# Direct FP64 multisize phase cache

The 2.0 integration ports the phase cache and memory guards from `290b112`,
retaining the current guarded vertex-phase arithmetic. The old streaming-edge,
adaptive-phi and mean-only FFT experiments are not enabled by this port.

A coherent, direct, full-only CUDA scan of 2–32 sizes uses the fused multi-k
path by default in FP64. FP32 keeps its previous opt-in selection. A failed
fused request falls back to the existing per-size GPU calculation.
`MBS_GPU_MULTI_K_FULL=0` selects the legacy per-size path;
`MBS_GPU_MULTI_K_PHASE_CACHE=0` disables only the phase cache.

Each beam/size has eight real Jones components. Their direction-independent
optical-path phase and absorption are evaluated once and reused at every
angular cell. The cache needs `8 * nSizes * nBeams * sizeof(GpuReal)` bytes.
It is refreshed for every chunk. Optional allocation is checked against the
available VRAM budget; if it cannot fit, the uncached calculation is retained.
The multi-size output allocation is checked separately. Repacking shared
beam data invalidates the single-size prepared cache before fallback.

Neither sampling, beam cutoffs nor coherent summation changes. Fullauto still
runs one particle size per controller; this optimization concerns the separate
legacy shared-size path.

```bash
CUDA_VISIBLE_DEVICES=0 python3 tests/regression_multik_phase_cache.py \
  --binary gpu/bin/mbs_po_gpu_double --particle examples/cube.particle
```

The regression compares cached, uncached and per-size results for three sizes,
multiple chunks, real and absorbing indices, and all 16 elements scaled by M11.
Its FP64 tolerance is 1e-9. Release-specific measurements belong in the release
notes; older branch timings are not measurements of the newly merged binary.
