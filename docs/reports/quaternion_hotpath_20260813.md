# Quaternion hot-path experiment, 2026-08-13

## Question

Can the SO(3) quaternion representation accelerate the production physical
optics calculation when moved into its hot path?

## Production profile

Completed GPU runs at wavelength 0.355 um spend a median 8.63% of wall time in
phase 1 (rotation, tracing, and beam preparation). For large target particles,
phase 1 falls to 3.43-4.0%; CUDA diffraction is the dominant phase.

With `MBS_ORIENTATION_TIMING=1`, 65,536 quaternion orientations of a small
hexagonal particle produced this summed CPU-time breakdown:

```text
rotation:         0.02 s (0.29%)
tracing:          3.41 s
beam preparation: 2.94 s
```

For a hollow column with `L=50 um`, `D/L=0.7`, `mu=1`, and 12 reflections,
the rotation share fell to approximately 0.01% of phase-1 CPU work.

## V100 FP64 probe

The command

```bash
CUDA_VISIBLE_DEVICES=4 ./bin/gpu_quaternion_rotation_probe 65536 64 50
```

measured 5.41 billion direct quaternion-vector rotations per second and 5.37
billion matrix-vector rotations per second. Direct quaternion rotation was
0.8% faster. Across 24, 64, and 256 vectors per orientation, its advantage was
0.8-2.7%. The maximum absolute difference was `6.66e-16`.

## Decision

Keep quaternions as the compact representation for uniform full-SO(3)
sampling. Do not move particle rotation alone to the GPU: current ray tracing
consumes the rotated facets on the CPU, so this would add device transfers to
optimize at most 0.01-0.29% of phase-1 work. Do not store a 3x3 matrix per
orientation either; it roughly doubles orientation storage and produced no
measurable end-to-end speedup.

A meaningful GPU architecture change would have to move ray tracing and beam
preparation together, not just quaternion arithmetic.
