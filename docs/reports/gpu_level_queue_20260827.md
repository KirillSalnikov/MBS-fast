# GPU level queue experiment (2026-08-27)

## Scope

The experimental `MBS_GPU_TRACE_LEVEL_QUEUE=1` scheduler processes active beam
apertures in increasing `Beam::nActs` order. All apertures at the current
reflection level are drained before the next level starts. Packed aperture and
candidate data use the existing reusable CUDA workspace; exact first hits are
enabled separately with `MBS_GPU_TRACE_EXACT_HIT=1`.

This is a hybrid implementation. CUDA performs candidate ordering, filtering,
and the first exact polygon intersection. CPU code still performs polygon
subtraction, Fresnel/Jones updates, optical-path updates, and construction of
the following reflection level.

## Validation case

- particle: concave hexagonal column, `L=31.5575`, `D=22.22807196634097`,
  cavity parameter `36.74404811468746`;
- wavelength `0.532`, refractive index `1.3116+0i`;
- 512 Sobol orientations, seed 42;
- 12 reflections, FP64 V100 build;
- cutoff parameters `J=0.001`, area `0.002`, importance `0.0001`;
- one host thread and one V100 for the deterministic comparison.

The level queue traversed 13 reflection levels in 19 CUDA batches. Its maximum
observed level width was 119 active apertures. With
`MBS_GPU_TRACE_VERIFY_EXACT=1`, all exact CUDA hits matched the CPU reference.

## Result

| mode | wall time | user time | peak RSS |
|---|---:|---:|---:|
| depth-first control | 5.56 s | 5.03 s | 200 MiB |
| reflection-level queue | 2.87 s | 2.46 s | 200 MiB |

The deterministic wall-time speedup is `1.94x`. The normalized differences are:

- `M11` relative L2: `1.95e-20`;
- largest L2 difference of any Mueller element, normalized by the control
  `M11` norm: `4.81e-20`;
- backscatter `M11`: identical at printed precision.

With eight simultaneous orientation workers, repeated 8192-orientation runs
showed a smaller `1.1-1.3x` wall-time improvement because independent
orientations already provide GPU occupancy. Those parallel runs differ at the
few-parts-per-million level because OpenMP accumulation order is not
deterministic; the single-thread comparison isolates the scheduler itself.

## Next boundary

The bounded GPU polygon-difference prerequisite is now available behind
`MBS_GPU_TRACE_PARTITION=1`. It returns the reached target aperture and one
residual polygon with up to 20 vertices. A complement requiring several
polygons has a non-unique convex decomposition, so it deliberately triggers a
per-beam CPU fallback. This compact representation avoids transferring empty
polygon slots for every hit and preserves the CPU branch topology in ambiguous
cases. Device and pinned-host buffers are allocated only when the mode is
enabled. Per-beam area validation, overflow detection, and CPU fallback remain
mandatory.

The initial combined exact-hit/partition kernel required a 2560-byte local
stack and 128 registers per thread. Splitting it into two ordered kernels
reduced the exact-hit kernel to 1600 bytes/118 registers and the partition
kernel to 2080 bytes/76 registers. The compact transfer format also reuses the
reached polygon as the exact-hit polygon, and levels narrower than
`MBS_GPU_TRACE_PARTITION_MIN_ITEMS` (64 by default) remain on the CPU.

The final deterministic 512-orientation validation used a threshold of 16.
Strict verification recomputed 380490 exact-hit decisions on the CPU with zero
mismatches. CUDA prepared 129233 unambiguous partitions; 87166 multi-part
complements used the local CPU fallback. Verification output was bitwise equal
to the CPU-partition control. With verification disabled, the differences were
`1.10e-7` in relative L2 for `M11`, at most `1.11e-7` of the control `M11` norm
for any Mueller element, and `0.0030%` for backscatter `M11`.

The isolated partition is not yet a speedup: interleaved runs measured
2.70-2.76 s with CPU subtraction and 2.98 s with GPU subtraction. An earlier
6.95/4.38 s pair was contaminated by a concurrent production job and must not
be interpreted as a `1.59x` partition speedup. The transfer and launch cost is
larger than the CPU subtraction saved while every reflection level still
returns to the host. Consequently the partition remains opt-in and is useful
as a validated prerequisite for a fully device-resident queue, not as a
standalone production optimization.

A fully device-resident `current`/`next` queue still requires CUDA versions of
polarization-basis rotation, Fresnel/Jones updates, optical paths, pruning, and
whole-orientation overflow recovery. Until those agree with the CPU
implementation, the level queue and GPU partition remain opt-in.
