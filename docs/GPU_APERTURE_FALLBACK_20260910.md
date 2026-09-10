# Unsupported-aperture CPU fallback

Base: `290b1123c4f6be9dcea38fac193084caf249c2b2`.

Reproducer: stretched Sahara particle p0175, r_eq=8 um, wavelength .532,
m=1.4, 276 beta/gamma orientations, 181 theta rows, latitude-phi cap 300,
8 internal reflections, production cutoffs unchanged.

`Handler::PrecomputeEdgeData` rejects polygons with 32 or more vertices
(`nVertices >= BeamEdgeData::MAX_EDGES`). It retains a full polygon in
`PreparedBeam::fallback` for the CPU calculation. The CUDA direct path used
to return `bad-edge-data`, which aborted the variable-phi run. Repeating
with orientation chunks 32, 16 and 8 could not fix this representation
limit. A debug replay (`MBS_FFT_DEBUG=1`, FFT still disabled) reproduced it
with ample free device memory; the repaired diagnostic identified a
32-vertex source contour.

The direct CUDA handler now splits the orientation batch around an
unsupported contour, evaluates the complete affected orientation using
the existing CPU coherent Jones path, then adds its weighted Mueller
matrix. Supported neighboring orientations still run on CUDA. Beams are
not dropped, clipped further, or independently added as Mueller matrices;
the latter would lose coherent cross terms. The GPU packing cache remains
keyed by token/start/count. The change is restricted to the default
scale=1, waveIndex=0 path; unsupported scaled/multisize calls still fail
explicitly rather than silently using incorrect CPU scale constants.

Build/validation artifacts on epyc2:
`campaigns/mbs_stretched_req6_30_xi2_div280_epyc2_20260909/repair_20260910/`.
The original source and binary in the 2..5 um campaign are untouched.
The production controller accepts both original and validated repaired
binary hashes, skips existing completed results, and records the actual
binary hash in every newly completed job. Physics grids/cutoffs are not
changed and partial ensembles are not averaged as complete.

Regression commands and comparison are recorded by `diagnose_failure.py`
and `validate_repair.py` in the extension campaign. They test the original
failing command, repaired GPU result against original full CPU result,
and p0174 at the same radius against its existing production GPU result.
The quantitative results are in `repair_20260910/validation.json`.

The repository-side regression checker can consume those retained runs:

```sh
python3 tests/check_cuda_aperture_fallback.py \
  --gpu-run /path/to/repair_20260910/repaired_gpu \
  --cpu-run /path/to/repair_20260910/original_cpu_reference
```

It requires a logged fallback involving at least 32 source vertices,
matching physical inputs and orientation/angle grids, finite complete
181-row output, and agreement of all 16 raw Mueller elements normalized
by reference M11. The default FP64 tolerance is 1e-6. The original broken
GPU run fails this check because it has no completed Mueller output.
For a new build, regenerate the paired runs before running the checker;
checking old output alone does not test a newly compiled binary.

## Validation on epyc2, 2026-09-10

FP64, precise phase trigonometry, CUDA 12.2 / sm_70, GCC 12.4,
two CPU threads per job, concurrent ADDA left running:

- Original p0175 CUDA command: failed at `bad-edge-data`, including an
  isolated debug replay (67.97 s).
- Repaired p0175 CUDA/CPU-hybrid: complete, 126.96 s; 32-vertex fallback
  exercised, all 276 orientations and all 181 output rows retained.
- Original p0175 full CPU reference: complete, 642.49 s.
- Maximum absolute difference of all 16 elements divided by CPU M11:
  **8.50188236050939e-9**; maximum relative M11 difference **4.165069178441172e-9**;
  backscatter relative M11 difference **9.89814230578645e-10**.
- p0174 unchanged-path control vs its original production CUDA result:
  maximum difference divided by M11 **2.021483497155442e-16**.
- `tests/check_cuda_aperture_fallback.py`: PASS on the paired p0175 outputs.

These timings are validation timings under concurrent CPU load, not an
isolated performance benchmark. Repair binary SHA256:
`7fd3f098a77839fb8755f5a89218c6f908f9665f4bae91ae9c6276cb7555a7ca`.
