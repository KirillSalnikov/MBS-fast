# Coherent CUDA fullauto benchmark

The initial commit of this branch (7f0f92f) preserves the CUDA source used for the RTX 3080 Ti shortwave timing experiment. The controller times candidate azimuth rules on matched scrambled Sobol pilots. Its working-grid cost model charges shared CPU tracing to the largest required orientation prefix and diffraction to the mean of the row requirements. The ordinary independent convergence checks still use the selected physical estimator.

The machine-readable record in `shortwave.json` contains the arguments, common and candidate-only environment, geometry checksum, timings, selected parameters, and worst differences over the complete angle grid. Replace `PARTICLE.dat` with the matching input geometry. Run the baseline from commit `4294c7f56`, the candidate from this branch, then the baseline again, using separate output directories. The reference used two persistent GPU workers with one dispatch thread.

For the 10 µm arbitrary particle at 0.355 µm the three complete runs took 758.1, 569.3, and 749.2 seconds. The faster of the two baseline runs divided by the candidate time gives 1.316×. All three selected depth 10 and 16 azimuth samples. The largest matrix difference was 0.904% of the reference M11 over 2454 angles. The largest difference of normalized elements was 0.519 percentage points.

Results depend on the workload and timing-based rule selection. A later four-core, single-worker variant selected 8 azimuth samples and exceeded the 1% matrix comparison threshold. This benchmark profile is intended for reproducing the two-worker result.

## Dense mean variant

The next commit adds bounded grouping of the original Gauss nodes and a conservative projected-facet test. The four-thread configuration in `dense_mean.json` uses one GPU worker. All three runs converged on the same 2454 angles with depth 10 and 16 azimuth samples. The baseline, candidate, and repeat baseline took 583.26, 267.59, and 606.95 seconds. Using the faster baseline gives 2.180× for the complete calculation. Compared with either baseline the candidate differed by at most 0.904% of M11 across the raw elements and 0.518 percentage points across normalized elements.

The same candidate source was tested on a second RTX 3080 Ti with the 31.62 µm droxtal at 1.550 µm. All three runs met the requested 4% tolerance over 1027 angles and selected depth 12 and 16 azimuths. The baseline took 1305.36 s, the candidate 871.48 s, and the repeat baseline 1280.86 s, giving **1.47×** against the faster baseline. The largest raw difference was 1.065% of reference M11 for M33 at 31.58°. The largest normalized difference was 0.750 percentage points. This case exceeded the 1% raw comparison threshold. [Arguments, binary checksums and both reference comparisons](dense_mean_droxtal.json).
