# Coherent Backscatter Audit, 2026-08-14

## Scope

This audit checks whether MBS-fast converts ray amplitudes to intensities too
early, whether CPU/GPU batching separates mutually coherent terms, and whether
the production cutoffs can bias exact backscatter.

The probe used the built-in concave hexagonal column at `beta=37 deg`,
`gamma=19 deg`, `k_eq=18.41289202624`, `m=1.3116`, and `lambda=0.532 um`.
Cutoff comparisons used exact `theta=180 deg` and `phi=0 deg`.

## Result

For one particle orientation, both CPU and CUDA paths retain Jones amplitudes
until every output beam assigned to that orientation has been summed. Mueller
conversion happens once per orientation. Multi-GPU partitioning assigns whole
orientations to devices and therefore does not split one coherent beam sum
between devices.

The per-beam audit at reflection depth 4 produced 416 outgoing contributions:

| Quantity | Value |
|---|---:|
| `0.5 * norm(sum_b J_b)^2` | 0.0137133849549 |
| `0.5 * sum_b norm(J_b)^2` | 0.0369448516738 |
| Interference cross term relative to the incoherent sum | -62.8815% |
| Relative error between `sum_b J_b` and written Jones output | 2.34e-17 |
| Relative error between summed per-beam Mueller and `--incoherent` output | 0 |

Interference is therefore present and strongly destructive in this probe.
There is no premature Jones-to-Mueller conversion within one orientation.

## Reflection-depth convergence

With all cutoffs disabled, exact-backscatter `M11` was:

| `--max-reflections` | `M11(180 deg)` |
|---:|---:|
| 2 | 0.00201826577 |
| 4 | 0.01371338495 |
| 6 | 0.06465435994 |
| 8 | 0.08294429804 |
| 12 | 0.09416216180 |
| 16 | 0.09280063544 |
| 20 | 0.09181165326 |
| 24 | 0.09182191258 |

Depth 12 is about 2.55% above the depth-24 value for this one orientation.
Reflection depth must therefore be converged separately from orientation and
scattering grids.

## Cutoff sensitivity

At depth 24, the production-like custom settings
`--beam-cutoff-jones 0.001 --beam-cutoff-area 0.002`
`--trace-cutoff-importance 0.0001 --trace-max-beams 20000` produced
`M11=0.08079797162`, about 12.0% below the cutoff-off reference. The tracer
evaluated 3009 branches and rejected 806 (26.8%).

The preset profiles were close to the cutoff-off result in this probe:

| Profile | `M11(180 deg)` |
|---|---:|
| off | 0.09182191258 |
| safe | 0.09181820867 |
| balanced | 0.09182698430 |
| fast | 0.09190307223 |

This is a single-orientation diagnostic, not a universal error bound. Custom
relative cutoffs can remove weak amplitudes whose cross terms are not weak.
Backscatter production runs should be repeated with `safe` or `off` at selected
orientations and sizes.

## Ensemble interpretation

MBS-fast models coherent interference among ray/aperture paths of one particle
at one orientation. Random independent orientations must be averaged as Mueller
matrices, because their relative Jones phase is unspecified.

The disabled `--coherent-orientations` option did not implement a coherent
ensemble: its legacy path cleared the Jones accumulator for every orientation
and produced the ordinary orientation average. A physically coherent sum over
distinct particles additionally requires their positions, incident phases, and
inter-particle multiple scattering. The current single-particle PO model does
not provide that coherent-backscattering-medium calculation.

## Regression

Run:

```bash
make test_coherence
```

The test reconstructs exact-backscatter Jones output from every retained beam,
checks the explicit incoherent branch, and requires a nonzero interference
cross term.
