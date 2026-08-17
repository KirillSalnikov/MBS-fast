# Coherent Backscatter Audit, updated 2026-08-17

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

The per-beam audit at reflection depth 4 now produces 607 outgoing
contributions. Before target facets were clipped at the source aperture plane,
the same case produced only 416 contributions:

| Quantity | Value |
|---|---:|
| `0.5 * norm(sum_b J_b)^2` | 0.0512922268562 |
| `0.5 * sum_b norm(J_b)^2` | 0.0579249895 |
| Interference cross term relative to the incoherent sum | -11.451% |
| Relative error between `sum_b J_b` and written Jones output | 2.13e-16 |
| Relative error between summed per-beam Mueller and `--incoherent` output | 0 |

Interference is therefore present and strongly destructive in this probe.
There is no premature Jones-to-Mueller conversion within one orientation.

## Reflection-depth convergence

With all cutoffs disabled, exact-backscatter `M11` was:

| `--max-reflections` | `M11(180 deg)` |
|---:|---:|
| 2 | 0.01100370283 |
| 4 | 0.05129222686 |
| 6 | 0.06511251830 |
| 8 | 0.08860667759 |
| 12 | 0.09273171485 |
| 16 | 0.09242936386 |
| 20 | 0.09181820477 |
| 24 | 0.09182458069 |

Depth 12 is about 0.99% above the depth-24 value for this one orientation.
Reflection depth must therefore be converged separately from orientation and
scattering grids.

## Cutoff sensitivity

At depth 24, the production-like custom settings
`--beam-cutoff-jones 0.001 --beam-cutoff-area 0.002`
`--trace-cutoff-importance 0.0001 --trace-max-beams 20000` produced
`M11=0.07956817292`, about 13.3% below the cutoff-off reference.

The preset profiles were close to the cutoff-off result in this probe:

| Profile | `M11(180 deg)` |
|---|---:|
| off | 0.09182458069 |
| safe | 0.09182079600 |
| balanced | 0.09183413793 |
| fast | 0.09187037688 |

This is a single-orientation diagnostic, not a universal error bound. Custom
relative cutoffs can remove weak amplitudes whose cross terms are not weak.
Backscatter production runs should be repeated with `safe` or `off` at selected
orientations and sizes.

For a 128-orientation `k_eq=92.06446` check over `170..180 deg`, depth 6 was
11.03% below depth 20 in weighted M11 L2, depth 12 differed by 0.69%, and depth
16 by 0.10%. The production-like manual cutoffs were 3.37% below the uncut
depth-12 result in the same metric, whereas the `safe` profile differed from
`off` by only 0.00069%.

## Forward-depth clipping regression

The exact intersection path clipped a target facet at the source aperture
plane, but `Difference`, used for parent subtraction and outgoing occlusion,
only checked whether any target vertex was forward. It then projected the
complete tilted facet, including its behind-plane portion. That portion could
shadow a physical outgoing aperture and postpone valid ray paths until a much
larger reflection limit.

`Difference` now clips the target polygon by `t <= 0` before projection. The
half-space predicate avoids division by a nearly parallel denominator and is
tested at geometry scales `1e-6`, `1`, and `1e6`. At depth 4 the regression
case changes from 416 beams and `M11=0.01371338495` to 607 beams and
`M11=0.05129222686`. Over 512 orientations at `k_eq=92.06446`, however, the
change is only 0.139% weighted M11 L2 over `160..180 deg` and 0.073% at exact
backscatter, so it does not by itself explain a large ensemble discrepancy.

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
