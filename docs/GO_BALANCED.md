# Balanced cutoffs for GO fullauto

`--method go --fullauto EPS --cutoff-profile balanced` supports a fixed
balanced preset. `--fullauto-cutoffs off` disables PO-specific preset search;
it does **not** disable the explicitly selected balanced preset.

The thresholds are dimensionless:

| Control | Balanced value | Reference |
|---|---:|---|
| Outgoing beam importance | 1e-8 | Strongest outgoing beam in this orientation |
| Trace branch importance | 1e-10 | Strongest incident branch |
| Small-fragment area ratio | 1e5 | Existing trace fragment simplification rule |
| Individual Jones/area cutoffs | 0 | Disabled |

GO outgoing importance is `||J||_F^2 * projected_cross_section`. References
and tests precede absorption and scattering-plane rotations. Positive
absorption can only reduce the absolute energy of a discarded beam.
Forward/backward poles retain their integrated-mass interpretation. Output
beam filtering is shared by the total and track-resolved GO handlers;
OpenMP handlers inherit thresholds and share atomic diagnostic counters.
An explicit `off` keeps the original accumulation path, without filtering.
Neither the preset nor this change enables mirror symmetry for unsupported
geometries or absorbing materials.

## Accuracy control

Balanced GO fullauto reserves `EPS/8` for cutoff bias. After statistical and
independent validation, it compares balanced and off on matched orientations
from eight fresh scrambles. The pointwise 95% paired bound, scaled by M11
(and polarization ratios when enabled), must pass this cutoff budget.
The existing depth+2 check is also performed. With automatic depth, EPS/4
is reserved for depth truncation; fixed-depth runs retain the existing
separate depth-check policy. No claim of simultaneous global confidence is
made. `--fullauto-final-checks off` or `--fullauto-physics-audit off` is rejected
with this preset. A failed audit returns status `physics_validation_failed`
(code 4), not a certified result. Hard beam limits and retry safeguards
remain independent controls.

The actual preset, reserved budget and measured cutoff residual are recorded
in the status/parameter files and `physics_validation.json`. A changed
profile or binary does not reuse an incompatible fullauto checkpoint.
Use a fresh output directory when changing an existing off campaign.

## Bounded checks

```
python3 tests/test_go_balanced.py --binary bin/mbs_po --output /tmp/go_balanced_new
python3 tests/test_go_mirror.py --binary bin/mbs_po --output /tmp/go_mirror_new
```

The balanced test compares all 16 elements and both pole masses against off
for six cavity depths, two trace depths, a fixed Sobol sample and repeated
wall-time measurements. It also checks threads, mirror reconstruction,
area scaling, fullauto auditing, resume identity and automatic-depth setup.
The report includes pointwise errors in weak bins as well as peak-scaled
and integrated differences: a small integrated difference alone is not a
production precision certificate. Timing includes process startup and I/O;
subsecond cases should not be interpreted as precise speed benchmarks.
