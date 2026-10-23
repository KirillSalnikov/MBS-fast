# GO fullauto mirror reduction

Explicit `--fullauto-mirror on` now supports the built-in regular hexagonal
column (type 1) and the column with identical cavities at both ends (type 10),
with a real refractive index (zero absorption).
With `--symmetry 2 6`, beta spans 0..90 degrees and the sampled gamma spans
0..30 degrees instead of 0..60. Alpha remains analytically averaged.

This is **not** the same as using `--symmetry 2 12`. The missing reflection
partner is reconstructed for every beam as `(M + P M P)/2`, with
`P=diag(1,1,-1,-1)`, in the scattering-plane basis. Reconstruction applies
before angular binning, including forward/backward delta-ray matrices, and
is preserved by OpenMP thread-local handlers. Normalization remains the
mean over the traced half-domain points; no additional factor of two is used.

Supported beta factors are 1 and 2, and gamma factors 1, 2, 3 and 6.
Particle files, absorbing materials and other built-in types are not silently assumed supported:
explicit `on` rejects them; explicit `auto` retains full gamma for them.
For GO, the default remains `off` to preserve existing calculations. For PO,
the existing geometry/pilot verification and default policy are unchanged.
Absorption is deliberately excluded: a paired test of the existing absorption
quadrature did not satisfy the reconstruction tolerance. Changing that
quadrature is outside this feature; the ice campaign uses zero imaginary index.

The controller applies the same mirror mode during automatic depth selection,
production, independent validation and cutoff/depth checks. Worker commands
and cache signatures distinguish mirror-on/off; old results cannot be resumed
under changed symmetry or a different binary without their existing explicit
compatibility checks. The geometry verification is a built-in-shape contract,
not a new per-run numerical mirror pilot; metadata states this explicitly.

## Independent bounded regression test

```bash
python3 tests/test_go_mirror.py --binary bin/mbs_po \
  --output /tmp/go_mirror_validation_new_directory
```

Optionally provide `--old-binary PATH` to compare mirror-off native results.
The test compares the reconstruction with **explicitly traced reflected
partners** using `--haar-mirror-audit`. For GO that audit needs no `--haar-alpha`:
alpha is analytic. An audit with N base Sobol points traces 2N orientations;
a mirrored calculation with the same seed and offset traces only N.

Tests cover all 16 elements and both delta rows, solid columns and plates,
several cavity depths including almost contacting cavities,
different reflection depths and rotational domains, one/four threads,
contiguous Sobol windows, fullauto worker propagation and resume identity,
and rejection/fallback for unsupported geometries/domains and absorbing inputs. These are short
software checks, not a claim of production convergence or GO model accuracy.
