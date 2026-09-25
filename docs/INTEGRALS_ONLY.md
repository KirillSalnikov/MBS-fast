# CPU forward-integral mode

`--integrals-only` computes orientation-averaged extinction from the coherent
forward amplitude (optical theorem), without calculating an angular Mueller
grid. It requires PO, explicit convex geometry, CPU execution and Hammersley
orientations. It always samples the full beta/gamma domain without symmetry
reduction; one MPI rank with OpenMP threads is supported.

```sh
make -j8
bin/mbs_po --method po --backend cpu --geometry convex \
  --particle-file examples/cube.particle --refractive-index 1.5 0.01 \
  --wavelength-um 1.064 --dmax-grid 6 100 16 \
  --hammersley 8192 --integrals-only --max-reflections 8 \
  --beam-cutoff-jones 0.001 --beam-cutoff-area 0.002 \
  --trace-cutoff-importance 0.0001 --trace-max-beams 20000 \
  --trace-limit-retries 0 --threads 8 --output integrals --close
```

Single sizes, `--dmax-grid`, `--k-eq-grid` and `--k-eq-list` are supported.
Ray geometry is traced once per orientation at the smallest requested size.
Areas, optical paths, phases and absorption are scaled for the other sizes.
Internal ray tracing is still necessary; this is not an assumption that
extinction equals twice the projected area. Validate cutoff, reflection-depth
and orientation convergence for each application. A failed orientation aborts
the calculation instead of silently contributing partial results.

The result is `<prefix>_fast_integrals.tsv`. `Cext_OT` is the forward optical
theorem result within the chosen PO model and numerical settings. `G` is the
mean projected area; `Qext_OT = Cext_OT/G`. Cross sections use the square of the
input length unit. Fixed-order reduction makes results independent of OpenMP
scheduling. Angular sampling flags do not determine the cost of this mode.

**No angular Csca integral is computed.** For albedo, combine Cext with an
adequately resolved absolute M11 for the same shape, optical constants and size:

```
Csca = 2*pi * integral(M11(theta)*sin(theta), theta=0..pi)
albedo = Csca/Cext
Cabs = Cext-Csca
```

Angles in this integral are radians. Different orientation grids remain a
source of mismatch; a coarse forward peak can make the Csca integral unreliable.
Do not clip nonphysical balances to hide these errors.

The other TSV columns are explicitly **approximate ray-energy diagnostics**.
In particular, `Cabs_ray_loss`, `Csca_OT_minus_ray_loss` and `albedo_OT_ray`
are not substitutes for coherent PO integrals. Finite ray depth, pruning and
wave interference can make their albedo disagree significantly with the full
PO result. `ray_closure_relative` records unclosed unabsorbed ray flux.

Run the self-contained CPU regression (no article data or external Python
packages required):

```sh
python3 tests/regression_integrals_only.py --binary bin/mbs_po
bash tests/run_cli_tests.sh
```

The regression compares fast/full forward Cext on identical orientations,
one/four-thread reproducibility, shared/independent sizes, and rejected options.
Existing solver modes are unchanged when the new flag is absent.
