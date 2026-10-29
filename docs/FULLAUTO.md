# Fullauto: practical and mathematical guide

This guide describes the integrated code prepared for MBS-fast 2.0. For the
extended Russian explanation and worked parameter definitions, see
[FULLAUTO_RU.md](FULLAUTO_RU.md) or [HTML](FULLAUTO_ADAPTIVE_RU.html).

## Start with the physical inputs

```bash
make -C cpu -j2
cpu/bin/mbs_po_mpi --method po --backend cpu \
  --particle 1 5 2 --ri 1.3116 0 --wavelength-um 0.532 \
  --fullauto 0.02 --output results/column_5um
```

The type-1 column has axial length 5 and basal diameter 2, in the same units as
the wavelength. `0.02` means 2%. For a file particle, replace `--particle 1 5 2`
with `--particle-file examples/cube.particle --resize-dmax-um 5`.
Dmax is the maximum vertex separation, not a volume-equivalent diameter.

No Python controller or hand-written orientation file is required. Do not add
an explicit Sobol selector to make fullauto use Sobol. Run one controller per
output directory, without a multi-rank MPI launch.

For CUDA, build `make -C gpu fp64 GPU_ARCH=70 -j2` and use
`gpu/bin/mbs_po_gpu_double --backend cuda` instead of the CPU executable/backend.
Architectures 70, 86 and 89 correspond to V100, RTX 3080 Ti and RTX 4090.
Add `--fullauto-theta-range 0 25` for the forward cone or
`--fullauto-theta-range 170 180` for the backward cone. The default GPU list is
`0`; `--fullauto-gpus 0,1` shares work across two visible devices. Indices are
relative to `CUDA_VISIBLE_DEVICES`.

A fullauto controller handles one size. Legacy multisize flags are a separate
workflow. Give independent sizes different output directories and GPU resources.

## Parameters: automatic versus fixed

| Quantity | Automatic behavior | Explicit option |
|---|---|---|
| Start per seed | Empirical size/wavelength/range model, at least 512 | `--fullauto-initial N` fixes only the start |
| Reflection target | Paired depth search from 6 through 18 | `--max-reflections N` fixes the target |
| Depth pilot | Starts at min(512, initial N) | `--fullauto-depth-pilot N` |
| Azimuth phi | Cost/variance calibration with sampled Haar alpha | `--phi-points N` fixes it |
| Output angles | Size-dependent reference grid | `--theta-grid-file FILE` fixes output angles |
| Interpolation | Local cubic support, matched residual checks | Allowance is EPS/4 |
| Pruning in PO | Balanced estimate plus coupled uncut correction, including incoherent PO | Explicit cutoff settings replace this automatic policy |
| Analytic controls | Independent fitting and validation | `--fullauto-analytic-controls off` disables them |
| Symmetry | Verify supported body rotations and mirror reduction | Explicit symmetry overrides are the user's physical assertion |
| Execution | Worker threads, chunks and cache/work limits | Explicit execution settings take precedence |

`--max-orientations` limits growth; exhausting it does not certify convergence.
Explicit depth disables automatic depth selection. A final depth check can
still reject that fixed target. `--fullauto-final-checks off` omits final
rechecks and is recorded as such; it does not make omitted checks pass.

Example with fixed depth and phi, without controls or pruning:

```bash
cpu/bin/mbs_po_mpi --method po --backend cpu \
  --particle 1 5 2 --ri 1.3116 0 --wavelength-um 0.532 \
  --fullauto 0.02 --fullauto-theta-range 170 180 \
  --max-reflections 8 --phi-points 4 --fullauto-initial 512 \
  --cutoff-profile off --fullauto-analytic-controls off \
  --output results/fixed_settings
```

## Controller sequence

The controller validates inputs and restores compatible saved work, chooses
initial budgets and the requested output grid, and selects reflection depth
unless fixed. CUDA uses a bounded kernel/phi pilot before the depth comparison;
complete production tuning follows at the selected depth.

Each new depth candidate restarts at the configured depth pilot; only resuming
the same unresolved candidate retains its expanded count. After two pilots,
the controller may try the next depth if the mean residual exceeds the budget,
or if its interval exceeds twice the remaining gap. The latter is a cost
heuristic based on an N^-1/2 forecast, not a Sobol convergence guarantee. It
may select a larger depth. The maximum depth retains the full sampling search.
Acceptance still requires the original upper confidence bound. Decisions are
logged as `depth_pilot_decision`.

Next it verifies supported symmetries, fits control coefficients on training
scrambles and validates them on different scrambles. It chooses phi and
execution settings from measured cost and variance. Unsupported controls are
disabled while the full PO calculation remains available.

Production evaluates support/check angles, depth levels and eight replicas.
Completed blocks are durable and reused. Each round estimates statistical
uncertainty, matched interpolation residuals and prefix stability. The next
round increases samples or angular support where needed. After working criteria
pass, enabled independent and physical audits run. An unresolved limit or failed
audit is not a converged result.

## Sobol and the orientation measure

An orientation is $R=R_z(\alpha)R_y(\beta)R_z(\gamma)$. The isotropic average is

$$\langle M\rangle=\frac1{8\pi^2}\int_0^{2\pi}\!d\alpha
\int_0^\pi\!\sin\beta\,d\beta\int_0^{2\pi}\!M(R,\theta)\,d\gamma.$$

This is Haar measure on SO(3), the group of rotations in three dimensions.
For full-domain Sobol points $(u,v,w)$, use $\cos\beta=1-2u$,
$\gamma=2\pi v$ and $\alpha=2\pi w$. Uniform beta without the sine weight is
incorrect. Verified symmetry can reduce the domain with the proper reconstruction.

Nested Owen scrambling randomizes a Sobol sequence while retaining its
space-filling structure. Eight independent seeds give eight complete estimates,
not eight independent points. Samples within a replica are not ordinary IID
Monte Carlo. Power-of-two prefixes allow incremental refinement and reuse.
`--sobol-offset` continues a native worker prefix without discarding prior work.

Let $Y_r$ be the complete level/control/cutoff estimate of replica $r$:

$$\bar Y=\tfrac18\sum_rY_r,\qquad
s^2=\tfrac17\sum_r(Y_r-\bar Y)^2,\qquad
h_{95}=t_{0.975,7}s/\sqrt8.$$

The Student-t factor is approximately 2.365. The combined replica preserves
covariances between coupled differences. Normalized ratios account for
numerator/denominator covariance and unresolved denominators; halved prefixes
supply an additional stability check.

Sobol is the current production default and a practical choice for the measured
PO campaigns. Its extensible prefixes and coupled differences suit irregular
integrands and local refinement. This does not establish a universal advantage
over Lebedev or Wigner rules on every shape or smooth integrand.

### Why both alpha and phi exist

Alpha rotates the particle around the illumination axis; phi changes the
observation azimuth. With full Haar alpha, extra phi samples retain the same
expectation and can reduce variance using the same traced beams. Polarized
quantities also require the correct Stokes-basis rotation. More phi can be
cheaper than retracing more orientations, so the controller chooses their cost
tradeoff. A sixfold body symmetry does not divide laboratory alpha by six.

Fullauto supplies `--haar-alpha` in PO worker recipes. `--phi-points 4` fixes
four observation azimuths. For a starting count of 32768 per seed use
`--fullauto-initial 32768`, not a competing manual orientation selector.

## Analytic work and numerical work

The diffraction of a planar polygonal beam contains
$F(\mathbf q)=\int_A e^{i\mathbf q\cdot\mathbf x}\,dA$. Constant-amplitude
polygon integration reduces to edge expressions, with stable limits/moments
near zero phase and small denominators. This is analytic diffraction at one
orientation, not a closed orientation average of the entire particle.

An analytic control $C$ with known mean $\mu_C$ uses

$$\widehat M=\overline{M-bC}+b\mu_C.$$

The coefficient $b$ is fitted independently. Internal paths, occlusion and
interference missing from $C$ remain in the full numerical residual. Averaging
Fresnel coefficients and substituting them at every reflection would not be
this estimator. In coherent PO, amplitudes are summed before forming Mueller
matrices. Incoherent PO sums beam intensities. Both support the same
intensity-valued facet control: its sampled value and known mean cancel in
expectation independently of the target's coherence. This does not add
interference between beams to the incoherent target. The separate experimental
`--analytic-backscatter` coherent-return control retains its coherent-only limit.

Facet/self-reflection and shadow components provide physical control means.
For incidence cosine $\mu$, the Fresnel coefficients $r_s,r_p$ give
$R(\mu)=(|r_s|^2+|r_p|^2)/2$. A Legendre moment is, for example,
$c_\ell=(2\ell+1)\int_0^1R(\mu)P_\ell(\mu)d\mu/2$; actual kernels also include
projection/geometric weights. Legendre addition and plane-wave expansions
integrate angular phase factors using spherical Bessel functions. Smooth
Fresnel quadrature, radial tables and parts of the edge preparation remain
numerical. Means are shared across angles and cached for compatible inputs.

`--analytic-facet-average` requests the facet control.
`--analytic-shadow-control circular` uses a circular-aperture control with a
geometry-derived area, with model differences retained in the residual.
`--analytic-mean-cache FILE` specifies its persistent mean cache.
Fullauto can choose these automatically; forcing a control is not a guarantee
that it reduces variance or setup time. Polarization controls undergo separate
validation; analytic M11 control does not analytically solve all 16 elements.

An analytic-control fallback does not reject the particle's mirror symmetry.
Fullauto rechecks `mirror-gamma` using ordinary PO at the selected phi and
reflection depth. A passing paired all-16 test retains the mirror shortcut
under the unchanged `epsilon/256` budget; automatic rejection uses full gamma.
Explicit `--fullauto-mirror off` remains off, while explicit `--mirror-gamma`
fails rather than silently changing modes if verification fails.
`analytic_fallback_mirror.json` records the recovery decision and
`mirror_verification.json` records the numerical comparison. Cached proofs
include the angular grid, phi, depth, executable, execution recipe and analytic
control coefficients. Checkpoints retain the selected mode; upgrading alone
does not switch an existing full-gamma calculation to mirror reconstruction.

## Depth and cutoff differences

For depth levels $n_0<\cdots<n_L$,

$$E[M_{n_L}]=E[M_{n_0}]+\sum_{\ell=1}^LE[M_{n_\ell}-M_{n_{\ell-1}}].$$

Both terms of a difference use identical orientations. Where supported, a trace
produces both depths. Each term is a complete coherent Mueller estimate at its
own depth, so interference terms are retained. A low-variance correction needs
fewer samples than the expensive full solution.

For incoherent PO, each term instead sums individual beam Mueller matrices at
its own depth. Convex shapes can also share the trace. Before enabling sharing,
two paired pilot seeds compare both depths against independent traces, including
all 16 entries and normalized ratios. The numerical allowance is
`min(1e-7, EPS/1024)` in double, `min(1e-6, EPS/4096)` in float.
`depth_sharing_verification.json` records the result; failed or unsupported
sharing uses separate traces. Facet intensity controls are available independently
of depth sharing, including incoherent PO.

Default PO pruning uses
$\widehat M=\overline B_N+\overline{(O-B)}_m$, where $B$ is balanced and $O$ is
uncut. The physical target remains uncut, with cutoff correction uncertainty
included in the complete replicas. Explicit `--cutoff-profile balanced` fixes a
preset and bypasses this automatic policy; it is not an alias for the corrected
estimator. `--cutoff-profile off` fixes no pruning.

### Incoherent startup and scheduling

Incoherent PO retains aperture diffraction, so its startup estimate is distinct
from the scale-invariant nonabsorbing GO estimate. Before the existing count
cap and power-of-two rounding it proposes
`clamp(max(pi*Dmax/lambda, 512*(1+sqrt(4*pi*abs(Im(m))*Dmax/lambda))), 512, 8192)`.
This is a bounded work heuristic, not an accuracy prediction. The measured
first-round work cap can reduce it; explicit `--fullauto-initial` takes priority.

Automatic GPU tracer comparison no longer skips an initial count of 512. It
benchmarks 512–2048 orientations per seed (bounded by `--max-orientations`) at
the selected phi/depth and compares CPU, resident GPU and hybrid tracing.
A different tracer requires numerical agreement and at least 5% measured time
saving. Unsupported geometry or absorption retains its existing path.

A correction-dominated error can trigger up to two short subrounds before a
costly base extension. Each selects up to 128 affected angles, including
different correction levels, and grows selected counts by 2–4 times. Its
forecast cost must remain at most 10% of the pending base extension. Eight
replicas, coupled orientations and acceptance tests remain unchanged; the
forecast cannot certify convergence. With automatic cutoffs, resumed incoherent
checkpoints retain completed direct-off blocks and use the corrected balanced
estimator for new eligible blocks. Their physical target is the same uncut
matrix; prefix assembly weights each block by its orientation count.
`cutoff_estimator_restore` logs the transition. Explicit cutoff settings remain
fixed. Saved analytic-control coefficients and depth-sharing decisions remain
unchanged on resume: this update enables facet controls for new calibrations,
not retroactively for an existing control-free sample.

The round limiter also ranks the worst interpolated output error and keeps its
support stencil together. Previously it ranked only direct-angle intervals;
an interpolated maximum could remain unchanged for many rounds while other
angles grew. This changes work ordering, not the accuracy criterion.

## Angular interpolation and acceptance

The output grid, interpolation supports and guard angles are distinct.
For fixed supports local cubic interpolation $I$ is linear in their matrix
values. The check is $\Delta=M(\theta)-I[M_{support}](\theta)$ on matched
orientations. A separate checkpointed residual budget prevents a difficult
neighbor's very large production count from automatically inflating this check.
Both sides use the same phi rule even if production phi varies by angle.

At a guard, confirmed residual structure can add support; a noisy residual
first grows its coupled sample. Newly inserted production angles use ordinary
startup budgets. Complete cached eight-seed prefixes remain reusable. The
round scheduler uses measured cost to limit new work; its 1800 worker-second
per-GPU target is advisory, not a wall-clock promise.

The first production round also uses that cost target. The empirical initial
count is an upper proposal: before reading production matrix values, the
controller chooses a power-of-two base count from measured timings, starting
at 512 when timings are missing. Complete cached eight-seed prefixes and the
correction minima are retained. `initial_production_work.json` records the
proposal, selected count and forecast. Explicit `--fullauto-initial`, a loaded
adaptive configuration and `--fullauto-local off` keep their requested start.
Later rounds use the ordinary error-driven refinement and unchanged acceptance
tests. A first-round estimate can therefore appear before the empirical count
has been reached; it is not a convergence declaration.

All 16 errors are scaled by local M11, not by a potentially zero Mij. Ratio
checks use absolute errors in Mij/M11. At each component/angle, statistical
and interpolation terms plus fixed reserves must fit EPS; take the maximum
**after** this componentwise sum. At unchecked output angles the full EPS/4
interpolation allowance is retained. Independent and physical checks have their
own recorded policies.

These are pointwise uncertainty estimates. They are not a global simultaneous
confidence bound, an error bound on PO versus Maxwell's equations, or proof
that every unsampled angular feature has been found.

## Results and recovery

Old output directories acquire a compact discovery journal on their first
resume with the current controller. The initial migration still reads the old
metadata once; later resumes read `completed_jobs.journal` and inspect only
directories missing from it. Existing matrices and completion markers are not
rewritten. The journal accelerates discovery; hashes and source dependencies
are still checked before a saved result enters an estimator.

A persistent worker that reports a recoverable RAM-budget failure is retried
once in a fresh process with exactly the same sampling. The failed attempt's
partial output is never committed, and its log is saved as
`stdout.memory-retry.log`. For the system RAM reserve guard, the controller
waits at most 60 seconds for the reserve to recover after releasing that
worker. It records `resource_wait.json` and memory-retry events. The reserve
is not lowered. Continued pressure or an unrelated numerical/CUDA failure
stops the task and preserves completed checkpoints for recovery.

Read `fullauto_status.json` for convergence and selected settings.
`intermediate_current.json` points to the last published matrix snapshot.
`execution_progress.json` reports active work, not a newly measured accuracy.
`fullauto_history.jsonl` records parameter changes and execution events.

`mueller_fullauto.dat` contains the requested output grid. Guard rows can contain
interpolation even when a separate direct check exists. Read `point_type` and
`value_source` in `mueller_fullauto_annotated.csv`; use `mueller_evaluated.dat`
for direct estimates. Empty error fields mean unavailable, not zero.
`parameters_accuracy.csv` and `interpolation_orientation_budget.csv` in each
snapshot distinguish production and residual sample budgets.

Resume the exact command with the same binary and output directory. Keep the
whole directory, including `jobs/`, configuration, source manifests and
checkpoints. An incomplete block can repeat; complete blocks are reused.
Schema 6 stores both budgets and residual phi; older checkpoint parsing does
not by itself prove cross-binary numerical compatibility.

A new binary requires an exact tested compatibility receipt before
`--fullauto-resume-compatible OLD_BINARY` can be used. Do not edit fingerprints.
`--fullauto-refine-window A B EPS` requests a documented child-window refinement
from a retained parent recipe rather than silently changing its physics.

`scripts/run_fullauto_service.sh UNIT BINARY LOG [arguments...]` installs and
starts a user systemd service for process recovery. Numerical decisions remain
native. It is an optional launcher, not a required Python audit.

## Tracing and preparation memory

`--fullauto-trace-backend auto` compares CPU, resident GPU and hybrid tracing
on identical orientations. A candidate must be at least 5% faster and satisfy
both the all-16/M11 and normalized-ratio gates: `min(1e-7, EPS/1024)` in FP64,
`min(1e-5, EPS/1024)` in FP32. See `trace_backend_calibration.json` for actual
times and errors. The benchmark uses 512–2048 orientations per seed, bounded
by `--max-orientations`, including jobs whose startup count is 512. An
explicit `cpu`, `resident` or `hybrid` choice fixes the backend.

Resident tracing supports nonabsorbing convex coherent and incoherent PO, up to 32 faces,
24 vertices per beam and depth 18. This includes 20-face droxtals. Geometry
and preparation phase use double even with FP32 diffraction. Shapes with
more than 16 faces use batches of at most 64 orientations; the total sample
count is unchanged. Queue capacity depends on face count and depth, with
live free-VRAM checks and a complete CPU retry on overflow. Incoherent mode
converts each beam to a Mueller matrix before accumulation, without inter-beam
interference. Absorbing and unsupported cases retain CPU preparation.

CPU convex tracing reuses its scratch child beam and constructs an outgoing
beam only after a facet intersection succeeds. Lossless Fresnel splitting
also reuses the incidence-independent effective squared index. These native
optimizations require no extra flags; they preserve clipping, ray order,
cutoffs, orientation samples and checkpoint estimators. They reduce CPU work
feeding CUDA, but do not guarantee a speedup for a GPU-diffraction-bound job.

Facet-edge geometry is prepared once per rotated convex particle. Clipping
classifies each vertex once per edge and avoids copying an aperture that is
wholly retained. The same geometric tolerances and intersection formulas apply;
the edge cache is rebuilt for every orientation, including size changes.

The controller reuses an execution-recipe digest only while every input is
unchanged; saved job identities retain their original format. Scheduling uses
numeric in-memory request keys and skips a leading native prefix only after
proving that every requested angle is saved for all eight seeds. This affects
cost forecasts, not matrix assembly or convergence checks. Coverage accelerators
share an 8 MiB bound across phi-specific planners and fall back to the exact
scan when full. These optimizations need no additional flags or new checkpoint
format. They shorten parts of the CPU work; they do not promise continuous
100% GPU utilization or eliminate checkpoint and validation time.

Persistent workers retain a conservative RAM estimate per physics recipe.
After 128 measured orientations, a later request may reuse a proven batch
size instead of restarting at 16. The warm start is bounded by 256, the
largest measured batch, explicit `--orientation-chunk`, and current RAM
with an additional 2x margin. `Prepared memory reuse` logs the estimate.
This stores no optical results and works independently of the beam cache.

Compatible angular packets with the same seed, prefix, depth and phi are
scheduled consecutively to reduce preparation-cache eviction. Sampling and
memory limits are unchanged; this helps only when a request's beams fit in
the cache. Short packets can remain after subtracting cached angles.
`execution_progress.json` also records phase-local committed block/work
counters and the last 16 durable blocks, with their angular bounds and sampling
parameters. Bounds may contain gaps. These are work counters, not a new error
estimate or a whole-round completion percentage.

## GO and other modes

```bash
cpu/bin/mbs_po_mpi --method go --backend cpu \
  --particle 1 100 50 --ri 1.3116 0 --wavelength-um 0.532 \
  --fullauto 0.02 --output results/column_go
```

GO analytically averages alpha and computes fixed finite solid-angle bins.
There is no PO cubic interpolation. Forward/backward delta matrices are separate
integrated masses in `mueller_forward.dat` and `mueller_backward.dat`, not
point values of a continuous phase function. GO mirror defaults off; supported
explicit reduction reconstructs $(M+PMP)/2$, $P=diag(1,1,-1,-1)$, for nonabsorbing
built-in types 1/10. See [GO mirror](GO_MIRROR.md) and [GO balanced](GO_BALANCED.md).

[Extinction-only](EXTINCTION_ONLY.md) is a separate convex CPU/Hammersley PO
workflow without an angular grid. It does not supply full coherent albedo.
