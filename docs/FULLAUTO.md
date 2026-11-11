# Fullauto: practical and mathematical guide

This guide describes the integrated code prepared for MBS-fast 2.0. For the
extended Russian explanation and worked parameter definitions, see
[FULLAUTO_RU.md](FULLAUTO_RU.md) or [HTML](FULLAUTO_ADAPTIVE_RU.html).

## Start with the physical inputs

```bash
make -C cpu -j2
cpu/bin/mbs_po_mpi --method po --backend cpu \
  --particle 1 5 2 --ri 1.3116 0 --wavelength-um 0.532 \
  --fullauto 0.02 --checkpoint --output results/column_5um
```

The type-1 column has axial length 5 and basal diameter 2, in the same units as
the wavelength. `0.02` means 2%. For a file particle, replace `--particle 1 5 2`
with `--particle-file examples/cube.particle --resize-dmax-um 5`.
Dmax is the maximum vertex separation, not a volume-equivalent diameter.

No Python controller or hand-written orientation file is required. Do not add
an explicit Sobol selector to make fullauto use Sobol. Run one controller per
output directory, without a multi-rank MPI launch.

New output directories use bounded RAM working storage by default.
`--checkpoint` enables recovery snapshots every 300 seconds; without it,
published matrices remain available but an interrupted adaptive run cannot
resume. Existing legacy outputs retain their previous storage policy.

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

CUDA layout calibration uses fresh timing requests after a warmup and repeats
the two fastest numerically valid candidates before changing the default.
The comparison includes the complete request time and checks all 16 elements
and enabled normalized ratios on identical orientations. Explicit layout
options remain fixed. Further probes are limited by a production-cost forecast
and a five-second soft budget; a running request is not interrupted. If the
budget is insufficient, the current layout remains selected. See
`kernel_calibration.json` for individual timings, numerical checks and the
fallback reason. Trial matrices are temporary and are not checkpoint sources.

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
Completed blocks are retained for reuse; their recovery lifetime depends on
the checkpoint policy described below. Each round estimates statistical
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

New runs default to `--fullauto-storage memory`. Native blocks travel over an
interprocess channel in chunks of at most 64 KiB. Requests with at most 64 MiB
of estimated inputs and outputs use temporary RAM staging when reserves allow;
larger jobs use disposable disk staging without per-file durability barriers.
Controller values
use a 512 MiB retained-cache budget; overflow values and, when needed, the
SQLite catalog use anonymous temporary files. These files are disposable,
not recovery checkpoints. The cache budget excludes other controller indexes,
workers, temporary read buffers, library code and the OS page cache. Disk spill
size is not bounded by the RAM budget. Catalog spill dynamically requires
`libsqlite3.so.0` or `libsqlite3.so`; the normal build links with `-ldl`.

`--checkpoint` enables periodic `fullauto_checkpoint.bin` snapshots, including
the separate eight-seed state and paired corrections. An explicit
`--fullauto-checkpoint-seconds N` wins, including `0` to disable snapshots.
Snapshots publish atomically at a consistent controller point; an unfinished
write cannot replace the last complete generation. They currently include all
live working records, so their size and duration grow with the stored state.
`fullauto_storage.json` reports the policy and actual memory/spill sizes.
See the [detailed storage options](FULLAUTO_RU.md#рабочее-состояние-и-периодический-чекпоинт).

Existing legacy output directories automatically keep `legacy` storage. Do not change
an existing legacy directory to `memory`; start a new output instead.
Legacy directories acquire a compact discovery journal on their first
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
stops the task. Recovery requires an enabled published memory checkpoint or
the retained legacy blocks; a DAT alone does not restore adaptive state.

Read `fullauto_status.json` for convergence and selected settings.
`intermediate_current.json` points to the last published matrix snapshot.
In memory mode this pointer changes atomically only after the complete immutable
`public_snapshots/<hash>/` generation has been written. The current and previous
generations are retained; follow the pointer instead of assuming a fixed folder.
`execution_progress.json` reports active work, not a newly measured accuracy.
`fullauto_history.jsonl` records parameter changes and execution events.

`mueller_fullauto.dat` contains the requested output grid. Guard rows can contain
interpolation even when a separate direct check exists. Read `point_type` and
`value_source` in `mueller_fullauto_annotated.csv`; use `mueller_evaluated.dat`
for direct estimates. Empty error fields mean unavailable, not zero.
`parameters_accuracy.csv` and `interpolation_orientation_budget.csv` in each
snapshot distinguish production and residual sample budgets.

With recovery enabled, resume the exact command with the same binary and
output directory. Keep the whole directory: memory mode needs its marker,
configuration and complete binary snapshot; legacy mode also needs `jobs/`
and source manifests. Work after the last memory snapshot can repeat.
Without a snapshot, restart in a new directory. Legacy schema 6 stores both
budgets and residual phi; older checkpoint parsing does
not by itself prove cross-binary numerical compatibility.
Keep the original particle file and its visibility sidecar accessible at their
original paths, including when they live inside the output directory. Inputs
are read before restoring the snapshot; it is not a portable input archive.

A new binary requires an exact tested compatibility receipt before
`--fullauto-resume-compatible OLD_BINARY` can be used. Do not edit fingerprints.
`--fullauto-refine-window A B EPS` requests a documented child-window refinement
from a retained parent recipe rather than silently changing its physics.

`scripts/run_fullauto_service.sh UNIT BINARY LOG [arguments...]` installs and
starts a user systemd service for process recovery. Numerical decisions remain
native. It is an optional launcher, not a required Python audit.
Pass `--checkpoint` to the calculation for recoverable new memory-mode runs;
restarting a service alone does not save adaptive state.

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

## Experimental: one worker with bounded resident buffers

The following opt-in combination is for separate validation runs. It is not the
simple default above and does not imply deployment to existing campaigns or a
GPU utilization guarantee of 50%. Use a CUDA binary built with these changes
(the path below is an example) and a new output directory:

```bash
env -u MBS_GPU_TRACE_CONVEX_WORKSPACE_MIB \
  MBS_GPU_TRACE_DEEP_CAPACITY=1 \
  MBS_GPU_TRACE_CAPACITY_HINT=1 \
  MBS_GPU_TRACE_BUFFER_POOL=1 \
  MBS_GPU_RESIDENT_REUSE=1 \
  MBS_GPU_RESIDENT_BOUNDED_REUSE=1 \
  MBS_GPU_RECYCLE_PREPARED=1 \
  MBS_HOST_MEM_BUDGET_MB=4096 MBS_GPU_MEM_FRACTION=.20 \
  gpu/bin/mbs_po_gpu_float_ds \
    --method po --backend cuda --incoherent \
    --particle 1 501.2 155.8 --wavelength-um 0.532 --ri 1.3116 0 \
    --symmetry 2 6 --fullauto 0.04 --fullauto-theta-range 150 160 \
    --fullauto-gpus 0 --fullauto-workers persistent --threads 8 \
    --po-trace-backend resident --allow-experimental-environment \
    --output results/column501_resident_experiment --close
```

This is one fullauto calculation: one controller and one persistent compute
worker on GPU 0, not several independent jobs. `0.04` requests 4% accuracy;
final checks remain enabled. The column dimensions are L=501.2 and D=155.8 µm.
Remove `--incoherent` for coherent PO. The explicit resident backend remains
subject to its support checks. The command above deliberately leaves resident
facet control disabled; ordinary facet-control stages retain CPU preparation.

| Setting | Effect |
|---|---|
| `MBS_GPU_TRACE_DEEP_CAPACITY=1` | Supports depth-22 reference tracing with bounded queues and at most two capacity retries. Does not fix the selected physical depth. |
| `MBS_GPU_TRACE_CAPACITY_HINT=1` | Reuses successful queue-size estimates for the same recipe; requires deep capacity. Does not reuse optical results. |
| `MBS_GPU_TRACE_BUFFER_POOL=1` | Reuses trace allocations, including scratch. Retention is at most 2 GiB and one quarter of current free VRAM plus this pool's owned buffers. |
| `MBS_GPU_RESIDENT_REUSE=1` and `MBS_GPU_RESIDENT_BOUNDED_REUSE=1` | Reuse resident preparation/diffraction allocations with a separate limit of at most 1 GiB and dynamic memory bounds. Explicit `MBS_GPU_RESIDENT_REUSE=0` wins. |
| Optional `MBS_GPU_RESIDENT_FACET_CONTROL=1` | For a separate large/deep-case experiment only: allows ordinary analytic facet control with resident tracing. Controls stay on the persistent consumer thread; GPU stages are ordered, including after a CPU fallback. |
| `MBS_GPU_RECYCLE_PREPARED=1` | Returns completed host blocks to their preparation thread for destruction, preserving the two-block limit. |

The 4096 MiB host budget and `.20` GPU fraction are test settings, not universal
values or hard limits on total RSS/VRAM. The pools bound their own allocations;
CUDA context, matrices, analytic scratch and other program storage also use memory.

For a separate large/deep-case comparison, add
`MBS_GPU_RESIDENT_FACET_CONTROL=1` to the `env` assignments above. Do not enable
it across all jobs: a shallow-depth fullauto test selected n=6 and became
**slower**, from 19.582 to 28.085 s in incoherent mode and 24.719 to 27.144 s in
coherent mode. Each comparison performed the same 160 native jobs and 236,448
orientations. The resident cap increased packet count from 569 to 2,400.
Numerical comparisons passed, with a maximum all-element difference of
`4.125e-7 * M11`. These measured regressions rule out recommending this option
as a blanket speedup; backend selection needs a representative timed comparison.

The facet estimator is unchanged: sampled reflection/shadow controls use the
same orientation weights, and the final M11 correction remains the known mean
minus that sample, with the existing theta-dependent coefficients. Paired depth
and prefix outputs retain their own matrix normalization. Backscatter controls,
per-orientation facet sample diagnostics, Gaussian azimuth controls, hybrid
tracing, absorption and nonconvex/unsupported geometry do not gain this resident
facet route; they retain existing CPU preparation and compatibility checks.
Check `accepted_GPU`, `CPU_fallback` and `GPU resident facet control` in the log.

With the trace pool enabled, the **implicit GPU raw-beam cache is zero**, while
persistent CPU `PreparedReuse` keeps its previous 256 MiB default and RAM clamps.
`Trace caches` prints both effective limits. Explicit `--trace-cache-mib N`,
including zero, still sets **both** caches subject to the existing RAM clamps.
`MBS_GPU_TRACE_CONVEX_WORKSPACE_MIB=0` disables trace-buffer retention; it does
not turn the raw-beam cache back on while `MBS_GPU_TRACE_BUFFER_POOL=1` is selected.
With that pool flag off, the implicit cache policy is unchanged.

New runs still default to memory working storage. Add `--checkpoint` for
300-second recovery snapshots; explicit `--fullauto-checkpoint-seconds N`
takes priority, including zero. Without a snapshot, an interrupted adaptive run
cannot resume from a published matrix alone. Recovery requires the same command,
binary, output and accessible particle inputs; work after the last snapshot may
repeat.

Validation so far includes four matched native GPU tests, coherent and
incoherent, with 512 orientations, five theta rows, 32 phi points, depth 22,
paired depth 20 and a 256-orientation prefix. All 16 matrix elements, normalized
ratios, controls and integrals passed comparison against CPU tracing plus GPU
diffraction. The largest matrix difference was `8.88e-7 * M11`; sampled-control
differences were about `1.5e-17 * M11`. These are numerical agreement tests for
that fixture, not proof of fullauto completion time or of speed on other shapes.
Whole-device `nvidia-smi` utilization on a shared GPU does not measure this
worker's own kernel activity.

### Depth reference tolerance

The paired deep-reference check uses the same tolerance as candidate depth
selection: `epsilon / 4`, rather than `epsilon / 16`. For `--fullauto 0.04`
both thresholds are 1%. The overall requested tolerance remains 4%, and the
default candidate range remains 6–18. The 20/22 reference depths do not become
production depths. Reports record the actual `reference_budget`. A candidate
or reference whose uncertainty exceeds its threshold remains unresolved or
fails; this change does not automatically mark previous failures converged.

### Continue at the maximum depth after an unresolved pilot

If automatic depth selection cannot accept a candidate, fullauto continues
at the configured `--fullauto-depth-max` (18 by default). It does not increase
that limit. The selected production depth and the last tested pilot depth
are recorded separately; `pilot_passed` remains false. Orientation and angular
refinement continue normally. Explicit `-n` still bypasses this selection.

Status, final parameters and the precision contract report
`depth_selection_fallback: true` and `accuracy_conditional_on_fixed_depth: true`.
A converged status in this case describes averaging at that fixed depth; it
does not certify the truncation error from omitted internal reflections.
Requested final physical checks remain active and can still reject a result.
A verified fallback receipt restores the same depth and warning on resume,
without repeating the unsuccessful selection. Missing or modified receipts
cannot be used to alter an existing production checkpoint.


Компактный каталог PO/GO: [`--fullauto-layout compact`](FULLAUTO_COMPACT_RESULTS_RU.md) — результаты через `latest`, рабочие данные в `work/`, восстановление и проверяемая упаковка старых снимков.


## Default completion policy

PO and GO fullauto now require one successful production pass and omit final
rechecks by default. Enable extended final audits explicitly with
`--fullauto-final-checks on`; this does not implicitly add a second production
pass. Explicit `--convergence-passes N` or `controller.stable_passes` still wins.
An adaptive configuration without that key leaves the fullauto default at one.
Explicit `--fullauto-final-audit on` or `--fullauto-physics-audit on` also opts in
unless the master switch is explicitly off. Startup depth selection and all
production convergence criteria remain enabled; omitted audits remain unverified.
Other adaptive modes retain their existing two-pass default.
