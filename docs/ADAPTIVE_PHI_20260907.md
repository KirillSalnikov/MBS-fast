# Adaptive direct phi quadrature for shared multisize total PO

Enable `MBS_PHI_ADAPTIVE=1` with the usual serial shared multisize CUDA command
and `--allow-experimental-environment`. Do **not** pass `--fft-factor`:
this mode improves the directly sampled azimuthal quadrature, rather than
interpolating an under-resolved phi field. Default off. Production unchanged.

The same prepared beams are reused for all sizes and all refinement levels.
For each size, theta and orientation chunk, compute means on 75,150,300,...
direct phi points. Test pointwise M11 relative change and the maximum absolute
change over all Mij/M11. Accept the finer result only after two consecutive
passing comparisons. Accepted theta rows leave the active set; the GPU then
evaluates only unresolved rows on the next grid. Endpoint handling uses actual
theta values, not the first/last indices of the active subset.

| Environment switch | Default | Meaning |
|---|---:|---|
| `MBS_PHI_ADAPTIVE` | 0 | Activate this shared multisize path |
| `MBS_PHI_ADAPT_MIN` | 75 | First direct phi count, integer >=16 |
| `MBS_PHI_ADAPT_MAX` | 9600 | Limit; MIN times a power of two, >=4*MIN, <=65536 |
| `MBS_PHI_ADAPT_M11_TOL` | 0.005 | Pointwise relative M11 change |
| `MBS_PHI_ADAPT_POL_TOL` | 0.0025 | Absolute normalized Mueller change |

No global forward-dominated norm or intensity masking is used. At poles the
same total-output Mueller convention is applied before comparisons. Positive
intensity denominators have only a 1e-300 numerical zero floor. Nonfinite
values, unsupported mode combinations, GPU failure, inability to write the
report, or unresolved rows at MAX cause an explicit error. There is no
unchecked fallthrough to a coarse solution. On fresh output paths a failed
scan does not write final size matrices; it retains diagnostics.

Only coherent full-only total CUDA shared multisize is currently supported.
Mirror gamma, explicit CPU scaling, disabled multi-k, FFT interpolation and
legacy global refinement are rejected. The original scattering grid and FFT
state are restored after each adaptive size calculation, including exceptions.
Nonuniform GPU grids are uploaded explicitly, avoiding signature collisions
between different active theta subsets.

`<result>_phi_convergence.csv` records size index, orientation chunk, GPU batch,
theta, selected phi count, the final two-grid changes and convergence status.
The default two-success rule also requires the preceding comparison to pass;
the report contains the last comparison, not an error bound against infinity.
The final Mueller output format, theta grid, normalization, optical-theorem
extinction and integral reporting remain the same.

## Limits of the claim

Two nested-grid comparisons are an error indicator, not a rigorous quadrature
error bound. Validation uses separate direct 4096 and 8192 grids (outside the
75*2^k sequence). Their agreement is checked independently. This work controls
phi sampling only; it does not establish convergence in theta, orientation
count, ensemble size, internal reflection count or beam cutoffs.

The first implementation recomputes each nested grid and evaluates sizes
separately to allow their theta sets to diverge. It therefore prioritizes
accuracy and explicit failure over maximum speed. Reusing old phi points and
grouping compatible size/theta sets can be investigated separately.

Regression: `tests/regression_adaptive_phi.py` on an idle GPU checks two
orientation chunks, multiple sizes, absorption, full and cropped theta ranges,
and a deliberately insufficient phi limit. Benchmark artifacts live in
`../adaptive_phi_bench_20260907/` and `../adaptive_phi_regression_20260907/`.

## Independent dense-grid benchmark (2026-09-07)

The benchmark holds all non-phi settings fixed: lambda=.532 um, m=1.4+0i,
max reflections=8, 8x8 orientations, 181 theta points, campaign beam cutoffs,
and four radii 4,35,67,100 um. Controls use direct FP64 grids with 4096 and
8192 phi points. The old 600-point reference uses the same settings and the
original production binary; prior regressions establish implementation
equivalence of its diffraction to the optimized direct kernel.

For p0000, maximum pointwise M11 error against direct 8192 was:

| Radius, um | Old direct 600 | Adaptive |
|---|---:|---:|
| 4 | 1.13e-10% | 7.17e-11% |
| 35 | 0.499% | 1.41e-7% |
| 67 | 11.80% | 1.42e-5% |
| 100 | 29.46% | 0.002663% |

At 100 um the maximum absolute Mij/M11 difference falls from 0.1489 to
0.00004722. Direct 4096 and 8192 agree to below 2.5e-10% in M11 across these
four sizes. Selected adaptive phi counts range from 300 to 9600.

Paired diffraction times for p0000: adaptive 31.95 s, direct 4096 24.11 s,
direct 8192 48.09 s. This adaptive implementation is slower than the sufficient
4096 grid for this particular shape; its benefit is automatic error checking
and avoiding a universally prescribed dense grid, not an unconditional speedup.

For p0070 (148 facets), the same comparisons give:

| Radius, um | Old direct 600 | Adaptive |
|---|---:|---:|
| 4 | 1.81e-10% | 8.48e-11% |
| 35 | 0.05141% | 0.003491% |
| 67 | 2.034% | 0.006070% |
| 100 | 10.564% | 0.003149% |

Maximum absolute normalized Mueller difference across all four sizes is
0.0002768. Direct 4096 vs 8192 maximum M11 difference is below 2.7e-10%.
Selected phi counts range from 300 to 4800. Diffraction times: adaptive
49.35 s, direct 4096 51.20 s, direct 8192 102.58 s. These are single paired
probes, not a full-ensemble throughput forecast. Both adaptive cases are
more expensive than the under-resolved 600-point direct calculation.

Additional regression passed on p0070 with 16 orientations in two chunks,
three sizes and 19 theta rows. For nonabsorbing full-angle output, maximum
relative M11 error was 1.76e-12. For m=1.4+0.01i and cropped theta 20..160 deg,
it was 2.48e-5 (0.00248%), with normalized Mueller error 4.80e-5. The intentional
MAX=300 failure cases returned nonzero without final size matrices. The local
`tests/adaptive_phi_config.cpp` test also passed invalid-parameter, minimum-limit
and inactive-mode checks. Artificial CUDA OOM was not forced on the shared GPU.

Validated build SHA256:
`e058249d0a6cd7c805fd9488dbce1a36d019b01d17664473d529f58b14461582`.
