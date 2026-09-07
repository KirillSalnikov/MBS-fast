# Compact phi moments for total shared multisize output

This is an opt-in implementation optimization, not a new phi-convergence rule.
Enable `MBS_FFT_PHI_AVERAGE_ONLY=1` alongside the existing `--fft-factor F`
and `--allow-experimental-environment` in the serial shared multisize total
CUDA calculation. Leave the flag absent (or set to 0) to retain interpolation.

## Why it preserves the FFT result

The CUDA phi FFT interpolates **Mueller matrices**, after coherent Jones
summation and Jones-to-Mueller conversion. `HandlerPOTotal::WriteMatricesToFile`
then averages `M(phi) * L(phi)`, except at the poles where it averages `M(phi)`
before applying the pole convention. The right Stokes rotation contains only
constant terms and sin/cos(2 phi). Thus only Fourier moments 0 and +/-2 are
observable in the final output. The existing FFT keeps these moments unchanged.

The new branch retains the same direct phi grid selected by the existing CUDA
FFT factor logic. It executes direct diffraction and accumulation on that grid
and lets the existing total writer perform the weighted moments. No Fourier
transform, dense inverse transform, host/device FFT transfers or dense phi
accumulation is needed. This is algebraically equivalent to averaging the
interpolant; it does not reconstruct or expose a reduced-harmonic phi field.
The smallest reduced grid is 16, keeping mode 2 away from Nyquist.

Examples: requested 600/factor 2 uses 300 direct phi points; requested
300/factor 4 uses 75, exactly as before. The latter is an odd grid and is tested.
All theta points, sizes, orientations, coherent beams, reflections and cutoffs
are unchanged. FP64 fused multi-k diffraction remains available.

## Scope and safety

Only `TracerPOTotal::TraceRandomMultiSize` activates this flag. It requires the
total averaged writer, coherent full-only CUDA FFT, no mirror gamma, no theta
FFT, no requested FFT tolerance/check, and no adaptive/global refinement.
Unsupported combinations in this path fail explicitly instead of silently
removing checks. Other tracing modes are unchanged and do not use this shortcut.
The original scattering sphere, arrays, and FFT state are restored after the
calculation, including an exception exit.

The flag is default-off. It has NOT been enabled in the production queue.
The standalone optimized binary is still isolated from the production binary.

## Validation

`tests/regression_fft_phi_average.py` compares every output Mueller component
on p0070 at three sizes (4,35,100 um), multiple orientation chunks, odd/even
direct phi grids, and refractive indices 1.4+0i and 1.4+0.01i. It tests factor 1
as a direct fallback edge case and verifies that four incompatible refinement
switches are rejected. Tests passed: for actual FFT interpolation, maximum
all-element error / M11 was 7.96e-15; factor-1 comparison against the existing
per-size fallback was 1.16e-12 or less.

Separate 8x8-orientation benchmarks use p0000 and p0070, radii 4,35,67,100 um,
181 theta points, the campaign cutoffs, and both 600/factor2 and 300/factor4.
All 16 paired matrices passed; maximum all-element error / M11 was 1.64e-14.
These are implementation equivalence tests, not physical convergence tests.

The production-shaped 32-size block (p0000, radii 4..35 um, 64 orientations,
181 theta points) gave the following paired diffraction times on a V100:

| Requested grid / factor | Fused FFT interpolation | Compact mean | Time reduction |
|---|---:|---:|---:|
| 600 / 2 | 16.63 s | 14.99 s | 9.9% |
| 300 / 4 | 5.23 s | 4.42 s | 15.5% |

Both controls use the previously optimized FP64 multi-k/streaming kernel.
This isolates the additional benefit of skipping interpolation. All 64 paired
size matrices passed, with maximum all-element error / M11 of 1.66e-14.
These are single paired timing probes, not a full-campaign throughput estimate.

Artifacts: `../fft_average_bench_20260907/`, `../fft_average32_bench_20260907/`,
and `../fft_average_regression_20260907/`.

Build SHA256:
`9187e1dfdbc3fb8f9668e2698a80a1b6e6ac3140c90256fb7f1f9aa840be7739`.

## Remaining accuracy problem

Skipping interpolation does not fix aliasing in the directly sampled phi grid.
Prior tests at large radii already showed significant side-angle differences
between coarse FFT sampling and direct 600, and between direct 600 and 1200.
Consequently no coarse FFT setting is recommended for the full campaign yet.
Adaptive per-size/per-theta sampling with nested grids and a dense independent
control is a separate next step; it is not implemented by this shortcut.
