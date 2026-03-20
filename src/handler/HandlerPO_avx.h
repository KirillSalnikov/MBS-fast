#pragma once
// AVX2 vectorized diffraction: process 4 directions simultaneously
#include <immintrin.h>
#include <cmath>
#include "compl.hpp"

// Vectorized sincos for 4 doubles using AVX2
// Each element computed via scalar sincos (no vector sincos in glibc)
// But arithmetic between sincos calls is vectorized
inline void sincos_4x(__m256d phases, __m256d &sin_out, __m256d &cos_out)
{
    alignas(32) double p[4], s[4], c[4];
    _mm256_store_pd(p, phases);
    sincos(p[0], &s[0], &c[0]);
    sincos(p[1], &s[1], &c[1]);
    sincos(p[2], &s[2], &c[2]);
    sincos(p[3], &s[3], &c[3]);
    sin_out = _mm256_load_pd(s);
    cos_out = _mm256_load_pd(c);
}

// Process 4 directions at once for one beam
// Returns 4 complex fresnel values (8 doubles: re0,im0,re1,im1,...)
inline void diffract_4dir(
    const double *vx, const double *vy, int nv,
    const double *edge_slope_yx, const bool *edge_valid_x,
    const double *edge_slope_xy, const bool *edge_valid_y,
    double horAx, double horAy, double horAz,
    double verAx, double verAy, double verAz,
    double beamDx, double beamDy, double beamDz,
    const double *dirx4, const double *diry4, const double *dirz4,
    double area, double waveIndex, double wi2,
    double eps1, double eps2,
    double complWaveRe, double complWaveIm,
    double invComplWaveRe, double invComplWaveIm,
    double *out_re, double *out_im)  // 4 complex outputs
{
    __m256d vdx = _mm256_load_pd(dirx4);
    __m256d vdy = _mm256_load_pd(diry4);
    __m256d vdz = _mm256_load_pd(dirz4);

    // k_k0 = -dir + beamDir (4 directions at once)
    __m256d vbd_x = _mm256_set1_pd(beamDx);
    __m256d vbd_y = _mm256_set1_pd(beamDy);
    __m256d vbd_z = _mm256_set1_pd(beamDz);
    __m256d kkx = _mm256_sub_pd(vbd_x, vdx);
    __m256d kky = _mm256_sub_pd(vbd_y, vdy);
    __m256d kkz = _mm256_sub_pd(vbd_z, vdz);

    // A = kk · horAxis, B = kk · verAxis
    __m256d vhx = _mm256_set1_pd(horAx), vhy = _mm256_set1_pd(horAy), vhz = _mm256_set1_pd(horAz);
    __m256d vvx = _mm256_set1_pd(verAx), vvy = _mm256_set1_pd(verAy), vvz = _mm256_set1_pd(verAz);

    __m256d A = _mm256_add_pd(_mm256_add_pd(_mm256_mul_pd(kkx, vhx), _mm256_mul_pd(kky, vhy)), _mm256_mul_pd(kkz, vhz));
    __m256d B = _mm256_add_pd(_mm256_add_pd(_mm256_mul_pd(kkx, vvx), _mm256_mul_pd(kky, vvy)), _mm256_mul_pd(kkz, vvz));

    __m256d vwi = _mm256_set1_pd(waveIndex);
    __m256d kA = _mm256_mul_pd(vwi, A);
    __m256d kB = _mm256_mul_pd(vwi, B);

    // Compute vertex phases and sincos for all 4 directions × all vertices
    // For each vertex j: phase[dir][j] = kA[dir]*vx[j] + kB[dir]*vy[j]
    // Store as vc[dir][j], vs[dir][j]

    // Using scalar fallback for sincos but vectorized phase computation
    __m256d sr = _mm256_setzero_pd(), si = _mm256_setzero_pd(); // accumulator

    __m256d absA = _mm256_andnot_pd(_mm256_set1_pd(-0.0), A);
    __m256d absB = _mm256_andnot_pd(_mm256_set1_pd(-0.0), B);
    __m256d veps2 = _mm256_set1_pd(eps2);

    // Check forward scattering for each direction
    __m256d fwd_mask = _mm256_and_pd(_mm256_cmp_pd(absA, veps2, _CMP_LT_OQ),
                                      _mm256_cmp_pd(absB, veps2, _CMP_LT_OQ));

    // For non-forward directions: compute edge contributions
    // Use absB > absA branch (simpler, handle all 4 dirs the same way)
    // This may be suboptimal for some dirs but avoids branch divergence

    // Precompute vertex sin/cos for all 4 directions
    alignas(32) double vc_arr[4*32], vs_arr[4*32]; // [dir*nv + vertex]

    for (int j = 0; j < nv; ++j)
    {
        __m256d vvxj = _mm256_set1_pd(vx[j]);
        __m256d vvyj = _mm256_set1_pd(vy[j]);
        __m256d phase = _mm256_add_pd(_mm256_mul_pd(kA, vvxj), _mm256_mul_pd(kB, vvyj));

        __m256d sv, cv;
        sincos_4x(phase, sv, cv);
        _mm256_store_pd(&vc_arr[j*4], cv);
        _mm256_store_pd(&vs_arr[j*4], sv);
    }

    // Process edges: for each edge, accumulate (exp[next]-exp[i])/Ci for all 4 directions
    __m256d useB = _mm256_cmp_pd(absB, absA, _CMP_GT_OQ); // mask: which dirs use B branch

    for (int i = 0; i < nv; ++i)
    {
        int inext = (i + 1 < nv) ? i + 1 : 0;

        // Load vertex exp_im for this edge (4 directions each)
        __m256d ci = _mm256_load_pd(&vc_arr[i*4]);
        __m256d si_v = _mm256_load_pd(&vs_arr[i*4]);
        __m256d cn = _mm256_load_pd(&vc_arr[inext*4]);
        __m256d sn = _mm256_load_pd(&vs_arr[inext*4]);

        // diff_c = cn - ci, diff_s = sn - si
        __m256d dc = _mm256_sub_pd(cn, ci);
        __m256d ds = _mm256_sub_pd(sn, si_v);

        // For B branch: Ci = A + slope_yx * B
        // For A branch: Ei = A * slope_xy + B
        double ai = edge_slope_yx[i];
        double ci_slope = edge_slope_xy[i];

        __m256d vai = _mm256_set1_pd(ai);
        __m256d vci = _mm256_set1_pd(ci_slope);

        __m256d Ci = _mm256_add_pd(A, _mm256_mul_pd(vai, B));  // B branch
        __m256d Ei = _mm256_add_pd(_mm256_mul_pd(vci, A), B);  // A branch

        // Select denominator based on branch
        __m256d denom = _mm256_blendv_pd(Ei, Ci, useB);

        // Avoid division by zero
        __m256d veps1 = _mm256_set1_pd(eps1);
        __m256d abs_denom = _mm256_andnot_pd(_mm256_set1_pd(-0.0), denom);
        __m256d safe_denom = _mm256_max_pd(abs_denom, veps1);
        __m256d sign_denom = _mm256_and_pd(denom, _mm256_set1_pd(-0.0)); // sign bit
        safe_denom = _mm256_or_pd(safe_denom, sign_denom); // restore sign

        // Skip edges with too-small dx/dy
        bool valid_x = edge_valid_x[i];
        bool valid_y = edge_valid_y[i];
        if (!valid_x && !valid_y) continue;

        // (diff_c / denom, diff_s / denom)
        __m256d inv_d = _mm256_div_pd(_mm256_set1_pd(1.0), safe_denom);
        __m256d contrib_r = _mm256_mul_pd(dc, inv_d);
        __m256d contrib_i = _mm256_mul_pd(ds, inv_d);

        // Zero out contribution where |denom| < eps1
        __m256d valid_mask = _mm256_cmp_pd(abs_denom, veps1, _CMP_GE_OQ);
        contrib_r = _mm256_and_pd(contrib_r, valid_mask);
        contrib_i = _mm256_and_pd(contrib_i, valid_mask);

        sr = _mm256_add_pd(sr, contrib_r);
        si = _mm256_add_pd(si, contrib_i);
    }

    // Divide by B (for B branch) or -A (for A branch)
    __m256d final_denom = _mm256_blendv_pd(_mm256_sub_pd(_mm256_setzero_pd(), A), B, useB);
    __m256d inv_final = _mm256_div_pd(_mm256_set1_pd(1.0), final_denom);
    sr = _mm256_mul_pd(sr, inv_final);
    si = _mm256_mul_pd(si, inv_final);

    // Multiply by complWave = (cwr + i*cwi)
    // result = (cwr + i*cwi) * (sr + i*si) = (cwr*sr - cwi*si) + i*(cwr*si + cwi*sr)
    __m256d vcwr = _mm256_set1_pd(complWaveRe);
    __m256d vcwi = _mm256_set1_pd(complWaveIm);
    __m256d rr = _mm256_sub_pd(_mm256_mul_pd(vcwr, sr), _mm256_mul_pd(vcwi, si));
    __m256d ri = _mm256_add_pd(_mm256_mul_pd(vcwr, si), _mm256_mul_pd(vcwi, sr));

    // Handle forward scattering: replace with invComplWave * area
    __m256d fwd_re = _mm256_set1_pd(invComplWaveRe * area);
    __m256d fwd_im = _mm256_set1_pd(invComplWaveIm * area);
    rr = _mm256_blendv_pd(rr, fwd_re, fwd_mask);
    ri = _mm256_blendv_pd(ri, fwd_im, fwd_mask);

    _mm256_store_pd(out_re, rr);
    _mm256_store_pd(out_im, ri);
}
