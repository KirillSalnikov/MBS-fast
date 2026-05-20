#include "HandlerPO.h"

#include <cuda_runtime.h>
#include <cufft.h>
#include <cmath>
#include <algorithm>
#include <cstdlib>
#include <cstdio>
#include <vector>

#ifdef MBS_GPU_FLOAT
using GpuReal = float;
using GpuComplex = cufftComplex;
static const cufftType GpuCufftType = CUFFT_C2C;
static inline cufftResult gpu_cufft_exec(cufftHandle plan, GpuComplex *in,
                                         GpuComplex *out, int direction)
{
    return cufftExecC2C(plan, in, out, direction);
}
#else
using GpuReal = double;
using GpuComplex = cufftDoubleComplex;
static const cufftType GpuCufftType = CUFFT_Z2Z;
static inline cufftResult gpu_cufft_exec(cufftHandle plan, GpuComplex *in,
                                         GpuComplex *out, int direction)
{
    return cufftExecZ2Z(plan, in, out, direction);
}
#endif

__device__ inline void gpu_sincos(GpuReal x, GpuReal *s, GpuReal *c)
{
#ifdef MBS_GPU_FLOAT
    sincosf(x, s, c);
#else
    sincos(x, s, c);
#endif
}

struct GpuBeam
{
    GpuReal x[32], y[32];
    GpuReal slope_yx[32], slope_xy[32];
    unsigned char edge_valid_x[32], edge_valid_y[32];
    int nVertices;
    int isExternal;

    GpuReal bdx, bdy, bdz;
    GpuReal horAx, horAy, horAz;
    GpuReal verAx, verAy, verAz;
    GpuReal cenx, ceny, cenz;
    GpuReal beam_area;
    GpuReal pNTx, pNTy, pNTz;
    GpuReal pNPx, pNPy, pNPz;
    GpuReal pnxDTx, pnxDTy, pnxDTz;
    GpuReal pnxDPx, pnxDPy, pnxDPz;
    GpuReal jp00r, jp00i, jp01r, jp01i;
    GpuReal jp10r, jp10i, jp11r, jp11i;
    int orientation;
};

struct GpuWorkspace
{
    GpuBeam *beams = nullptr;
    GpuReal *sinTheta = nullptr, *cosTheta = nullptr;
    GpuReal *sinPhi = nullptr, *cosPhi = nullptr;
    GpuReal *vf = nullptr;
    GpuReal *j = nullptr, *jNoShadow = nullptr;
    GpuReal *weights = nullptr;
    int *beamOffsets = nullptr;
    GpuReal *m = nullptr, *mNoShadow = nullptr;
    size_t beamCap = 0;
    size_t sinThetaCap = 0, cosThetaCap = 0;
    size_t sinPhiCap = 0, cosPhiCap = 0;
    size_t vfCap = 0, jCap = 0, jNoShadowCap = 0;
    size_t weightsCap = 0, beamOffsetsCap = 0, mCap = 0, mNoShadowCap = 0;
    int gridNAz = -1, gridNZen = -1;
    std::vector<GpuBeam> hBeams;
    std::vector<GpuReal> hWeights;
    std::vector<int> hBeamOffsets;
    std::vector<GpuReal> hM;
    std::vector<GpuReal> hMNoShadow;
};

static GpuWorkspace g_gpuWorkspace;

static GpuReal gpu_memory_fraction()
{
    const char *value = std::getenv("MBS_GPU_MEM_FRACTION");
    if (!value || !*value)
        return 0.88;
    char *end = nullptr;
    GpuReal parsed = std::strtod(value, &end);
    if (!end || *end != '\0' || parsed <= 0.1 || parsed > 0.98)
        return 0.88;
    return parsed;
}

static bool gpu_no_atomics_enabled()
{
    const char *value = std::getenv("MBS_GPU_NO_ATOMICS");
    return value && value[0] == '1' && value[1] == '\0';
}

static int gpu_fft_phi_factor()
{
    const char *value = std::getenv("MBS_FFT_PHI_FACTOR");
    if (!value || !*value)
        return 0;
    if ((value[0] == 'a' || value[0] == 'A') &&
        (value[1] == 'u' || value[1] == 'U') &&
        (value[2] == 't' || value[2] == 'T') &&
        (value[3] == 'o' || value[3] == 'O') &&
        value[4] == '\0')
        return 0;
    char *end = nullptr;
    long parsed = std::strtol(value, &end, 10);
    if (!end || *end != '\0' || parsed < 1 || parsed > 64)
        return 0;
    return (int)parsed;
}

static bool gpu_fft_check_enabled()
{
    const char *value = std::getenv("MBS_FFT_CHECK");
    return value && value[0] == '1' && value[1] == '\0';
}

static bool gpu_report_cuda_error(cudaError_t err, const char *where)
{
    if (err == cudaSuccess)
        return true;
    std::fprintf(stderr, "CUDA backend error at %s: %s\n",
                 where, cudaGetErrorString(err));
    return false;
}

static int choose_fft_phi_factor(int nFull)
{
    int requested = gpu_fft_phi_factor();
    if (requested > 0)
        return requested;
    if (nFull >= 2400)
        return 8;
    if (nFull >= 600)
        return 10;
    if (nFull >= 240)
        return 3;
    return 2;
}

template <typename T>
static bool ensure_device_capacity(T *&ptr, size_t &cap, size_t count)
{
    if (cap >= count)
        return true;
    cudaFree(ptr);
    ptr = nullptr;
    cap = 0;
    if (count == 0)
        return true;
    if (cudaMalloc(&ptr, count * sizeof(T)) != cudaSuccess)
        return false;
    cap = count;
    return true;
}

__device__ inline void cmul(GpuReal ar, GpuReal ai, GpuReal br, GpuReal bi,
                            GpuReal &cr, GpuReal &ci)
{
    cr = ar * br - ai * bi;
    ci = ar * bi + ai * br;
}

__global__ void fft_phi_pad_kernel(const GpuComplex *low,
                                   GpuComplex *full,
                                   int nLow,
                                   int nFull,
                                   int batch,
                                   GpuReal scale)
{
    int idx = (int)(blockIdx.x * blockDim.x + threadIdx.x);
    int total = nLow * batch;
    if (idx >= total) return;

    int k = idx % nLow;
    int b = idx / nLow;
    int freq = (k <= nLow / 2) ? k : k - nLow;
    int dstK = (freq >= 0) ? freq : nFull + freq;
    if (dstK < 0 || dstK >= nFull) return;

    GpuComplex v = low[(size_t)b * nLow + k];
    full[(size_t)b * nFull + dstK].x = v.x * scale;
    full[(size_t)b * nFull + dstK].y = v.y * scale;
}

__device__ inline void rotate_jones_gpu(
    GpuReal NTx, GpuReal NTy, GpuReal NTz,
    GpuReal NPx, GpuReal NPy, GpuReal NPz,
    GpuReal nxDTx, GpuReal nxDTy, GpuReal nxDTz,
    GpuReal nxDPx, GpuReal nxDPy, GpuReal nxDPz,
    GpuReal vfx, GpuReal vfy, GpuReal vfz,
    GpuReal dirx, GpuReal diry, GpuReal dirz,
    GpuReal &r00, GpuReal &r01, GpuReal &r10, GpuReal &r11)
{
    GpuReal vtx = vfy * dirz - vfz * diry;
    GpuReal vty = vfz * dirx - vfx * dirz;
    GpuReal vtz = vfx * diry - vfy * dirx;
    GpuReal vtLen2 = vtx * vtx + vty * vty + vtz * vtz;
    if (vtLen2 > 1e-30)
    {
        GpuReal invLen = rsqrt(vtLen2);
        vtx *= invLen; vty *= invLen; vtz *= invLen;
    }

    GpuReal dot_dir_nxDT = dirx * nxDTx + diry * nxDTy + dirz * nxDTz;
    GpuReal cpTx = (diry * NTz - dirz * NTy) + nxDTx - dirx * dot_dir_nxDT;
    GpuReal cpTy = (dirz * NTx - dirx * NTz) + nxDTy - diry * dot_dir_nxDT;
    GpuReal cpTz = (dirx * NTy - diry * NTx) + nxDTz - dirz * dot_dir_nxDT;

    GpuReal dot_dir_nxDP = dirx * nxDPx + diry * nxDPy + dirz * nxDPz;
    GpuReal cpPx = (diry * NPz - dirz * NPy) + nxDPx - dirx * dot_dir_nxDP;
    GpuReal cpPy = (dirz * NPx - dirx * NPz) + nxDPy - diry * dot_dir_nxDP;
    GpuReal cpPz = (dirx * NPy - diry * NPx) + nxDPz - dirz * dot_dir_nxDP;

    r00 = cpTx * vtx + cpTy * vty + cpTz * vtz;
    r01 = cpPx * vtx + cpPy * vty + cpPz * vtz;
    r10 = cpTx * vfx + cpTy * vfy + cpTz * vfz;
    r11 = cpPx * vfx + cpPy * vfy + cpPz * vfz;
}

__device__ inline bool compute_beam_jones_gpu(const GpuBeam &b,
                                              int grid,
                                              int nZen,
                                              const GpuReal *sinTheta,
                                              const GpuReal *cosTheta,
                                              const GpuReal *sinPhi,
                                              const GpuReal *cosPhi,
                                              const GpuReal *vf,
                                              GpuReal waveIndex,
                                              GpuReal wi2,
                                              GpuReal eps1,
                                              GpuReal eps2,
                                              GpuReal complWaveR,
                                              GpuReal complWaveI,
                                              GpuReal invComplWaveR,
                                              GpuReal invComplWaveI,
                                              int legacySign,
                                              int singularCorrection,
                                              GpuReal &d00r, GpuReal &d00i,
                                              GpuReal &d01r, GpuReal &d01i,
                                              GpuReal &d10r, GpuReal &d10i,
                                              GpuReal &d11r, GpuReal &d11i)
{
    int nv = b.nVertices;
    if (nv <= 0) return false;

    int p = grid / (nZen + 1);
    int t = grid - p * (nZen + 1);
    GpuReal cp = cosPhi[p], sp = sinPhi[p];
    GpuReal sin_t = sinTheta[t], cos_t = cosTheta[t];
    GpuReal dx = sin_t * cp;
    GpuReal dy = sin_t * sp;
    GpuReal dz = -cos_t;

    GpuReal neg_cp = -cp, neg_sp = -sp;
    GpuReal a_sin = neg_cp * b.horAx + neg_sp * b.horAy;
    GpuReal a_cos = b.horAz;
    GpuReal a0 = b.bdx * b.horAx + b.bdy * b.horAy + b.bdz * b.horAz;
    GpuReal b_sin = neg_cp * b.verAx + neg_sp * b.verAy;
    GpuReal b_cos = b.verAz;
    GpuReal b0 = b.bdx * b.verAx + b.bdy * b.verAy + b.bdz * b.verAz;

    GpuReal A = sin_t * a_sin + cos_t * a_cos + a0;
    GpuReal B = sin_t * b_sin + cos_t * b_cos + b0;
    GpuReal absA = fabs(A);
    GpuReal absB = fabs(B);

    GpuReal fr, fi;
    if (absA < eps2 && absB < eps2)
    {
        GpuReal sign = legacySign ? 1.0 : -1.0;
        fr = sign * invComplWaveR * b.beam_area;
        fi = sign * invComplWaveI * b.beam_area;
    }
    else
    {
        GpuReal vc[32], vs[32];
        for (int v = 0; v < nv; ++v)
        {
            GpuReal psin = waveIndex * (a_sin * b.x[v] + b_sin * b.y[v]);
            GpuReal pcos = waveIndex * (a_cos * b.x[v] + b_cos * b.y[v]);
            GpuReal p0 = waveIndex * (a0 * b.x[v] + b0 * b.y[v]);
            gpu_sincos(sin_t * psin + cos_t * pcos + p0, &vs[v], &vc[v]);
        }

        GpuReal sr = 0.0, si = 0.0;
        if (absB > absA)
        {
            for (int e = 0; e < nv; ++e)
            {
                if (!b.edge_valid_x[e]) continue;
                int en = (e + 1 < nv) ? e + 1 : 0;
                GpuReal Ci = A + b.slope_yx[e] * B;
                GpuReal absCi = fabs(Ci);
                GpuReal inv = (absCi > eps1) ? (1.0 / Ci) : 0.0;
                sr += (vc[en] - vc[e]) * inv;
                si += (vs[en] - vs[e]) * inv;
                if (singularCorrection && absCi <= eps1)
                {
                    GpuReal p1x = b.x[e], p2x = b.x[en];
                    GpuReal tr = -wi2 * Ci * (p2x * p2x - p1x * p1x) * 0.5;
                    GpuReal ti = waveIndex * (p2x - p1x);
                    sr += vc[e] * tr - vs[e] * ti;
                    si += vc[e] * ti + vs[e] * tr;
                }
            }
            sr /= B; si /= B;
        }
        else
        {
            for (int e = 0; e < nv; ++e)
            {
                if (!b.edge_valid_y[e]) continue;
                int en = (e + 1 < nv) ? e + 1 : 0;
                GpuReal Ei = A * b.slope_xy[e] + B;
                GpuReal absEi = fabs(Ei);
                GpuReal inv = (absEi > eps1) ? (1.0 / Ei) : 0.0;
                sr += (vc[en] - vc[e]) * inv;
                si += (vs[en] - vs[e]) * inv;
                if (singularCorrection && absEi <= eps1)
                {
                    GpuReal p1y = b.y[e], p2y = b.y[en];
                    GpuReal tr = -wi2 * Ei * (p2y * p2y - p1y * p1y) * 0.5;
                    GpuReal ti = waveIndex * (p2y - p1y);
                    sr += vc[e] * tr - vs[e] * ti;
                    si += vc[e] * ti + vs[e] * tr;
                }
            }
            GpuReal inv_nA = -1.0 / A;
            sr *= inv_nA; si *= inv_nA;
        }
        cmul(complWaveR, complWaveI, sr, si, fr, fi);
    }
    if (isnan(fr)) return false;

    GpuReal dpr = 1.0, dpi = 0.0;
    if (!b.isExternal)
    {
        GpuReal dp_sin = cp * b.cenx + sp * b.ceny;
        GpuReal dp_cos = -b.cenz;
        gpu_sincos(-waveIndex * (sin_t * dp_sin + cos_t * dp_cos), &dpi, &dpr);
    }

    GpuReal vfx = vf[(grid * 3) + 0];
    GpuReal vfy = vf[(grid * 3) + 1];
    GpuReal vfz = vf[(grid * 3) + 2];
    GpuReal r00, r01, r10, r11;
    rotate_jones_gpu(b.pNTx, b.pNTy, b.pNTz, b.pNPx, b.pNPy, b.pNPz,
                     b.pnxDTx, b.pnxDTy, b.pnxDTz,
                     b.pnxDPx, b.pnxDPy, b.pnxDPz,
                     vfx, vfy, vfz, dx, dy, dz,
                     r00, r01, r10, r11);

    GpuReal cpr, cpi;
    cmul(fr, fi, dpr, dpi, cpr, cpi);
    GpuReal sr00r = cpr * r00, sr00i = cpi * r00;
    GpuReal sr01r = cpr * r01, sr01i = cpi * r01;
    GpuReal sr10r = cpr * r10, sr10i = cpi * r10;
    GpuReal sr11r = cpr * r11, sr11i = cpi * r11;

    d00r = sr00r * b.jp00r - sr00i * b.jp00i + sr01r * b.jp10r - sr01i * b.jp10i;
    d00i = sr00r * b.jp00i + sr00i * b.jp00r + sr01r * b.jp10i + sr01i * b.jp10r;
    d01r = sr00r * b.jp01r - sr00i * b.jp01i + sr01r * b.jp11r - sr01i * b.jp11i;
    d01i = sr00r * b.jp01i + sr00i * b.jp01r + sr01r * b.jp11i + sr01i * b.jp11r;
    d10r = sr10r * b.jp00r - sr10i * b.jp00i + sr11r * b.jp10r - sr11i * b.jp10i;
    d10i = sr10r * b.jp00i + sr10i * b.jp00r + sr11r * b.jp10i + sr11i * b.jp10r;
    d11r = sr10r * b.jp01r - sr10i * b.jp01i + sr11r * b.jp11r - sr11i * b.jp11i;
    d11i = sr10r * b.jp01i + sr10i * b.jp01r + sr11r * b.jp11i + sr11i * b.jp11r;
    return true;
}

__global__ void diffraction_kernel(const GpuBeam *beams, int nBeams,
                                   const GpuReal *sinTheta,
                                   const GpuReal *cosTheta,
                                   const GpuReal *sinPhi,
                                   const GpuReal *cosPhi,
                                   const GpuReal *vf,
                                   int nAz, int nZen,
                                   GpuReal waveIndex,
                                   GpuReal wi2,
                                   GpuReal eps1,
                                   GpuReal eps2,
                                   GpuReal complWaveR,
                                   GpuReal complWaveI,
                                   GpuReal invComplWaveR,
                                   GpuReal invComplWaveI,
                                   int legacySign,
                                   GpuReal *jFull,
                                   GpuReal *jNoShadow)
{
    int gridCount = nAz * (nZen + 1);
    long long total = (long long)nBeams * gridCount;
    long long idx = (long long)blockIdx.x * blockDim.x + threadIdx.x;
    if (idx >= total) return;

    int grid = (int)(idx % gridCount);
    int beamIdx = (int)(idx / gridCount);
    int p = grid / (nZen + 1);
    int t = grid - p * (nZen + 1);
    const GpuBeam &b = beams[beamIdx];
    int nv = b.nVertices;
    if (nv <= 0) return;

    GpuReal cp = cosPhi[p], sp = sinPhi[p];
    GpuReal sin_t = sinTheta[t], cos_t = cosTheta[t];
    GpuReal dx = sin_t * cp;
    GpuReal dy = sin_t * sp;
    GpuReal dz = -cos_t;

    GpuReal neg_cp = -cp, neg_sp = -sp;
    GpuReal a_sin = neg_cp * b.horAx + neg_sp * b.horAy;
    GpuReal a_cos = b.horAz;
    GpuReal a0 = b.bdx * b.horAx + b.bdy * b.horAy + b.bdz * b.horAz;
    GpuReal b_sin = neg_cp * b.verAx + neg_sp * b.verAy;
    GpuReal b_cos = b.verAz;
    GpuReal b0 = b.bdx * b.verAx + b.bdy * b.verAy + b.bdz * b.verAz;

    GpuReal A = sin_t * a_sin + cos_t * a_cos + a0;
    GpuReal B = sin_t * b_sin + cos_t * b_cos + b0;
    GpuReal absA = fabs(A);
    GpuReal absB = fabs(B);

    GpuReal fr, fi;
    if (absA < eps2 && absB < eps2)
    {
        GpuReal sign = legacySign ? 1.0 : -1.0;
        fr = sign * invComplWaveR * b.beam_area;
        fi = sign * invComplWaveI * b.beam_area;
    }
    else
    {
        GpuReal vc[32], vs[32];
        for (int v = 0; v < nv; ++v)
        {
            GpuReal psin = waveIndex * (a_sin * b.x[v] + b_sin * b.y[v]);
            GpuReal pcos = waveIndex * (a_cos * b.x[v] + b_cos * b.y[v]);
            GpuReal p0 = waveIndex * (a0 * b.x[v] + b0 * b.y[v]);
            gpu_sincos(sin_t * psin + cos_t * pcos + p0, &vs[v], &vc[v]);
        }

        GpuReal sr = 0.0, si = 0.0;
        if (absB > absA)
        {
            for (int e = 0; e < nv; ++e)
            {
                if (!b.edge_valid_x[e]) continue;
                int en = (e + 1 < nv) ? e + 1 : 0;
                GpuReal Ci = A + b.slope_yx[e] * B;
                GpuReal absCi = fabs(Ci);
                GpuReal inv = (absCi > eps1) ? (1.0 / Ci) : 0.0;
                sr += (vc[en] - vc[e]) * inv;
                si += (vs[en] - vs[e]) * inv;
                if (absCi <= eps1)
                {
                    GpuReal p1x = b.x[e], p2x = b.x[en];
                    GpuReal tr = -wi2 * Ci * (p2x * p2x - p1x * p1x) * 0.5;
                    GpuReal ti = waveIndex * (p2x - p1x);
                    sr += vc[e] * tr - vs[e] * ti;
                    si += vc[e] * ti + vs[e] * tr;
                }
            }
            sr /= B; si /= B;
        }
        else
        {
            for (int e = 0; e < nv; ++e)
            {
                if (!b.edge_valid_y[e]) continue;
                int en = (e + 1 < nv) ? e + 1 : 0;
                GpuReal Ei = A * b.slope_xy[e] + B;
                GpuReal absEi = fabs(Ei);
                GpuReal inv = (absEi > eps1) ? (1.0 / Ei) : 0.0;
                sr += (vc[en] - vc[e]) * inv;
                si += (vs[en] - vs[e]) * inv;
                if (absEi <= eps1)
                {
                    GpuReal p1y = b.y[e], p2y = b.y[en];
                    GpuReal tr = -wi2 * Ei * (p2y * p2y - p1y * p1y) * 0.5;
                    GpuReal ti = waveIndex * (p2y - p1y);
                    sr += vc[e] * tr - vs[e] * ti;
                    si += vc[e] * ti + vs[e] * tr;
                }
            }
            GpuReal inv_nA = -1.0 / A;
            sr *= inv_nA; si *= inv_nA;
        }
        cmul(complWaveR, complWaveI, sr, si, fr, fi);
    }
    if (isnan(fr)) return;

    GpuReal dpr = 1.0, dpi = 0.0;
    if (!b.isExternal)
    {
        GpuReal dp_sin = cp * b.cenx + sp * b.ceny;
        GpuReal dp_cos = -b.cenz;
        gpu_sincos(-waveIndex * (sin_t * dp_sin + cos_t * dp_cos), &dpi, &dpr);
    }

    GpuReal vfx = vf[(grid * 3) + 0];
    GpuReal vfy = vf[(grid * 3) + 1];
    GpuReal vfz = vf[(grid * 3) + 2];
    GpuReal r00, r01, r10, r11;
    rotate_jones_gpu(b.pNTx, b.pNTy, b.pNTz, b.pNPx, b.pNPy, b.pNPz,
                     b.pnxDTx, b.pnxDTy, b.pnxDTz,
                     b.pnxDPx, b.pnxDPy, b.pnxDPz,
                     vfx, vfy, vfz, dx, dy, dz,
                     r00, r01, r10, r11);

    GpuReal cpr, cpi;
    cmul(fr, fi, dpr, dpi, cpr, cpi);
    GpuReal sr00r = cpr * r00, sr00i = cpi * r00;
    GpuReal sr01r = cpr * r01, sr01i = cpi * r01;
    GpuReal sr10r = cpr * r10, sr10i = cpi * r10;
    GpuReal sr11r = cpr * r11, sr11i = cpi * r11;

    GpuReal d00r = sr00r * b.jp00r - sr00i * b.jp00i + sr01r * b.jp10r - sr01i * b.jp10i;
    GpuReal d00i = sr00r * b.jp00i + sr00i * b.jp00r + sr01r * b.jp10i + sr01i * b.jp10r;
    GpuReal d01r = sr00r * b.jp01r - sr00i * b.jp01i + sr01r * b.jp11r - sr01i * b.jp11i;
    GpuReal d01i = sr00r * b.jp01i + sr00i * b.jp01r + sr01r * b.jp11i + sr01i * b.jp11r;
    GpuReal d10r = sr10r * b.jp00r - sr10i * b.jp00i + sr11r * b.jp10r - sr11i * b.jp10i;
    GpuReal d10i = sr10r * b.jp00i + sr10i * b.jp00r + sr11r * b.jp10i + sr11i * b.jp10r;
    GpuReal d11r = sr10r * b.jp01r - sr10i * b.jp01i + sr11r * b.jp11r - sr11i * b.jp11i;
    GpuReal d11i = sr10r * b.jp01i + sr10i * b.jp01r + sr11r * b.jp11i + sr11i * b.jp11r;

    int off = (b.orientation * gridCount + grid) * 8;
    atomicAdd(&jFull[off + 0], d00r); atomicAdd(&jFull[off + 1], d00i);
    atomicAdd(&jFull[off + 2], d01r); atomicAdd(&jFull[off + 3], d01i);
    atomicAdd(&jFull[off + 4], d10r); atomicAdd(&jFull[off + 5], d10i);
    atomicAdd(&jFull[off + 6], d11r); atomicAdd(&jFull[off + 7], d11i);
    if (!b.isExternal)
    {
        atomicAdd(&jNoShadow[off + 0], d00r); atomicAdd(&jNoShadow[off + 1], d00i);
        atomicAdd(&jNoShadow[off + 2], d01r); atomicAdd(&jNoShadow[off + 3], d01i);
        atomicAdd(&jNoShadow[off + 4], d10r); atomicAdd(&jNoShadow[off + 5], d10i);
        atomicAdd(&jNoShadow[off + 6], d11r); atomicAdd(&jNoShadow[off + 7], d11i);
    }
}

__global__ void theta_m11_kernel(const GpuBeam *beams,
                                 int nBeams,
                                 const GpuReal *weights,
                                 const GpuReal *sinTheta,
                                 const GpuReal *cosTheta,
                                 const GpuReal *sinPhi,
                                 const GpuReal *cosPhi,
                                 const GpuReal *vf,
                                 int nAz,
                                 int nTheta,
                                 GpuReal waveIndex,
                                 GpuReal wi2,
                                 GpuReal eps1,
                                 GpuReal eps2,
                                 GpuReal complWaveR,
                                 GpuReal complWaveI,
                                 GpuReal invComplWaveR,
                                 GpuReal invComplWaveI,
                                 int legacySign,
                                 GpuReal *m11)
{
    int gridCount = nAz * nTheta;
    long long total = (long long)nBeams * gridCount;
    long long idx = (long long)blockIdx.x * blockDim.x + threadIdx.x;
    if (idx >= total) return;

    int grid = (int)(idx % gridCount);
    int thetaIdx = grid % nTheta;
    int beamIdx = (int)(idx / gridCount);

    GpuReal d00r, d00i, d01r, d01i, d10r, d10i, d11r, d11i;
    const GpuBeam &b = beams[beamIdx];
    if (!compute_beam_jones_gpu(b, grid, nTheta - 1, sinTheta, cosTheta,
                                sinPhi, cosPhi, vf, waveIndex, wi2,
                                eps1, eps2, complWaveR, complWaveI,
                                invComplWaveR, invComplWaveI, legacySign,
                                0,
                                d00r, d00i, d01r, d01i, d10r, d10i,
                                d11r, d11i))
        return;

    GpuReal value = 0.5 * (d00r*d00r + d00i*d00i
                        + d01r*d01r + d01i*d01i
                        + d10r*d10r + d10i*d10i
                        + d11r*d11r + d11i*d11i)
                 * weights[b.orientation] / (nAz + 1);
    atomicAdd(&m11[thetaIdx], value);
}

__global__ void diffraction_grid_kernel(const GpuBeam *beams,
                                        const int *beamOffsets,
                                        const GpuReal *sinTheta,
                                        const GpuReal *cosTheta,
                                        const GpuReal *sinPhi,
                                        const GpuReal *cosPhi,
                                        const GpuReal *vf,
                                        int nAz, int nZen,
                                        int nOrient,
                                        GpuReal waveIndex,
                                        GpuReal wi2,
                                        GpuReal eps1,
                                        GpuReal eps2,
                                        GpuReal complWaveR,
                                        GpuReal complWaveI,
                                        GpuReal invComplWaveR,
                                        GpuReal invComplWaveI,
                                        int legacySign,
                                        GpuReal *jFull,
                                        GpuReal *jNoShadow)
{
    int gridCount = nAz * (nZen + 1);
    long long total = (long long)nOrient * gridCount;
    long long idx = (long long)blockIdx.x * blockDim.x + threadIdx.x;
    if (idx >= total) return;

    int orient = (int)(idx / gridCount);
    int grid = (int)(idx - (long long)orient * gridCount);
    int begin = beamOffsets[orient];
    int end = beamOffsets[orient + 1];

    GpuReal j00r = 0.0, j00i = 0.0;
    GpuReal j01r = 0.0, j01i = 0.0;
    GpuReal j10r = 0.0, j10i = 0.0;
    GpuReal j11r = 0.0, j11i = 0.0;
    GpuReal n00r = 0.0, n00i = 0.0;
    GpuReal n01r = 0.0, n01i = 0.0;
    GpuReal n10r = 0.0, n10i = 0.0;
    GpuReal n11r = 0.0, n11i = 0.0;

    for (int bi = begin; bi < end; ++bi)
    {
        GpuReal d00r, d00i, d01r, d01i, d10r, d10i, d11r, d11i;
        const GpuBeam &b = beams[bi];
        if (!compute_beam_jones_gpu(b, grid, nZen, sinTheta, cosTheta,
                                    sinPhi, cosPhi, vf, waveIndex, wi2,
                                    eps1, eps2, complWaveR, complWaveI,
                                    invComplWaveR, invComplWaveI, legacySign,
                                    1,
                                    d00r, d00i, d01r, d01i, d10r, d10i,
                                    d11r, d11i))
            continue;

        j00r += d00r; j00i += d00i;
        j01r += d01r; j01i += d01i;
        j10r += d10r; j10i += d10i;
        j11r += d11r; j11i += d11i;
        if (!b.isExternal)
        {
            n00r += d00r; n00i += d00i;
            n01r += d01r; n01i += d01i;
            n10r += d10r; n10i += d10i;
            n11r += d11r; n11i += d11i;
        }
    }

    int off = ((orient * gridCount) + grid) * 8;
    jFull[off + 0] = j00r; jFull[off + 1] = j00i;
    jFull[off + 2] = j01r; jFull[off + 3] = j01i;
    jFull[off + 4] = j10r; jFull[off + 5] = j10i;
    jFull[off + 6] = j11r; jFull[off + 7] = j11i;
    jNoShadow[off + 0] = n00r; jNoShadow[off + 1] = n00i;
    jNoShadow[off + 2] = n01r; jNoShadow[off + 3] = n01i;
    jNoShadow[off + 4] = n10r; jNoShadow[off + 5] = n10i;
    jNoShadow[off + 6] = n11r; jNoShadow[off + 7] = n11i;
}

__global__ void mueller_batch_kernel(const GpuReal *jFull,
                                     const GpuReal *jNoShadow,
                                     const GpuReal *weights,
                                     int nOrient,
                                     int gridCount,
                                     GpuReal *mFull,
                                     GpuReal *mNoShadow)
{
    long long total = (long long)nOrient * gridCount;
    long long idx = (long long)blockIdx.x * blockDim.x + threadIdx.x;
    if (idx >= total) return;

    int orient = (int)(idx / gridCount);
    int grid = (int)(idx - (long long)orient * gridCount);
    GpuReal weight = weights[orient] * 0.25;
    const GpuReal *j = &jFull[idx * 8];
    GpuReal *out = &mFull[grid * 16];

    GpuReal j00r = j[0], j00i = j[1];
    GpuReal j01r = j[2], j01i = j[3];
    GpuReal j10r = j[4], j10i = j[5];
    GpuReal j11r = j[6], j11i = j[7];

    GpuReal a11 = j00r*j00r + j00i*j00i;
    GpuReal a12 = j01r*j01r + j01i*j01i;
    GpuReal a21 = j10r*j10r + j10i*j10i;
    GpuReal a22 = j11r*j11r + j11i*j11i;

    GpuReal A1 = a11 + a21, A2 = a12 + a22;
    atomicAdd(&out[0], ((A1 + A2) * 0.5) * weight);
    atomicAdd(&out[1], ((A1 - A2) * 0.5) * weight);
    A1 = a11 - a21; A2 = a12 - a22;
    atomicAdd(&out[4], ((A1 + A2) * 0.5) * weight);
    atomicAdd(&out[5], ((A1 - A2) * 0.5) * weight);

    GpuReal c1r = j00r*j01r + j00i*j01i;
    GpuReal c1i = j00i*j01r - j00r*j01i;
    GpuReal c2r = j11r*j10r + j11i*j10i;
    GpuReal c2i = j11i*j10r - j11r*j10i;
    atomicAdd(&out[2], (-c1r - c2r) * weight);
    atomicAdd(&out[3], ( c2i - c1i) * weight);
    atomicAdd(&out[6], ( c2r - c1r) * weight);
    atomicAdd(&out[7], (-c1i - c2i) * weight);

    c1r = j00r*j10r + j00i*j10i;
    c1i = j00i*j10r - j00r*j10i;
    c2r = j11r*j01r + j11i*j01i;
    c2i = j11i*j01r - j11r*j01i;
    atomicAdd(&out[8], (-c1r - c2r) * weight);
    atomicAdd(&out[9], ( c2r - c1r) * weight);
    atomicAdd(&out[12], ( c1i - c2i) * weight);
    atomicAdd(&out[13], ( c2i + c1i) * weight);

    c1r = j00r*j11r + j00i*j11i;
    c1i = j00i*j11r - j00r*j11i;
    c2r = j01r*j10r + j01i*j10i;
    c2i = j01i*j10r - j01r*j10i;
    atomicAdd(&out[10], ( c1r + c2r) * weight);
    atomicAdd(&out[11], ( c1i - c2i) * weight);
    atomicAdd(&out[14], (-c1i - c2i) * weight);
    atomicAdd(&out[15], ( c1r - c2r) * weight);

    j = &jNoShadow[idx * 8];
    out = &mNoShadow[grid * 16];
    j00r = j[0]; j00i = j[1];
    j01r = j[2]; j01i = j[3];
    j10r = j[4]; j10i = j[5];
    j11r = j[6]; j11i = j[7];

    a11 = j00r*j00r + j00i*j00i;
    a12 = j01r*j01r + j01i*j01i;
    a21 = j10r*j10r + j10i*j10i;
    a22 = j11r*j11r + j11i*j11i;

    A1 = a11 + a21; A2 = a12 + a22;
    atomicAdd(&out[0], ((A1 + A2) * 0.5) * weight);
    atomicAdd(&out[1], ((A1 - A2) * 0.5) * weight);
    A1 = a11 - a21; A2 = a12 - a22;
    atomicAdd(&out[4], ((A1 + A2) * 0.5) * weight);
    atomicAdd(&out[5], ((A1 - A2) * 0.5) * weight);

    c1r = j00r*j01r + j00i*j01i;
    c1i = j00i*j01r - j00r*j01i;
    c2r = j11r*j10r + j11i*j10i;
    c2i = j11i*j10r - j11r*j10i;
    atomicAdd(&out[2], (-c1r - c2r) * weight);
    atomicAdd(&out[3], ( c2i - c1i) * weight);
    atomicAdd(&out[6], ( c2r - c1r) * weight);
    atomicAdd(&out[7], (-c1i - c2i) * weight);

    c1r = j00r*j10r + j00i*j10i;
    c1i = j00i*j10r - j00r*j10i;
    c2r = j11r*j01r + j11i*j01i;
    c2i = j11i*j01r - j11r*j01i;
    atomicAdd(&out[8], (-c1r - c2r) * weight);
    atomicAdd(&out[9], ( c2r - c1r) * weight);
    atomicAdd(&out[12], ( c1i - c2i) * weight);
    atomicAdd(&out[13], ( c2i + c1i) * weight);

    c1r = j00r*j11r + j00i*j11i;
    c1i = j00i*j11r - j00r*j11i;
    c2r = j01r*j10r + j01i*j10i;
    c2i = j01i*j10r - j01r*j10i;
    atomicAdd(&out[10], ( c1r + c2r) * weight);
    atomicAdd(&out[11], ( c1i - c2i) * weight);
    atomicAdd(&out[14], (-c1i - c2i) * weight);
    atomicAdd(&out[15], ( c1r - c2r) * weight);
}

static void add_mueller_from_jones(const std::vector<GpuReal> &jones,
                                   GpuReal weight,
                                   int nAz, int nZen,
                                   Arr2D &out)
{
    // Match the existing CPU coherent accumulator normalization.
    weight *= 0.25;

    for (int p = 0; p < nAz; ++p)
    {
        for (int t = 0; t <= nZen; ++t)
        {
            int grid = p * (nZen + 1) + t;
            const GpuReal *j = &jones[grid * 8];
            GpuReal j00r = j[0], j00i = j[1];
            GpuReal j01r = j[2], j01i = j[3];
            GpuReal j10r = j[4], j10i = j[5];
            GpuReal j11r = j[6], j11i = j[7];

            GpuReal a11 = j00r*j00r + j00i*j00i;
            GpuReal a12 = j01r*j01r + j01i*j01i;
            GpuReal a21 = j10r*j10r + j10i*j10i;
            GpuReal a22 = j11r*j11r + j11i*j11i;

            GpuReal A1 = a11 + a21, A2 = a12 + a22;
            out(p, t, 0, 0) += ((A1 + A2) * 0.5) * weight;
            out(p, t, 0, 1) += ((A1 - A2) * 0.5) * weight;
            A1 = a11 - a21; A2 = a12 - a22;
            out(p, t, 1, 0) += ((A1 + A2) * 0.5) * weight;
            out(p, t, 1, 1) += ((A1 - A2) * 0.5) * weight;

            GpuReal c1r = j00r*j01r + j00i*j01i;
            GpuReal c1i = j00i*j01r - j00r*j01i;
            GpuReal c2r = j11r*j10r + j11i*j10i;
            GpuReal c2i = j11i*j10r - j11r*j10i;
            out(p, t, 0, 2) += (-c1r - c2r) * weight;
            out(p, t, 0, 3) += ( c2i - c1i) * weight;
            out(p, t, 1, 2) += ( c2r - c1r) * weight;
            out(p, t, 1, 3) += (-c1i - c2i) * weight;

            c1r = j00r*j10r + j00i*j10i;
            c1i = j00i*j10r - j00r*j10i;
            c2r = j11r*j01r + j11i*j01i;
            c2i = j11i*j01r - j11r*j01i;
            out(p, t, 2, 0) += (-c1r - c2r) * weight;
            out(p, t, 2, 1) += ( c2r - c1r) * weight;
            out(p, t, 3, 0) += ( c1i - c2i) * weight;
            out(p, t, 3, 1) += ( c2i + c1i) * weight;

            c1r = j00r*j11r + j00i*j11i;
            c1i = j00i*j11r - j00r*j11i;
            c2r = j01r*j10r + j01i*j10i;
            c2i = j01i*j10r - j01r*j10i;
            out(p, t, 2, 2) += ( c1r + c2r) * weight;
            out(p, t, 2, 3) += ( c1i - c2i) * weight;
            out(p, t, 3, 2) += (-c1i - c2i) * weight;
            out(p, t, 3, 3) += ( c1r - c2r) * weight;
        }
    }
}

bool HandlerPO::HandleBeamsToLocalGpu(const PreparedOrientation &prepared,
                                       Arr2D &localM,
                                       Arr2D &localM_noshadow)
{
    if (!isCoh)
        return false;

    const int nBeams = (int)prepared.beams.size();
    const int nAz = m_sphere.nAzimuth;
    const int nZen = m_sphere.nZenith;
    const int gridCount = nAz * (nZen + 1);
    if (nBeams == 0 || gridCount == 0)
        return true;

    std::vector<GpuBeam> hBeams(nBeams);
    for (int i = 0; i < nBeams; ++i)
    {
        const PreparedBeam &pb = prepared.beams[i];
        if (!pb.edgeData.valid || pb.edgeData.nVertices <= 0 || pb.edgeData.nVertices > 32)
            return false;
        GpuBeam &b = hBeams[i];
        b.nVertices = pb.edgeData.nVertices;
        b.isExternal = pb.isExternal ? 1 : 0;
        b.orientation = 0;
        for (int e = 0; e < 32; ++e)
        {
            b.x[e] = pb.edgeData.x[e];
            b.y[e] = pb.edgeData.y[e];
            b.slope_yx[e] = pb.edgeData.slope_yx[e];
            b.slope_xy[e] = pb.edgeData.slope_xy[e];
            b.edge_valid_x[e] = pb.edgeData.edge_valid_x[e] ? 1 : 0;
            b.edge_valid_y[e] = pb.edgeData.edge_valid_y[e] ? 1 : 0;
        }
        b.bdx = pb.bdx; b.bdy = pb.bdy; b.bdz = pb.bdz;
        b.horAx = pb.horAx; b.horAy = pb.horAy; b.horAz = pb.horAz;
        b.verAx = pb.verAx; b.verAy = pb.verAy; b.verAz = pb.verAz;
        b.cenx = pb.cenx; b.ceny = pb.ceny; b.cenz = pb.cenz;
        b.beam_area = pb.beam_area;
        b.pNTx = pb.pNTx; b.pNTy = pb.pNTy; b.pNTz = pb.pNTz;
        b.pNPx = pb.pNPx; b.pNPy = pb.pNPy; b.pNPz = pb.pNPz;
        b.pnxDTx = pb.pnxDTx; b.pnxDTy = pb.pnxDTy; b.pnxDTz = pb.pnxDTz;
        b.pnxDPx = pb.pnxDPx; b.pnxDPy = pb.pnxDPy; b.pnxDPz = pb.pnxDPz;
        b.jp00r = pb.jp00r; b.jp00i = pb.jp00i;
        b.jp01r = pb.jp01r; b.jp01i = pb.jp01i;
        b.jp10r = pb.jp10r; b.jp10i = pb.jp10i;
        b.jp11r = pb.jp11r; b.jp11i = pb.jp11i;
    }

    std::vector<GpuReal> hJ(gridCount * 8, 0.0), hJns(gridCount * 8, 0.0);
    GpuWorkspace &ws = g_gpuWorkspace;

    if (!ensure_device_capacity(ws.beams, ws.beamCap, hBeams.size())) return false;
    if (!ensure_device_capacity(ws.sinTheta, ws.sinThetaCap, (size_t)nZen + 1)) return false;
    if (!ensure_device_capacity(ws.cosTheta, ws.cosThetaCap, (size_t)nZen + 1)) return false;
    if (!ensure_device_capacity(ws.sinPhi, ws.sinPhiCap, (size_t)nAz)) return false;
    if (!ensure_device_capacity(ws.cosPhi, ws.cosPhiCap, (size_t)nAz)) return false;
    if (!ensure_device_capacity(ws.vf, ws.vfCap, (size_t)gridCount * 3)) return false;
    if (!ensure_device_capacity(ws.j, ws.jCap, hJ.size())) return false;
    if (!ensure_device_capacity(ws.jNoShadow, ws.jNoShadowCap, hJns.size())) return false;

    if (ws.gridNAz != nAz || ws.gridNZen != nZen)
    {
        std::vector<GpuReal> hSinTheta(nZen + 1), hCosTheta(nZen + 1);
        for (int t = 0; t <= nZen; ++t)
        {
            GpuReal theta = m_sphere.GetZenith(t);
            hSinTheta[t] = sin(theta);
            hCosTheta[t] = cos(theta);
        }
        std::vector<GpuReal> hSinPhi(nAz), hCosPhi(nAz);
        for (int p = 0; p < nAz; ++p)
        {
            GpuReal phi = p * m_sphere.azinuthStep;
            hSinPhi[p] = sin(phi);
            hCosPhi[p] = cos(phi);
        }
        std::vector<GpuReal> hVf(gridCount * 3);
        for (int p = 0; p < nAz; ++p)
            for (int t = 0; t <= nZen; ++t)
            {
                int grid = p * (nZen + 1) + t;
                Point3d &v = m_sphere.vf[p][t];
                hVf[grid * 3 + 0] = v.x;
                hVf[grid * 3 + 1] = v.y;
                hVf[grid * 3 + 2] = v.z;
            }

        if (cudaMemcpy(ws.sinTheta, hSinTheta.data(), hSinTheta.size() * sizeof(GpuReal), cudaMemcpyHostToDevice) != cudaSuccess) return false;
        if (cudaMemcpy(ws.cosTheta, hCosTheta.data(), hCosTheta.size() * sizeof(GpuReal), cudaMemcpyHostToDevice) != cudaSuccess) return false;
        if (cudaMemcpy(ws.sinPhi, hSinPhi.data(), hSinPhi.size() * sizeof(GpuReal), cudaMemcpyHostToDevice) != cudaSuccess) return false;
        if (cudaMemcpy(ws.cosPhi, hCosPhi.data(), hCosPhi.size() * sizeof(GpuReal), cudaMemcpyHostToDevice) != cudaSuccess) return false;
        if (cudaMemcpy(ws.vf, hVf.data(), hVf.size() * sizeof(GpuReal), cudaMemcpyHostToDevice) != cudaSuccess) return false;
        ws.gridNAz = nAz;
        ws.gridNZen = nZen;
    }

    if (cudaMemcpy(ws.beams, hBeams.data(), hBeams.size() * sizeof(GpuBeam), cudaMemcpyHostToDevice) != cudaSuccess) return false;
    if (cudaMemset(ws.j, 0, hJ.size() * sizeof(GpuReal)) != cudaSuccess) return false;
    if (cudaMemset(ws.jNoShadow, 0, hJns.size() * sizeof(GpuReal)) != cudaSuccess) return false;

    long long total = (long long)nBeams * gridCount;
    int block = 256;
    int grid = (int)((total + block - 1) / block);
    diffraction_kernel<<<grid, block>>>(
        ws.beams, nBeams, ws.sinTheta, ws.cosTheta, ws.sinPhi, ws.cosPhi, ws.vf,
        nAz, nZen, m_waveIndex, m_wi2, m_eps1, m_eps2,
        real(m_complWave), imag(m_complWave),
        real(m_invComplWave), imag(m_invComplWave),
        m_legacySign ? 1 : 0, ws.j, ws.jNoShadow);

    if (cudaGetLastError() != cudaSuccess || cudaDeviceSynchronize() != cudaSuccess)
        return false;

    if (cudaMemcpy(hJ.data(), ws.j, hJ.size() * sizeof(GpuReal), cudaMemcpyDeviceToHost) != cudaSuccess) return false;
    if (cudaMemcpy(hJns.data(), ws.jNoShadow, hJns.size() * sizeof(GpuReal), cudaMemcpyDeviceToHost) != cudaSuccess) return false;

    add_mueller_from_jones(hJ, prepared.sinZenith, nAz, nZen, localM);
    add_mueller_from_jones(hJns, prepared.sinZenith, nAz, nZen, localM_noshadow);
    return true;
}

int HandlerPO::SelectGpuOrientationBatchSize(const std::vector<PreparedOrientation> &prepared,
                                             int start,
                                             int maxCount) const
{
    if (maxCount <= 1)
        return std::max(1, maxCount);

    size_t freeBytes = 0, totalBytes = 0;
    if (cudaMemGetInfo(&freeBytes, &totalBytes) != cudaSuccess || freeBytes == 0)
        return std::min(64, maxCount);

    const int nAz = m_sphere.nAzimuth;
    const int nZen = m_sphere.nZenith;
    const size_t gridCount = (size_t)nAz * ((size_t)nZen + 1);

    const size_t gridVectorsBytes =
        ((size_t)nZen + 1) * 2 * sizeof(GpuReal) +
        (size_t)nAz * 2 * sizeof(GpuReal) +
        gridCount * 3 * sizeof(GpuReal);
    const size_t muellerBytes = gridCount * 16 * 2 * sizeof(GpuReal);
    const size_t fixedBytes = gridVectorsBytes + muellerBytes + (256ULL << 20);
    const size_t jonesPerOrientation = gridCount * 8 * 2 * sizeof(GpuReal);

    size_t usable = (size_t)(freeBytes * gpu_memory_fraction());
    if (usable <= fixedBytes + jonesPerOrientation)
        return 1;
    usable -= fixedBytes;

    int count = 0;
    size_t used = 0;
    for (int i = start; i < (int)prepared.size() && count < maxCount; ++i)
    {
        size_t candidate = jonesPerOrientation +
                           prepared[i].beams.size() * sizeof(GpuBeam);
        if (count > 0 && used + candidate > usable)
            break;
        used += candidate;
        ++count;
    }

    return std::max(1, count);
}

static void add_arr2d_inplace(const Arr2D &src, int nAz, int nZen, Arr2D &dst)
{
    for (int p = 0; p < nAz; ++p)
        for (int t = 0; t <= nZen; ++t)
            for (int r = 0; r < 4; ++r)
                for (int c = 0; c < 4; ++c)
                    dst(p, t, r, c) += src(p, t, r, c);
}

static void report_fft_check(const Arr2D &fftM, const Arr2D &refM,
                             int nAz, int nZen, const char *name)
{
    double refMax = 0.0;
    for (int p = 0; p < nAz; ++p)
        for (int t = 0; t <= nZen; ++t)
            refMax = std::max(refMax, fabs(refM(p, t, 0, 0)));
    double significant = 1e-3 * refMax;
    double maxAbs = 0.0;
    double maxRel = 0.0;
    double sigMaxRel = 0.0;
    double rmsRel = 0.0;
    int count = 0;
    int sigCount = 0;
    for (int p = 0; p < nAz; ++p)
        for (int t = 0; t <= nZen; ++t)
        {
            double a = refM(p, t, 0, 0);
            double b = fftM(p, t, 0, 0);
            double absErr = fabs(b - a);
            double relErr = absErr / std::max(fabs(a), 1e-12);
            maxAbs = std::max(maxAbs, absErr);
            maxRel = std::max(maxRel, relErr);
            if (fabs(a) >= significant)
            {
                sigMaxRel = std::max(sigMaxRel, relErr);
                ++sigCount;
            }
            rmsRel += relErr * relErr;
            ++count;
        }
    rmsRel = sqrt(rmsRel / std::max(1, count));
    std::fprintf(stderr,
                 "GPU FFT check %s M11: max_abs=%.6g max_rel=%.6g sig_max_rel=%.6g rms_rel=%.6g samples=%d sig_samples=%d\n",
                 name, maxAbs, maxRel, sigMaxRel, rmsRel, count, sigCount);
}

static bool fft_upsample_phi_arr2d(const Arr2D &low,
                                   int nLow,
                                   int nFull,
                                   int nZen,
                                   Arr2D &dst)
{
    if (nLow <= 0 || nFull <= 0 || nZen < 0)
        return false;

    const int nElem = 16;
    const int batch = (nZen + 1) * nElem;
    const size_t lowCount = (size_t)batch * nLow;
    const size_t fullCount = (size_t)batch * nFull;

    std::vector<GpuComplex> hLow(lowCount);
    for (int t = 0; t <= nZen; ++t)
        for (int r = 0; r < 4; ++r)
            for (int c = 0; c < 4; ++c)
            {
                int b = (t * nElem) + r * 4 + c;
                for (int p = 0; p < nLow; ++p)
                {
                    GpuComplex z;
                    z.x = (GpuReal)low(p, t, r, c);
                    z.y = (GpuReal)0.0;
                    hLow[(size_t)b * nLow + p] = z;
                }
            }

    GpuComplex *dLow = nullptr;
    GpuComplex *dFull = nullptr;
    cufftHandle planLow = 0;
    cufftHandle planFull = 0;
    bool ok = false;

    do
    {
        if (cudaMalloc(&dLow, lowCount * sizeof(GpuComplex)) != cudaSuccess) break;
        if (cudaMalloc(&dFull, fullCount * sizeof(GpuComplex)) != cudaSuccess) break;
        if (cudaMemcpy(dLow, hLow.data(), lowCount * sizeof(GpuComplex),
                       cudaMemcpyHostToDevice) != cudaSuccess) break;

        int nLowArr[1] = { nLow };
        int nFullArr[1] = { nFull };
        if (cufftPlanMany(&planLow, 1, nLowArr,
                          nullptr, 1, nLow,
                          nullptr, 1, nLow,
                          GpuCufftType, batch) != CUFFT_SUCCESS) break;
        if (cufftPlanMany(&planFull, 1, nFullArr,
                          nullptr, 1, nFull,
                          nullptr, 1, nFull,
                          GpuCufftType, batch) != CUFFT_SUCCESS) break;
        if (gpu_cufft_exec(planLow, dLow, dLow, CUFFT_FORWARD) != CUFFT_SUCCESS) break;
        if (cudaMemset(dFull, 0, fullCount * sizeof(GpuComplex)) != cudaSuccess) break;

        int block = 256;
        int grid = (int)((lowCount + block - 1) / block);
        GpuReal scale = (GpuReal)(1.0 / (double)nLow);
        fft_phi_pad_kernel<<<grid, block>>>(dLow, dFull, nLow, nFull, batch, scale);
        if (cudaGetLastError() != cudaSuccess) break;
        if (gpu_cufft_exec(planFull, dFull, dFull, CUFFT_INVERSE) != CUFFT_SUCCESS) break;

        std::vector<GpuComplex> hFull(fullCount);
        if (cudaMemcpy(hFull.data(), dFull, fullCount * sizeof(GpuComplex),
                       cudaMemcpyDeviceToHost) != cudaSuccess) break;

        for (int t = 0; t <= nZen; ++t)
            for (int r = 0; r < 4; ++r)
                for (int c = 0; c < 4; ++c)
                {
                    int b = (t * nElem) + r * 4 + c;
                    for (int p = 0; p < nFull; ++p)
                        dst(p, t, r, c) += hFull[(size_t)b * nFull + p].x;
                }

        ok = true;
    } while (false);

    if (planLow) cufftDestroy(planLow);
    if (planFull) cufftDestroy(planFull);
    cudaFree(dLow);
    cudaFree(dFull);
    return ok;
}

bool HandlerPO::HandleOrientationsToLocalGpu(const std::vector<PreparedOrientation> &prepared,
                                             Arr2D &localM,
                                             Arr2D &localM_noshadow)
{
    return HandleOrientationsToLocalGpu(prepared, 0, (int)prepared.size(),
                                        localM, localM_noshadow);
}

bool HandlerPO::HandleOrientationsToLocalGpuFftPhi(const std::vector<PreparedOrientation> &prepared,
                                                   int start,
                                                   int count,
                                                   Arr2D &localM,
                                                   Arr2D &localM_noshadow)
{
    if (!isCoh)
        return false;

    const int nFull = m_sphere.nAzimuth;
    const int nZen = m_sphere.nZenith;
    int factor = choose_fft_phi_factor(nFull);
    if (factor <= 1 || nFull < 32)
        return HandleOrientationsToLocalGpu(prepared, start, count,
                                            localM, localM_noshadow);

    int nLow = nFull / factor;
    if (nLow < 16) nLow = std::min(nFull, 16);
    if (nLow >= nFull)
        return HandleOrientationsToLocalGpu(prepared, start, count,
                                            localM, localM_noshadow);

    static bool printed = false;
    if (!printed)
    {
        std::fprintf(stderr,
                     "GPU FFT phi interpolation: direct Nphi=%d, output Nphi=%d (factor=%d). "
                     "This is angular Fourier interpolation, not aperture pFFT/FMM.\n",
                     nLow, nFull, factor);
        printed = true;
    }

    ScatteringRange fullSphere = m_sphere;
    ScatteringRange lowSphere = fullSphere;
    lowSphere.nAzimuth = nLow;
    lowSphere.azinuthStep = M_2PI / nLow;
    lowSphere.ComputeSphereDirections(*m_incidentLight);

    Arr2D lowM(nLow + 1, nZen + 1, 4, 4); lowM.ClearArr();
    Arr2D lowMns(nLow + 1, nZen + 1, 4, 4); lowMns.ClearArr();

    bool savedFft = m_fftEnabled;
    m_fftEnabled = false;
    m_sphere = lowSphere;
    bool ok = HandleOrientationsToLocalGpu(prepared, start, count, lowM, lowMns);
    m_sphere = fullSphere;
    m_fftEnabled = savedFft;
    if (!ok)
        return false;

    bool doCheck = gpu_fft_check_enabled() && factor > 2;
    if (!doCheck)
    {
        if (!fft_upsample_phi_arr2d(lowM, nLow, nFull, nZen, localM))
            return false;
        if (!fft_upsample_phi_arr2d(lowMns, nLow, nFull, nZen, localM_noshadow))
            return false;
        return true;
    }

    Arr2D fftM(nFull + 1, nZen + 1, 4, 4); fftM.ClearArr();
    Arr2D fftMns(nFull + 1, nZen + 1, 4, 4); fftMns.ClearArr();
    if (!fft_upsample_phi_arr2d(lowM, nLow, nFull, nZen, fftM))
        return false;
    if (!fft_upsample_phi_arr2d(lowMns, nLow, nFull, nZen, fftMns))
        return false;

    if (doCheck)
    {
        int nCheck = std::min(nFull, nLow * 2);
        ScatteringRange checkSphere = fullSphere;
        checkSphere.nAzimuth = nCheck;
        checkSphere.azinuthStep = M_2PI / nCheck;
        checkSphere.ComputeSphereDirections(*m_incidentLight);

        Arr2D checkLowM(nCheck + 1, nZen + 1, 4, 4); checkLowM.ClearArr();
        Arr2D checkLowMns(nCheck + 1, nZen + 1, 4, 4); checkLowMns.ClearArr();
        m_sphere = checkSphere;
        bool checkOk = HandleOrientationsToLocalGpu(prepared, start, count,
                                                    checkLowM, checkLowMns);
        m_sphere = fullSphere;
        if (!checkOk)
            return false;

        Arr2D checkM(nFull + 1, nZen + 1, 4, 4); checkM.ClearArr();
        Arr2D checkMns(nFull + 1, nZen + 1, 4, 4); checkMns.ClearArr();
        if (!fft_upsample_phi_arr2d(checkLowM, nCheck, nFull, nZen, checkM))
            return false;
        if (!fft_upsample_phi_arr2d(checkLowMns, nCheck, nFull, nZen, checkMns))
            return false;
        report_fft_check(fftM, checkM, nFull, nZen, "full");
        report_fft_check(fftMns, checkMns, nFull, nZen, "no-shadow");
    }

    add_arr2d_inplace(fftM, nFull, nZen, localM);
    add_arr2d_inplace(fftMns, nFull, nZen, localM_noshadow);

    return true;
}

bool HandlerPO::HandleOrientationsToLocalGpu(const std::vector<PreparedOrientation> &prepared,
                                             int start,
                                             int count,
                                             Arr2D &localM,
                                             Arr2D &localM_noshadow)
{
    if (!isCoh)
        return false;

    if (start < 0 || count < 0 || start > (int)prepared.size())
        return false;
    const int nOrient = std::min(count, (int)prepared.size() - start);
    const int nAz = m_sphere.nAzimuth;
    const int nZen = m_sphere.nZenith;
    const int gridCount = nAz * (nZen + 1);
    if (nOrient == 0 || gridCount == 0)
        return true;

    size_t nBeams = 0;
    for (int oi = 0; oi < nOrient; ++oi)
    {
        const PreparedOrientation &po = prepared[start + oi];
        for (const PreparedBeam &pb : po.beams)
        {
            if (!pb.edgeData.valid || pb.edgeData.nVertices <= 0 || pb.edgeData.nVertices > 32)
                return false;
            ++nBeams;
        }
    }
    if (nBeams == 0)
        return true;

    GpuWorkspace &ws = g_gpuWorkspace;
    ws.hBeams.resize(nBeams);
    ws.hWeights.assign(nOrient, 0.0);
    ws.hBeamOffsets.assign(nOrient + 1, 0);
    std::vector<GpuBeam> &hBeams = ws.hBeams;
    std::vector<GpuReal> &hWeights = ws.hWeights;
    std::vector<int> &hBeamOffsets = ws.hBeamOffsets;
    size_t bi = 0;
    for (int oi = 0; oi < nOrient; ++oi)
    {
        hBeamOffsets[oi] = (int)bi;
        const PreparedOrientation &po = prepared[start + oi];
        hWeights[oi] = po.sinZenith;
        for (const PreparedBeam &pb : po.beams)
        {
            GpuBeam &b = hBeams[bi++];
            b.nVertices = pb.edgeData.nVertices;
            b.isExternal = pb.isExternal ? 1 : 0;
            b.orientation = oi;
            for (int e = 0; e < 32; ++e)
            {
                b.x[e] = pb.edgeData.x[e];
                b.y[e] = pb.edgeData.y[e];
                b.slope_yx[e] = pb.edgeData.slope_yx[e];
                b.slope_xy[e] = pb.edgeData.slope_xy[e];
                b.edge_valid_x[e] = pb.edgeData.edge_valid_x[e] ? 1 : 0;
                b.edge_valid_y[e] = pb.edgeData.edge_valid_y[e] ? 1 : 0;
            }
            b.bdx = pb.bdx; b.bdy = pb.bdy; b.bdz = pb.bdz;
            b.horAx = pb.horAx; b.horAy = pb.horAy; b.horAz = pb.horAz;
            b.verAx = pb.verAx; b.verAy = pb.verAy; b.verAz = pb.verAz;
            b.cenx = pb.cenx; b.ceny = pb.ceny; b.cenz = pb.cenz;
            b.beam_area = pb.beam_area;
            b.pNTx = pb.pNTx; b.pNTy = pb.pNTy; b.pNTz = pb.pNTz;
            b.pNPx = pb.pNPx; b.pNPy = pb.pNPy; b.pNPz = pb.pNPz;
            b.pnxDTx = pb.pnxDTx; b.pnxDTy = pb.pnxDTy; b.pnxDTz = pb.pnxDTz;
            b.pnxDPx = pb.pnxDPx; b.pnxDPy = pb.pnxDPy; b.pnxDPz = pb.pnxDPz;
            b.jp00r = pb.jp00r; b.jp00i = pb.jp00i;
            b.jp01r = pb.jp01r; b.jp01i = pb.jp01i;
            b.jp10r = pb.jp10r; b.jp10i = pb.jp10i;
            b.jp11r = pb.jp11r; b.jp11i = pb.jp11i;
        }
    }
    hBeamOffsets[nOrient] = (int)bi;

    const size_t jCount = (size_t)nOrient * gridCount * 8;
    const size_t mCount = (size_t)gridCount * 16;

    if (!ensure_device_capacity(ws.beams, ws.beamCap, hBeams.size())) return false;
    if (!ensure_device_capacity(ws.weights, ws.weightsCap, hWeights.size())) return false;
    if (!ensure_device_capacity(ws.beamOffsets, ws.beamOffsetsCap, hBeamOffsets.size())) return false;
    if (!ensure_device_capacity(ws.sinTheta, ws.sinThetaCap, (size_t)nZen + 1)) return false;
    if (!ensure_device_capacity(ws.cosTheta, ws.cosThetaCap, (size_t)nZen + 1)) return false;
    if (!ensure_device_capacity(ws.sinPhi, ws.sinPhiCap, (size_t)nAz)) return false;
    if (!ensure_device_capacity(ws.cosPhi, ws.cosPhiCap, (size_t)nAz)) return false;
    if (!ensure_device_capacity(ws.vf, ws.vfCap, (size_t)gridCount * 3)) return false;
    if (!ensure_device_capacity(ws.j, ws.jCap, jCount)) return false;
    if (!ensure_device_capacity(ws.jNoShadow, ws.jNoShadowCap, jCount)) return false;
    if (!ensure_device_capacity(ws.m, ws.mCap, mCount)) return false;
    if (!ensure_device_capacity(ws.mNoShadow, ws.mNoShadowCap, mCount)) return false;

    if (ws.gridNAz != nAz || ws.gridNZen != nZen)
    {
        std::vector<GpuReal> hSinTheta(nZen + 1), hCosTheta(nZen + 1);
        for (int t = 0; t <= nZen; ++t)
        {
            GpuReal theta = m_sphere.GetZenith(t);
            hSinTheta[t] = sin(theta);
            hCosTheta[t] = cos(theta);
        }
        std::vector<GpuReal> hSinPhi(nAz), hCosPhi(nAz);
        for (int p = 0; p < nAz; ++p)
        {
            GpuReal phi = p * m_sphere.azinuthStep;
            hSinPhi[p] = sin(phi);
            hCosPhi[p] = cos(phi);
        }
        std::vector<GpuReal> hVf(gridCount * 3);
        for (int p = 0; p < nAz; ++p)
            for (int t = 0; t <= nZen; ++t)
            {
                int grid = p * (nZen + 1) + t;
                Point3d &v = m_sphere.vf[p][t];
                hVf[grid * 3 + 0] = v.x;
                hVf[grid * 3 + 1] = v.y;
                hVf[grid * 3 + 2] = v.z;
            }

        if (cudaMemcpy(ws.sinTheta, hSinTheta.data(), hSinTheta.size() * sizeof(GpuReal), cudaMemcpyHostToDevice) != cudaSuccess) return false;
        if (cudaMemcpy(ws.cosTheta, hCosTheta.data(), hCosTheta.size() * sizeof(GpuReal), cudaMemcpyHostToDevice) != cudaSuccess) return false;
        if (cudaMemcpy(ws.sinPhi, hSinPhi.data(), hSinPhi.size() * sizeof(GpuReal), cudaMemcpyHostToDevice) != cudaSuccess) return false;
        if (cudaMemcpy(ws.cosPhi, hCosPhi.data(), hCosPhi.size() * sizeof(GpuReal), cudaMemcpyHostToDevice) != cudaSuccess) return false;
        if (cudaMemcpy(ws.vf, hVf.data(), hVf.size() * sizeof(GpuReal), cudaMemcpyHostToDevice) != cudaSuccess) return false;
        ws.gridNAz = nAz;
        ws.gridNZen = nZen;
    }

    if (cudaMemcpy(ws.beams, hBeams.data(), hBeams.size() * sizeof(GpuBeam), cudaMemcpyHostToDevice) != cudaSuccess) return false;
    if (cudaMemcpy(ws.weights, hWeights.data(), hWeights.size() * sizeof(GpuReal), cudaMemcpyHostToDevice) != cudaSuccess) return false;
    if (cudaMemcpy(ws.beamOffsets, hBeamOffsets.data(), hBeamOffsets.size() * sizeof(int), cudaMemcpyHostToDevice) != cudaSuccess) return false;
    if (cudaMemset(ws.j, 0, jCount * sizeof(GpuReal)) != cudaSuccess) return false;
    if (cudaMemset(ws.jNoShadow, 0, jCount * sizeof(GpuReal)) != cudaSuccess) return false;
    if (cudaMemset(ws.m, 0, mCount * sizeof(GpuReal)) != cudaSuccess) return false;
    if (cudaMemset(ws.mNoShadow, 0, mCount * sizeof(GpuReal)) != cudaSuccess) return false;

    int block = 256;
    bool noAtomics = gpu_no_atomics_enabled();
    long long diffractionTotal = (long long)(noAtomics ? nOrient : (int)hBeams.size()) * gridCount;
    int diffractionGrid = (int)((diffractionTotal + block - 1) / block);
    if (noAtomics)
    {
        diffraction_grid_kernel<<<diffractionGrid, block>>>(
            ws.beams, ws.beamOffsets, ws.sinTheta, ws.cosTheta, ws.sinPhi, ws.cosPhi, ws.vf,
            nAz, nZen, nOrient, m_waveIndex, m_wi2, m_eps1, m_eps2,
            real(m_complWave), imag(m_complWave),
            real(m_invComplWave), imag(m_invComplWave),
            m_legacySign ? 1 : 0, ws.j, ws.jNoShadow);
    }
    else
    {
        diffraction_kernel<<<diffractionGrid, block>>>(
            ws.beams, (int)hBeams.size(), ws.sinTheta, ws.cosTheta, ws.sinPhi, ws.cosPhi, ws.vf,
            nAz, nZen, m_waveIndex, m_wi2, m_eps1, m_eps2,
            real(m_complWave), imag(m_complWave),
            real(m_invComplWave), imag(m_invComplWave),
            m_legacySign ? 1 : 0, ws.j, ws.jNoShadow);
    }
    cudaError_t err = cudaGetLastError();
    if (!gpu_report_cuda_error(err, "diffraction kernel launch"))
        return false;

    long long muellerTotal = (long long)nOrient * gridCount;
    int muellerGrid = (int)((muellerTotal + block - 1) / block);
    mueller_batch_kernel<<<muellerGrid, block>>>(ws.j, ws.jNoShadow, ws.weights,
                                                 nOrient, gridCount,
                                                 ws.m, ws.mNoShadow);
    err = cudaGetLastError();
    if (!gpu_report_cuda_error(err, "mueller kernel launch"))
        return false;
    err = cudaDeviceSynchronize();
    if (!gpu_report_cuda_error(err, "mueller kernel sync"))
        return false;

    ws.hM.resize(mCount);
    ws.hMNoShadow.resize(mCount);
    std::vector<GpuReal> &hM = ws.hM;
    std::vector<GpuReal> &hMns = ws.hMNoShadow;
    if (cudaMemcpy(hM.data(), ws.m, mCount * sizeof(GpuReal), cudaMemcpyDeviceToHost) != cudaSuccess) return false;
    if (cudaMemcpy(hMns.data(), ws.mNoShadow, mCount * sizeof(GpuReal), cudaMemcpyDeviceToHost) != cudaSuccess) return false;

    for (int p = 0; p < nAz; ++p)
        for (int t = 0; t <= nZen; ++t)
        {
            int grid = p * (nZen + 1) + t;
            const GpuReal *m = &hM[grid * 16];
            const GpuReal *mn = &hMns[grid * 16];
            for (int r = 0; r < 4; ++r)
                for (int c = 0; c < 4; ++c)
                {
                    localM(p, t, r, c) += m[r * 4 + c];
                    localM_noshadow(p, t, r, c) += mn[r * 4 + c];
                }
        }

    return true;
}

bool HandlerPO::DiffractThetasGpu(const std::vector<PreparedOrientation> &prepared,
                                  const double *theta_rads,
                                  int nPoints,
                                  std::vector<double> &m11_out)
{
    if (nPoints <= 0)
        return true;

    const int nOrient = (int)prepared.size();
    const int nAz = m_sphere.nAzimuth;
    if (nOrient == 0 || nAz <= 0)
        return true;

    size_t nBeams = 0;
    for (int oi = 0; oi < nOrient; ++oi)
    {
        for (const PreparedBeam &pb : prepared[oi].beams)
        {
            if (!pb.edgeData.valid || pb.edgeData.nVertices <= 0 || pb.edgeData.nVertices > 32)
                return false;
            ++nBeams;
        }
    }

    m11_out.assign(nPoints, 0.0);
    if (nBeams == 0)
        return true;

    GpuWorkspace &ws = g_gpuWorkspace;
    ws.hBeams.resize(nBeams);
    ws.hWeights.assign(nOrient, 0.0);

    size_t bi = 0;
    for (int oi = 0; oi < nOrient; ++oi)
    {
        ws.hWeights[oi] = prepared[oi].sinZenith;
        for (const PreparedBeam &pb : prepared[oi].beams)
        {
            GpuBeam &b = ws.hBeams[bi++];
            b.nVertices = pb.edgeData.nVertices;
            b.isExternal = pb.isExternal ? 1 : 0;
            b.orientation = oi;
            for (int e = 0; e < 32; ++e)
            {
                b.x[e] = pb.edgeData.x[e];
                b.y[e] = pb.edgeData.y[e];
                b.slope_yx[e] = pb.edgeData.slope_yx[e];
                b.slope_xy[e] = pb.edgeData.slope_xy[e];
                b.edge_valid_x[e] = pb.edgeData.edge_valid_x[e] ? 1 : 0;
                b.edge_valid_y[e] = pb.edgeData.edge_valid_y[e] ? 1 : 0;
            }
            b.bdx = pb.bdx; b.bdy = pb.bdy; b.bdz = pb.bdz;
            b.horAx = pb.horAx; b.horAy = pb.horAy; b.horAz = pb.horAz;
            b.verAx = pb.verAx; b.verAy = pb.verAy; b.verAz = pb.verAz;
            b.cenx = pb.cenx; b.ceny = pb.ceny; b.cenz = pb.cenz;
            b.beam_area = pb.beam_area;
            b.pNTx = pb.pNTx; b.pNTy = pb.pNTy; b.pNTz = pb.pNTz;
            b.pNPx = pb.pNPx; b.pNPy = pb.pNPy; b.pNPz = pb.pNPz;
            b.pnxDTx = pb.pnxDTx; b.pnxDTy = pb.pnxDTy; b.pnxDTz = pb.pnxDTz;
            b.pnxDPx = pb.pnxDPx; b.pnxDPy = pb.pnxDPy; b.pnxDPz = pb.pnxDPz;
            b.jp00r = pb.jp00r; b.jp00i = pb.jp00i;
            b.jp01r = pb.jp01r; b.jp01i = pb.jp01i;
            b.jp10r = pb.jp10r; b.jp10i = pb.jp10i;
            b.jp11r = pb.jp11r; b.jp11i = pb.jp11i;
        }
    }

    std::vector<GpuReal> hSinTheta(nPoints), hCosTheta(nPoints);
    for (int t = 0; t < nPoints; ++t)
    {
        hSinTheta[t] = (GpuReal)sin(theta_rads[t]);
        hCosTheta[t] = (GpuReal)cos(theta_rads[t]);
    }

    std::vector<GpuReal> hSinPhi(nAz), hCosPhi(nAz);
    for (int p = 0; p < nAz; ++p)
    {
        GpuReal phi = p * m_sphere.azinuthStep;
        hSinPhi[p] = sin(phi);
        hCosPhi[p] = cos(phi);
    }

    const int gridCount = nAz * nPoints;
    std::vector<GpuReal> hVf((size_t)gridCount * 3);
    for (int p = 0; p < nAz; ++p)
    {
        GpuReal cp = hCosPhi[p], sp = hSinPhi[p];
        for (int t = 0; t < nPoints; ++t)
        {
            GpuReal dz = -hCosTheta[t];
            GpuReal vfx, vfy, vfz;
            if (dz >= 1.0 - 1e-15) { vfx = 0.0; vfy = -1.0; vfz = 0.0; }
            else if (dz <= -1.0 + 1e-15) { vfx = 0.0; vfy = 1.0; vfz = 0.0; }
            else { vfx = -sp; vfy = cp; vfz = 0.0; }
            int grid = p * nPoints + t;
            hVf[grid * 3 + 0] = vfx;
            hVf[grid * 3 + 1] = vfy;
            hVf[grid * 3 + 2] = vfz;
        }
    }

    if (!ensure_device_capacity(ws.beams, ws.beamCap, ws.hBeams.size())) return false;
    if (!ensure_device_capacity(ws.weights, ws.weightsCap, ws.hWeights.size())) return false;
    if (!ensure_device_capacity(ws.sinTheta, ws.sinThetaCap, (size_t)nPoints)) return false;
    if (!ensure_device_capacity(ws.cosTheta, ws.cosThetaCap, (size_t)nPoints)) return false;
    if (!ensure_device_capacity(ws.sinPhi, ws.sinPhiCap, (size_t)nAz)) return false;
    if (!ensure_device_capacity(ws.cosPhi, ws.cosPhiCap, (size_t)nAz)) return false;
    if (!ensure_device_capacity(ws.vf, ws.vfCap, hVf.size())) return false;
    if (!ensure_device_capacity(ws.m, ws.mCap, (size_t)nPoints)) return false;

    if (cudaMemcpy(ws.beams, ws.hBeams.data(), ws.hBeams.size() * sizeof(GpuBeam), cudaMemcpyHostToDevice) != cudaSuccess) return false;
    if (cudaMemcpy(ws.weights, ws.hWeights.data(), ws.hWeights.size() * sizeof(GpuReal), cudaMemcpyHostToDevice) != cudaSuccess) return false;
    if (cudaMemcpy(ws.sinTheta, hSinTheta.data(), hSinTheta.size() * sizeof(GpuReal), cudaMemcpyHostToDevice) != cudaSuccess) return false;
    if (cudaMemcpy(ws.cosTheta, hCosTheta.data(), hCosTheta.size() * sizeof(GpuReal), cudaMemcpyHostToDevice) != cudaSuccess) return false;
    if (cudaMemcpy(ws.sinPhi, hSinPhi.data(), hSinPhi.size() * sizeof(GpuReal), cudaMemcpyHostToDevice) != cudaSuccess) return false;
    if (cudaMemcpy(ws.cosPhi, hCosPhi.data(), hCosPhi.size() * sizeof(GpuReal), cudaMemcpyHostToDevice) != cudaSuccess) return false;
    if (cudaMemcpy(ws.vf, hVf.data(), hVf.size() * sizeof(GpuReal), cudaMemcpyHostToDevice) != cudaSuccess) return false;
    if (cudaMemset(ws.m, 0, (size_t)nPoints * sizeof(GpuReal)) != cudaSuccess) return false;

    int block = 256;
    long long total = (long long)nBeams * gridCount;
    int launchGrid = (int)((total + block - 1) / block);
    theta_m11_kernel<<<launchGrid, block>>>(
        ws.beams, (int)nBeams, ws.weights,
        ws.sinTheta, ws.cosTheta, ws.sinPhi, ws.cosPhi, ws.vf,
        nAz, nPoints, m_waveIndex, m_wi2, m_eps1, m_eps2,
        real(m_complWave), imag(m_complWave),
        real(m_invComplWave), imag(m_invComplWave),
        m_legacySign ? 1 : 0, ws.m);

    if (cudaGetLastError() != cudaSuccess || cudaDeviceSynchronize() != cudaSuccess)
        return false;

    std::vector<GpuReal> hM11(nPoints);
    if (cudaMemcpy(hM11.data(), ws.m, (size_t)nPoints * sizeof(GpuReal), cudaMemcpyDeviceToHost) != cudaSuccess)
        return false;
    for (int i = 0; i < nPoints; ++i)
        m11_out[i] = hM11[i];

    ws.gridNAz = -1;
    ws.gridNZen = -1;
    return true;
}
