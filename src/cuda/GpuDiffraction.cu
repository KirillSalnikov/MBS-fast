#include "HandlerPO.h"

#include <cuda_runtime.h>
#include <cmath>
#include <algorithm>
#include <cstdlib>
#include <vector>

struct GpuBeam
{
    double x[32], y[32];
    double slope_yx[32], slope_xy[32];
    unsigned char edge_valid_x[32], edge_valid_y[32];
    int nVertices;
    int isExternal;

    double bdx, bdy, bdz;
    double horAx, horAy, horAz;
    double verAx, verAy, verAz;
    double cenx, ceny, cenz;
    double beam_area;
    double pNTx, pNTy, pNTz;
    double pNPx, pNPy, pNPz;
    double pnxDTx, pnxDTy, pnxDTz;
    double pnxDPx, pnxDPy, pnxDPz;
    double jp00r, jp00i, jp01r, jp01i;
    double jp10r, jp10i, jp11r, jp11i;
    int orientation;
};

struct GpuWorkspace
{
    GpuBeam *beams = nullptr;
    double *sinTheta = nullptr, *cosTheta = nullptr;
    double *sinPhi = nullptr, *cosPhi = nullptr;
    double *vf = nullptr;
    double *j = nullptr, *jNoShadow = nullptr;
    double *weights = nullptr;
    int *beamOffsets = nullptr;
    double *m = nullptr, *mNoShadow = nullptr;
    size_t beamCap = 0;
    size_t sinThetaCap = 0, cosThetaCap = 0;
    size_t sinPhiCap = 0, cosPhiCap = 0;
    size_t vfCap = 0, jCap = 0, jNoShadowCap = 0;
    size_t weightsCap = 0, beamOffsetsCap = 0, mCap = 0, mNoShadowCap = 0;
    int gridNAz = -1, gridNZen = -1;
};

static GpuWorkspace g_gpuWorkspace;

static double gpu_memory_fraction()
{
    const char *value = std::getenv("MBS_GPU_MEM_FRACTION");
    if (!value || !*value)
        return 0.88;
    char *end = nullptr;
    double parsed = std::strtod(value, &end);
    if (!end || *end != '\0' || parsed <= 0.1 || parsed > 0.98)
        return 0.88;
    return parsed;
}

static bool gpu_no_atomics_enabled()
{
    const char *value = std::getenv("MBS_GPU_NO_ATOMICS");
    return value && value[0] == '1' && value[1] == '\0';
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

__device__ inline void cmul(double ar, double ai, double br, double bi,
                            double &cr, double &ci)
{
    cr = ar * br - ai * bi;
    ci = ar * bi + ai * br;
}

__device__ inline void rotate_jones_gpu(
    double NTx, double NTy, double NTz,
    double NPx, double NPy, double NPz,
    double nxDTx, double nxDTy, double nxDTz,
    double nxDPx, double nxDPy, double nxDPz,
    double vfx, double vfy, double vfz,
    double dirx, double diry, double dirz,
    double &r00, double &r01, double &r10, double &r11)
{
    double vtx = vfy * dirz - vfz * diry;
    double vty = vfz * dirx - vfx * dirz;
    double vtz = vfx * diry - vfy * dirx;
    double vtLen2 = vtx * vtx + vty * vty + vtz * vtz;
    if (vtLen2 > 1e-30)
    {
        double invLen = rsqrt(vtLen2);
        vtx *= invLen; vty *= invLen; vtz *= invLen;
    }

    double dot_dir_nxDT = dirx * nxDTx + diry * nxDTy + dirz * nxDTz;
    double cpTx = (diry * NTz - dirz * NTy) + nxDTx - dirx * dot_dir_nxDT;
    double cpTy = (dirz * NTx - dirx * NTz) + nxDTy - diry * dot_dir_nxDT;
    double cpTz = (dirx * NTy - diry * NTx) + nxDTz - dirz * dot_dir_nxDT;

    double dot_dir_nxDP = dirx * nxDPx + diry * nxDPy + dirz * nxDPz;
    double cpPx = (diry * NPz - dirz * NPy) + nxDPx - dirx * dot_dir_nxDP;
    double cpPy = (dirz * NPx - dirx * NPz) + nxDPy - diry * dot_dir_nxDP;
    double cpPz = (dirx * NPy - diry * NPx) + nxDPz - dirz * dot_dir_nxDP;

    r00 = cpTx * vtx + cpTy * vty + cpTz * vtz;
    r01 = cpPx * vtx + cpPy * vty + cpPz * vtz;
    r10 = cpTx * vfx + cpTy * vfy + cpTz * vfz;
    r11 = cpPx * vfx + cpPy * vfy + cpPz * vfz;
}

__device__ inline bool compute_beam_jones_gpu(const GpuBeam &b,
                                              int grid,
                                              int nZen,
                                              const double *sinTheta,
                                              const double *cosTheta,
                                              const double *sinPhi,
                                              const double *cosPhi,
                                              const double *vf,
                                              double waveIndex,
                                              double wi2,
                                              double eps1,
                                              double eps2,
                                              double complWaveR,
                                              double complWaveI,
                                              double invComplWaveR,
                                              double invComplWaveI,
                                              int legacySign,
                                              double &d00r, double &d00i,
                                              double &d01r, double &d01i,
                                              double &d10r, double &d10i,
                                              double &d11r, double &d11i)
{
    int nv = b.nVertices;
    if (nv <= 0) return false;

    int p = grid / (nZen + 1);
    int t = grid - p * (nZen + 1);
    double cp = cosPhi[p], sp = sinPhi[p];
    double sin_t = sinTheta[t], cos_t = cosTheta[t];
    double dx = sin_t * cp;
    double dy = sin_t * sp;
    double dz = -cos_t;

    double neg_cp = -cp, neg_sp = -sp;
    double a_sin = neg_cp * b.horAx + neg_sp * b.horAy;
    double a_cos = b.horAz;
    double a0 = b.bdx * b.horAx + b.bdy * b.horAy + b.bdz * b.horAz;
    double b_sin = neg_cp * b.verAx + neg_sp * b.verAy;
    double b_cos = b.verAz;
    double b0 = b.bdx * b.verAx + b.bdy * b.verAy + b.bdz * b.verAz;

    double A = sin_t * a_sin + cos_t * a_cos + a0;
    double B = sin_t * b_sin + cos_t * b_cos + b0;
    double absA = fabs(A);
    double absB = fabs(B);

    double fr, fi;
    if (absA < eps2 && absB < eps2)
    {
        double sign = legacySign ? 1.0 : -1.0;
        fr = sign * invComplWaveR * b.beam_area;
        fi = sign * invComplWaveI * b.beam_area;
    }
    else
    {
        double vc[32], vs[32];
        for (int v = 0; v < nv; ++v)
        {
            double psin = waveIndex * (a_sin * b.x[v] + b_sin * b.y[v]);
            double pcos = waveIndex * (a_cos * b.x[v] + b_cos * b.y[v]);
            double p0 = waveIndex * (a0 * b.x[v] + b0 * b.y[v]);
            sincos(sin_t * psin + cos_t * pcos + p0, &vs[v], &vc[v]);
        }

        double sr = 0.0, si = 0.0;
        if (absB > absA)
        {
            for (int e = 0; e < nv; ++e)
            {
                if (!b.edge_valid_x[e]) continue;
                int en = (e + 1 < nv) ? e + 1 : 0;
                double Ci = A + b.slope_yx[e] * B;
                double absCi = fabs(Ci);
                double inv = (absCi > eps1) ? (1.0 / Ci) : 0.0;
                sr += (vc[en] - vc[e]) * inv;
                si += (vs[en] - vs[e]) * inv;
                if (absCi <= eps1)
                {
                    double p1x = b.x[e], p2x = b.x[en];
                    double tr = -wi2 * Ci * (p2x * p2x - p1x * p1x) * 0.5;
                    double ti = waveIndex * (p2x - p1x);
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
                double Ei = A * b.slope_xy[e] + B;
                double absEi = fabs(Ei);
                double inv = (absEi > eps1) ? (1.0 / Ei) : 0.0;
                sr += (vc[en] - vc[e]) * inv;
                si += (vs[en] - vs[e]) * inv;
                if (absEi <= eps1)
                {
                    double p1y = b.y[e], p2y = b.y[en];
                    double tr = -wi2 * Ei * (p2y * p2y - p1y * p1y) * 0.5;
                    double ti = waveIndex * (p2y - p1y);
                    sr += vc[e] * tr - vs[e] * ti;
                    si += vc[e] * ti + vs[e] * tr;
                }
            }
            double inv_nA = -1.0 / A;
            sr *= inv_nA; si *= inv_nA;
        }
        cmul(complWaveR, complWaveI, sr, si, fr, fi);
    }
    if (isnan(fr)) return false;

    double dpr = 1.0, dpi = 0.0;
    if (!b.isExternal)
    {
        double dp_sin = cp * b.cenx + sp * b.ceny;
        double dp_cos = -b.cenz;
        sincos(-waveIndex * (sin_t * dp_sin + cos_t * dp_cos), &dpi, &dpr);
    }

    double vfx = vf[(grid * 3) + 0];
    double vfy = vf[(grid * 3) + 1];
    double vfz = vf[(grid * 3) + 2];
    double r00, r01, r10, r11;
    rotate_jones_gpu(b.pNTx, b.pNTy, b.pNTz, b.pNPx, b.pNPy, b.pNPz,
                     b.pnxDTx, b.pnxDTy, b.pnxDTz,
                     b.pnxDPx, b.pnxDPy, b.pnxDPz,
                     vfx, vfy, vfz, dx, dy, dz,
                     r00, r01, r10, r11);

    double cpr, cpi;
    cmul(fr, fi, dpr, dpi, cpr, cpi);
    double sr00r = cpr * r00, sr00i = cpi * r00;
    double sr01r = cpr * r01, sr01i = cpi * r01;
    double sr10r = cpr * r10, sr10i = cpi * r10;
    double sr11r = cpr * r11, sr11i = cpi * r11;

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
                                   const double *sinTheta,
                                   const double *cosTheta,
                                   const double *sinPhi,
                                   const double *cosPhi,
                                   const double *vf,
                                   int nAz, int nZen,
                                   double waveIndex,
                                   double wi2,
                                   double eps1,
                                   double eps2,
                                   double complWaveR,
                                   double complWaveI,
                                   double invComplWaveR,
                                   double invComplWaveI,
                                   int legacySign,
                                   double *jFull,
                                   double *jNoShadow)
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

    double cp = cosPhi[p], sp = sinPhi[p];
    double sin_t = sinTheta[t], cos_t = cosTheta[t];
    double dx = sin_t * cp;
    double dy = sin_t * sp;
    double dz = -cos_t;

    double neg_cp = -cp, neg_sp = -sp;
    double a_sin = neg_cp * b.horAx + neg_sp * b.horAy;
    double a_cos = b.horAz;
    double a0 = b.bdx * b.horAx + b.bdy * b.horAy + b.bdz * b.horAz;
    double b_sin = neg_cp * b.verAx + neg_sp * b.verAy;
    double b_cos = b.verAz;
    double b0 = b.bdx * b.verAx + b.bdy * b.verAy + b.bdz * b.verAz;

    double A = sin_t * a_sin + cos_t * a_cos + a0;
    double B = sin_t * b_sin + cos_t * b_cos + b0;
    double absA = fabs(A);
    double absB = fabs(B);

    double fr, fi;
    if (absA < eps2 && absB < eps2)
    {
        double sign = legacySign ? 1.0 : -1.0;
        fr = sign * invComplWaveR * b.beam_area;
        fi = sign * invComplWaveI * b.beam_area;
    }
    else
    {
        double vc[32], vs[32];
        for (int v = 0; v < nv; ++v)
        {
            double psin = waveIndex * (a_sin * b.x[v] + b_sin * b.y[v]);
            double pcos = waveIndex * (a_cos * b.x[v] + b_cos * b.y[v]);
            double p0 = waveIndex * (a0 * b.x[v] + b0 * b.y[v]);
            sincos(sin_t * psin + cos_t * pcos + p0, &vs[v], &vc[v]);
        }

        double sr = 0.0, si = 0.0;
        if (absB > absA)
        {
            for (int e = 0; e < nv; ++e)
            {
                if (!b.edge_valid_x[e]) continue;
                int en = (e + 1 < nv) ? e + 1 : 0;
                double Ci = A + b.slope_yx[e] * B;
                double absCi = fabs(Ci);
                double inv = (absCi > eps1) ? (1.0 / Ci) : 0.0;
                sr += (vc[en] - vc[e]) * inv;
                si += (vs[en] - vs[e]) * inv;
                if (absCi <= eps1)
                {
                    double p1x = b.x[e], p2x = b.x[en];
                    double tr = -wi2 * Ci * (p2x * p2x - p1x * p1x) * 0.5;
                    double ti = waveIndex * (p2x - p1x);
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
                double Ei = A * b.slope_xy[e] + B;
                double absEi = fabs(Ei);
                double inv = (absEi > eps1) ? (1.0 / Ei) : 0.0;
                sr += (vc[en] - vc[e]) * inv;
                si += (vs[en] - vs[e]) * inv;
                if (absEi <= eps1)
                {
                    double p1y = b.y[e], p2y = b.y[en];
                    double tr = -wi2 * Ei * (p2y * p2y - p1y * p1y) * 0.5;
                    double ti = waveIndex * (p2y - p1y);
                    sr += vc[e] * tr - vs[e] * ti;
                    si += vc[e] * ti + vs[e] * tr;
                }
            }
            double inv_nA = -1.0 / A;
            sr *= inv_nA; si *= inv_nA;
        }
        cmul(complWaveR, complWaveI, sr, si, fr, fi);
    }
    if (isnan(fr)) return;

    double dpr = 1.0, dpi = 0.0;
    if (!b.isExternal)
    {
        double dp_sin = cp * b.cenx + sp * b.ceny;
        double dp_cos = -b.cenz;
        sincos(-waveIndex * (sin_t * dp_sin + cos_t * dp_cos), &dpi, &dpr);
    }

    double vfx = vf[(grid * 3) + 0];
    double vfy = vf[(grid * 3) + 1];
    double vfz = vf[(grid * 3) + 2];
    double r00, r01, r10, r11;
    rotate_jones_gpu(b.pNTx, b.pNTy, b.pNTz, b.pNPx, b.pNPy, b.pNPz,
                     b.pnxDTx, b.pnxDTy, b.pnxDTz,
                     b.pnxDPx, b.pnxDPy, b.pnxDPz,
                     vfx, vfy, vfz, dx, dy, dz,
                     r00, r01, r10, r11);

    double cpr, cpi;
    cmul(fr, fi, dpr, dpi, cpr, cpi);
    double sr00r = cpr * r00, sr00i = cpi * r00;
    double sr01r = cpr * r01, sr01i = cpi * r01;
    double sr10r = cpr * r10, sr10i = cpi * r10;
    double sr11r = cpr * r11, sr11i = cpi * r11;

    double d00r = sr00r * b.jp00r - sr00i * b.jp00i + sr01r * b.jp10r - sr01i * b.jp10i;
    double d00i = sr00r * b.jp00i + sr00i * b.jp00r + sr01r * b.jp10i + sr01i * b.jp10r;
    double d01r = sr00r * b.jp01r - sr00i * b.jp01i + sr01r * b.jp11r - sr01i * b.jp11i;
    double d01i = sr00r * b.jp01i + sr00i * b.jp01r + sr01r * b.jp11i + sr01i * b.jp11r;
    double d10r = sr10r * b.jp00r - sr10i * b.jp00i + sr11r * b.jp10r - sr11i * b.jp10i;
    double d10i = sr10r * b.jp00i + sr10i * b.jp00r + sr11r * b.jp10i + sr11i * b.jp10r;
    double d11r = sr10r * b.jp01r - sr10i * b.jp01i + sr11r * b.jp11r - sr11i * b.jp11i;
    double d11i = sr10r * b.jp01i + sr10i * b.jp01r + sr11r * b.jp11i + sr11i * b.jp11r;

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

__global__ void diffraction_grid_kernel(const GpuBeam *beams,
                                        const int *beamOffsets,
                                        const double *sinTheta,
                                        const double *cosTheta,
                                        const double *sinPhi,
                                        const double *cosPhi,
                                        const double *vf,
                                        int nAz, int nZen,
                                        int nOrient,
                                        double waveIndex,
                                        double wi2,
                                        double eps1,
                                        double eps2,
                                        double complWaveR,
                                        double complWaveI,
                                        double invComplWaveR,
                                        double invComplWaveI,
                                        int legacySign,
                                        double *jFull,
                                        double *jNoShadow)
{
    int gridCount = nAz * (nZen + 1);
    long long total = (long long)nOrient * gridCount;
    long long idx = (long long)blockIdx.x * blockDim.x + threadIdx.x;
    if (idx >= total) return;

    int orient = (int)(idx / gridCount);
    int grid = (int)(idx - (long long)orient * gridCount);
    int begin = beamOffsets[orient];
    int end = beamOffsets[orient + 1];

    double j00r = 0.0, j00i = 0.0;
    double j01r = 0.0, j01i = 0.0;
    double j10r = 0.0, j10i = 0.0;
    double j11r = 0.0, j11i = 0.0;
    double n00r = 0.0, n00i = 0.0;
    double n01r = 0.0, n01i = 0.0;
    double n10r = 0.0, n10i = 0.0;
    double n11r = 0.0, n11i = 0.0;

    for (int bi = begin; bi < end; ++bi)
    {
        double d00r, d00i, d01r, d01i, d10r, d10i, d11r, d11i;
        const GpuBeam &b = beams[bi];
        if (!compute_beam_jones_gpu(b, grid, nZen, sinTheta, cosTheta,
                                    sinPhi, cosPhi, vf, waveIndex, wi2,
                                    eps1, eps2, complWaveR, complWaveI,
                                    invComplWaveR, invComplWaveI, legacySign,
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

__global__ void mueller_batch_kernel(const double *jFull,
                                     const double *jNoShadow,
                                     const double *weights,
                                     int nOrient,
                                     int gridCount,
                                     double *mFull,
                                     double *mNoShadow)
{
    long long total = (long long)nOrient * gridCount;
    long long idx = (long long)blockIdx.x * blockDim.x + threadIdx.x;
    if (idx >= total) return;

    int orient = (int)(idx / gridCount);
    int grid = (int)(idx - (long long)orient * gridCount);
    double weight = weights[orient] * 0.25;
    const double *j = &jFull[idx * 8];
    double *out = &mFull[grid * 16];

    double j00r = j[0], j00i = j[1];
    double j01r = j[2], j01i = j[3];
    double j10r = j[4], j10i = j[5];
    double j11r = j[6], j11i = j[7];

    double a11 = j00r*j00r + j00i*j00i;
    double a12 = j01r*j01r + j01i*j01i;
    double a21 = j10r*j10r + j10i*j10i;
    double a22 = j11r*j11r + j11i*j11i;

    double A1 = a11 + a21, A2 = a12 + a22;
    atomicAdd(&out[0], ((A1 + A2) * 0.5) * weight);
    atomicAdd(&out[1], ((A1 - A2) * 0.5) * weight);
    A1 = a11 - a21; A2 = a12 - a22;
    atomicAdd(&out[4], ((A1 + A2) * 0.5) * weight);
    atomicAdd(&out[5], ((A1 - A2) * 0.5) * weight);

    double c1r = j00r*j01r + j00i*j01i;
    double c1i = j00i*j01r - j00r*j01i;
    double c2r = j11r*j10r + j11i*j10i;
    double c2i = j11i*j10r - j11r*j10i;
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

static void add_mueller_from_jones(const std::vector<double> &jones,
                                   double weight,
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
            const double *j = &jones[grid * 8];
            double j00r = j[0], j00i = j[1];
            double j01r = j[2], j01i = j[3];
            double j10r = j[4], j10i = j[5];
            double j11r = j[6], j11i = j[7];

            double a11 = j00r*j00r + j00i*j00i;
            double a12 = j01r*j01r + j01i*j01i;
            double a21 = j10r*j10r + j10i*j10i;
            double a22 = j11r*j11r + j11i*j11i;

            double A1 = a11 + a21, A2 = a12 + a22;
            out(p, t, 0, 0) += ((A1 + A2) * 0.5) * weight;
            out(p, t, 0, 1) += ((A1 - A2) * 0.5) * weight;
            A1 = a11 - a21; A2 = a12 - a22;
            out(p, t, 1, 0) += ((A1 + A2) * 0.5) * weight;
            out(p, t, 1, 1) += ((A1 - A2) * 0.5) * weight;

            double c1r = j00r*j01r + j00i*j01i;
            double c1i = j00i*j01r - j00r*j01i;
            double c2r = j11r*j10r + j11i*j10i;
            double c2i = j11i*j10r - j11r*j10i;
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

    std::vector<double> hJ(gridCount * 8, 0.0), hJns(gridCount * 8, 0.0);
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
        std::vector<double> hSinTheta(nZen + 1), hCosTheta(nZen + 1);
        for (int t = 0; t <= nZen; ++t)
        {
            double theta = m_sphere.GetZenith(t);
            hSinTheta[t] = sin(theta);
            hCosTheta[t] = cos(theta);
        }
        std::vector<double> hSinPhi(nAz), hCosPhi(nAz);
        for (int p = 0; p < nAz; ++p)
        {
            double phi = p * m_sphere.azinuthStep;
            hSinPhi[p] = sin(phi);
            hCosPhi[p] = cos(phi);
        }
        std::vector<double> hVf(gridCount * 3);
        for (int p = 0; p < nAz; ++p)
            for (int t = 0; t <= nZen; ++t)
            {
                int grid = p * (nZen + 1) + t;
                Point3d &v = m_sphere.vf[p][t];
                hVf[grid * 3 + 0] = v.x;
                hVf[grid * 3 + 1] = v.y;
                hVf[grid * 3 + 2] = v.z;
            }

        if (cudaMemcpy(ws.sinTheta, hSinTheta.data(), hSinTheta.size() * sizeof(double), cudaMemcpyHostToDevice) != cudaSuccess) return false;
        if (cudaMemcpy(ws.cosTheta, hCosTheta.data(), hCosTheta.size() * sizeof(double), cudaMemcpyHostToDevice) != cudaSuccess) return false;
        if (cudaMemcpy(ws.sinPhi, hSinPhi.data(), hSinPhi.size() * sizeof(double), cudaMemcpyHostToDevice) != cudaSuccess) return false;
        if (cudaMemcpy(ws.cosPhi, hCosPhi.data(), hCosPhi.size() * sizeof(double), cudaMemcpyHostToDevice) != cudaSuccess) return false;
        if (cudaMemcpy(ws.vf, hVf.data(), hVf.size() * sizeof(double), cudaMemcpyHostToDevice) != cudaSuccess) return false;
        ws.gridNAz = nAz;
        ws.gridNZen = nZen;
    }

    if (cudaMemcpy(ws.beams, hBeams.data(), hBeams.size() * sizeof(GpuBeam), cudaMemcpyHostToDevice) != cudaSuccess) return false;
    if (cudaMemset(ws.j, 0, hJ.size() * sizeof(double)) != cudaSuccess) return false;
    if (cudaMemset(ws.jNoShadow, 0, hJns.size() * sizeof(double)) != cudaSuccess) return false;

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

    if (cudaMemcpy(hJ.data(), ws.j, hJ.size() * sizeof(double), cudaMemcpyDeviceToHost) != cudaSuccess) return false;
    if (cudaMemcpy(hJns.data(), ws.jNoShadow, hJns.size() * sizeof(double), cudaMemcpyDeviceToHost) != cudaSuccess) return false;

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
        ((size_t)nZen + 1) * 2 * sizeof(double) +
        (size_t)nAz * 2 * sizeof(double) +
        gridCount * 3 * sizeof(double);
    const size_t muellerBytes = gridCount * 16 * 2 * sizeof(double);
    const size_t fixedBytes = gridVectorsBytes + muellerBytes + (256ULL << 20);
    const size_t jonesPerOrientation = gridCount * 8 * 2 * sizeof(double);

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

bool HandlerPO::HandleOrientationsToLocalGpu(const std::vector<PreparedOrientation> &prepared,
                                             Arr2D &localM,
                                             Arr2D &localM_noshadow)
{
    if (!isCoh)
        return false;

    const int nOrient = (int)prepared.size();
    const int nAz = m_sphere.nAzimuth;
    const int nZen = m_sphere.nZenith;
    const int gridCount = nAz * (nZen + 1);
    if (nOrient == 0 || gridCount == 0)
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
    if (nBeams == 0)
        return true;

    std::vector<GpuBeam> hBeams(nBeams);
    std::vector<double> hWeights(nOrient, 0.0);
    std::vector<int> hBeamOffsets(nOrient + 1, 0);
    size_t bi = 0;
    for (int oi = 0; oi < nOrient; ++oi)
    {
        hBeamOffsets[oi] = (int)bi;
        hWeights[oi] = prepared[oi].sinZenith;
        for (const PreparedBeam &pb : prepared[oi].beams)
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

    GpuWorkspace &ws = g_gpuWorkspace;
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
        std::vector<double> hSinTheta(nZen + 1), hCosTheta(nZen + 1);
        for (int t = 0; t <= nZen; ++t)
        {
            double theta = m_sphere.GetZenith(t);
            hSinTheta[t] = sin(theta);
            hCosTheta[t] = cos(theta);
        }
        std::vector<double> hSinPhi(nAz), hCosPhi(nAz);
        for (int p = 0; p < nAz; ++p)
        {
            double phi = p * m_sphere.azinuthStep;
            hSinPhi[p] = sin(phi);
            hCosPhi[p] = cos(phi);
        }
        std::vector<double> hVf(gridCount * 3);
        for (int p = 0; p < nAz; ++p)
            for (int t = 0; t <= nZen; ++t)
            {
                int grid = p * (nZen + 1) + t;
                Point3d &v = m_sphere.vf[p][t];
                hVf[grid * 3 + 0] = v.x;
                hVf[grid * 3 + 1] = v.y;
                hVf[grid * 3 + 2] = v.z;
            }

        if (cudaMemcpy(ws.sinTheta, hSinTheta.data(), hSinTheta.size() * sizeof(double), cudaMemcpyHostToDevice) != cudaSuccess) return false;
        if (cudaMemcpy(ws.cosTheta, hCosTheta.data(), hCosTheta.size() * sizeof(double), cudaMemcpyHostToDevice) != cudaSuccess) return false;
        if (cudaMemcpy(ws.sinPhi, hSinPhi.data(), hSinPhi.size() * sizeof(double), cudaMemcpyHostToDevice) != cudaSuccess) return false;
        if (cudaMemcpy(ws.cosPhi, hCosPhi.data(), hCosPhi.size() * sizeof(double), cudaMemcpyHostToDevice) != cudaSuccess) return false;
        if (cudaMemcpy(ws.vf, hVf.data(), hVf.size() * sizeof(double), cudaMemcpyHostToDevice) != cudaSuccess) return false;
        ws.gridNAz = nAz;
        ws.gridNZen = nZen;
    }

    if (cudaMemcpy(ws.beams, hBeams.data(), hBeams.size() * sizeof(GpuBeam), cudaMemcpyHostToDevice) != cudaSuccess) return false;
    if (cudaMemcpy(ws.weights, hWeights.data(), hWeights.size() * sizeof(double), cudaMemcpyHostToDevice) != cudaSuccess) return false;
    if (cudaMemcpy(ws.beamOffsets, hBeamOffsets.data(), hBeamOffsets.size() * sizeof(int), cudaMemcpyHostToDevice) != cudaSuccess) return false;
    if (cudaMemset(ws.j, 0, jCount * sizeof(double)) != cudaSuccess) return false;
    if (cudaMemset(ws.jNoShadow, 0, jCount * sizeof(double)) != cudaSuccess) return false;
    if (cudaMemset(ws.m, 0, mCount * sizeof(double)) != cudaSuccess) return false;
    if (cudaMemset(ws.mNoShadow, 0, mCount * sizeof(double)) != cudaSuccess) return false;

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
    if (cudaGetLastError() != cudaSuccess)
        return false;

    long long muellerTotal = (long long)nOrient * gridCount;
    int muellerGrid = (int)((muellerTotal + block - 1) / block);
    mueller_batch_kernel<<<muellerGrid, block>>>(ws.j, ws.jNoShadow, ws.weights,
                                                 nOrient, gridCount,
                                                 ws.m, ws.mNoShadow);
    if (cudaGetLastError() != cudaSuccess || cudaDeviceSynchronize() != cudaSuccess)
        return false;

    std::vector<double> hM(mCount), hMns(mCount);
    if (cudaMemcpy(hM.data(), ws.m, mCount * sizeof(double), cudaMemcpyDeviceToHost) != cudaSuccess) return false;
    if (cudaMemcpy(hMns.data(), ws.mNoShadow, mCount * sizeof(double), cudaMemcpyDeviceToHost) != cudaSuccess) return false;

    for (int p = 0; p < nAz; ++p)
        for (int t = 0; t <= nZen; ++t)
        {
            int grid = p * (nZen + 1) + t;
            const double *m = &hM[grid * 16];
            const double *mn = &hMns[grid * 16];
            for (int r = 0; r < 4; ++r)
                for (int c = 0; c < 4; ++c)
                {
                    localM(p, t, r, c) += m[r * 4 + c];
                    localM_noshadow(p, t, r, c) += mn[r * 4 + c];
                }
        }

    return true;
}
