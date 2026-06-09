#include "HandlerPO.h"

#include <cuda_runtime.h>
#include <cufft.h>
#include <cmath>
#include <algorithm>
#include <chrono>
#include <cstdlib>
#include <cfloat>
#include <cstdio>
#include <vector>
#include <utility>
#include <future>

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

__device__ inline void gpu_sincos_phase(double x, GpuReal *s, GpuReal *c)
{
#ifdef MBS_GPU_FLOAT
    double sd, cd;
    sincos(x, &sd, &cd);
    *s = (GpuReal)sd;
    *c = (GpuReal)cd;
#else
    sincos(x, s, c);
#endif
}

__device__ inline void gpu_atomic_add(GpuReal *address, GpuReal value)
{
    atomicAdd(address, value);
}

static inline double gpu_theta(const ScatteringRange &sphere, int t)
{
    return sphere.GetZenith(t);
}

static inline double gpu_grid_signature(const ScatteringRange &sphere, int nZen)
{
    double sig = sphere.azinuthStep * 0.12582917;
    sig += sphere.isNonUniform ? 0.31415926 : 0.27182818;
    for (int t = 0; t <= nZen; ++t)
    {
        const double theta = gpu_theta(sphere, t);
        sig += theta * (double)(t + 1) * 1.6180339887498948;
    }
    return sig;
}

struct GpuBeam
{
    GpuReal x[32], y[32];
    GpuReal slope_yx[32], slope_xy[32];
    // Compact edge index lists after host packing.
    unsigned char edge_valid_x[32], edge_valid_y[32];
    int nVertices;
    int nEdgeX, nEdgeY;
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
    double raw00r, raw00i, raw01r, raw01i;
    double raw10r, raw10i, raw11r, raw11i;
    double phasePath;
    double phaseSign;
    int absOffset;
    int absCount;
    int orientation;
};

struct GpuBeam8
{
    GpuReal x[8], y[8];
    GpuReal slope_yx[8], slope_xy[8];
    unsigned char edge_valid_x[8], edge_valid_y[8];
    int nVertices;
    int nEdgeX, nEdgeY;
    int isExternal;

    GpuReal bdx, bdy;
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
};

struct GpuWorkspace
{
    GpuBeam *beams = nullptr;
    GpuBeam8 *beams8 = nullptr;
    GpuReal *sinTheta = nullptr, *cosTheta = nullptr;
    GpuReal *sinPhi = nullptr, *cosPhi = nullptr;
    GpuReal *vf = nullptr;
    GpuReal *j = nullptr, *jNoShadow = nullptr;
    GpuReal *weights = nullptr;
    GpuReal *scales = nullptr;
    double *absPaths = nullptr;
    int *beamOffsets = nullptr;
    int *beamOffsets8 = nullptr;
    GpuReal *m = nullptr, *mNoShadow = nullptr, *mOrient = nullptr;
    size_t beamCap = 0;
    size_t beam8Cap = 0;
    size_t sinThetaCap = 0, cosThetaCap = 0;
    size_t sinPhiCap = 0, cosPhiCap = 0;
    size_t vfCap = 0, jCap = 0, jNoShadowCap = 0;
    size_t weightsCap = 0, scalesCap = 0, absPathsCap = 0, beamOffsetsCap = 0, mCap = 0, mNoShadowCap = 0, mOrientCap = 0;
    size_t beamOffsets8Cap = 0;
    int gridNAz = -1, gridNZen = -1;
    double gridSignature = 0.0;
    std::vector<GpuBeam> hBeams;
    std::vector<GpuBeam8> hBeams8;
    std::vector<GpuReal> hWeights;
    std::vector<GpuReal> hScales;
    std::vector<double> hAbsPaths;
    std::vector<int> hBeamOffsets;
    std::vector<int> hBeamOffsets8;
    std::vector<GpuReal> hM;
    std::vector<GpuReal> hMNoShadow;
    GpuComplex *fftLow = nullptr;
    GpuComplex *fftFull = nullptr;
    size_t fftLowCap = 0;
    size_t fftFullCap = 0;
    cufftHandle fftPlanLow = 0;
    cufftHandle fftPlanFull = 0;
    int fftNLow = 0;
    int fftNFull = 0;
    int fftBatch = 0;
    std::vector<GpuComplex> hFftLow;
    std::vector<GpuComplex> hFftFull;
};

static thread_local GpuWorkspace g_gpuWorkspace;
static thread_local bool g_gpuMultiWorker = false;

static inline double prepared_absorption_factor(const PreparedBeam &pb,
                                                double scale,
                                                double cAbs)
{
    if (pb.absorptionPaths.empty())
        return 1.0;
    double sum = 0.0;
    int count = 0;
    for (double path : pb.absorptionPaths)
    {
        sum += (path > DBL_EPSILON) ? std::exp(cAbs * path * scale) : 1.0;
        ++count;
    }
    return count > 0 ? sum / count : 1.0;
}

template <typename BeamT, int MaxEdges>
static inline void pack_prepared_gpu_beam(const PreparedBeam &pb,
                                          int orientation,
                                          double scale,
                                          double scale2,
                                          bool scaleOnPack,
                                          double waveIndex,
                                          double cAbs,
                                          BeamT &b)
{
    b.nVertices = pb.edgeData.nVertices;
    b.nEdgeX = 0;
    b.nEdgeY = 0;
    b.isExternal = pb.isExternal ? 1 : 0;
    b.orientation = orientation;
    for (int e = 0; e < MaxEdges; ++e)
    {
        b.x[e] = pb.edgeData.x[e] * scale;
        b.y[e] = pb.edgeData.y[e] * scale;
        b.slope_yx[e] = pb.edgeData.slope_yx[e];
        b.slope_xy[e] = pb.edgeData.slope_xy[e];
        const bool validX = e < b.nVertices && pb.edgeData.edge_valid_x[e];
        const bool validY = e < b.nVertices && pb.edgeData.edge_valid_y[e];
        if (validX)
            b.edge_valid_x[b.nEdgeX++] = (unsigned char)e;
        if (validY)
            b.edge_valid_y[b.nEdgeY++] = (unsigned char)e;
    }
    b.bdx = (GpuReal)(pb.bdx * pb.horAx + pb.bdy * pb.horAy + pb.bdz * pb.horAz);
    b.bdy = (GpuReal)(pb.bdx * pb.verAx + pb.bdy * pb.verAy + pb.bdz * pb.verAz);
    b.bdz = 0.0;
    b.horAx = pb.horAx; b.horAy = pb.horAy; b.horAz = pb.horAz;
    b.verAx = pb.verAx; b.verAy = pb.verAy; b.verAz = pb.verAz;
    b.cenx = pb.cenx * scale; b.ceny = pb.ceny * scale; b.cenz = pb.cenz * scale;
    b.beam_area = pb.beam_area * scale2;
    b.pNTx = pb.pNTx; b.pNTy = pb.pNTy; b.pNTz = pb.pNTz;
    b.pNPx = pb.pNPx; b.pNPy = pb.pNPy; b.pNPz = pb.pNPz;
    b.pnxDTx = pb.pnxDTx; b.pnxDTy = pb.pnxDTy; b.pnxDTz = pb.pnxDTz;
    b.pnxDPx = pb.pnxDPx; b.pnxDPy = pb.pnxDPy; b.pnxDPz = pb.pnxDPz;
    const double phasePath = pb.isExternal
        ? pb.origBeam.opticalPath
        : pb.info.projLenght;
    double phaseSign = 1.0;
    if (pb.isExternal)
        phaseSign = -phaseSign;
    if (!pb.isExternal && (pb.origBeam.nActs & 1))
        phaseSign = -phaseSign;
    b.raw00r = real(pb.origBeam.J.m11);
    b.raw00i = imag(pb.origBeam.J.m11);
    b.raw01r = real(pb.origBeam.J.m12);
    b.raw01i = imag(pb.origBeam.J.m12);
    b.raw10r = real(pb.origBeam.J.m21);
    b.raw10i = imag(pb.origBeam.J.m21);
    b.raw11r = real(pb.origBeam.J.m22);
    b.raw11i = imag(pb.origBeam.J.m22);
    b.phasePath = phasePath;
    b.phaseSign = phaseSign;
    b.absOffset = 0;
    b.absCount = 0;
    if (scaleOnPack)
    {
        const double path = phasePath * scale;
        const double arg = waveIndex * path;
        double sn = std::sin(arg);
        double cs = std::cos(arg);
        const double absorption = prepared_absorption_factor(pb, scale, cAbs);
        auto set_j = [&](const ::complex &jc, GpuReal &reOut, GpuReal &imOut)
        {
            const double jr = real(jc) * absorption;
            const double ji = imag(jc) * absorption;
            reOut = (GpuReal)(phaseSign * (jr * cs - ji * sn));
            imOut = (GpuReal)(phaseSign * (jr * sn + ji * cs));
        };
        set_j(pb.origBeam.J.m11, b.jp00r, b.jp00i);
        set_j(pb.origBeam.J.m12, b.jp01r, b.jp01i);
        set_j(pb.origBeam.J.m21, b.jp10r, b.jp10i);
        set_j(pb.origBeam.J.m22, b.jp11r, b.jp11i);
    }
    else
    {
        b.jp00r = pb.jp00r; b.jp00i = pb.jp00i;
        b.jp01r = pb.jp01r; b.jp01i = pb.jp01i;
        b.jp10r = pb.jp10r; b.jp10i = pb.jp10i;
        b.jp11r = pb.jp11r; b.jp11i = pb.jp11i;
    }
}

static inline void pack_prepared_gpu_beam8(const PreparedBeam &pb,
                                           int orientation,
                                           double scale,
                                           double scale2,
                                           bool scaleOnPack,
                                           double waveIndex,
                                           double cAbs,
                                           GpuBeam8 &b)
{
    b.nVertices = pb.edgeData.nVertices;
    b.nEdgeX = 0;
    b.nEdgeY = 0;
    b.isExternal = pb.isExternal ? 1 : 0;
    for (int e = 0; e < 8; ++e)
    {
        b.x[e] = pb.edgeData.x[e] * scale;
        b.y[e] = pb.edgeData.y[e] * scale;
        b.slope_yx[e] = pb.edgeData.slope_yx[e];
        b.slope_xy[e] = pb.edgeData.slope_xy[e];
        const bool validX = e < b.nVertices && pb.edgeData.edge_valid_x[e];
        const bool validY = e < b.nVertices && pb.edgeData.edge_valid_y[e];
        if (validX)
            b.edge_valid_x[b.nEdgeX++] = (unsigned char)e;
        if (validY)
            b.edge_valid_y[b.nEdgeY++] = (unsigned char)e;
    }
    b.bdx = (GpuReal)(pb.bdx * pb.horAx + pb.bdy * pb.horAy + pb.bdz * pb.horAz);
    b.bdy = (GpuReal)(pb.bdx * pb.verAx + pb.bdy * pb.verAy + pb.bdz * pb.verAz);
    b.horAx = pb.horAx; b.horAy = pb.horAy; b.horAz = pb.horAz;
    b.verAx = pb.verAx; b.verAy = pb.verAy; b.verAz = pb.verAz;
    b.cenx = pb.cenx * scale; b.ceny = pb.ceny * scale; b.cenz = pb.cenz * scale;
    b.beam_area = pb.beam_area * scale2;
    b.pNTx = pb.pNTx; b.pNTy = pb.pNTy; b.pNTz = pb.pNTz;
    b.pNPx = pb.pNPx; b.pNPy = pb.pNPy; b.pNPz = pb.pNPz;
    b.pnxDTx = pb.pnxDTx; b.pnxDTy = pb.pnxDTy; b.pnxDTz = pb.pnxDTz;
    b.pnxDPx = pb.pnxDPx; b.pnxDPy = pb.pnxDPy; b.pnxDPz = pb.pnxDPz;
    if (scaleOnPack)
    {
        const double path = pb.isExternal
            ? pb.origBeam.opticalPath * scale
            : pb.info.projLenght * scale;
        const double arg = waveIndex * path;
        double sn = std::sin(arg);
        double cs = std::cos(arg);
        double sign = 1.0;
        if (pb.isExternal)
            sign = -sign;
        if (!pb.isExternal && (pb.origBeam.nActs & 1))
            sign = -sign;
        const double absorption = prepared_absorption_factor(pb, scale, cAbs);
        auto set_j = [&](const ::complex &jc, GpuReal &reOut, GpuReal &imOut)
        {
            const double jr = real(jc) * absorption;
            const double ji = imag(jc) * absorption;
            reOut = (GpuReal)(sign * (jr * cs - ji * sn));
            imOut = (GpuReal)(sign * (jr * sn + ji * cs));
        };
        set_j(pb.origBeam.J.m11, b.jp00r, b.jp00i);
        set_j(pb.origBeam.J.m12, b.jp01r, b.jp01i);
        set_j(pb.origBeam.J.m21, b.jp10r, b.jp10i);
        set_j(pb.origBeam.J.m22, b.jp11r, b.jp11i);
    }
    else
    {
        b.jp00r = pb.jp00r; b.jp00i = pb.jp00i;
        b.jp01r = pb.jp01r; b.jp01i = pb.jp01i;
        b.jp10r = pb.jp10r; b.jp10i = pb.jp10i;
        b.jp11r = pb.jp11r; b.jp11i = pb.jp11i;
    }
    (void)orientation;
}

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

static int gpu_no_atomics_mode()
{
    const char *value = std::getenv("MBS_GPU_NO_ATOMICS");
    if (!value || !*value)
        return -1;
    if (value[0] == '1' && value[1] == '\0')
        return 1;
    if (value[0] == '0' && value[1] == '\0')
        return 0;
    return -1;
}

static int gpu_fused_mueller_mode()
{
    const char *value = std::getenv("MBS_GPU_FUSED_MUELLER");
    if (!value || !*value)
        return -1;
    if (value[0] == '1' && value[1] == '\0')
        return 1;
    if (value[0] == '0' && value[1] == '\0')
        return 0;
    return -1;
}

static int gpu_no_vertex_cache_mode()
{
    const char *value = std::getenv("MBS_GPU_NO_VERTEX_CACHE");
    if (!value || !*value)
        return -1;
    if (value[0] == '1' && value[1] == '\0')
        return 1;
    if (value[0] == '0' && value[1] == '\0')
        return 0;
    return -1;
}

static int gpu_stage_mueller_mode()
{
    const char *value = std::getenv("MBS_GPU_STAGE_MUELLER");
    if (!value || !*value)
        return -1;
    if (value[0] == '1' && value[1] == '\0')
        return 1;
    if (value[0] == '0' && value[1] == '\0')
        return 0;
    return -1;
}

static bool gpu_beam_stats_enabled()
{
    const char *value = std::getenv("MBS_GPU_BEAM_STATS");
    return value && value[0] == '1' && value[1] == '\0';
}

static GpuReal gpu_effective_eps1(double eps1)
{
    double value = eps1;
    if (const char *env = std::getenv("MBS_GPU_EPS1"))
    {
        char *end = nullptr;
        double parsed = std::strtod(env, &end);
        if (end && *end == '\0' && parsed > 0.0)
            value = parsed;
    }
#ifdef MBS_GPU_FLOAT
    value = std::max(value, 3e-5);
#endif
    return (GpuReal)value;
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

static int gpu_fft_theta_factor()
{
    const char *value = std::getenv("MBS_FFT_THETA_FACTOR");
    if (!value || !*value)
        return 1;
    char *end = nullptr;
    long parsed = std::strtol(value, &end, 10);
    if (!end || *end != '\0' || parsed < 1 || parsed > 16)
        return 1;
    return (int)parsed;
}

static bool gpu_fft_check_enabled()
{
    const char *value = std::getenv("MBS_FFT_CHECK");
    return value && value[0] == '1' && value[1] == '\0';
}

static bool gpu_timing_enabled()
{
    const char *value = std::getenv("MBS_GPU_TIMING");
    return value && value[0] == '1' && value[1] == '\0';
}

static bool gpu_multik_full_enabled()
{
    const char *value = std::getenv("MBS_GPU_MULTI_K_FULL");
    return value && value[0] == '1' && value[1] == '\0';
}

static bool gpu_fft_pair_enabled()
{
    const char *value = std::getenv("MBS_GPU_FFT_PAIR");
    return value && value[0] == '1' && value[1] == '\0';
}

static bool gpu_fft_debug_enabled()
{
    const char *value = std::getenv("MBS_FFT_DEBUG");
    return value && value[0] == '1' && value[1] == '\0';
}

static int gpu_multi_device_count(int workCount)
{
    if (workCount <= 1)
        return 1;
    const char *value = std::getenv("MBS_GPU_MULTI");
    if (value && value[0] == '0' && value[1] == '\0')
        return 1;

    int count = 0;
    if (cudaGetDeviceCount(&count) != cudaSuccess || count <= 1)
        return 1;

    int limit = count;
    if (value && *value)
    {
        char *end = nullptr;
        long parsed = std::strtol(value, &end, 10);
        if (end && *end == '\0' && parsed > 0)
            limit = (int)parsed;
    }
    if (const char *capEnv = std::getenv("MBS_GPU_MULTI_MAX"))
    {
        char *end = nullptr;
        long parsed = std::strtol(capEnv, &end, 10);
        if (end && *end == '\0' && parsed > 0)
            limit = std::min(limit, (int)parsed);
    }

    return std::max(1, std::min(std::min(count, limit), workCount));
}

static bool gpu_multi_debug_enabled()
{
    const char *value = std::getenv("MBS_GPU_MULTI_DEBUG");
    return value && value[0] == '1' && value[1] == '\0';
}

static bool gpu_fft_adaptive_phi_enabled()
{
    const char *value = std::getenv("MBS_FFT_ADAPTIVE_PHI");
    if (!value || !*value)
        return false;
    return !(value[0] == '0' && value[1] == '\0');
}

static double gpu_fft_adaptive_phi_threshold()
{
    const char *value = std::getenv("MBS_FFT_ADAPTIVE_PHI_THRESHOLD");
    if (!value || !*value)
        return 0.04;
    char *end = nullptr;
    double parsed = std::strtod(value, &end);
    if (!end || *end != '\0' || parsed <= 0.0)
        return 0.04;
    return parsed;
}

static int gpu_fft_adaptive_phi_max_rows()
{
    const char *value = std::getenv("MBS_FFT_ADAPTIVE_PHI_MAX_ROWS");
    if (!value || !*value)
        return 8;
    char *end = nullptr;
    long parsed = std::strtol(value, &end, 10);
    if (!end || *end != '\0' || parsed < 0 || parsed > 128)
        return 8;
    return (int)parsed;
}


static int gpu_block_size()
{
    const char *value = std::getenv("MBS_GPU_BLOCK");
    if (!value || !*value)
        return 128;
    char *end = nullptr;
    long parsed = std::strtol(value, &end, 10);
    if (!end || *end != '\0')
        return 128;
    if (parsed != 128 && parsed != 256 && parsed != 512)
        return 128;
    return (int)parsed;
}

static bool gpu_block_size_overridden()
{
    const char *value = std::getenv("MBS_GPU_BLOCK");
    return value && *value;
}

static double gpu_now_ms()
{
    using clock = std::chrono::steady_clock;
    return std::chrono::duration<double, std::milli>(
        clock::now().time_since_epoch()).count();
}

static bool gpu_report_cuda_error(cudaError_t err, const char *where)
{
    if (err == cudaSuccess)
        return true;
    std::fprintf(stderr, "CUDA backend error at %s: %s\n",
                 where, cudaGetErrorString(err));
    const char *allowFallback = std::getenv("MBS_GPU_ALLOW_FALLBACK");
    if (!(allowFallback && allowFallback[0] == '1' && allowFallback[1] == '\0'))
    {
        std::fprintf(stderr,
                     "FATAL: CUDA backend failed. Set MBS_GPU_ALLOW_FALLBACK=1 "
                     "to allow the legacy CPU fallback after GPU errors.\n");
        std::fflush(stderr);
        std::exit(2);
    }
    return false;
}

static int choose_fft_phi_factor(int nFull, int configured)
{
    int requested = gpu_fft_phi_factor();
    if (requested > 0)
        return requested;
    if (configured > 0)
        return configured;
    if (nFull >= 2400)
        return 8;
    if (nFull >= 600)
        return 6;
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

static bool ensure_fft_workspace(GpuWorkspace &ws, int nLow, int nFull, int batch)
{
    const size_t lowCount = (size_t)batch * nLow;
    const size_t fullCount = (size_t)batch * nFull;
    if (!ensure_device_capacity(ws.fftLow, ws.fftLowCap, lowCount))
        return false;
    if (!ensure_device_capacity(ws.fftFull, ws.fftFullCap, fullCount))
        return false;

    if (ws.fftPlanLow && ws.fftNLow == nLow && ws.fftNFull == nFull
        && ws.fftBatch == batch)
        return true;

    if (ws.fftPlanLow)
    {
        cufftDestroy(ws.fftPlanLow);
        ws.fftPlanLow = 0;
    }
    if (ws.fftPlanFull)
    {
        cufftDestroy(ws.fftPlanFull);
        ws.fftPlanFull = 0;
    }

    int nLowArr[1] = { nLow };
    int nFullArr[1] = { nFull };
    if (cufftPlanMany(&ws.fftPlanLow, 1, nLowArr,
                      nullptr, 1, nLow,
                      nullptr, 1, nLow,
                      GpuCufftType, batch) != CUFFT_SUCCESS)
        return false;
    if (cufftPlanMany(&ws.fftPlanFull, 1, nFullArr,
                      nullptr, 1, nFull,
                      nullptr, 1, nFull,
                      GpuCufftType, batch) != CUFFT_SUCCESS)
        return false;

    ws.fftNLow = nLow;
    ws.fftNFull = nFull;
    ws.fftBatch = batch;
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

template <int MaxVertices, typename BeamT>
__device__ inline bool compute_beam_integral_cached_gpu(const BeamT &b,
                                                              int nv,
                                                              GpuReal sin_t,
                                                              GpuReal cos_t,
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
                                                              GpuReal a_sin,
                                                              GpuReal a_cos,
                                                              GpuReal a0,
                                                              GpuReal b_sin,
                                                              GpuReal b_cos,
                                                              GpuReal b0,
                                                              GpuReal A,
                                                              GpuReal B,
                                                              GpuReal absA,
                                                              GpuReal absB,
                                                              GpuReal &fr,
                                                              GpuReal &fi)
{
    if (nv > MaxVertices)
        return false;
    if (absA < eps2 && absB < eps2)
    {
        GpuReal sign = legacySign ? 1.0 : -1.0;
        fr = sign * invComplWaveR * b.beam_area;
        fi = sign * invComplWaveI * b.beam_area;
        return true;
    }

    GpuReal vc[MaxVertices], vs[MaxVertices];
    for (int v = 0; v < nv; ++v)
    {
        const double x = (double)b.x[v];
        const double y = (double)b.y[v];
        const double psin = (double)waveIndex * ((double)a_sin * x + (double)b_sin * y);
        const double pcos = (double)waveIndex * ((double)a_cos * x + (double)b_cos * y);
        const double p0 = (double)waveIndex * ((double)a0 * x + (double)b0 * y);
        gpu_sincos_phase((double)sin_t * psin + (double)cos_t * pcos + p0, &vs[v], &vc[v]);
    }

    GpuReal sr = 0.0, si = 0.0;
    if (absB > absA)
    {
        const int nEdge = b.nEdgeX;
        for (int ii = 0; ii < nEdge; ++ii)
        {
            const int e = b.edge_valid_x[ii];
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
        const int nEdge = b.nEdgeY;
        for (int ii = 0; ii < nEdge; ++ii)
        {
            const int e = b.edge_valid_y[ii];
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
    return true;
}

template <typename BeamT, int MaxVertices>
__device__ inline bool compute_beam_jones_context_limited_gpu(const BeamT &b,
                                                              GpuReal cp,
                                                              GpuReal sp,
                                                              GpuReal sin_t,
                                                              GpuReal cos_t,
                                                              GpuReal dx,
                                                              GpuReal dy,
                                                              GpuReal dz,
                                                              GpuReal vfx,
                                                              GpuReal vfy,
                                                              GpuReal vfz,
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

    GpuReal neg_cp = -cp, neg_sp = -sp;
    GpuReal a_sin = neg_cp * b.horAx + neg_sp * b.horAy;
    GpuReal a_cos = b.horAz;
    GpuReal a0 = b.bdx;
    GpuReal b_sin = neg_cp * b.verAx + neg_sp * b.verAy;
    GpuReal b_cos = b.verAz;
    GpuReal b0 = b.bdy;

    GpuReal A = sin_t * a_sin + cos_t * a_cos + a0;
    GpuReal B = sin_t * b_sin + cos_t * b_cos + b0;
    GpuReal absA = fabs(A);
    GpuReal absB = fabs(B);

    GpuReal fr, fi;
    bool integralOk = compute_beam_integral_cached_gpu<MaxVertices>(
        b, nv, sin_t, cos_t, waveIndex, wi2, eps1, eps2,
        complWaveR, complWaveI, invComplWaveR, invComplWaveI,
        legacySign, singularCorrection,
        a_sin, a_cos, a0, b_sin, b_cos, b0,
        A, B, absA, absB, fr, fi);
    if (!integralOk)
        return false;
    if (isnan(fr)) return false;

    GpuReal dpr = 1.0, dpi = 0.0;
    if (!b.isExternal)
    {
        GpuReal dp_sin = cp * b.cenx + sp * b.ceny;
        GpuReal dp_cos = -b.cenz;
        gpu_sincos(-waveIndex * (sin_t * dp_sin + cos_t * dp_cos), &dpi, &dpr);
    }

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

__device__ inline bool compute_beam_jones_context_gpu(const GpuBeam &b,
                                                      GpuReal cp,
                                                      GpuReal sp,
                                                      GpuReal sin_t,
                                                      GpuReal cos_t,
                                                      GpuReal dx,
                                                      GpuReal dy,
                                                      GpuReal dz,
                                                      GpuReal vfx,
                                                      GpuReal vfy,
                                                      GpuReal vfz,
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
    return compute_beam_jones_context_limited_gpu<GpuBeam, 32>(
        b, cp, sp, sin_t, cos_t, dx, dy, dz, vfx, vfy, vfz,
        waveIndex, wi2, eps1, eps2, complWaveR, complWaveI,
        invComplWaveR, invComplWaveI, legacySign, singularCorrection,
        d00r, d00i, d01r, d01i, d10r, d10i, d11r, d11i);
}

__device__ inline void phase_raw_jones_scaled(const GpuBeam &b,
                                              GpuReal waveIndex,
                                              GpuReal scale,
                                              const double *__restrict__ absPaths,
                                              GpuReal cAbs,
                                              GpuReal &jp00r, GpuReal &jp00i,
                                              GpuReal &jp01r, GpuReal &jp01i,
                                              GpuReal &jp10r, GpuReal &jp10i,
                                              GpuReal &jp11r, GpuReal &jp11i)
{
    double sn, cs;
    sincos((double)waveIndex * b.phasePath * (double)scale, &sn, &cs);
    const double sign = b.phaseSign;
    double absorption = 1.0;
    if (absPaths && b.absCount > 0)
    {
        double sum = 0.0;
        for (int i = 0; i < b.absCount; ++i)
        {
            const double path = absPaths[b.absOffset + i];
            sum += (path > DBL_EPSILON)
                ? exp((double)cAbs * path * (double)scale)
                : 1.0;
        }
        absorption = sum / (double)b.absCount;
    }
    jp00r = (GpuReal)(sign * absorption * (b.raw00r * cs - b.raw00i * sn));
    jp00i = (GpuReal)(sign * absorption * (b.raw00r * sn + b.raw00i * cs));
    jp01r = (GpuReal)(sign * absorption * (b.raw01r * cs - b.raw01i * sn));
    jp01i = (GpuReal)(sign * absorption * (b.raw01r * sn + b.raw01i * cs));
    jp10r = (GpuReal)(sign * absorption * (b.raw10r * cs - b.raw10i * sn));
    jp10i = (GpuReal)(sign * absorption * (b.raw10r * sn + b.raw10i * cs));
    jp11r = (GpuReal)(sign * absorption * (b.raw11r * cs - b.raw11i * sn));
    jp11i = (GpuReal)(sign * absorption * (b.raw11r * sn + b.raw11i * cs));
}

__device__ inline bool compute_beam_jones_context_gpu_multik(
                                                      const GpuBeam &b,
                                                      GpuReal cp,
                                                      GpuReal sp,
                                                      GpuReal sin_t,
                                                      GpuReal cos_t,
                                                      GpuReal dx,
                                                      GpuReal dy,
                                                      GpuReal dz,
                                                      GpuReal vfx,
                                                      GpuReal vfy,
                                                      GpuReal vfz,
                                                      GpuReal waveIndex,
                                                      GpuReal wi2,
                                                      GpuReal eps1,
                                                      GpuReal eps2,
                                                      GpuReal complWaveR,
                                                      GpuReal complWaveI,
                                                      GpuReal invComplWaveR,
                                                      GpuReal invComplWaveI,
                                                      int legacySign,
                                                      GpuReal scale,
                                                      const double *__restrict__ absPaths,
                                                      GpuReal cAbs,
                                                      GpuReal &d00r, GpuReal &d00i,
                                                      GpuReal &d01r, GpuReal &d01i,
                                                      GpuReal &d10r, GpuReal &d10i,
                                                      GpuReal &d11r, GpuReal &d11i)
{
    int nv = b.nVertices;
    if (nv <= 0) return false;

    GpuReal neg_cp = -cp, neg_sp = -sp;
    GpuReal a_sin = neg_cp * b.horAx + neg_sp * b.horAy;
    GpuReal a_cos = b.horAz;
    GpuReal a0 = b.bdx;
    GpuReal b_sin = neg_cp * b.verAx + neg_sp * b.verAy;
    GpuReal b_cos = b.verAz;
    GpuReal b0 = b.bdy;

    GpuReal A = sin_t * a_sin + cos_t * a_cos + a0;
    GpuReal B = sin_t * b_sin + cos_t * b_cos + b0;
    GpuReal absA = fabs(A);
    GpuReal absB = fabs(B);
    const GpuReal scale2 = scale * scale;
    const GpuReal waveIndexEff = waveIndex * scale;
    const GpuReal wi2Eff = wi2 * scale2;

    GpuReal fr, fi;
    bool integralOk = compute_beam_integral_cached_gpu<32>(
        b, nv, sin_t, cos_t, waveIndexEff, wi2Eff, eps1, eps2,
        complWaveR, complWaveI, invComplWaveR, invComplWaveI,
        legacySign, 1,
        a_sin, a_cos, a0, b_sin, b_cos, b0,
        A, B, absA, absB, fr, fi);
    if (!integralOk)
        return false;
    if (absA < eps2 && absB < eps2)
    {
        fr *= scale2;
        fi *= scale2;
    }
    if (isnan(fr)) return false;

    GpuReal dpr = 1.0, dpi = 0.0;
    if (!b.isExternal)
    {
        GpuReal dp_sin = cp * b.cenx + sp * b.ceny;
        GpuReal dp_cos = -b.cenz;
        gpu_sincos(-waveIndexEff * (sin_t * dp_sin + cos_t * dp_cos), &dpi, &dpr);
    }

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

    GpuReal jp00r, jp00i, jp01r, jp01i, jp10r, jp10i, jp11r, jp11i;
    phase_raw_jones_scaled(b, waveIndex, scale, absPaths, cAbs,
                           jp00r, jp00i, jp01r, jp01i,
                           jp10r, jp10i, jp11r, jp11i);

    d00r = sr00r * jp00r - sr00i * jp00i + sr01r * jp10r - sr01i * jp10i;
    d00i = sr00r * jp00i + sr00i * jp00r + sr01r * jp10i + sr01i * jp10r;
    d01r = sr00r * jp01r - sr00i * jp01i + sr01r * jp11r - sr01i * jp11i;
    d01i = sr00r * jp01i + sr00i * jp01r + sr01r * jp11i + sr01i * jp11r;
    d10r = sr10r * jp00r - sr10i * jp00i + sr11r * jp10r - sr11i * jp10i;
    d10i = sr10r * jp00i + sr10i * jp00r + sr11r * jp10i + sr11i * jp10r;
    d11r = sr10r * jp01r - sr10i * jp01i + sr11r * jp11r - sr11i * jp11i;
    d11i = sr10r * jp01i + sr10i * jp01r + sr11r * jp11i + sr11i * jp11r;
    return true;
}

__device__ inline bool compute_beam_jones_context_gpu8_auto(const GpuBeam8 &b,
                                                            GpuReal cp,
                                                            GpuReal sp,
                                                            GpuReal sin_t,
                                                            GpuReal cos_t,
                                                            GpuReal dx,
                                                            GpuReal dy,
                                                            GpuReal dz,
                                                            GpuReal vfx,
                                                            GpuReal vfy,
                                                            GpuReal vfz,
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
    if (b.nVertices <= 4)
    {
        return compute_beam_jones_context_limited_gpu<GpuBeam8, 4>(
            b, cp, sp, sin_t, cos_t, dx, dy, dz, vfx, vfy, vfz,
            waveIndex, wi2, eps1, eps2, complWaveR, complWaveI,
            invComplWaveR, invComplWaveI, legacySign, singularCorrection,
            d00r, d00i, d01r, d01i, d10r, d10i, d11r, d11i);
    }
    return compute_beam_jones_context_limited_gpu<GpuBeam8, 8>(
        b, cp, sp, sin_t, cos_t, dx, dy, dz, vfx, vfy, vfz,
        waveIndex, wi2, eps1, eps2, complWaveR, complWaveI,
        invComplWaveR, invComplWaveI, legacySign, singularCorrection,
        d00r, d00i, d01r, d01i, d10r, d10i, d11r, d11i);
}

__device__ inline bool compute_beam_jones_gpu(const GpuBeam &b,
                                              int grid,
                                              int nZen,
                                              const GpuReal *__restrict__ sinTheta,
                                              const GpuReal *__restrict__ cosTheta,
                                              const GpuReal *__restrict__ sinPhi,
                                              const GpuReal *__restrict__ cosPhi,
                                              const GpuReal *__restrict__ vf,
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
    int p = grid / (nZen + 1);
    int t = grid - p * (nZen + 1);
    GpuReal cp = cosPhi[p], sp = sinPhi[p];
    GpuReal sin_t = sinTheta[t], cos_t = cosTheta[t];
    GpuReal dx = sin_t * cp;
    GpuReal dy = sin_t * sp;
    GpuReal dz = -cos_t;
    GpuReal vfx = vf[(grid * 3) + 0];
    GpuReal vfy = vf[(grid * 3) + 1];
    GpuReal vfz = vf[(grid * 3) + 2];
    return compute_beam_jones_context_gpu(
        b, cp, sp, sin_t, cos_t, dx, dy, dz, vfx, vfy, vfz,
        waveIndex, wi2, eps1, eps2, complWaveR, complWaveI,
        invComplWaveR, invComplWaveI, legacySign, singularCorrection,
        d00r, d00i, d01r, d01i, d10r, d10i, d11r, d11i);
}

__device__ inline void beam_vertex_sincos_gpu(const GpuBeam &b,
                                              int v,
                                              GpuReal sin_t,
                                              GpuReal cos_t,
                                              GpuReal waveIndex,
                                              GpuReal a_sin,
                                              GpuReal a_cos,
                                              GpuReal a0,
                                              GpuReal b_sin,
                                              GpuReal b_cos,
                                              GpuReal b0,
                                              GpuReal *vs,
                                              GpuReal *vc)
{
    const double x = (double)b.x[v];
    const double y = (double)b.y[v];
    const double psin = (double)waveIndex * ((double)a_sin * x + (double)b_sin * y);
    const double pcos = (double)waveIndex * ((double)a_cos * x + (double)b_cos * y);
    const double p0 = (double)waveIndex * ((double)a0 * x + (double)b0 * y);
    gpu_sincos_phase((double)sin_t * psin + (double)cos_t * pcos + p0, vs, vc);
}

__device__ inline bool compute_beam_jones_nocache_gpu(const GpuBeam &b,
                                                      int grid,
                                                      int nZen,
                                                      const GpuReal *__restrict__ sinTheta,
                                                      const GpuReal *__restrict__ cosTheta,
                                                      const GpuReal *__restrict__ sinPhi,
                                                      const GpuReal *__restrict__ cosPhi,
                                                      const GpuReal *__restrict__ vf,
                                                      GpuReal waveIndex,
                                                      GpuReal wi2,
                                                      GpuReal eps1,
                                                      GpuReal eps2,
                                                      GpuReal complWaveR,
                                                      GpuReal complWaveI,
                                                      GpuReal invComplWaveR,
                                                      GpuReal invComplWaveI,
                                                      int legacySign,
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
    GpuReal a0 = b.bdx;
    GpuReal b_sin = neg_cp * b.verAx + neg_sp * b.verAy;
    GpuReal b_cos = b.verAz;
    GpuReal b0 = b.bdy;

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
        GpuReal sr = 0.0, si = 0.0;
        if (absB > absA)
        {
            const int nEdge = b.nEdgeX;
            int cachedVertex = -1;
            GpuReal cachedS = 0.0, cachedC = 0.0;
            for (int ii = 0; ii < nEdge; ++ii)
            {
                const int e = b.edge_valid_x[ii];
                int en = (e + 1 < nv) ? e + 1 : 0;
                GpuReal vce, vse, vcn, vsn;
                if (e == cachedVertex)
                {
                    vse = cachedS;
                    vce = cachedC;
                }
                else
                {
                    beam_vertex_sincos_gpu(b, e, sin_t, cos_t, waveIndex,
                                           a_sin, a_cos, a0, b_sin, b_cos, b0,
                                           &vse, &vce);
                }
                beam_vertex_sincos_gpu(b, en, sin_t, cos_t, waveIndex,
                                       a_sin, a_cos, a0, b_sin, b_cos, b0,
                                       &vsn, &vcn);
                cachedVertex = en;
                cachedS = vsn;
                cachedC = vcn;
                GpuReal Ci = A + b.slope_yx[e] * B;
                GpuReal absCi = fabs(Ci);
                GpuReal inv = (absCi > eps1) ? (1.0 / Ci) : 0.0;
                sr += (vcn - vce) * inv;
                si += (vsn - vse) * inv;
                if (absCi <= eps1)
                {
                    GpuReal p1x = b.x[e], p2x = b.x[en];
                    GpuReal tr = -wi2 * Ci * (p2x * p2x - p1x * p1x) * 0.5;
                    GpuReal ti = waveIndex * (p2x - p1x);
                    sr += vce * tr - vse * ti;
                    si += vce * ti + vse * tr;
                }
            }
            sr /= B; si /= B;
        }
        else
        {
            const int nEdge = b.nEdgeY;
            int cachedVertex = -1;
            GpuReal cachedS = 0.0, cachedC = 0.0;
            for (int ii = 0; ii < nEdge; ++ii)
            {
                const int e = b.edge_valid_y[ii];
                int en = (e + 1 < nv) ? e + 1 : 0;
                GpuReal vce, vse, vcn, vsn;
                if (e == cachedVertex)
                {
                    vse = cachedS;
                    vce = cachedC;
                }
                else
                {
                    beam_vertex_sincos_gpu(b, e, sin_t, cos_t, waveIndex,
                                           a_sin, a_cos, a0, b_sin, b_cos, b0,
                                           &vse, &vce);
                }
                beam_vertex_sincos_gpu(b, en, sin_t, cos_t, waveIndex,
                                       a_sin, a_cos, a0, b_sin, b_cos, b0,
                                       &vsn, &vcn);
                cachedVertex = en;
                cachedS = vsn;
                cachedC = vcn;
                GpuReal Ei = A * b.slope_xy[e] + B;
                GpuReal absEi = fabs(Ei);
                GpuReal inv = (absEi > eps1) ? (1.0 / Ei) : 0.0;
                sr += (vcn - vce) * inv;
                si += (vsn - vse) * inv;
                if (absEi <= eps1)
                {
                    GpuReal p1y = b.y[e], p2y = b.y[en];
                    GpuReal tr = -wi2 * Ei * (p2y * p2y - p1y * p1y) * 0.5;
                    GpuReal ti = waveIndex * (p2y - p1y);
                    sr += vce * tr - vse * ti;
                    si += vce * ti + vse * tr;
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

__device__ inline void mueller_accum_from_jones(GpuReal j00r, GpuReal j00i,
                                                GpuReal j01r, GpuReal j01i,
                                                GpuReal j10r, GpuReal j10i,
                                                GpuReal j11r, GpuReal j11i,
                                                GpuReal weight,
                                                GpuReal *__restrict__ out)
{
    GpuReal a11 = j00r*j00r + j00i*j00i;
    GpuReal a12 = j01r*j01r + j01i*j01i;
    GpuReal a21 = j10r*j10r + j10i*j10i;
    GpuReal a22 = j11r*j11r + j11i*j11i;

    const GpuReal half = (GpuReal)0.5;
    GpuReal A1 = a11 + a21, A2 = a12 + a22;
    gpu_atomic_add(&out[0], ((A1 + A2) * half) * weight);
    gpu_atomic_add(&out[1], ((A1 - A2) * half) * weight);
    A1 = a11 - a21; A2 = a12 - a22;
    gpu_atomic_add(&out[4], ((A1 + A2) * half) * weight);
    gpu_atomic_add(&out[5], ((A1 - A2) * half) * weight);

    GpuReal c1r = j00r*j01r + j00i*j01i;
    GpuReal c1i = j00i*j01r - j00r*j01i;
    GpuReal c2r = j11r*j10r + j11i*j10i;
    GpuReal c2i = j11i*j10r - j11r*j10i;
    gpu_atomic_add(&out[2], (-c1r - c2r) * weight);
    gpu_atomic_add(&out[3], ( c2i - c1i) * weight);
    gpu_atomic_add(&out[6], ( c2r - c1r) * weight);
    gpu_atomic_add(&out[7], (-c1i - c2i) * weight);

    c1r = j00r*j10r + j00i*j10i;
    c1i = j00i*j10r - j00r*j10i;
    c2r = j11r*j01r + j11i*j01i;
    c2i = j11i*j01r - j11r*j01i;
    gpu_atomic_add(&out[8], (-c1r - c2r) * weight);
    gpu_atomic_add(&out[9], ( c2r - c1r) * weight);
    gpu_atomic_add(&out[12], ( c1i - c2i) * weight);
    gpu_atomic_add(&out[13], ( c2i + c1i) * weight);

    c1r = j00r*j11r + j00i*j11i;
    c1i = j00i*j11r - j00r*j11i;
    c2r = j01r*j10r + j01i*j10i;
    c2i = j01i*j10r - j01r*j10i;
    gpu_atomic_add(&out[10], ( c1r + c2r) * weight);
    gpu_atomic_add(&out[11], ( c1i - c2i) * weight);
    gpu_atomic_add(&out[14], (-c1i - c2i) * weight);
    gpu_atomic_add(&out[15], ( c1r - c2r) * weight);
}

__device__ inline void mueller_store_from_jones(GpuReal j00r, GpuReal j00i,
                                                GpuReal j01r, GpuReal j01i,
                                                GpuReal j10r, GpuReal j10i,
                                                GpuReal j11r, GpuReal j11i,
                                                GpuReal weight,
                                                GpuReal *__restrict__ out)
{
    GpuReal a11 = j00r*j00r + j00i*j00i;
    GpuReal a12 = j01r*j01r + j01i*j01i;
    GpuReal a21 = j10r*j10r + j10i*j10i;
    GpuReal a22 = j11r*j11r + j11i*j11i;

    const GpuReal half = (GpuReal)0.5;
    GpuReal A1 = a11 + a21, A2 = a12 + a22;
    out[0] = ((A1 + A2) * half) * weight;
    out[1] = ((A1 - A2) * half) * weight;
    A1 = a11 - a21; A2 = a12 - a22;
    out[4] = ((A1 + A2) * half) * weight;
    out[5] = ((A1 - A2) * half) * weight;

    GpuReal c1r = j00r*j01r + j00i*j01i;
    GpuReal c1i = j00i*j01r - j00r*j01i;
    GpuReal c2r = j11r*j10r + j11i*j10i;
    GpuReal c2i = j11i*j10r - j11r*j10i;
    out[2] = (-c1r - c2r) * weight;
    out[3] = ( c2i - c1i) * weight;
    out[6] = ( c2r - c1r) * weight;
    out[7] = (-c1i - c2i) * weight;

    c1r = j00r*j10r + j00i*j10i;
    c1i = j00i*j10r - j00r*j10i;
    c2r = j11r*j01r + j11i*j01i;
    c2i = j11i*j01r - j11r*j01i;
    out[8] = (-c1r - c2r) * weight;
    out[9] = ( c2r - c1r) * weight;
    out[12] = ( c1i - c2i) * weight;
    out[13] = ( c2i + c1i) * weight;

    c1r = j00r*j11r + j00i*j11i;
    c1i = j00i*j11r - j00r*j11i;
    c2r = j01r*j10r + j01i*j10i;
    c2i = j01i*j10r - j01r*j10i;
    out[10] = ( c1r + c2r) * weight;
    out[11] = ( c1i - c2i) * weight;
    out[14] = (-c1i - c2i) * weight;
    out[15] = ( c1r - c2r) * weight;
}

__global__ void diffraction_kernel(const GpuBeam *__restrict__ beams, int nBeams,
                                   const GpuReal *__restrict__ sinTheta,
                                   const GpuReal *__restrict__ cosTheta,
                                   const GpuReal *__restrict__ sinPhi,
                                   const GpuReal *__restrict__ cosPhi,
                                   const GpuReal *__restrict__ vf,
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
                                   GpuReal *__restrict__ jFull,
                                   GpuReal *__restrict__ jNoShadow)
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
    GpuReal a0 = b.bdx;
    GpuReal b_sin = neg_cp * b.verAx + neg_sp * b.verAy;
    GpuReal b_cos = b.verAz;
    GpuReal b0 = b.bdy;

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
            const double x = (double)b.x[v];
            const double y = (double)b.y[v];
            const double psin = (double)waveIndex * ((double)a_sin * x + (double)b_sin * y);
            const double pcos = (double)waveIndex * ((double)a_cos * x + (double)b_cos * y);
            const double p0 = (double)waveIndex * ((double)a0 * x + (double)b0 * y);
            gpu_sincos_phase((double)sin_t * psin + (double)cos_t * pcos + p0, &vs[v], &vc[v]);
        }

        GpuReal sr = 0.0, si = 0.0;
        if (absB > absA)
        {
            const int nEdge = b.nEdgeX;
            for (int ii = 0; ii < nEdge; ++ii)
            {
                const int e = b.edge_valid_x[ii];
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
            const int nEdge = b.nEdgeY;
            for (int ii = 0; ii < nEdge; ++ii)
            {
                const int e = b.edge_valid_y[ii];
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
    if (jNoShadow && !b.isExternal)
    {
        atomicAdd(&jNoShadow[off + 0], d00r); atomicAdd(&jNoShadow[off + 1], d00i);
        atomicAdd(&jNoShadow[off + 2], d01r); atomicAdd(&jNoShadow[off + 3], d01i);
        atomicAdd(&jNoShadow[off + 4], d10r); atomicAdd(&jNoShadow[off + 5], d10i);
        atomicAdd(&jNoShadow[off + 6], d11r); atomicAdd(&jNoShadow[off + 7], d11i);
    }
}

__global__ void theta_m11_kernel(const GpuBeam *__restrict__ beams,
                                 int nBeams,
                                 const GpuReal *__restrict__ weights,
                                 const GpuReal *__restrict__ sinTheta,
                                 const GpuReal *__restrict__ cosTheta,
                                 const GpuReal *__restrict__ sinPhi,
                                 const GpuReal *__restrict__ cosPhi,
                                 const GpuReal *__restrict__ vf,
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
                                 long long startIdx,
                                 GpuReal *__restrict__ m11)
{
    int gridCount = nAz * nTheta;
    long long total = (long long)nBeams * gridCount;
    long long idx = startIdx + (long long)blockIdx.x * blockDim.x + threadIdx.x;
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

__global__ void diffraction_grid_kernel(const GpuBeam *__restrict__ beams,
                                        const int *__restrict__ beamOffsets,
                                        const GpuReal *__restrict__ sinTheta,
                                        const GpuReal *__restrict__ cosTheta,
                                        const GpuReal *__restrict__ sinPhi,
                                        const GpuReal *__restrict__ cosPhi,
                                        const GpuReal *__restrict__ vf,
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
                                        GpuReal *__restrict__ jFull,
                                        GpuReal *__restrict__ jNoShadow)
{
    int gridCount = nAz * (nZen + 1);
    long long total = (long long)nOrient * gridCount;
    long long idx = (long long)blockIdx.x * blockDim.x + threadIdx.x;
    if (idx >= total) return;

    int orient = (int)(idx / gridCount);
    int grid = (int)(idx - (long long)orient * gridCount);
    int begin = beamOffsets[orient];
    int end = beamOffsets[orient + 1];
    int p = grid / (nZen + 1);
    int t = grid - p * (nZen + 1);
    GpuReal cp = cosPhi[p], sp = sinPhi[p];
    GpuReal sin_t = sinTheta[t], cos_t = cosTheta[t];
    GpuReal dx = sin_t * cp;
    GpuReal dy = sin_t * sp;
    GpuReal dz = -cos_t;
    GpuReal vfx = vf[(grid * 3) + 0];
    GpuReal vfy = vf[(grid * 3) + 1];
    GpuReal vfz = vf[(grid * 3) + 2];

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
        if (!compute_beam_jones_context_gpu(
                b, cp, sp, sin_t, cos_t, dx, dy, dz, vfx, vfy, vfz,
                waveIndex, wi2, eps1, eps2, complWaveR, complWaveI,
                invComplWaveR, invComplWaveI, legacySign, 1,
                d00r, d00i, d01r, d01i, d10r, d10i, d11r, d11i))
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
    if (jNoShadow)
    {
        jNoShadow[off + 0] = n00r; jNoShadow[off + 1] = n00i;
        jNoShadow[off + 2] = n01r; jNoShadow[off + 3] = n01i;
        jNoShadow[off + 4] = n10r; jNoShadow[off + 5] = n10i;
        jNoShadow[off + 6] = n11r; jNoShadow[off + 7] = n11i;
    }
}

__global__ void mueller_batch_kernel(const GpuReal *__restrict__ jFull,
                                     const GpuReal *__restrict__ jNoShadow,
                                     const GpuReal *__restrict__ weights,
                                     int nOrient,
                                     int gridCount,
                                     GpuReal *__restrict__ mFull,
                                     GpuReal *__restrict__ mNoShadow)
{
    long long total = (long long)nOrient * gridCount;
    long long idx = (long long)blockIdx.x * blockDim.x + threadIdx.x;
    if (idx >= total) return;

    int orient = (int)(idx / gridCount);
    int grid = (int)(idx - (long long)orient * gridCount);
    const GpuReal half = (GpuReal)0.5;
    GpuReal weight = weights[orient] * (GpuReal)0.25;
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
    gpu_atomic_add(&out[0], ((A1 + A2) * half) * weight);
    gpu_atomic_add(&out[1], ((A1 - A2) * half) * weight);
    A1 = a11 - a21; A2 = a12 - a22;
    gpu_atomic_add(&out[4], ((A1 + A2) * half) * weight);
    gpu_atomic_add(&out[5], ((A1 - A2) * half) * weight);

    GpuReal c1r = j00r*j01r + j00i*j01i;
    GpuReal c1i = j00i*j01r - j00r*j01i;
    GpuReal c2r = j11r*j10r + j11i*j10i;
    GpuReal c2i = j11i*j10r - j11r*j10i;
    gpu_atomic_add(&out[2], (-c1r - c2r) * weight);
    gpu_atomic_add(&out[3], ( c2i - c1i) * weight);
    gpu_atomic_add(&out[6], ( c2r - c1r) * weight);
    gpu_atomic_add(&out[7], (-c1i - c2i) * weight);

    c1r = j00r*j10r + j00i*j10i;
    c1i = j00i*j10r - j00r*j10i;
    c2r = j11r*j01r + j11i*j01i;
    c2i = j11i*j01r - j11r*j01i;
    gpu_atomic_add(&out[8], (-c1r - c2r) * weight);
    gpu_atomic_add(&out[9], ( c2r - c1r) * weight);
    gpu_atomic_add(&out[12], ( c1i - c2i) * weight);
    gpu_atomic_add(&out[13], ( c2i + c1i) * weight);

    c1r = j00r*j11r + j00i*j11i;
    c1i = j00i*j11r - j00r*j11i;
    c2r = j01r*j10r + j01i*j10i;
    c2i = j01i*j10r - j01r*j10i;
    gpu_atomic_add(&out[10], ( c1r + c2r) * weight);
    gpu_atomic_add(&out[11], ( c1i - c2i) * weight);
    gpu_atomic_add(&out[14], (-c1i - c2i) * weight);
    gpu_atomic_add(&out[15], ( c1r - c2r) * weight);

    if (!jNoShadow || !mNoShadow)
        return;

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
    gpu_atomic_add(&out[0], ((A1 + A2) * half) * weight);
    gpu_atomic_add(&out[1], ((A1 - A2) * half) * weight);
    A1 = a11 - a21; A2 = a12 - a22;
    gpu_atomic_add(&out[4], ((A1 + A2) * half) * weight);
    gpu_atomic_add(&out[5], ((A1 - A2) * half) * weight);

    c1r = j00r*j01r + j00i*j01i;
    c1i = j00i*j01r - j00r*j01i;
    c2r = j11r*j10r + j11i*j10i;
    c2i = j11i*j10r - j11r*j10i;
    gpu_atomic_add(&out[2], (-c1r - c2r) * weight);
    gpu_atomic_add(&out[3], ( c2i - c1i) * weight);
    gpu_atomic_add(&out[6], ( c2r - c1r) * weight);
    gpu_atomic_add(&out[7], (-c1i - c2i) * weight);

    c1r = j00r*j10r + j00i*j10i;
    c1i = j00i*j10r - j00r*j10i;
    c2r = j11r*j01r + j11i*j01i;
    c2i = j11i*j01r - j11r*j01i;
    gpu_atomic_add(&out[8], (-c1r - c2r) * weight);
    gpu_atomic_add(&out[9], ( c2r - c1r) * weight);
    gpu_atomic_add(&out[12], ( c1i - c2i) * weight);
    gpu_atomic_add(&out[13], ( c2i + c1i) * weight);

    c1r = j00r*j11r + j00i*j11i;
    c1i = j00i*j11r - j00r*j11i;
    c2r = j01r*j10r + j01i*j10i;
    c2i = j01i*j10r - j01r*j10i;
    gpu_atomic_add(&out[10], ( c1r + c2r) * weight);
    gpu_atomic_add(&out[11], ( c1i - c2i) * weight);
    gpu_atomic_add(&out[14], (-c1i - c2i) * weight);
    gpu_atomic_add(&out[15], ( c1r - c2r) * weight);
}

__global__ void diffraction_grid_mueller_kernel(const GpuBeam *__restrict__ beams,
                                                const int *__restrict__ beamOffsets,
                                                const GpuReal *__restrict__ sinTheta,
                                                const GpuReal *__restrict__ cosTheta,
                                                const GpuReal *__restrict__ sinPhi,
                                                const GpuReal *__restrict__ cosPhi,
                                                const GpuReal *__restrict__ vf,
                                                const GpuReal *__restrict__ weights,
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
                                                GpuReal *__restrict__ mFull,
                                                GpuReal *__restrict__ mNoShadow)
{
    const int gridCount = nAz * (nZen + 1);
    const long long total = (long long)nOrient * gridCount;
    const long long idx = (long long)blockIdx.x * blockDim.x + threadIdx.x;
    if (idx >= total) return;

    const int orient = (int)(idx / gridCount);
    const int grid = (int)(idx - (long long)orient * gridCount);
    const int begin = beamOffsets[orient];
    const int end = beamOffsets[orient + 1];
    const int p = grid / (nZen + 1);
    const int t = grid - p * (nZen + 1);
    const GpuReal cp = cosPhi[p], sp = sinPhi[p];
    const GpuReal sin_t = sinTheta[t], cos_t = cosTheta[t];
    const GpuReal dx = sin_t * cp;
    const GpuReal dy = sin_t * sp;
    const GpuReal dz = -cos_t;
    const GpuReal vfx = vf[(grid * 3) + 0];
    const GpuReal vfy = vf[(grid * 3) + 1];
    const GpuReal vfz = vf[(grid * 3) + 2];

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
        if (!compute_beam_jones_context_gpu(
                b, cp, sp, sin_t, cos_t, dx, dy, dz, vfx, vfy, vfz,
                waveIndex, wi2, eps1, eps2, complWaveR, complWaveI,
                invComplWaveR, invComplWaveI, legacySign, 1,
                d00r, d00i, d01r, d01i, d10r, d10i, d11r, d11i))
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

    const GpuReal weight = weights[orient] * (GpuReal)0.25;
    mueller_accum_from_jones(j00r, j00i, j01r, j01i, j10r, j10i, j11r, j11i,
                             weight, &mFull[grid * 16]);
    if (mNoShadow)
        mueller_accum_from_jones(n00r, n00i, n01r, n01i, n10r, n10i, n11r, n11i,
                                 weight, &mNoShadow[grid * 16]);
}

__global__ void diffraction_grid_mueller_full_kernel(const GpuBeam *__restrict__ beams,
                                                     const int *__restrict__ beamOffsets,
                                                     const GpuReal *__restrict__ sinTheta,
                                                     const GpuReal *__restrict__ cosTheta,
                                                     const GpuReal *__restrict__ sinPhi,
                                                     const GpuReal *__restrict__ cosPhi,
                                                     const GpuReal *__restrict__ vf,
                                                     const GpuReal *__restrict__ weights,
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
                                                     GpuReal *__restrict__ mFull)
{
    const int gridCount = nAz * (nZen + 1);
    const long long total = (long long)nOrient * gridCount;
    const long long idx = (long long)blockIdx.x * blockDim.x + threadIdx.x;
    if (idx >= total) return;

    const int orient = (int)(idx / gridCount);
    const int grid = (int)(idx - (long long)orient * gridCount);
    const int begin = beamOffsets[orient];
    const int end = beamOffsets[orient + 1];
    const int p = grid / (nZen + 1);
    const int t = grid - p * (nZen + 1);
    const GpuReal cp = cosPhi[p], sp = sinPhi[p];
    const GpuReal sin_t = sinTheta[t], cos_t = cosTheta[t];
    const GpuReal dx = sin_t * cp;
    const GpuReal dy = sin_t * sp;
    const GpuReal dz = -cos_t;
    const GpuReal vfx = vf[(grid * 3) + 0];
    const GpuReal vfy = vf[(grid * 3) + 1];
    const GpuReal vfz = vf[(grid * 3) + 2];

    GpuReal j00r = 0.0, j00i = 0.0;
    GpuReal j01r = 0.0, j01i = 0.0;
    GpuReal j10r = 0.0, j10i = 0.0;
    GpuReal j11r = 0.0, j11i = 0.0;

    for (int bi = begin; bi < end; ++bi)
    {
        GpuReal d00r, d00i, d01r, d01i, d10r, d10i, d11r, d11i;
        const GpuBeam &b = beams[bi];
        if (!compute_beam_jones_context_gpu(
                b, cp, sp, sin_t, cos_t, dx, dy, dz, vfx, vfy, vfz,
                waveIndex, wi2, eps1, eps2, complWaveR, complWaveI,
                invComplWaveR, invComplWaveI, legacySign, 1,
                d00r, d00i, d01r, d01i, d10r, d10i, d11r, d11i))
            continue;

        j00r += d00r; j00i += d00i;
        j01r += d01r; j01i += d01i;
        j10r += d10r; j10i += d10i;
        j11r += d11r; j11i += d11i;
    }

    mueller_accum_from_jones(j00r, j00i, j01r, j01i, j10r, j10i, j11r, j11i,
                             weights[orient] * (GpuReal)0.25,
                             &mFull[grid * 16]);
}

__global__ void diffraction_grid_mueller_multik_kernel(const GpuBeam *__restrict__ beams,
                                                       const int *__restrict__ beamOffsets,
                                                       const GpuReal *__restrict__ sinTheta,
                                                       const GpuReal *__restrict__ cosTheta,
                                                       const GpuReal *__restrict__ sinPhi,
                                                       const GpuReal *__restrict__ cosPhi,
                                                       const GpuReal *__restrict__ vf,
                                                       const GpuReal *__restrict__ weights,
                                                       const GpuReal *__restrict__ scales,
                                                       const double *__restrict__ absPaths,
                                                       int nAz, int nZen,
                                                       int nOrient,
                                                       int nSizes,
                                                       GpuReal waveIndex,
                                                       GpuReal wi2,
                                                       GpuReal eps1,
                                                       GpuReal eps2,
                                                       GpuReal complWaveR,
                                                       GpuReal complWaveI,
                                                       GpuReal invComplWaveR,
                                                       GpuReal invComplWaveI,
                                                       int legacySign,
                                                       GpuReal cAbs,
                                                       GpuReal *__restrict__ mFull)
{
    const int gridCount = nAz * (nZen + 1);
    const long long perSize = (long long)nOrient * gridCount;
    const long long total = (long long)nSizes * perSize;
    const long long idx = (long long)blockIdx.x * blockDim.x + threadIdx.x;
    if (idx >= total) return;

    const int sizeIdx = (int)(idx / perSize);
    const long long localIdx = idx - (long long)sizeIdx * perSize;
    const int orient = (int)(localIdx / gridCount);
    const int grid = (int)(localIdx - (long long)orient * gridCount);
    const int begin = beamOffsets[orient];
    const int end = beamOffsets[orient + 1];
    const int p = grid / (nZen + 1);
    const int t = grid - p * (nZen + 1);
    const GpuReal cp = cosPhi[p], sp = sinPhi[p];
    const GpuReal sin_t = sinTheta[t], cos_t = cosTheta[t];
    const GpuReal dx = sin_t * cp;
    const GpuReal dy = sin_t * sp;
    const GpuReal dz = -cos_t;
    const GpuReal vfx = vf[(grid * 3) + 0];
    const GpuReal vfy = vf[(grid * 3) + 1];
    const GpuReal vfz = vf[(grid * 3) + 2];
    const GpuReal scale = scales[sizeIdx];

    GpuReal j00r = 0.0, j00i = 0.0;
    GpuReal j01r = 0.0, j01i = 0.0;
    GpuReal j10r = 0.0, j10i = 0.0;
    GpuReal j11r = 0.0, j11i = 0.0;

    for (int bi = begin; bi < end; ++bi)
    {
        GpuReal d00r, d00i, d01r, d01i, d10r, d10i, d11r, d11i;
        const GpuBeam &b = beams[bi];
        if (!compute_beam_jones_context_gpu_multik(
                b, cp, sp, sin_t, cos_t, dx, dy, dz, vfx, vfy, vfz,
                waveIndex, wi2, eps1, eps2, complWaveR, complWaveI,
                invComplWaveR, invComplWaveI, legacySign, scale,
                absPaths, cAbs,
                d00r, d00i, d01r, d01i, d10r, d10i, d11r, d11i))
            continue;

        j00r += d00r; j00i += d00i;
        j01r += d01r; j01i += d01i;
        j10r += d10r; j10i += d10i;
        j11r += d11r; j11i += d11i;
    }

    mueller_accum_from_jones(j00r, j00i, j01r, j01i, j10r, j10i, j11r, j11i,
                             weights[orient] * (GpuReal)0.25,
                             &mFull[((size_t)sizeIdx * gridCount + grid) * 16]);
}

__global__ void diffraction_grid_mueller_full8_kernel(const GpuBeam8 *__restrict__ beams,
                                                      const int *__restrict__ beamOffsets,
                                                      const GpuReal *__restrict__ sinTheta,
                                                      const GpuReal *__restrict__ cosTheta,
                                                      const GpuReal *__restrict__ sinPhi,
                                                      const GpuReal *__restrict__ cosPhi,
                                                      const GpuReal *__restrict__ vf,
                                                      const GpuReal *__restrict__ weights,
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
                                                      GpuReal *__restrict__ mFull)
{
    const int gridCount = nAz * (nZen + 1);
    const long long total = (long long)nOrient * gridCount;
    const long long idx = (long long)blockIdx.x * blockDim.x + threadIdx.x;
    if (idx >= total) return;

    const int orient = (int)(idx / gridCount);
    const int grid = (int)(idx - (long long)orient * gridCount);
    const int begin = beamOffsets[orient];
    const int end = beamOffsets[orient + 1];
    const int p = grid / (nZen + 1);
    const int t = grid - p * (nZen + 1);
    const GpuReal cp = cosPhi[p], sp = sinPhi[p];
    const GpuReal sin_t = sinTheta[t], cos_t = cosTheta[t];
    const GpuReal dx = sin_t * cp;
    const GpuReal dy = sin_t * sp;
    const GpuReal dz = -cos_t;
    const GpuReal vfx = vf[(grid * 3) + 0];
    const GpuReal vfy = vf[(grid * 3) + 1];
    const GpuReal vfz = vf[(grid * 3) + 2];

    GpuReal j00r = 0.0, j00i = 0.0;
    GpuReal j01r = 0.0, j01i = 0.0;
    GpuReal j10r = 0.0, j10i = 0.0;
    GpuReal j11r = 0.0, j11i = 0.0;

    for (int bi = begin; bi < end; ++bi)
    {
        GpuReal d00r, d00i, d01r, d01i, d10r, d10i, d11r, d11i;
        const GpuBeam8 &b = beams[bi];
        if (!compute_beam_jones_context_gpu8_auto(
                b, cp, sp, sin_t, cos_t, dx, dy, dz, vfx, vfy, vfz,
                waveIndex, wi2, eps1, eps2, complWaveR, complWaveI,
                invComplWaveR, invComplWaveI, legacySign, 1,
                d00r, d00i, d01r, d01i, d10r, d10i, d11r, d11i))
            continue;

        j00r += d00r; j00i += d00i;
        j01r += d01r; j01i += d01i;
        j10r += d10r; j10i += d10i;
        j11r += d11r; j11i += d11i;
    }

    mueller_accum_from_jones(j00r, j00i, j01r, j01i, j10r, j10i, j11r, j11i,
                             weights[orient] * (GpuReal)0.25,
                             &mFull[grid * 16]);
}

__global__ void diffraction_grid_mueller_mixed8_kernel(const GpuBeam8 *__restrict__ beams8,
                                                       const int *__restrict__ beamOffsets8,
                                                       const GpuBeam *__restrict__ beams,
                                                       const int *__restrict__ beamOffsets,
                                                       const GpuReal *__restrict__ sinTheta,
                                                       const GpuReal *__restrict__ cosTheta,
                                                       const GpuReal *__restrict__ sinPhi,
                                                       const GpuReal *__restrict__ cosPhi,
                                                       const GpuReal *__restrict__ vf,
                                                       const GpuReal *__restrict__ weights,
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
                                                       GpuReal *__restrict__ mFull)
{
    const int gridCount = nAz * (nZen + 1);
    const long long total = (long long)nOrient * gridCount;
    const long long idx = (long long)blockIdx.x * blockDim.x + threadIdx.x;
    if (idx >= total) return;

    const int orient = (int)(idx / gridCount);
    const int grid = (int)(idx - (long long)orient * gridCount);
    const int p = grid / (nZen + 1);
    const int t = grid - p * (nZen + 1);
    const GpuReal cp = cosPhi[p], sp = sinPhi[p];
    const GpuReal sin_t = sinTheta[t], cos_t = cosTheta[t];
    const GpuReal dx = sin_t * cp;
    const GpuReal dy = sin_t * sp;
    const GpuReal dz = -cos_t;
    const GpuReal vfx = vf[(grid * 3) + 0];
    const GpuReal vfy = vf[(grid * 3) + 1];
    const GpuReal vfz = vf[(grid * 3) + 2];

    GpuReal j00r = 0.0, j00i = 0.0;
    GpuReal j01r = 0.0, j01i = 0.0;
    GpuReal j10r = 0.0, j10i = 0.0;
    GpuReal j11r = 0.0, j11i = 0.0;

    for (int bi = beamOffsets8[orient]; bi < beamOffsets8[orient + 1]; ++bi)
    {
        GpuReal d00r, d00i, d01r, d01i, d10r, d10i, d11r, d11i;
        const GpuBeam8 &b = beams8[bi];
        if (!compute_beam_jones_context_gpu8_auto(
                b, cp, sp, sin_t, cos_t, dx, dy, dz, vfx, vfy, vfz,
                waveIndex, wi2, eps1, eps2, complWaveR, complWaveI,
                invComplWaveR, invComplWaveI, legacySign, 1,
                d00r, d00i, d01r, d01i, d10r, d10i, d11r, d11i))
            continue;

        j00r += d00r; j00i += d00i;
        j01r += d01r; j01i += d01i;
        j10r += d10r; j10i += d10i;
        j11r += d11r; j11i += d11i;
    }

    for (int bi = beamOffsets[orient]; bi < beamOffsets[orient + 1]; ++bi)
    {
        GpuReal d00r, d00i, d01r, d01i, d10r, d10i, d11r, d11i;
        const GpuBeam &b = beams[bi];
        if (!compute_beam_jones_context_gpu(
                b, cp, sp, sin_t, cos_t, dx, dy, dz, vfx, vfy, vfz,
                waveIndex, wi2, eps1, eps2, complWaveR, complWaveI,
                invComplWaveR, invComplWaveI, legacySign, 1,
                d00r, d00i, d01r, d01i, d10r, d10i, d11r, d11i))
            continue;

        j00r += d00r; j00i += d00i;
        j01r += d01r; j01i += d01i;
        j10r += d10r; j10i += d10i;
        j11r += d11r; j11i += d11i;
    }

    mueller_accum_from_jones(j00r, j00i, j01r, j01i, j10r, j10i, j11r, j11i,
                             weights[orient] * (GpuReal)0.25,
                             &mFull[grid * 16]);
}

__global__ void diffraction_grid_mueller_orient_kernel(const GpuBeam *__restrict__ beams,
                                                       const int *__restrict__ beamOffsets,
                                                       const GpuReal *__restrict__ sinTheta,
                                                       const GpuReal *__restrict__ cosTheta,
                                                       const GpuReal *__restrict__ sinPhi,
                                                       const GpuReal *__restrict__ cosPhi,
                                                       const GpuReal *__restrict__ vf,
                                                       const GpuReal *__restrict__ weights,
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
                                                       GpuReal *__restrict__ mOrient)
{
    const int gridCount = nAz * (nZen + 1);
    const long long total = (long long)nOrient * gridCount;
    const long long idx = (long long)blockIdx.x * blockDim.x + threadIdx.x;
    if (idx >= total) return;

    const int orient = (int)(idx / gridCount);
    const int grid = (int)(idx - (long long)orient * gridCount);
    const int begin = beamOffsets[orient];
    const int end = beamOffsets[orient + 1];
    const int p = grid / (nZen + 1);
    const int t = grid - p * (nZen + 1);
    const GpuReal cp = cosPhi[p], sp = sinPhi[p];
    const GpuReal sin_t = sinTheta[t], cos_t = cosTheta[t];
    const GpuReal dx = sin_t * cp;
    const GpuReal dy = sin_t * sp;
    const GpuReal dz = -cos_t;
    const GpuReal vfx = vf[(grid * 3) + 0];
    const GpuReal vfy = vf[(grid * 3) + 1];
    const GpuReal vfz = vf[(grid * 3) + 2];

    GpuReal j00r = 0.0, j00i = 0.0;
    GpuReal j01r = 0.0, j01i = 0.0;
    GpuReal j10r = 0.0, j10i = 0.0;
    GpuReal j11r = 0.0, j11i = 0.0;

    for (int bi = begin; bi < end; ++bi)
    {
        GpuReal d00r, d00i, d01r, d01i, d10r, d10i, d11r, d11i;
        const GpuBeam &b = beams[bi];
        if (!compute_beam_jones_context_gpu(
                b, cp, sp, sin_t, cos_t, dx, dy, dz, vfx, vfy, vfz,
                waveIndex, wi2, eps1, eps2, complWaveR, complWaveI,
                invComplWaveR, invComplWaveI, legacySign, 1,
                d00r, d00i, d01r, d01i, d10r, d10i, d11r, d11i))
            continue;

        j00r += d00r; j00i += d00i;
        j01r += d01r; j01i += d01i;
        j10r += d10r; j10i += d10i;
        j11r += d11r; j11i += d11i;
    }

    mueller_store_from_jones(j00r, j00i, j01r, j01i, j10r, j10i, j11r, j11i,
                             weights[orient] * (GpuReal)0.25,
                             &mOrient[idx * 16]);
}

__global__ void reduce_mueller_orient_kernel(const GpuReal *__restrict__ mOrient,
                                             int nOrient,
                                             int gridCount,
                                             GpuReal *__restrict__ mFull)
{
    const int grid = (int)((long long)blockIdx.x * blockDim.x + threadIdx.x);
    if (grid >= gridCount) return;

    GpuReal sum[16];
#pragma unroll
    for (int k = 0; k < 16; ++k)
        sum[k] = 0.0;

    for (int orient = 0; orient < nOrient; ++orient)
    {
        const GpuReal *src = &mOrient[((long long)orient * gridCount + grid) * 16];
#pragma unroll
        for (int k = 0; k < 16; ++k)
            sum[k] += src[k];
    }

    GpuReal *dst = &mFull[grid * 16];
#pragma unroll
    for (int k = 0; k < 16; ++k)
        dst[k] = sum[k];
}

__global__ void diffraction_grid_mueller_full_nocache_kernel(const GpuBeam *__restrict__ beams,
                                                            const int *__restrict__ beamOffsets,
                                                            const GpuReal *__restrict__ sinTheta,
                                                            const GpuReal *__restrict__ cosTheta,
                                                            const GpuReal *__restrict__ sinPhi,
                                                            const GpuReal *__restrict__ cosPhi,
                                                            const GpuReal *__restrict__ vf,
                                                            const GpuReal *__restrict__ weights,
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
                                                            GpuReal *__restrict__ mFull)
{
    const int gridCount = nAz * (nZen + 1);
    const long long total = (long long)nOrient * gridCount;
    const long long idx = (long long)blockIdx.x * blockDim.x + threadIdx.x;
    if (idx >= total) return;

    const int orient = (int)(idx / gridCount);
    const int grid = (int)(idx - (long long)orient * gridCount);
    const int begin = beamOffsets[orient];
    const int end = beamOffsets[orient + 1];

    GpuReal j00r = 0.0, j00i = 0.0;
    GpuReal j01r = 0.0, j01i = 0.0;
    GpuReal j10r = 0.0, j10i = 0.0;
    GpuReal j11r = 0.0, j11i = 0.0;

    for (int bi = begin; bi < end; ++bi)
    {
        GpuReal d00r, d00i, d01r, d01i, d10r, d10i, d11r, d11i;
        const GpuBeam &b = beams[bi];
        if (!compute_beam_jones_nocache_gpu(b, grid, nZen, sinTheta, cosTheta,
                                            sinPhi, cosPhi, vf, waveIndex, wi2,
                                            eps1, eps2, complWaveR, complWaveI,
                                            invComplWaveR, invComplWaveI, legacySign,
                                            d00r, d00i, d01r, d01i,
                                            d10r, d10i, d11r, d11i))
            continue;

        j00r += d00r; j00i += d00i;
        j01r += d01r; j01i += d01i;
        j10r += d10r; j10i += d10i;
        j11r += d11r; j11i += d11i;
    }

    mueller_accum_from_jones(j00r, j00i, j01r, j01i, j10r, j10i, j11r, j11i,
                             weights[orient] * (GpuReal)0.25,
                             &mFull[grid * 16]);
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
    auto failSingle = [](const char *where) {
        if (gpu_fft_debug_enabled())
            std::fprintf(stderr, "GPU single-orientation failed at %s\n", where);
        return false;
    };
    if (!isCoh)
        return failSingle("not-coherent");

    const int nBeams = (int)prepared.beams.size();
    const int nAz = m_sphere.nAzimuth;
    const int nZen = m_sphere.nZenith;
    const int gridCount = nAz * (nZen + 1);
    const bool computeNoShadow = ComputeNoShadow();
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
        b.nEdgeX = 0;
        b.nEdgeY = 0;
        b.isExternal = pb.isExternal ? 1 : 0;
        b.orientation = 0;
        for (int e = 0; e < 32; ++e)
        {
            b.x[e] = pb.edgeData.x[e];
            b.y[e] = pb.edgeData.y[e];
            b.slope_yx[e] = pb.edgeData.slope_yx[e];
            b.slope_xy[e] = pb.edgeData.slope_xy[e];
            const bool validX = e < b.nVertices && pb.edgeData.edge_valid_x[e];
            const bool validY = e < b.nVertices && pb.edgeData.edge_valid_y[e];
            if (validX)
                b.edge_valid_x[b.nEdgeX++] = (unsigned char)e;
            if (validY)
                b.edge_valid_y[b.nEdgeY++] = (unsigned char)e;
        }
        b.bdx = (GpuReal)(pb.bdx * pb.horAx + pb.bdy * pb.horAy + pb.bdz * pb.horAz);
        b.bdy = (GpuReal)(pb.bdx * pb.verAx + pb.bdy * pb.verAy + pb.bdz * pb.verAz);
        b.bdz = 0.0;
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

    const double gridSignature = gpu_grid_signature(m_sphere, nZen);
    if (ws.gridNAz != nAz || ws.gridNZen != nZen || fabs(ws.gridSignature - gridSignature) > 1e-12)
    {
        std::vector<GpuReal> hSinTheta(nZen + 1), hCosTheta(nZen + 1);
        for (int t = 0; t <= nZen; ++t)
        {
            GpuReal theta = (GpuReal)gpu_theta(m_sphere, t);
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

        if (cudaMemcpy(ws.sinTheta, hSinTheta.data(), hSinTheta.size() * sizeof(GpuReal), cudaMemcpyHostToDevice) != cudaSuccess) return failSingle("copy-sintheta");
        if (cudaMemcpy(ws.cosTheta, hCosTheta.data(), hCosTheta.size() * sizeof(GpuReal), cudaMemcpyHostToDevice) != cudaSuccess) return failSingle("copy-costheta");
        if (cudaMemcpy(ws.sinPhi, hSinPhi.data(), hSinPhi.size() * sizeof(GpuReal), cudaMemcpyHostToDevice) != cudaSuccess) return failSingle("copy-sinphi");
        if (cudaMemcpy(ws.cosPhi, hCosPhi.data(), hCosPhi.size() * sizeof(GpuReal), cudaMemcpyHostToDevice) != cudaSuccess) return failSingle("copy-cosphi");
        if (cudaMemcpy(ws.vf, hVf.data(), hVf.size() * sizeof(GpuReal), cudaMemcpyHostToDevice) != cudaSuccess) return failSingle("copy-vf");
        ws.gridNAz = nAz;
        ws.gridNZen = nZen;
        ws.gridSignature = gridSignature;
    }

    cudaError_t errSingle = cudaMemcpy(ws.beams, hBeams.data(), hBeams.size() * sizeof(GpuBeam), cudaMemcpyHostToDevice);
    if (!gpu_report_cuda_error(errSingle, "single copy-beams")) return failSingle("copy-beams");
    errSingle = cudaMemset(ws.j, 0, hJ.size() * sizeof(GpuReal));
    if (!gpu_report_cuda_error(errSingle, "single memset-j")) return failSingle("memset-j");
    errSingle = cudaMemset(ws.jNoShadow, 0, hJns.size() * sizeof(GpuReal));
    if (!gpu_report_cuda_error(errSingle, "single memset-j-ns")) return failSingle("memset-j-ns");

    long long total = (long long)nBeams * gridCount;
    int block = gpu_block_size();
    int grid = (int)((total + block - 1) / block);
    diffraction_kernel<<<grid, block>>>(
        ws.beams, nBeams, ws.sinTheta, ws.cosTheta, ws.sinPhi, ws.cosPhi, ws.vf,
        nAz, nZen, m_waveIndex, m_wi2, gpu_effective_eps1(m_eps1), m_eps2,
        real(m_complWave), imag(m_complWave),
        real(m_invComplWave), imag(m_invComplWave),
        m_legacySign ? 1 : 0, ws.j, ws.jNoShadow);

    errSingle = cudaGetLastError();
    if (!gpu_report_cuda_error(errSingle, "single diffraction kernel launch"))
        return failSingle("kernel-launch");
    errSingle = cudaDeviceSynchronize();
    if (!gpu_report_cuda_error(errSingle, "single diffraction kernel sync"))
        return failSingle("kernel-sync");

    errSingle = cudaMemcpy(hJ.data(), ws.j, hJ.size() * sizeof(GpuReal), cudaMemcpyDeviceToHost);
    if (!gpu_report_cuda_error(errSingle, "single copy-j-host")) return failSingle("copy-j-host");
    if (computeNoShadow)
    {
        errSingle = cudaMemcpy(hJns.data(), ws.jNoShadow, hJns.size() * sizeof(GpuReal), cudaMemcpyDeviceToHost);
        if (!gpu_report_cuda_error(errSingle, "single copy-j-ns-host")) return failSingle("copy-j-ns-host");
    }

    add_mueller_from_jones(hJ, prepared.sinZenith, nAz, nZen, localM);
    if (computeNoShadow)
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

    const bool computeNoShadow = ComputeNoShadow();
    const int noAtomicsMode = gpu_no_atomics_mode();
    const bool definitelyNoAtomics = (noAtomicsMode > 0)
        || (noAtomicsMode < 0 && !computeNoShadow);
    const int fusedMuellerMode = gpu_fused_mueller_mode();
    const bool fusedMueller = definitelyNoAtomics && ((fusedMuellerMode >= 0)
        ? (fusedMuellerMode != 0)
        : true);
    GpuWorkspace &ws = g_gpuWorkspace;

    size_t usable = (size_t)(freeBytes * gpu_memory_fraction());
    const size_t safetyBytes = 256ULL << 20;
    if (usable <= safetyBytes)
        return 1;
    usable -= safetyBytes;

    auto grow_bytes = [](size_t required, size_t capacity, size_t elemSize) -> size_t {
        return required > capacity ? (required - capacity) * elemSize : 0;
    };

    const size_t fixedGrowth =
        grow_bytes((size_t)nZen + 1, ws.sinThetaCap, sizeof(GpuReal)) +
        grow_bytes((size_t)nZen + 1, ws.cosThetaCap, sizeof(GpuReal)) +
        grow_bytes((size_t)nAz, ws.sinPhiCap, sizeof(GpuReal)) +
        grow_bytes((size_t)nAz, ws.cosPhiCap, sizeof(GpuReal)) +
        grow_bytes(gridCount * 3, ws.vfCap, sizeof(GpuReal)) +
        grow_bytes(gridCount * 16, ws.mCap, sizeof(GpuReal)) +
        (computeNoShadow
            ? grow_bytes(gridCount * 16, ws.mNoShadowCap, sizeof(GpuReal))
            : 0);
    if (usable <= fixedGrowth)
        return 1;
    usable -= fixedGrowth;

    int count = 0;
    size_t beamCount = 0;
    for (int i = start; i < (int)prepared.size() && count < maxCount; ++i)
    {
        const int nextCount = count + 1;
        beamCount += prepared[i].beams.size();
        size_t used =
            grow_bytes(beamCount, ws.beamCap, sizeof(GpuBeam)) +
            grow_bytes((size_t)nextCount, ws.weightsCap, sizeof(GpuReal)) +
            grow_bytes((size_t)nextCount + 1, ws.beamOffsetsCap, sizeof(int));
        if (!fusedMueller)
        {
            used += grow_bytes((size_t)nextCount * gridCount * 8,
                               ws.jCap, sizeof(GpuReal));
        }
        if (!fusedMueller && computeNoShadow)
            used += grow_bytes((size_t)nextCount * gridCount * 8,
                               ws.jNoShadowCap, sizeof(GpuReal));

        if (count > 0 && used > usable)
            break;
        ++count;
    }

    int selected = std::max(1, count);
    if (!g_gpuMultiWorker)
    {
        const int nDevices = gpu_multi_device_count(maxCount);
        if (nDevices > 1)
            selected = std::min(maxCount, selected * nDevices);
    }
    return selected;
}

static void add_arr2d_inplace(const Arr2D &src, int nAz, int nZen, Arr2D &dst)
{
    for (int p = 0; p < nAz; ++p)
        for (int t = 0; t <= nZen; ++t)
            for (int r = 0; r < 4; ++r)
                for (int c = 0; c < 4; ++c)
                    dst(p, t, r, c) += src(p, t, r, c);
}

static void copy_arr2d_theta_row(const Arr2D &src, int srcT,
                                 int nAz, Arr2D &dst, int dstT)
{
    for (int p = 0; p < nAz; ++p)
    {
        const double *s = src.RawCell(p, srcT);
        double *d = dst.RawCell(p, dstT);
        for (int k = 0; k < 16; ++k)
            d[k] = s[k];
    }
}

static void azimuth_average_output_cell(const Arr2D &src, int nAz, int nZen,
                                        int t, double out[16])
{
    for (int k = 0; k < 16; ++k)
        out[k] = 0.0;

    const bool isPole = (t == 0 || t == nZen);
    for (int p = 0; p < nAz; ++p)
    {
        const double *m = src.RawCell(p, t);
        if (isPole)
        {
            for (int k = 0; k < 16; ++k)
                out[k] += m[k];
            continue;
        }

        const double radAz = -p * (M_2PI / (double)nAz);
        const double cs = cos(2.0 * radAz);
        const double sn = sin(2.0 * radAz);
        for (int r = 0; r < 4; ++r)
        {
            const int base = r * 4;
            out[base + 0] += m[base + 0];
            out[base + 1] += m[base + 1] * cs - m[base + 2] * sn;
            out[base + 2] += m[base + 1] * sn + m[base + 2] * cs;
            out[base + 3] += m[base + 3];
        }
    }

    const double inv = 1.0 / (double)nAz;
    for (int k = 0; k < 16; ++k)
        out[k] *= inv;
}

static std::vector<int> detect_fft_phi_spike_rows(const Arr2D &m,
                                                  int nAz,
                                                  int nZen,
                                                  double threshold,
                                                  int maxRows)
{
    std::vector<int> rows;
    if (maxRows <= 0 || nAz <= 0 || nZen < 2)
        return rows;

    const int elems[] = {
        1 * 4 + 2, 1 * 4 + 3,
        2 * 4 + 1, 2 * 4 + 3,
        3 * 4 + 1, 3 * 4 + 2,
        0 * 4 + 2, 0 * 4 + 3
    };
    const int nElems = (int)(sizeof(elems) / sizeof(elems[0]));

    double prev[16], cur[16], next[16];
    std::vector<std::pair<double, int>> scored;
    for (int t = 1; t < nZen; ++t)
    {
        azimuth_average_output_cell(m, nAz, nZen, t - 1, prev);
        azimuth_average_output_cell(m, nAz, nZen, t, cur);
        azimuth_average_output_cell(m, nAz, nZen, t + 1, next);
        const double c11 = std::max(fabs(cur[0]), 1e-30);
        const double p11 = std::max(fabs(prev[0]), 1e-30);
        const double n11 = std::max(fabs(next[0]), 1e-30);
        double score = 0.0;
        for (int i = 0; i < nElems; ++i)
        {
            const int e = elems[i];
            const double yc = cur[e] / c11;
            const double yn = 0.5 * (prev[e] / p11 + next[e] / n11);
            score = std::max(score, fabs(yc - yn));
        }
        if (score >= threshold)
            scored.push_back(std::make_pair(score, t));
    }

    std::sort(scored.begin(), scored.end(),
              [](const std::pair<double, int> &a,
                 const std::pair<double, int> &b) {
                  return a.first > b.first;
              });
    for (const auto &entry : scored)
    {
        rows.push_back(entry.second);
        if ((int)rows.size() >= maxRows)
            break;
    }
    std::sort(rows.begin(), rows.end());
    return rows;
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

static bool fft_upsample_phi_arr2d_impl(const Arr2D &low,
                                        const Arr2D *low2,
                                        int nLow,
                                        int nFull,
                                        int nZen,
                                        Arr2D &dst,
                                        Arr2D *dst2)
{
    if (nLow <= 0 || nFull <= 0 || nZen < 0)
        return false;

    const int nElem = 16;
    const int sideCount = (low2 && dst2) ? 2 : 1;
    const int sideStride = (nZen + 1) * nElem;
    const int batch = sideCount * sideStride;
    const size_t lowCount = (size_t)batch * nLow;
    const size_t fullCount = (size_t)batch * nFull;
    const bool timing = gpu_timing_enabled();
    const double tStart = timing ? gpu_now_ms() : 0.0;
    double tPack = 0.0, tEnsure = 0.0, tH2d = 0.0, tForward = 0.0;
    double tMemset = 0.0, tPad = 0.0, tInverse = 0.0, tD2h = 0.0, tAdd = 0.0;

    GpuWorkspace &ws = g_gpuWorkspace;
    ws.hFftLow.resize(lowCount);
    std::vector<GpuComplex> &hLow = ws.hFftLow;
    for (int side = 0; side < sideCount; ++side)
    {
        const Arr2D &src = (side == 0) ? low : *low2;
        for (int t = 0; t <= nZen; ++t)
            for (int r = 0; r < 4; ++r)
                for (int c = 0; c < 4; ++c)
                {
                    int b = side * sideStride + (t * nElem) + r * 4 + c;
                    for (int p = 0; p < nLow; ++p)
                    {
                        GpuComplex z;
                        z.x = (GpuReal)src(p, t, r, c);
                        z.y = (GpuReal)0.0;
                        hLow[(size_t)b * nLow + p] = z;
                    }
                }
    }
    if (timing) tPack = gpu_now_ms() - tStart;

    bool ok = false;

    do
    {
        double t0 = timing ? gpu_now_ms() : 0.0;
        if (!ensure_fft_workspace(ws, nLow, nFull, batch)) break;
        if (timing) tEnsure = gpu_now_ms() - t0;
        t0 = timing ? gpu_now_ms() : 0.0;
        if (cudaMemcpy(ws.fftLow, hLow.data(), lowCount * sizeof(GpuComplex),
                       cudaMemcpyHostToDevice) != cudaSuccess) break;
        if (timing) tH2d = gpu_now_ms() - t0;

        t0 = timing ? gpu_now_ms() : 0.0;
        if (gpu_cufft_exec(ws.fftPlanLow, ws.fftLow, ws.fftLow, CUFFT_FORWARD) != CUFFT_SUCCESS) break;
        if (timing && cudaDeviceSynchronize() != cudaSuccess) break;
        if (timing) tForward = gpu_now_ms() - t0;
        t0 = timing ? gpu_now_ms() : 0.0;
        if (cudaMemset(ws.fftFull, 0, fullCount * sizeof(GpuComplex)) != cudaSuccess) break;
        if (timing) tMemset = gpu_now_ms() - t0;

        t0 = timing ? gpu_now_ms() : 0.0;
        int block = gpu_block_size();
        int grid = (int)((lowCount + block - 1) / block);
        GpuReal scale = (GpuReal)(1.0 / (double)nLow);
        fft_phi_pad_kernel<<<grid, block>>>(ws.fftLow, ws.fftFull, nLow, nFull, batch, scale);
        if (cudaGetLastError() != cudaSuccess) break;
        if (timing && cudaDeviceSynchronize() != cudaSuccess) break;
        if (timing) tPad = gpu_now_ms() - t0;
        t0 = timing ? gpu_now_ms() : 0.0;
        if (gpu_cufft_exec(ws.fftPlanFull, ws.fftFull, ws.fftFull, CUFFT_INVERSE) != CUFFT_SUCCESS) break;
        if (timing && cudaDeviceSynchronize() != cudaSuccess) break;
        if (timing) tInverse = gpu_now_ms() - t0;

        t0 = timing ? gpu_now_ms() : 0.0;
        ws.hFftFull.resize(fullCount);
        std::vector<GpuComplex> &hFull = ws.hFftFull;
        if (cudaMemcpy(hFull.data(), ws.fftFull, fullCount * sizeof(GpuComplex),
                       cudaMemcpyDeviceToHost) != cudaSuccess) break;
        if (timing) tD2h = gpu_now_ms() - t0;

        t0 = timing ? gpu_now_ms() : 0.0;
        const int addTasks = sideCount * nFull;
    #pragma omp parallel for schedule(static) if(addTasks * (nZen + 1) >= 8192)
        for (int task = 0; task < addTasks; ++task)
        {
            const int side = task / nFull;
            const int p = task - side * nFull;
            Arr2D &out = (side == 0) ? dst : *dst2;
            for (int t = 0; t <= nZen; ++t)
            {
                double *cell = out.RawCell(p, t);
                const int base = side * sideStride + t * nElem;
                for (int k = 0; k < nElem; ++k)
                    cell[k] += hFull[(size_t)(base + k) * nFull + p].x;
            }
        }
        if (timing) tAdd = gpu_now_ms() - t0;

        ok = true;
    } while (false);

    if (timing)
        std::fprintf(stderr,
                     "GPU timing fft_phi nLow=%d nFull=%d batch=%d sides=%d pack=%.3fms ensure=%.3fms h2d=%.3fms fwd=%.3fms memset=%.3fms pad=%.3fms inv=%.3fms d2h=%.3fms add=%.3fms total=%.3fms ok=%d\n",
                     nLow, nFull, batch, sideCount, tPack, tEnsure, tH2d, tForward,
                     tMemset, tPad, tInverse, tD2h, tAdd,
                     gpu_now_ms() - tStart, ok ? 1 : 0);

    return ok;
}

static bool fft_upsample_phi_arr2d(const Arr2D &low,
                                   int nLow,
                                   int nFull,
                                   int nZen,
                                   Arr2D &dst)
{
    return fft_upsample_phi_arr2d_impl(low, nullptr, nLow, nFull, nZen,
                                      dst, nullptr);
}

static bool fft_upsample_phi_arr2d_pair(const Arr2D &low,
                                        const Arr2D &low2,
                                        int nLow,
                                        int nFull,
                                        int nZen,
                                        Arr2D &dst,
                                        Arr2D &dst2)
{
    return fft_upsample_phi_arr2d_impl(low, &low2, nLow, nFull, nZen,
                                      dst, &dst2);
}

static bool fft_upsample_theta_arr2d_impl(const Arr2D &low,
                                          const Arr2D *low2,
                                          int nAz,
                                          int nThetaLow,
                                          int nThetaFull,
                                          Arr2D &dst,
                                          Arr2D *dst2)
{
    if (nAz <= 0 || nThetaLow <= 0 || nThetaFull <= 0)
        return false;

    const int nElem = 16;
    const int sideCount = (low2 && dst2) ? 2 : 1;
    const int nThetaLowPts = nThetaLow + 1;
    const int nThetaFullPts = nThetaFull + 1;
    const int nLow = 2 * nThetaLow;
    const int nFull = 2 * nThetaFull;
    const int sideStride = nAz * nElem;
    const int batch = sideCount * sideStride;
    const size_t lowCount = (size_t)batch * nLow;
    const size_t fullCount = (size_t)batch * nFull;
    const bool timing = gpu_timing_enabled();
    const double tStart = timing ? gpu_now_ms() : 0.0;

    GpuWorkspace &ws = g_gpuWorkspace;
    ws.hFftLow.resize(lowCount);
    std::vector<GpuComplex> &hLow = ws.hFftLow;
    for (int side = 0; side < sideCount; ++side)
    {
        const Arr2D &src = (side == 0) ? low : *low2;
        for (int p = 0; p < nAz; ++p)
            for (int r = 0; r < 4; ++r)
                for (int c = 0; c < 4; ++c)
                {
                    int b = side * sideStride + p * nElem + r * 4 + c;
                    for (int t = 0; t < nLow; ++t)
                    {
                        int srcT = (t <= nThetaLow) ? t : (nLow - t);
                        GpuComplex z;
                        z.x = (GpuReal)src(p, srcT, r, c);
                        z.y = (GpuReal)0.0;
                        hLow[(size_t)b * nLow + t] = z;
                    }
                }
    }

    bool ok = false;
    do
    {
        if (!ensure_fft_workspace(ws, nLow, nFull, batch)) break;
        if (cudaMemcpy(ws.fftLow, hLow.data(), lowCount * sizeof(GpuComplex),
                       cudaMemcpyHostToDevice) != cudaSuccess) break;
        if (gpu_cufft_exec(ws.fftPlanLow, ws.fftLow, ws.fftLow, CUFFT_FORWARD) != CUFFT_SUCCESS) break;
        if (cudaMemset(ws.fftFull, 0, fullCount * sizeof(GpuComplex)) != cudaSuccess) break;

        int block = gpu_block_size();
        int grid = (int)((lowCount + block - 1) / block);
        GpuReal scale = (GpuReal)(1.0 / (double)nLow);
        fft_phi_pad_kernel<<<grid, block>>>(ws.fftLow, ws.fftFull, nLow, nFull, batch, scale);
        if (cudaGetLastError() != cudaSuccess) break;
        if (gpu_cufft_exec(ws.fftPlanFull, ws.fftFull, ws.fftFull, CUFFT_INVERSE) != CUFFT_SUCCESS) break;

        ws.hFftFull.resize(fullCount);
        std::vector<GpuComplex> &hFull = ws.hFftFull;
        if (cudaMemcpy(hFull.data(), ws.fftFull, fullCount * sizeof(GpuComplex),
                       cudaMemcpyDeviceToHost) != cudaSuccess) break;

        for (int side = 0; side < sideCount; ++side)
        {
            Arr2D &out = (side == 0) ? dst : *dst2;
            for (int p = 0; p < nAz; ++p)
                for (int t = 0; t < nThetaFullPts; ++t)
                {
                    double *cell = out.RawCell(p, t);
                    const int base = side * sideStride + p * nElem;
                    for (int k = 0; k < nElem; ++k)
                        cell[k] += hFull[(size_t)(base + k) * nFull + t].x;
                }
        }
        ok = true;
    } while (false);

    if (timing)
        std::fprintf(stderr,
                     "GPU timing fft_theta_even nLowPts=%d nFullPts=%d extLow=%d extFull=%d batch=%d sides=%d total=%.3fms ok=%d\n",
                     nThetaLowPts, nThetaFullPts, nLow, nFull, batch, sideCount,
                     gpu_now_ms() - tStart, ok ? 1 : 0);

    return ok;
}

static bool fft_upsample_theta_arr2d(const Arr2D &low,
                                     int nAz,
                                     int nThetaLow,
                                     int nThetaFull,
                                     Arr2D &dst)
{
    return fft_upsample_theta_arr2d_impl(low, nullptr, nAz, nThetaLow,
                                        nThetaFull, dst, nullptr);
}

static bool fft_upsample_theta_arr2d_pair(const Arr2D &low,
                                          const Arr2D &low2,
                                          int nAz,
                                          int nThetaLow,
                                          int nThetaFull,
                                          Arr2D &dst,
                                          Arr2D &dst2)
{
    return fft_upsample_theta_arr2d_impl(low, &low2, nAz, nThetaLow,
                                        nThetaFull, dst, &dst2);
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
                                                   Arr2D &localM_noshadow,
                                                   double scale,
                                                   double waveIndex)
{
    auto failFft = [](const char *where) {
        if (gpu_fft_debug_enabled())
            std::fprintf(stderr, "GPU FFT failed at %s\n", where);
        return false;
    };
    if (!isCoh)
        return failFft("not-coherent");

    const int nFull = m_sphere.nAzimuth;
    const int nZen = m_sphere.nZenith;
    const bool computeNoShadow = ComputeNoShadow();
    int factor = choose_fft_phi_factor(nFull, m_fftPhiFactor);
    if (factor <= 1 || nFull < 32)
        return HandleOrientationsToLocalGpu(prepared, start, count,
                                            localM, localM_noshadow,
                                            scale, waveIndex);

    int nLow = nFull / factor;
    if (nLow < 16) nLow = std::min(nFull, 16);
    if (nLow >= nFull)
        return HandleOrientationsToLocalGpu(prepared, start, count,
                                            localM, localM_noshadow,
                                            scale, waveIndex);

    int thetaFactor = gpu_fft_theta_factor();
    int nZenLow = nZen;
    if (thetaFactor > 1)
    {
        if (m_sphere.isNonUniform || nZen < 32 || (nZen % thetaFactor) != 0)
            thetaFactor = 1;
        else
            nZenLow = std::max(2, nZen / thetaFactor);
    }

    static bool printed = false;
    if (!printed)
    {
        std::fprintf(stderr,
                     "GPU FFT angular interpolation: direct Nphi=%d, output Nphi=%d (factor=%d), "
                     "direct Ntheta=%d, output Ntheta=%d (factor=%d). "
                     "This is angular Fourier interpolation, not aperture pFFT/FMM.\n",
                     nLow, nFull, factor, nZenLow + 1, nZen + 1, thetaFactor);
        printed = true;
    }

    ScatteringRange fullSphere = m_sphere;
    ScatteringRange lowSphere = fullSphere;
    lowSphere.nAzimuth = nLow;
    lowSphere.azinuthStep = M_2PI / nLow;
    lowSphere.nZenith = nZenLow;
    lowSphere.zenithStep = (lowSphere.zenithEnd - lowSphere.zenithStart) / nZenLow;
    lowSphere.ComputeSphereDirections(*m_incidentLight);

    Arr2D lowM(nLow + 1, nZenLow + 1, 4, 4); lowM.ClearArr();
    Arr2D lowMns(computeNoShadow ? nLow + 1 : 0,
                 computeNoShadow ? nZenLow + 1 : 0,
                 computeNoShadow ? 4 : 0,
                 computeNoShadow ? 4 : 0);
    if (computeNoShadow)
        lowMns.ClearArr();

    bool savedFft = m_fftEnabled;
    m_fftEnabled = false;
    m_sphere = lowSphere;
    bool ok = HandleOrientationsToLocalGpu(prepared, start, count, lowM, lowMns,
                                           scale, waveIndex);
    m_sphere = fullSphere;
    m_fftEnabled = savedFft;
    if (!ok)
        return failFft("low-direct");

    bool doCheck = gpu_fft_check_enabled() && factor > 2;
    if (!doCheck)
    {
        if (!computeNoShadow)
        {
            const Arr2D *phiSource = &lowM;
            Arr2D thetaM(nLow + 1, nZen + 1, 4, 4); thetaM.ClearArr();
            if (thetaFactor > 1)
            {
                if (!fft_upsample_theta_arr2d(lowM, nLow, nZenLow, nZen, thetaM))
                    return failFft("theta-upsample-full-only");
                phiSource = &thetaM;
            }
            if (!fft_upsample_phi_arr2d(*phiSource, nLow, nFull, nZen, localM))
                return failFft("phi-upsample-full-only");
        }
        else if (gpu_fft_pair_enabled())
        {
            const Arr2D *phiSource = &lowM;
            const Arr2D *phiSourceNs = &lowMns;
            Arr2D thetaM(nLow + 1, nZen + 1, 4, 4); thetaM.ClearArr();
            Arr2D thetaMns(nLow + 1, nZen + 1, 4, 4); thetaMns.ClearArr();
            if (thetaFactor > 1)
            {
                if (!fft_upsample_theta_arr2d_pair(lowM, lowMns, nLow, nZenLow,
                                                   nZen, thetaM, thetaMns))
                    return failFft("theta-upsample-pair");
                phiSource = &thetaM;
                phiSourceNs = &thetaMns;
            }
            if (!fft_upsample_phi_arr2d_pair(*phiSource, *phiSourceNs, nLow, nFull, nZen,
                                             localM, localM_noshadow))
                return failFft("phi-upsample-pair");
        }
        else
        {
            const Arr2D *phiSource = &lowM;
            const Arr2D *phiSourceNs = &lowMns;
            Arr2D thetaM(nLow + 1, nZen + 1, 4, 4); thetaM.ClearArr();
            Arr2D thetaMns(nLow + 1, nZen + 1, 4, 4); thetaMns.ClearArr();
            if (thetaFactor > 1)
            {
                if (!fft_upsample_theta_arr2d_pair(lowM, lowMns, nLow, nZenLow,
                                                   nZen, thetaM, thetaMns))
                    return failFft("theta-upsample-split");
                phiSource = &thetaM;
                phiSourceNs = &thetaMns;
            }
            if (!fft_upsample_phi_arr2d(*phiSource, nLow, nFull, nZen, localM))
                return failFft("phi-upsample-full-split");
            if (!fft_upsample_phi_arr2d(*phiSourceNs, nLow, nFull, nZen, localM_noshadow))
                return failFft("phi-upsample-noshadow-split");
        }
        return true;
    }

    Arr2D fftM(nFull + 1, nZen + 1, 4, 4); fftM.ClearArr();
    Arr2D fftMns(nFull + 1, nZen + 1, 4, 4); fftMns.ClearArr();
    if (!computeNoShadow)
    {
        const Arr2D *phiSource = &lowM;
        Arr2D thetaM(nLow + 1, nZen + 1, 4, 4); thetaM.ClearArr();
        if (thetaFactor > 1)
        {
            if (!fft_upsample_theta_arr2d(lowM, nLow, nZenLow, nZen, thetaM))
                return false;
            phiSource = &thetaM;
        }
        if (!fft_upsample_phi_arr2d(*phiSource, nLow, nFull, nZen, fftM))
            return false;
    }
    else if (gpu_fft_pair_enabled())
    {
        const Arr2D *phiSource = &lowM;
        const Arr2D *phiSourceNs = &lowMns;
        Arr2D thetaM(nLow + 1, nZen + 1, 4, 4); thetaM.ClearArr();
        Arr2D thetaMns(nLow + 1, nZen + 1, 4, 4); thetaMns.ClearArr();
        if (thetaFactor > 1)
        {
            if (!fft_upsample_theta_arr2d_pair(lowM, lowMns, nLow, nZenLow,
                                               nZen, thetaM, thetaMns))
                return false;
            phiSource = &thetaM;
            phiSourceNs = &thetaMns;
        }
        if (!fft_upsample_phi_arr2d_pair(*phiSource, *phiSourceNs, nLow, nFull, nZen,
                                         fftM, fftMns))
            return false;
    }
    else
    {
        const Arr2D *phiSource = &lowM;
        const Arr2D *phiSourceNs = &lowMns;
        Arr2D thetaM(nLow + 1, nZen + 1, 4, 4); thetaM.ClearArr();
        Arr2D thetaMns(nLow + 1, nZen + 1, 4, 4); thetaMns.ClearArr();
        if (thetaFactor > 1)
        {
            if (!fft_upsample_theta_arr2d_pair(lowM, lowMns, nLow, nZenLow,
                                               nZen, thetaM, thetaMns))
                return false;
            phiSource = &thetaM;
            phiSourceNs = &thetaMns;
        }
        if (!fft_upsample_phi_arr2d(*phiSource, nLow, nFull, nZen, fftM))
            return false;
        if (!fft_upsample_phi_arr2d(*phiSourceNs, nLow, nFull, nZen, fftMns))
            return false;
    }

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
                                                    checkLowM, checkLowMns,
                                                    scale, waveIndex);
        m_sphere = fullSphere;
        if (!checkOk)
            return false;

        Arr2D checkM(nFull + 1, nZen + 1, 4, 4); checkM.ClearArr();
        Arr2D checkMns(nFull + 1, nZen + 1, 4, 4); checkMns.ClearArr();
        if (!computeNoShadow)
        {
            if (!fft_upsample_phi_arr2d(checkLowM, nCheck, nFull, nZen, checkM))
                return false;
        }
        else if (gpu_fft_pair_enabled())
        {
            if (!fft_upsample_phi_arr2d_pair(checkLowM, checkLowMns, nCheck,
                                             nFull, nZen, checkM, checkMns))
                return false;
        }
        else
        {
            if (!fft_upsample_phi_arr2d(checkLowM, nCheck, nFull, nZen, checkM))
                return false;
            if (!fft_upsample_phi_arr2d(checkLowMns, nCheck, nFull, nZen, checkMns))
                return false;
        }
        report_fft_check(fftM, checkM, nFull, nZen, "full");
        if (computeNoShadow)
            report_fft_check(fftMns, checkMns, nFull, nZen, "no-shadow");
    }

    if (gpu_fft_adaptive_phi_enabled())
    {
        const double threshold = gpu_fft_adaptive_phi_threshold();
        const int maxRows = gpu_fft_adaptive_phi_max_rows();
        std::vector<int> rows = detect_fft_phi_spike_rows(fftM, nFull, nZen,
                                                          threshold, maxRows);
        if (!rows.empty())
        {
            ScatteringRange savedSphere = m_sphere;
            ScatteringRange rowSphere = fullSphere;
            rowSphere.isNonUniform = true;
            rowSphere.thetaValues.clear();
            rowSphere.thetaValues.reserve(rows.size());
            for (int t : rows)
                rowSphere.thetaValues.push_back(fullSphere.GetZenith(t));
            rowSphere.nAzimuth = nFull;
            rowSphere.azinuthStep = M_2PI / nFull;
            rowSphere.nZenith = (int)rowSphere.thetaValues.size() - 1;
            rowSphere.zenithStart = rowSphere.thetaValues.front();
            rowSphere.zenithEnd = rowSphere.thetaValues.back();
            rowSphere.zenithStep = rowSphere.nZenith > 0
                ? (rowSphere.zenithEnd - rowSphere.zenithStart) / rowSphere.nZenith
                : 0.0;
            rowSphere.ComputeSphereDirections(*m_incidentLight);

            Arr2D directM(nFull + 1, rows.size(), 4, 4); directM.ClearArr();
            Arr2D directMns(nFull + 1, rows.size(), 4, 4); directMns.ClearArr();
            m_sphere = rowSphere;
            bool refineOk = HandleOrientationsToLocalGpu(prepared, start, count,
                                                         directM, directMns,
                                                         scale, waveIndex);
            m_sphere = savedSphere;
            if (!refineOk)
                return false;

            for (size_t i = 0; i < rows.size(); ++i)
            {
                copy_arr2d_theta_row(directM, (int)i, nFull, fftM, rows[i]);
                if (computeNoShadow)
                    copy_arr2d_theta_row(directMns, (int)i, nFull, fftMns, rows[i]);
            }

            std::fprintf(stderr,
                         "GPU FFT adaptive phi: refined %zu/%d theta rows at full Nphi=%d (threshold=%.3g):",
                         rows.size(), nZen + 1, nFull, threshold);
            for (int t : rows)
                std::fprintf(stderr, " %.6gdeg", RadToDeg(fullSphere.GetZenith(t)));
            std::fprintf(stderr, "\n");
        }
    }

    add_arr2d_inplace(fftM, nFull, nZen, localM);
    if (computeNoShadow)
        add_arr2d_inplace(fftMns, nFull, nZen, localM_noshadow);

    return true;
}

bool HandlerPO::HandleOrientationsToLocalGpuFftPhiMultiK(
                                                   const std::vector<PreparedOrientation> &prepared,
                                                   int start,
                                                   int count,
                                                   const std::vector<double> &scales,
                                                   double waveIndex,
                                                   std::vector<Arr2D> &localMs)
{
    auto failFft = [](const char *where) {
        if (gpu_fft_debug_enabled())
            std::fprintf(stderr, "GPU FFT multik failed at %s\n", where);
        return false;
    };
    if (!gpu_multik_full_enabled())
        return failFft("disabled");
    if (!isCoh || ComputeNoShadow())
        return failFft("unsupported-mode");
    if (gpu_fft_check_enabled() || gpu_fft_adaptive_phi_enabled())
        return failFft("check-or-adaptive-enabled");

    const int nFull = m_sphere.nAzimuth;
    const int nZen = m_sphere.nZenith;
    int factor = choose_fft_phi_factor(nFull, m_fftPhiFactor);
    if (factor <= 1 || nFull < 32)
        return HandleOrientationsToLocalGpuMultiK(prepared, start, count,
                                                 scales, waveIndex, localMs);

    int nLow = nFull / factor;
    if (nLow < 16) nLow = std::min(nFull, 16);
    if (nLow >= nFull)
        return HandleOrientationsToLocalGpuMultiK(prepared, start, count,
                                                 scales, waveIndex, localMs);

    int thetaFactor = gpu_fft_theta_factor();
    int nZenLow = nZen;
    if (thetaFactor > 1)
    {
        if (m_sphere.isNonUniform || nZen < 32 || (nZen % thetaFactor) != 0)
            thetaFactor = 1;
        else
            nZenLow = std::max(2, nZen / thetaFactor);
    }

    static bool printed = false;
    if (!printed)
    {
        std::fprintf(stderr,
                     "GPU FFT multik angular interpolation: direct Nphi=%d, output Nphi=%d (factor=%d), "
                     "direct Ntheta=%d, output Ntheta=%d (factor=%d).\n",
                     nLow, nFull, factor, nZenLow + 1, nZen + 1,
                     thetaFactor);
        printed = true;
    }

    ScatteringRange fullSphere = m_sphere;
    ScatteringRange lowSphere = fullSphere;
    lowSphere.nAzimuth = nLow;
    lowSphere.azinuthStep = M_2PI / nLow;
    lowSphere.nZenith = nZenLow;
    lowSphere.zenithStep = (lowSphere.zenithEnd - lowSphere.zenithStart) / nZenLow;
    lowSphere.ComputeSphereDirections(*m_incidentLight);

    std::vector<Arr2D> lowMs;
    lowMs.reserve(scales.size());
    for (size_t s = 0; s < scales.size(); ++s)
    {
        lowMs.push_back(Arr2D(nLow + 1, nZenLow + 1, 4, 4));
        lowMs.back().ClearArr();
    }

    bool savedFft = m_fftEnabled;
    m_fftEnabled = false;
    m_sphere = lowSphere;
    bool ok = HandleOrientationsToLocalGpuMultiK(prepared, start, count,
                                                 scales, waveIndex, lowMs);
    m_sphere = fullSphere;
    m_fftEnabled = savedFft;
    if (!ok)
        return failFft("low-direct-multik");

    const bool timing = gpu_timing_enabled();
    const double tStart = timing ? gpu_now_ms() : 0.0;
    for (size_t s = 0; s < scales.size(); ++s)
    {
        const Arr2D *phiSource = &lowMs[s];
        Arr2D thetaM(nLow + 1, nZen + 1, 4, 4);
        if (thetaFactor > 1)
        {
            thetaM.ClearArr();
            if (!fft_upsample_theta_arr2d(lowMs[s], nLow, nZenLow, nZen,
                                          thetaM))
                return failFft("theta-upsample");
            phiSource = &thetaM;
        }
        if (!fft_upsample_phi_arr2d(*phiSource, nLow, nFull, nZen,
                                    localMs[s]))
            return failFft("phi-upsample");
    }
    if (timing)
        std::fprintf(stderr,
                     "GPU timing fft_multik sizes=%zu nLow=%d nFull=%d nZenLow=%d nZen=%d total=%.3fms\n",
                     scales.size(), nLow, nFull, nZenLow, nZen,
                     gpu_now_ms() - tStart);

    return true;
}

bool HandlerPO::HandleOrientationsToLocalGpuMultiK(const std::vector<PreparedOrientation> &prepared,
                                                   int start,
                                                   int count,
                                                   const std::vector<double> &scales,
                                                   double waveIndex,
                                                   std::vector<Arr2D> &localMs)
{
    auto failMultiK = [](const char *where) {
        if (gpu_fft_debug_enabled())
            std::fprintf(stderr, "GPU multik failed at %s\n", where);
        return false;
    };
    if (!gpu_multik_full_enabled())
        return failMultiK("disabled");
    if (!isCoh || ComputeNoShadow() || m_fftEnabled)
        return failMultiK("unsupported-mode");
    if (count <= 0 || scales.empty())
        return true;
    if (localMs.size() < scales.size())
        return failMultiK("local-output-size");

    const int nOrient = count;
    const int nSizes = (int)scales.size();
    const int nAz = m_sphere.nAzimuth;
    const int nZen = m_sphere.nZenith;
    const int gridCount = nAz * (nZen + 1);
    size_t nBeams = 0;
    for (int oi = 0; oi < nOrient; ++oi)
    {
        const PreparedOrientation &po = prepared[start + oi];
        for (const PreparedBeam &pb : po.beams)
        {
            if (!pb.edgeData.valid || pb.edgeData.nVertices <= 0 || pb.edgeData.nVertices > 32)
                return failMultiK("bad-edge-data");
            ++nBeams;
        }
    }
    if (nBeams == 0)
        return true;

    const bool timing = gpu_timing_enabled();
    const double tStart = timing ? gpu_now_ms() : 0.0;
    double tPack = 0.0, tEnsure = 0.0, tGrid = 0.0;
    double tCopy = 0.0, tKernel = 0.0, tD2h = 0.0, tAdd = 0.0;

    GpuWorkspace &ws = g_gpuWorkspace;
    double t0 = timing ? gpu_now_ms() : 0.0;
    ws.hBeams.resize(nBeams);
    ws.hWeights.assign(nOrient, 0.0);
    ws.hScales.resize(nSizes);
    ws.hBeamOffsets.assign(nOrient + 1, 0);
    size_t nAbsPaths = 0;
    for (int oi = 0; oi < nOrient; ++oi)
        for (const PreparedBeam &pb : prepared[start + oi].beams)
            nAbsPaths += pb.absorptionPaths.size();
    ws.hAbsPaths.clear();
    ws.hAbsPaths.reserve(nAbsPaths);
    for (int s = 0; s < nSizes; ++s)
        ws.hScales[s] = (GpuReal)scales[s];

    size_t bi = 0;
    for (int oi = 0; oi < nOrient; ++oi)
    {
        ws.hBeamOffsets[oi] = (int)bi;
        const PreparedOrientation &po = prepared[start + oi];
        ws.hWeights[oi] = po.sinZenith;
        for (const PreparedBeam &pb : po.beams)
        {
            GpuBeam &gb = ws.hBeams[bi++];
            pack_prepared_gpu_beam<GpuBeam, 32>(
                pb, oi, 1.0, 1.0, false, waveIndex,
                AbsorptionCoefficient(), gb);
            gb.absOffset = (int)ws.hAbsPaths.size();
            gb.absCount = (int)pb.absorptionPaths.size();
            for (double path : pb.absorptionPaths)
                ws.hAbsPaths.push_back(path);
        }
    }
    ws.hBeamOffsets[nOrient] = (int)bi;
    if (timing) tPack = gpu_now_ms() - t0;

    const size_t mCount = (size_t)nSizes * gridCount * 16;
    t0 = timing ? gpu_now_ms() : 0.0;
    if (!ensure_device_capacity(ws.beams, ws.beamCap, nBeams)) return failMultiK("alloc-beams");
    if (!ensure_device_capacity(ws.weights, ws.weightsCap, ws.hWeights.size())) return failMultiK("alloc-weights");
    if (!ensure_device_capacity(ws.scales, ws.scalesCap, ws.hScales.size())) return failMultiK("alloc-scales");
    if (!ensure_device_capacity(ws.absPaths, ws.absPathsCap, ws.hAbsPaths.size())) return failMultiK("alloc-abs-paths");
    if (!ensure_device_capacity(ws.beamOffsets, ws.beamOffsetsCap, ws.hBeamOffsets.size())) return failMultiK("alloc-offsets");
    if (!ensure_device_capacity(ws.sinTheta, ws.sinThetaCap, (size_t)nZen + 1)) return failMultiK("alloc-sintheta");
    if (!ensure_device_capacity(ws.cosTheta, ws.cosThetaCap, (size_t)nZen + 1)) return failMultiK("alloc-costheta");
    if (!ensure_device_capacity(ws.sinPhi, ws.sinPhiCap, (size_t)nAz)) return failMultiK("alloc-sinphi");
    if (!ensure_device_capacity(ws.cosPhi, ws.cosPhiCap, (size_t)nAz)) return failMultiK("alloc-cosphi");
    if (!ensure_device_capacity(ws.vf, ws.vfCap, (size_t)gridCount * 3)) return failMultiK("alloc-vf");
    if (!ensure_device_capacity(ws.m, ws.mCap, mCount)) return failMultiK("alloc-mueller");
    if (timing) tEnsure = gpu_now_ms() - t0;

    t0 = timing ? gpu_now_ms() : 0.0;
    const double gridSignature = gpu_grid_signature(m_sphere, nZen);
    if (ws.gridNAz != nAz || ws.gridNZen != nZen || fabs(ws.gridSignature - gridSignature) > 1e-12)
    {
        std::vector<GpuReal> hSinTheta(nZen + 1), hCosTheta(nZen + 1);
        for (int t = 0; t <= nZen; ++t)
        {
            GpuReal theta = (GpuReal)gpu_theta(m_sphere, t);
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
        std::vector<GpuReal> hVf((size_t)gridCount * 3);
        for (int p = 0; p < nAz; ++p)
            for (int t = 0; t <= nZen; ++t)
            {
                int grid = p * (nZen + 1) + t;
                Point3d &v = m_sphere.vf[p][t];
                hVf[grid * 3 + 0] = v.x;
                hVf[grid * 3 + 1] = v.y;
                hVf[grid * 3 + 2] = v.z;
            }
        if (cudaMemcpy(ws.sinTheta, hSinTheta.data(), hSinTheta.size() * sizeof(GpuReal), cudaMemcpyHostToDevice) != cudaSuccess) return failMultiK("copy-sintheta");
        if (cudaMemcpy(ws.cosTheta, hCosTheta.data(), hCosTheta.size() * sizeof(GpuReal), cudaMemcpyHostToDevice) != cudaSuccess) return failMultiK("copy-costheta");
        if (cudaMemcpy(ws.sinPhi, hSinPhi.data(), hSinPhi.size() * sizeof(GpuReal), cudaMemcpyHostToDevice) != cudaSuccess) return failMultiK("copy-sinphi");
        if (cudaMemcpy(ws.cosPhi, hCosPhi.data(), hCosPhi.size() * sizeof(GpuReal), cudaMemcpyHostToDevice) != cudaSuccess) return failMultiK("copy-cosphi");
        if (cudaMemcpy(ws.vf, hVf.data(), hVf.size() * sizeof(GpuReal), cudaMemcpyHostToDevice) != cudaSuccess) return failMultiK("copy-vf");
        ws.gridNAz = nAz;
        ws.gridNZen = nZen;
        ws.gridSignature = gridSignature;
    }
    if (timing) tGrid = gpu_now_ms() - t0;

    t0 = timing ? gpu_now_ms() : 0.0;
    if (cudaMemcpy(ws.beams, ws.hBeams.data(), ws.hBeams.size() * sizeof(GpuBeam), cudaMemcpyHostToDevice) != cudaSuccess) return failMultiK("copy-beams");
    if (cudaMemcpy(ws.weights, ws.hWeights.data(), ws.hWeights.size() * sizeof(GpuReal), cudaMemcpyHostToDevice) != cudaSuccess) return failMultiK("copy-weights");
    if (cudaMemcpy(ws.scales, ws.hScales.data(), ws.hScales.size() * sizeof(GpuReal), cudaMemcpyHostToDevice) != cudaSuccess) return failMultiK("copy-scales");
    if (!ws.hAbsPaths.empty()
        && cudaMemcpy(ws.absPaths, ws.hAbsPaths.data(),
                      ws.hAbsPaths.size() * sizeof(double),
                      cudaMemcpyHostToDevice) != cudaSuccess)
        return failMultiK("copy-abs-paths");
    if (cudaMemcpy(ws.beamOffsets, ws.hBeamOffsets.data(), ws.hBeamOffsets.size() * sizeof(int), cudaMemcpyHostToDevice) != cudaSuccess) return failMultiK("copy-offsets");
    if (cudaMemset(ws.m, 0, mCount * sizeof(GpuReal)) != cudaSuccess) return failMultiK("memset-m");
    if (timing) tCopy = gpu_now_ms() - t0;

    t0 = timing ? gpu_now_ms() : 0.0;
    int block = gpu_block_size();
    long long total = (long long)nSizes * nOrient * gridCount;
    int kernelGrid = (int)((total + block - 1) / block);
    diffraction_grid_mueller_multik_kernel<<<kernelGrid, block>>>(
        ws.beams, ws.beamOffsets, ws.sinTheta, ws.cosTheta, ws.sinPhi,
        ws.cosPhi, ws.vf, ws.weights, ws.scales,
        ws.hAbsPaths.empty() ? nullptr : ws.absPaths,
        nAz, nZen, nOrient, nSizes,
        m_waveIndex, m_wi2, gpu_effective_eps1(m_eps1), m_eps2,
        real(m_complWave), imag(m_complWave),
        real(m_invComplWave), imag(m_invComplWave),
        m_legacySign ? 1 : 0, (GpuReal)AbsorptionCoefficient(), ws.m);
    cudaError_t err = cudaGetLastError();
    if (!gpu_report_cuda_error(err, "multik diffraction kernel launch"))
        return failMultiK("kernel-launch");
    if (timing)
    {
        err = cudaDeviceSynchronize();
        if (!gpu_report_cuda_error(err, "multik diffraction kernel sync"))
            return failMultiK("kernel-sync");
        tKernel = gpu_now_ms() - t0;
    }

    t0 = timing ? gpu_now_ms() : 0.0;
    ws.hM.resize(mCount);
    if (cudaMemcpy(ws.hM.data(), ws.m, mCount * sizeof(GpuReal), cudaMemcpyDeviceToHost) != cudaSuccess)
        return failMultiK("copy-m-host");
    if (timing) tD2h = gpu_now_ms() - t0;

    t0 = timing ? gpu_now_ms() : 0.0;
    for (int s = 0; s < nSizes; ++s)
    {
        const size_t sizeBase = (size_t)s * gridCount * 16;
        for (int p = 0; p < nAz; ++p)
            for (int t = 0; t <= nZen; ++t)
            {
                int grid = p * (nZen + 1) + t;
                const GpuReal *m = &ws.hM[sizeBase + (size_t)grid * 16];
                double *cell = localMs[s].RawCell(p, t);
                for (int k = 0; k < 16; ++k)
                    cell[k] += m[k];
            }
    }
    if (timing)
    {
        tAdd = gpu_now_ms() - t0;
        std::fprintf(stderr,
                     "GPU timing multik orient=%d sizes=%d beams=%zu nAz=%d nZen=%d pack=%.3fms ensure=%.3fms grid=%.3fms copy=%.3fms kernels=%.3fms d2h=%.3fms add=%.3fms total=%.3fms\n",
                     nOrient, nSizes, nBeams, nAz, nZen, tPack, tEnsure,
                     tGrid, tCopy, tKernel, tD2h, tAdd, gpu_now_ms() - tStart);
    }
    return true;
}

bool HandlerPO::HandleOrientationsToLocalGpu(const std::vector<PreparedOrientation> &prepared,
                                             int start,
                                             int count,
                                             Arr2D &localM,
                                             Arr2D &localM_noshadow,
                                             double scale,
                                             double waveIndex)
{
    auto failDirect = [](const char *where) {
        if (gpu_fft_debug_enabled())
            std::fprintf(stderr, "GPU direct failed at %s\n", where);
        return false;
    };
    if (!isCoh)
        return failDirect("not-coherent");

    if (start < 0 || count < 0 || start > (int)prepared.size())
        return failDirect("bad-range");
    const int nOrient = std::min(count, (int)prepared.size() - start);
    const int nAz = m_sphere.nAzimuth;
    const int nZen = m_sphere.nZenith;
    const int gridCount = nAz * (nZen + 1);
    if (nOrient == 0 || gridCount == 0)
        return true;
    const bool computeNoShadow = ComputeNoShadow();

    if (!g_gpuMultiWorker)
    {
        const int nDevices = gpu_multi_device_count(nOrient);
        if (nDevices > 1)
        {
            int savedDevice = 0;
            cudaGetDevice(&savedDevice);
            static bool printedMulti = false;
            if (!printedMulti || gpu_multi_debug_enabled())
            {
                std::fprintf(stderr,
                             "GPU multi-device: using %d visible CUDA devices for %d orientations in this batch "
                             "(disable with MBS_GPU_MULTI=0, cap with MBS_GPU_MULTI_MAX=N).\n",
                             nDevices, nOrient);
                printedMulti = true;
            }

            std::vector<Arr2D> partialM;
            std::vector<Arr2D> partialMns;
            partialM.reserve(nDevices);
            partialMns.reserve(nDevices);
            for (int i = 0; i < nDevices; ++i)
            {
                partialM.emplace_back(nAz + 1, nZen + 1, 4, 4);
                partialM.back().ClearArr();
                if (computeNoShadow)
                {
                    partialMns.emplace_back(nAz + 1, nZen + 1, 4, 4);
                    partialMns.back().ClearArr();
                }
                else
                {
                    partialMns.emplace_back(0, 0, 0, 0);
                }
            }

            std::vector<std::future<bool>> jobs;
            jobs.reserve(nDevices);
            for (int dev = 0; dev < nDevices; ++dev)
            {
                const int begin = start + (nOrient * dev) / nDevices;
                const int end = start + (nOrient * (dev + 1)) / nDevices;
                const int subCount = end - begin;
                jobs.push_back(std::async(std::launch::async,
                    [this, &prepared, &partialM, &partialMns, dev, begin, subCount,
                     scale, waveIndex]() -> bool {
                        if (subCount <= 0)
                            return true;
                        if (cudaSetDevice(dev) != cudaSuccess)
                            return false;
                        g_gpuMultiWorker = true;
                        bool ok = this->HandleOrientationsToLocalGpu(
                            prepared, begin, subCount, partialM[dev],
                            partialMns[dev], scale, waveIndex);
                        g_gpuMultiWorker = false;
                        return ok;
                    }));
            }

            bool ok = true;
            for (auto &job : jobs)
                ok = job.get() && ok;
            cudaSetDevice(savedDevice);
            if (!ok)
                return false;

            for (int dev = 0; dev < nDevices; ++dev)
            {
                add_arr2d_inplace(partialM[dev], nAz, nZen, localM);
                if (computeNoShadow)
                    add_arr2d_inplace(partialMns[dev], nAz, nZen, localM_noshadow);
            }
            return true;
        }
    }

    const bool timing = gpu_timing_enabled();
    const double tStart = timing ? gpu_now_ms() : 0.0;
    double tCount = 0.0, tPack = 0.0, tEnsure = 0.0, tGrid = 0.0;
    double tCopy = 0.0, tKernel = 0.0, tD2h = 0.0, tAdd = 0.0;

    size_t nBeams = 0;
    size_t nBeams8 = 0;
    size_t nBeamsLarge = 0;
    bool allBeam8 = true;
    for (int oi = 0; oi < nOrient; ++oi)
    {
        const PreparedOrientation &po = prepared[start + oi];
        for (const PreparedBeam &pb : po.beams)
        {
            if (!pb.edgeData.valid || pb.edgeData.nVertices <= 0 || pb.edgeData.nVertices > 32)
                return failDirect("bad-edge-data");
            if (pb.edgeData.nVertices > 8)
            {
                allBeam8 = false;
                ++nBeamsLarge;
            }
            else
            {
                ++nBeams8;
            }
            ++nBeams;
        }
    }
    if (nBeams == 0)
        return true;
    if (timing) tCount = gpu_now_ms() - tStart;

    GpuWorkspace &ws = g_gpuWorkspace;
    double t0 = timing ? gpu_now_ms() : 0.0;
    const bool scaleOnPack = (fabs(scale - 1.0) > 1e-15);
    const double scale2 = scale * scale;
    const int packNoAtomicsMode = gpu_no_atomics_mode();
    const bool packNoAtomics = (packNoAtomicsMode >= 0)
        ? (packNoAtomicsMode != 0)
        : !computeNoShadow;
    const int packFusedMuellerMode = gpu_fused_mueller_mode();
    const bool packFusedMueller = packNoAtomics && ((packFusedMuellerMode >= 0)
        ? (packFusedMuellerMode != 0)
        : true);
    const int packStageMuellerMode = gpu_stage_mueller_mode();
    const bool packStageMueller = packFusedMueller && !computeNoShadow
        && ((packStageMuellerMode >= 0) ? (packStageMuellerMode != 0) : false);
    const int packNoVertexCacheMode = gpu_no_vertex_cache_mode();
    const bool packNoVertexCache = (packNoVertexCacheMode >= 0)
        ? (packNoVertexCacheMode != 0)
        : false;
    const bool canUseBeam8 = packFusedMueller && !computeNoShadow
        && !packStageMueller && !packNoVertexCache;
    const bool packBeam8 = allBeam8 && canUseBeam8;
    const bool packMixedBeam8 = !allBeam8 && canUseBeam8
                             && nBeams8 > 0 && nBeamsLarge > 0;
    if (packBeam8)
        ws.hBeams8.resize(nBeams);
    else if (packMixedBeam8)
    {
        ws.hBeams8.resize(nBeams8);
        ws.hBeams.resize(nBeamsLarge);
    }
    else
        ws.hBeams.resize(nBeams);
    ws.hWeights.assign(nOrient, 0.0);
    ws.hBeamOffsets.assign(nOrient + 1, 0);
    if (packMixedBeam8)
        ws.hBeamOffsets8.assign(nOrient + 1, 0);
    std::vector<GpuBeam> &hBeams = ws.hBeams;
    std::vector<GpuBeam8> &hBeams8 = ws.hBeams8;
    std::vector<GpuReal> &hWeights = ws.hWeights;
    std::vector<int> &hBeamOffsets = ws.hBeamOffsets;
    std::vector<int> &hBeamOffsets8 = ws.hBeamOffsets8;
    const bool beamStats = gpu_beam_stats_enabled();
    size_t statVertexHist[33] = {};
    size_t statExternal = 0, statInternal = 0;
    size_t statEdgesX = 0, statEdgesY = 0;
    size_t statMinBeams = (size_t)-1, statMaxBeams = 0;
    size_t bi = 0;
    size_t bi8 = 0;
    size_t biLarge = 0;
    for (int oi = 0; oi < nOrient; ++oi)
    {
        hBeamOffsets[oi] = (int)(packMixedBeam8 ? biLarge : bi);
        if (packMixedBeam8)
            hBeamOffsets8[oi] = (int)bi8;
        const PreparedOrientation &po = prepared[start + oi];
        hWeights[oi] = po.sinZenith;
        if (beamStats)
        {
            const size_t beamsInOrient = po.beams.size();
            statMinBeams = std::min(statMinBeams, beamsInOrient);
            statMaxBeams = std::max(statMaxBeams, beamsInOrient);
        }
        for (const PreparedBeam &pb : po.beams)
        {
            if (beamStats)
            {
                ++statVertexHist[pb.edgeData.nVertices];
                if (pb.isExternal)
                    ++statExternal;
                else
                    ++statInternal;
                for (int e = 0; e < pb.edgeData.nVertices; ++e)
                {
                    statEdgesX += pb.edgeData.edge_valid_x[e] ? 1 : 0;
                    statEdgesY += pb.edgeData.edge_valid_y[e] ? 1 : 0;
                }
            }
        }
        if (packMixedBeam8)
        {
            for (const PreparedBeam &pb : po.beams)
            {
                if (pb.edgeData.nVertices <= 8)
                    ++bi8;
                else
                    ++biLarge;
            }
        }
        bi += po.beams.size();
    }
    hBeamOffsets[nOrient] = (int)(packMixedBeam8 ? biLarge : bi);
    if (packMixedBeam8)
        hBeamOffsets8[nOrient] = (int)bi8;

#pragma omp parallel for schedule(static) if(nOrient >= 256 && !g_gpuMultiWorker)
    for (int oi = 0; oi < nOrient; ++oi)
    {
        size_t out = (size_t)hBeamOffsets[oi];
        size_t out8 = packMixedBeam8 ? (size_t)hBeamOffsets8[oi] : 0;
        const PreparedOrientation &po = prepared[start + oi];
        for (const PreparedBeam &pb : po.beams)
        {
            if (packBeam8)
                pack_prepared_gpu_beam8(
                    pb, oi, scale, scale2, scaleOnPack, waveIndex,
                    AbsorptionCoefficient(), hBeams8[out++]);
            else if (packMixedBeam8 && pb.edgeData.nVertices <= 8)
                pack_prepared_gpu_beam8(
                    pb, oi, scale, scale2, scaleOnPack, waveIndex,
                    AbsorptionCoefficient(), hBeams8[out8++]);
            else
                pack_prepared_gpu_beam<GpuBeam, 32>(
                    pb, oi, scale, scale2, scaleOnPack, waveIndex,
                    AbsorptionCoefficient(), hBeams[out++]);
        }
    }
    if (beamStats)
    {
        std::fprintf(stderr,
                     "GPU beam stats orient=%d beams=%zu beams/orient min=%zu mean=%.3f max=%zu external=%zu internal=%zu valid_edges_x/beams=%.3f valid_edges_y/beams=%.3f vertices:",
                     nOrient, nBeams,
                     statMinBeams == (size_t)-1 ? 0 : statMinBeams,
                     nBeams == 0 ? 0.0 : (double)nBeams / (double)nOrient,
                     statMaxBeams, statExternal, statInternal,
                     nBeams == 0 ? 0.0 : (double)statEdgesX / (double)nBeams,
                     nBeams == 0 ? 0.0 : (double)statEdgesY / (double)nBeams);
        for (int nv = 1; nv <= 32; ++nv)
            if (statVertexHist[nv] != 0)
                std::fprintf(stderr, " %d:%zu", nv, statVertexHist[nv]);
        std::fprintf(stderr, "\n");
    }
    if (timing) tPack = gpu_now_ms() - t0;

    const size_t jCount = (size_t)nOrient * gridCount * 8;
    const size_t mCount = (size_t)gridCount * 16;
    const int noAtomicsMode = gpu_no_atomics_mode();
    bool noAtomics = false;
    if (noAtomicsMode >= 0)
    {
        noAtomics = (noAtomicsMode != 0);
    }
    else if (!computeNoShadow)
    {
        // For full-only output the fused orientation-grid kernel avoids
        // large Jones buffers and a second Mueller pass. Short production
        // probes are consistently faster than the beam-grid atomics path.
        noAtomics = true;
    }
    else
    {
        size_t freeBytes = 0, totalBytes = 0;
        bool preferAtomics = false;
        if (cudaMemGetInfo(&freeBytes, &totalBytes) == cudaSuccess && freeBytes > 0)
        {
            const size_t jonesBytes =
                jCount * sizeof(GpuReal) * (computeNoShadow ? 2 : 1);
            const size_t muellerBytes =
                mCount * sizeof(GpuReal) * (computeNoShadow ? 2 : 1);
            const size_t beamBytes =
                (packMixedBeam8
                    ? nBeams8 * sizeof(GpuBeam8) + nBeamsLarge * sizeof(GpuBeam)
                    : nBeams * (packBeam8 ? sizeof(GpuBeam8) : sizeof(GpuBeam)))
                + hWeights.size() * sizeof(GpuReal)
                + hBeamOffsets.size() * sizeof(int)
                + (packMixedBeam8 ? hBeamOffsets8.size() * sizeof(int) : 0);
            const size_t gridBytes =
                ((size_t)nZen + 1) * 2 * sizeof(GpuReal)
                + (size_t)nAz * 2 * sizeof(GpuReal)
                + (size_t)gridCount * 3 * sizeof(GpuReal);
            const size_t needBytes =
                jonesBytes + muellerBytes + beamBytes + gridBytes + (128ULL << 20);
            preferAtomics = needBytes < (size_t)(freeBytes * gpu_memory_fraction());
        }
        noAtomics = !preferAtomics;
    }
    const int fusedMuellerMode = gpu_fused_mueller_mode();
    const bool fusedMueller = noAtomics && ((fusedMuellerMode >= 0)
        ? (fusedMuellerMode != 0)
        : true);
    const int stageMuellerMode = gpu_stage_mueller_mode();
    const bool stageMueller = fusedMueller && !computeNoShadow && ((stageMuellerMode >= 0)
        ? (stageMuellerMode != 0)
        : false);
    const int noVertexCacheMode = gpu_no_vertex_cache_mode();
    const bool noVertexCache = (noVertexCacheMode >= 0)
        ? (noVertexCacheMode != 0)
        : false;
    const bool useBeam8 = packBeam8 && fusedMueller && !computeNoShadow
        && !stageMueller && !noVertexCache;
    const bool useMixedBeam8 = packMixedBeam8 && fusedMueller && !computeNoShadow
        && !stageMueller && !noVertexCache;
    t0 = timing ? gpu_now_ms() : 0.0;
    if (useBeam8)
    {
        if (!ensure_device_capacity(ws.beams8, ws.beam8Cap, nBeams)) return failDirect("alloc-beams8");
    }
    else if (useMixedBeam8)
    {
        if (!ensure_device_capacity(ws.beams8, ws.beam8Cap, nBeams8)) return failDirect("alloc-beams8");
        if (!ensure_device_capacity(ws.beams, ws.beamCap, nBeamsLarge)) return failDirect("alloc-beams");
    }
    else if (!ensure_device_capacity(ws.beams, ws.beamCap, nBeams)) return failDirect("alloc-beams");
    if (!ensure_device_capacity(ws.weights, ws.weightsCap, hWeights.size())) return failDirect("alloc-weights");
    if (!ensure_device_capacity(ws.beamOffsets, ws.beamOffsetsCap, hBeamOffsets.size())) return failDirect("alloc-offsets");
    if (useMixedBeam8 && !ensure_device_capacity(ws.beamOffsets8, ws.beamOffsets8Cap, hBeamOffsets8.size())) return failDirect("alloc-offsets8");
    if (!ensure_device_capacity(ws.sinTheta, ws.sinThetaCap, (size_t)nZen + 1)) return failDirect("alloc-sintheta");
    if (!ensure_device_capacity(ws.cosTheta, ws.cosThetaCap, (size_t)nZen + 1)) return failDirect("alloc-costheta");
    if (!ensure_device_capacity(ws.sinPhi, ws.sinPhiCap, (size_t)nAz)) return failDirect("alloc-sinphi");
    if (!ensure_device_capacity(ws.cosPhi, ws.cosPhiCap, (size_t)nAz)) return failDirect("alloc-cosphi");
    if (!ensure_device_capacity(ws.vf, ws.vfCap, (size_t)gridCount * 3)) return failDirect("alloc-vf");
    if (!fusedMueller && !ensure_device_capacity(ws.j, ws.jCap, jCount)) return failDirect("alloc-jones");
    if (!fusedMueller && computeNoShadow && !ensure_device_capacity(ws.jNoShadow, ws.jNoShadowCap, jCount)) return failDirect("alloc-jones-ns");
    if (!ensure_device_capacity(ws.m, ws.mCap, mCount)) return failDirect("alloc-mueller");
    if (stageMueller && !ensure_device_capacity(ws.mOrient, ws.mOrientCap, (size_t)nOrient * gridCount * 16)) return failDirect("alloc-orient-mueller");
    if (computeNoShadow && !ensure_device_capacity(ws.mNoShadow, ws.mNoShadowCap, mCount)) return failDirect("alloc-mueller-ns");
    if (timing) tEnsure = gpu_now_ms() - t0;

    t0 = timing ? gpu_now_ms() : 0.0;
    const double gridSignature = gpu_grid_signature(m_sphere, nZen);
    if (ws.gridNAz != nAz || ws.gridNZen != nZen || fabs(ws.gridSignature - gridSignature) > 1e-12)
    {
        std::vector<GpuReal> hSinTheta(nZen + 1), hCosTheta(nZen + 1);
        for (int t = 0; t <= nZen; ++t)
        {
            GpuReal theta = (GpuReal)gpu_theta(m_sphere, t);
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
        ws.gridSignature = gridSignature;
    }
    if (timing) tGrid = gpu_now_ms() - t0;

    t0 = timing ? gpu_now_ms() : 0.0;
    if (useBeam8)
    {
        if (cudaMemcpy(ws.beams8, ws.hBeams8.data(), ws.hBeams8.size() * sizeof(GpuBeam8), cudaMemcpyHostToDevice) != cudaSuccess) return false;
    }
    else if (useMixedBeam8)
    {
        if (cudaMemcpy(ws.beams8, ws.hBeams8.data(), ws.hBeams8.size() * sizeof(GpuBeam8), cudaMemcpyHostToDevice) != cudaSuccess) return false;
        if (cudaMemcpy(ws.beams, hBeams.data(), hBeams.size() * sizeof(GpuBeam), cudaMemcpyHostToDevice) != cudaSuccess) return false;
    }
    else if (cudaMemcpy(ws.beams, hBeams.data(), hBeams.size() * sizeof(GpuBeam), cudaMemcpyHostToDevice) != cudaSuccess) return false;
    if (cudaMemcpy(ws.weights, hWeights.data(), hWeights.size() * sizeof(GpuReal), cudaMemcpyHostToDevice) != cudaSuccess) return failDirect("copy-weights");
    if (cudaMemcpy(ws.beamOffsets, hBeamOffsets.data(), hBeamOffsets.size() * sizeof(int), cudaMemcpyHostToDevice) != cudaSuccess) return failDirect("copy-offsets");
    if (useMixedBeam8 && cudaMemcpy(ws.beamOffsets8, hBeamOffsets8.data(), hBeamOffsets8.size() * sizeof(int), cudaMemcpyHostToDevice) != cudaSuccess) return failDirect("copy-offsets8");
    if (!fusedMueller && cudaMemset(ws.j, 0, jCount * sizeof(GpuReal)) != cudaSuccess) return failDirect("memset-j");
    if (!fusedMueller && computeNoShadow && cudaMemset(ws.jNoShadow, 0, jCount * sizeof(GpuReal)) != cudaSuccess) return failDirect("memset-j-ns");
    if (cudaMemset(ws.m, 0, mCount * sizeof(GpuReal)) != cudaSuccess) return failDirect("memset-m");
    if (computeNoShadow && cudaMemset(ws.mNoShadow, 0, mCount * sizeof(GpuReal)) != cudaSuccess) return failDirect("memset-m-ns");
    if (timing) tCopy = gpu_now_ms() - t0;

    t0 = timing ? gpu_now_ms() : 0.0;
    int block = gpu_block_size();
    if (!gpu_block_size_overridden() && fusedMueller && !computeNoShadow && nOrient >= 512)
        block = 128;
    long long diffractionTotal = (long long)(noAtomics ? nOrient : (int)nBeams) * gridCount;
    int diffractionGrid = (int)((diffractionTotal + block - 1) / block);
    if (fusedMueller)
    {
        if (computeNoShadow)
        {
            diffraction_grid_mueller_kernel<<<diffractionGrid, block>>>(
                ws.beams, ws.beamOffsets, ws.sinTheta, ws.cosTheta, ws.sinPhi, ws.cosPhi,
                ws.vf, ws.weights, nAz, nZen, nOrient, m_waveIndex, m_wi2,
                gpu_effective_eps1(m_eps1), m_eps2, real(m_complWave), imag(m_complWave),
                real(m_invComplWave), imag(m_invComplWave), m_legacySign ? 1 : 0,
                ws.m, ws.mNoShadow);
        }
        else
        {
            if (stageMueller)
            {
                diffraction_grid_mueller_orient_kernel<<<diffractionGrid, block>>>(
                    ws.beams, ws.beamOffsets, ws.sinTheta, ws.cosTheta, ws.sinPhi, ws.cosPhi,
                    ws.vf, ws.weights, nAz, nZen, nOrient, m_waveIndex, m_wi2,
                    gpu_effective_eps1(m_eps1), m_eps2, real(m_complWave), imag(m_complWave),
                    real(m_invComplWave), imag(m_invComplWave), m_legacySign ? 1 : 0,
                    ws.mOrient);
                int reduceGrid = (gridCount + block - 1) / block;
                reduce_mueller_orient_kernel<<<reduceGrid, block>>>(ws.mOrient, nOrient, gridCount, ws.m);
            }
            else if (useBeam8)
            {
                diffraction_grid_mueller_full8_kernel<<<diffractionGrid, block>>>(
                    ws.beams8, ws.beamOffsets, ws.sinTheta, ws.cosTheta, ws.sinPhi, ws.cosPhi,
                    ws.vf, ws.weights, nAz, nZen, nOrient, m_waveIndex, m_wi2,
                    gpu_effective_eps1(m_eps1), m_eps2, real(m_complWave), imag(m_complWave),
                    real(m_invComplWave), imag(m_invComplWave), m_legacySign ? 1 : 0,
                    ws.m);
            }
            else if (useMixedBeam8)
            {
                diffraction_grid_mueller_mixed8_kernel<<<diffractionGrid, block>>>(
                    ws.beams8, ws.beamOffsets8, ws.beams, ws.beamOffsets,
                    ws.sinTheta, ws.cosTheta, ws.sinPhi, ws.cosPhi,
                    ws.vf, ws.weights, nAz, nZen, nOrient, m_waveIndex, m_wi2,
                    gpu_effective_eps1(m_eps1), m_eps2, real(m_complWave), imag(m_complWave),
                    real(m_invComplWave), imag(m_invComplWave), m_legacySign ? 1 : 0,
                    ws.m);
            }
            else if (noVertexCache)
            {
                diffraction_grid_mueller_full_nocache_kernel<<<diffractionGrid, block>>>(
                    ws.beams, ws.beamOffsets, ws.sinTheta, ws.cosTheta, ws.sinPhi, ws.cosPhi,
                    ws.vf, ws.weights, nAz, nZen, nOrient, m_waveIndex, m_wi2,
                    gpu_effective_eps1(m_eps1), m_eps2, real(m_complWave), imag(m_complWave),
                    real(m_invComplWave), imag(m_invComplWave), m_legacySign ? 1 : 0,
                    ws.m);
            }
            else
            {
                diffraction_grid_mueller_full_kernel<<<diffractionGrid, block>>>(
                    ws.beams, ws.beamOffsets, ws.sinTheta, ws.cosTheta, ws.sinPhi, ws.cosPhi,
                    ws.vf, ws.weights, nAz, nZen, nOrient, m_waveIndex, m_wi2,
                    gpu_effective_eps1(m_eps1), m_eps2, real(m_complWave), imag(m_complWave),
                    real(m_invComplWave), imag(m_invComplWave), m_legacySign ? 1 : 0,
                    ws.m);
            }
        }
    }
    else if (noAtomics)
    {
        diffraction_grid_kernel<<<diffractionGrid, block>>>(
            ws.beams, ws.beamOffsets, ws.sinTheta, ws.cosTheta, ws.sinPhi, ws.cosPhi, ws.vf,
            nAz, nZen, nOrient, m_waveIndex, m_wi2, gpu_effective_eps1(m_eps1), m_eps2,
            real(m_complWave), imag(m_complWave),
            real(m_invComplWave), imag(m_invComplWave),
            m_legacySign ? 1 : 0, ws.j, computeNoShadow ? ws.jNoShadow : nullptr);
    }
    else
    {
        diffraction_kernel<<<diffractionGrid, block>>>(
            ws.beams, (int)hBeams.size(), ws.sinTheta, ws.cosTheta, ws.sinPhi, ws.cosPhi, ws.vf,
            nAz, nZen, m_waveIndex, m_wi2, gpu_effective_eps1(m_eps1), m_eps2,
            real(m_complWave), imag(m_complWave),
            real(m_invComplWave), imag(m_invComplWave),
            m_legacySign ? 1 : 0, ws.j, computeNoShadow ? ws.jNoShadow : nullptr);
    }
    cudaError_t err = cudaGetLastError();
    if (!gpu_report_cuda_error(err, "diffraction kernel launch"))
        return failDirect("diffraction-launch");

    if (!fusedMueller)
    {
        long long muellerTotal = (long long)nOrient * gridCount;
        int muellerGrid = (int)((muellerTotal + block - 1) / block);
        mueller_batch_kernel<<<muellerGrid, block>>>(ws.j, computeNoShadow ? ws.jNoShadow : nullptr, ws.weights,
                                                     nOrient, gridCount,
                                                     ws.m, computeNoShadow ? ws.mNoShadow : nullptr);
        err = cudaGetLastError();
        if (!gpu_report_cuda_error(err, "mueller kernel launch"))
            return failDirect("mueller-launch");
    }
    if (timing)
    {
        err = cudaDeviceSynchronize();
        if (!gpu_report_cuda_error(err, "mueller kernel sync"))
            return failDirect("kernel-sync");
        tKernel = gpu_now_ms() - t0;
    }

    t0 = timing ? gpu_now_ms() : 0.0;
    ws.hM.resize(mCount);
    if (computeNoShadow)
        ws.hMNoShadow.resize(mCount);
    std::vector<GpuReal> &hM = ws.hM;
    std::vector<GpuReal> &hMns = ws.hMNoShadow;
    if (cudaMemcpy(hM.data(), ws.m, mCount * sizeof(GpuReal), cudaMemcpyDeviceToHost) != cudaSuccess) return failDirect("copy-m-host");
    if (computeNoShadow && cudaMemcpy(hMns.data(), ws.mNoShadow, mCount * sizeof(GpuReal), cudaMemcpyDeviceToHost) != cudaSuccess) return failDirect("copy-m-ns-host");
    if (timing) tD2h = gpu_now_ms() - t0;

    t0 = timing ? gpu_now_ms() : 0.0;
    for (int p = 0; p < nAz; ++p)
        for (int t = 0; t <= nZen; ++t)
        {
            int grid = p * (nZen + 1) + t;
            const GpuReal *m = &hM[grid * 16];
            double *cell = localM.RawCell(p, t);
            for (int k = 0; k < 16; ++k)
                cell[k] += m[k];
            if (computeNoShadow)
            {
                const GpuReal *mn = &hMns[grid * 16];
                double *cellNs = localM_noshadow.RawCell(p, t);
                for (int k = 0; k < 16; ++k)
                    cellNs[k] += mn[k];
            }
        }
    if (timing)
    {
        tAdd = gpu_now_ms() - t0;
        std::fprintf(stderr,
                     "GPU timing diffract orient=%d beams=%zu nAz=%d nZen=%d count=%.3fms pack=%.3fms ensure=%.3fms grid=%.3fms copy=%.3fms kernels=%.3fms d2h=%.3fms add=%.3fms total=%.3fms\n",
                     nOrient, nBeams, nAz, nZen, tCount, tPack, tEnsure,
                     tGrid, tCopy, tKernel, tD2h, tAdd, gpu_now_ms() - tStart);
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
            b.nEdgeX = 0;
            b.nEdgeY = 0;
            b.isExternal = pb.isExternal ? 1 : 0;
            b.orientation = oi;
            for (int e = 0; e < 32; ++e)
            {
                b.x[e] = pb.edgeData.x[e];
                b.y[e] = pb.edgeData.y[e];
                b.slope_yx[e] = pb.edgeData.slope_yx[e];
                b.slope_xy[e] = pb.edgeData.slope_xy[e];
                const bool validX = e < b.nVertices && pb.edgeData.edge_valid_x[e];
                const bool validY = e < b.nVertices && pb.edgeData.edge_valid_y[e];
                if (validX)
                    b.edge_valid_x[b.nEdgeX++] = (unsigned char)e;
                if (validY)
                    b.edge_valid_y[b.nEdgeY++] = (unsigned char)e;
            }
            b.bdx = (GpuReal)(pb.bdx * pb.horAx + pb.bdy * pb.horAy + pb.bdz * pb.horAz);
            b.bdy = (GpuReal)(pb.bdx * pb.verAx + pb.bdy * pb.verAy + pb.bdz * pb.verAz);
            b.bdz = 0.0;
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
        double theta = theta_rads[t];
        if (theta <= 1e-14)
            theta = 1e-6;
        else if (theta >= M_PI - 1e-14)
            theta = M_PI - 1e-6;
        hSinTheta[t] = (GpuReal)sin(theta);
        hCosTheta[t] = (GpuReal)cos(theta);
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

    cudaError_t errTheta = cudaMemcpy(ws.beams, ws.hBeams.data(), ws.hBeams.size() * sizeof(GpuBeam), cudaMemcpyHostToDevice);
    if (!gpu_report_cuda_error(errTheta, "theta copy-beams")) return false;
    errTheta = cudaMemcpy(ws.weights, ws.hWeights.data(), ws.hWeights.size() * sizeof(GpuReal), cudaMemcpyHostToDevice);
    if (!gpu_report_cuda_error(errTheta, "theta copy-weights")) return false;
    errTheta = cudaMemcpy(ws.sinTheta, hSinTheta.data(), hSinTheta.size() * sizeof(GpuReal), cudaMemcpyHostToDevice);
    if (!gpu_report_cuda_error(errTheta, "theta copy-sintheta")) return false;
    errTheta = cudaMemcpy(ws.cosTheta, hCosTheta.data(), hCosTheta.size() * sizeof(GpuReal), cudaMemcpyHostToDevice);
    if (!gpu_report_cuda_error(errTheta, "theta copy-costheta")) return false;
    errTheta = cudaMemcpy(ws.sinPhi, hSinPhi.data(), hSinPhi.size() * sizeof(GpuReal), cudaMemcpyHostToDevice);
    if (!gpu_report_cuda_error(errTheta, "theta copy-sinphi")) return false;
    errTheta = cudaMemcpy(ws.cosPhi, hCosPhi.data(), hCosPhi.size() * sizeof(GpuReal), cudaMemcpyHostToDevice);
    if (!gpu_report_cuda_error(errTheta, "theta copy-cosphi")) return false;
    errTheta = cudaMemcpy(ws.vf, hVf.data(), hVf.size() * sizeof(GpuReal), cudaMemcpyHostToDevice);
    if (!gpu_report_cuda_error(errTheta, "theta copy-vf")) return false;
    errTheta = cudaMemset(ws.m, 0, (size_t)nPoints * sizeof(GpuReal));
    if (!gpu_report_cuda_error(errTheta, "theta memset-m")) return false;

    int block = gpu_block_size();
    long long total = (long long)nBeams * gridCount;
    cudaDeviceProp prop;
    int maxGridX = 2147483647;
    int dev = 0;
    if (cudaGetDevice(&dev) == cudaSuccess
        && cudaGetDeviceProperties(&prop, dev) == cudaSuccess
        && prop.maxGridSize[0] > 0)
        maxGridX = prop.maxGridSize[0];
    const long long maxBlocksPerLaunch =
        std::min<long long>(maxGridX, 2147483000LL);
    const long long maxWorkPerLaunch = maxBlocksPerLaunch * block;
    for (long long startIdx = 0; startIdx < total; startIdx += maxWorkPerLaunch)
    {
        const long long remaining = total - startIdx;
        const int launchGrid =
            (int)std::min(maxBlocksPerLaunch,
                          (remaining + block - 1) / block);
        theta_m11_kernel<<<launchGrid, block>>>(
            ws.beams, (int)nBeams, ws.weights,
            ws.sinTheta, ws.cosTheta, ws.sinPhi, ws.cosPhi, ws.vf,
            nAz, nPoints, m_waveIndex, m_wi2, gpu_effective_eps1(m_eps1), m_eps2,
            real(m_complWave), imag(m_complWave),
            real(m_invComplWave), imag(m_invComplWave),
            m_legacySign ? 1 : 0, startIdx, ws.m);

        errTheta = cudaGetLastError();
        if (!gpu_report_cuda_error(errTheta, "theta m11 kernel launch"))
            return false;
        errTheta = cudaDeviceSynchronize();
        if (!gpu_report_cuda_error(errTheta, "theta m11 kernel sync"))
            return false;
    }

    std::vector<GpuReal> hM11(nPoints);
    errTheta = cudaMemcpy(hM11.data(), ws.m, (size_t)nPoints * sizeof(GpuReal), cudaMemcpyDeviceToHost);
    if (!gpu_report_cuda_error(errTheta, "theta copy-m11-host"))
        return false;
    for (int i = 0; i < nPoints; ++i)
        m11_out[i] = hM11[i];

    ws.gridNAz = -1;
    ws.gridNZen = -1;
    ws.gridSignature = 0.0;
    return true;
}
