#include <cufft.h>
#include <cuda_runtime.h>

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <random>
#include <queue>
#include <vector>

struct RectBeam
{
    double width, height;
    double angleRad;
    double centerX, centerY;
    double phaseFx, phaseFy;
    double ampR, ampI;
};

__device__ static cufftDoubleComplex d_cmul(cufftDoubleComplex a, cufftDoubleComplex b)
{
    cufftDoubleComplex out;
    out.x = a.x * b.x - a.y * b.y;
    out.y = a.x * b.y + a.y * b.x;
    return out;
}

__device__ static int d_freq_index(int k, int n)
{
    return (k >= 0) ? k : n + k;
}

__device__ static double d_sinc_pi(double x)
{
    if (fabs(x) < 1e-12)
        return 1.0;
    return sin(M_PI * x) / (M_PI * x);
}

__device__ static double d_lanczos_weight(double x, int radius)
{
    double ax = fabs(x);
    if (ax >= radius)
        return 0.0;
    return d_sinc_pi(x) * d_sinc_pi(x / radius);
}

__device__ static cufftDoubleComplex d_phase_corrected(const cufftDoubleComplex *fft,
                                                       int nx, int ny,
                                                       int kx, int ky,
                                                       double dx, double dy)
{
    int ix = d_freq_index(kx, nx);
    int iy = d_freq_index(ky, ny);
    cufftDoubleComplex z = fft[(size_t)iy * nx + ix];
    double phase = 2.0 * M_PI * (
        ((double)kx / nx) * (0.5 - nx / 2.0) +
        ((double)ky / ny) * (0.5 - ny / 2.0));
    double s, c;
    sincos(phase, &s, &c);
    cufftDoubleComplex out;
    out.x = (z.x * c - z.y * s) * dx * dy;
    out.y = (z.x * s + z.y * c) * dx * dy;
    return out;
}

__device__ static cufftDoubleComplex d_sample_fft_lanczos(const cufftDoubleComplex *fft,
                                                          int nx, int ny,
                                                          double dx, double dy,
                                                          double fx, double fy,
                                                          int radius)
{
    double ux = fx * nx * dx;
    double uy = fy * ny * dy;
    int cx = (int)floor(ux);
    int cy = (int)floor(uy);
    cufftDoubleComplex out{0.0, 0.0};
    double wsum = 0.0;
    for (int jy = cy - radius + 1; jy <= cy + radius; ++jy)
    {
        double wy = d_lanczos_weight(uy - jy, radius);
        if (wy == 0.0) continue;
        for (int ix = cx - radius + 1; ix <= cx + radius; ++ix)
        {
            double wx = d_lanczos_weight(ux - ix, radius);
            double w = wx * wy;
            if (w == 0.0) continue;
            cufftDoubleComplex z = d_phase_corrected(fft, nx, ny, ix, jy, dx, dy);
            out.x += w * z.x;
            out.y += w * z.y;
            wsum += w;
        }
    }
    if (fabs(wsum) > 1e-30)
    {
        out.x /= wsum;
        out.y /= wsum;
    }
    return out;
}

__global__ void angular_sample_abs_kernel(const cufftDoubleComplex *fft,
                                          double *outAbs,
                                          int nx,
                                          int ny,
                                          double dx,
                                          double dy,
                                          double lambda,
                                          int nPhi,
                                          int nTheta,
                                          int interpRadius)
{
    int total = (nTheta + 1) * nPhi;
    int idx = (int)(blockIdx.x * blockDim.x + threadIdx.x);
    if (idx >= total) return;
    int ip = idx % nPhi;
    int it = idx / nPhi;
    double theta = M_PI * it / nTheta;
    double phi = 2.0 * M_PI * ip / nPhi;
    double st = sin(theta);
    double fx = st * cos(phi) / lambda;
    double fy = st * sin(phi) / lambda;
    cufftDoubleComplex z = d_sample_fft_lanczos(fft, nx, ny, dx, dy,
                                                fx, fy, interpRadius);
    outAbs[idx] = hypot(z.x, z.y);
}

__global__ void rasterize_rect_beams_kernel(cufftDoubleComplex *grid,
                                            int nx,
                                            int ny,
                                            double dx,
                                            double dy,
                                            const RectBeam *beams,
                                            int nBeams)
{
    size_t total = (size_t)nx * (size_t)ny;
    size_t idx = (size_t)blockIdx.x * blockDim.x + threadIdx.x;
    if (idx >= total) return;

    int x = (int)(idx % nx);
    int y = (int)(idx / nx);
    double xx = (x + 0.5 - nx / 2.0) * dx;
    double yy = (y + 0.5 - ny / 2.0) * dy;

    cufftDoubleComplex value{0.0, 0.0};
    for (int i = 0; i < nBeams; ++i)
    {
        const RectBeam &beam = beams[i];
        double ca = cos(beam.angleRad), sa = sin(beam.angleRad);
        double rx = xx - beam.centerX;
        double ry = yy - beam.centerY;
        double localX = rx * ca + ry * sa;
        double localY = -rx * sa + ry * ca;
        if (fabs(localX) >= 0.5 * beam.width || fabs(localY) >= 0.5 * beam.height)
            continue;

        double phase = 2.0 * M_PI * (beam.phaseFx * xx + beam.phaseFy * yy);
        double s, c;
        sincos(phase, &s, &c);
        cufftDoubleComplex carrier{c, s};
        cufftDoubleComplex amp{beam.ampR, beam.ampI};
        cufftDoubleComplex add = d_cmul(amp, carrier);
        value.x += add.x;
        value.y += add.y;
    }
    grid[idx] = value;
}

static void cadd(cufftDoubleComplex &a, const cufftDoubleComplex &b)
{
    a.x += b.x;
    a.y += b.y;
}

static cufftDoubleComplex cmul(const cufftDoubleComplex &a, const cufftDoubleComplex &b)
{
    cufftDoubleComplex out;
    out.x = a.x * b.x - a.y * b.y;
    out.y = a.x * b.y + a.y * b.x;
    return out;
}

static void die(const char *msg)
{
    std::fprintf(stderr, "ERROR: %s\n", msg);
    std::exit(1);
}

static void checkCuda(cudaError_t err, const char *where)
{
    if (err != cudaSuccess)
    {
        std::fprintf(stderr, "CUDA error at %s: %s\n", where, cudaGetErrorString(err));
        std::exit(1);
    }
}

static void checkCufft(cufftResult err, const char *where)
{
    if (err != CUFFT_SUCCESS)
    {
        std::fprintf(stderr, "cuFFT error at %s: %d\n", where, (int)err);
        std::exit(1);
    }
}

static double sinc_pi(double x)
{
    if (std::fabs(x) < 1e-12)
        return 1.0;
    return std::sin(M_PI * x) / (M_PI * x);
}

static cufftDoubleComplex analytic_rect_complex(double fx, double fy,
                                                const RectBeam &beam)
{
    double ca = std::cos(beam.angleRad), sa = std::sin(beam.angleRad);
    double dfx = fx - beam.phaseFx;
    double dfy = fy - beam.phaseFy;
    double localFx = dfx * ca + dfy * sa;
    double localFy = -dfx * sa + dfy * ca;
    double amp = beam.width * beam.height
        * sinc_pi(beam.width * localFx)
        * sinc_pi(beam.height * localFy);
    double phase = 2.0 * M_PI * ((beam.phaseFx - fx) * beam.centerX
                               + (beam.phaseFy - fy) * beam.centerY);
    cufftDoubleComplex carrier{amp * std::cos(phase), amp * std::sin(phase)};
    cufftDoubleComplex beamAmp{beam.ampR, beam.ampI};
    return cmul(beamAmp, carrier);
}

static double analytic_rect_abs(double fx, double fy, const std::vector<RectBeam> &beams)
{
    cufftDoubleComplex sum{0.0, 0.0};
    for (const RectBeam &beam : beams)
        cadd(sum, analytic_rect_complex(fx, fy, beam));
    return std::hypot(sum.x, sum.y);
}

static int freq_index(int k, int n)
{
    int idx = (k >= 0) ? k : n + k;
    if (idx < 0 || idx >= n)
        die("frequency index out of range");
    return idx;
}

static cufftDoubleComplex phase_corrected(const std::vector<cufftDoubleComplex> &fft,
                                          int nx, int ny,
                                          int kx, int ky,
                                          double dx, double dy)
{
    int ix = freq_index(kx, nx);
    int iy = freq_index(ky, ny);
    cufftDoubleComplex z = fft[(size_t)iy * nx + ix];

    double phase = 2.0 * M_PI * (
        ((double)kx / nx) * (0.5 - nx / 2.0) +
        ((double)ky / ny) * (0.5 - ny / 2.0));
    double c = std::cos(phase), s = std::sin(phase);
    cufftDoubleComplex out;
    out.x = (z.x * c - z.y * s) * dx * dy;
    out.y = (z.x * s + z.y * c) * dx * dy;
    return out;
}

static cufftDoubleComplex sample_fft_bilinear(const std::vector<cufftDoubleComplex> &fft,
                                              int nx, int ny,
                                              double dx, double dy,
                                              double fx, double fy)
{
    double ux = fx * nx * dx;
    double uy = fy * ny * dy;
    int kx0 = (int)std::floor(ux);
    int ky0 = (int)std::floor(uy);
    double tx = ux - kx0;
    double ty = uy - ky0;

    cufftDoubleComplex z00 = phase_corrected(fft, nx, ny, kx0,     ky0,     dx, dy);
    cufftDoubleComplex z10 = phase_corrected(fft, nx, ny, kx0 + 1, ky0,     dx, dy);
    cufftDoubleComplex z01 = phase_corrected(fft, nx, ny, kx0,     ky0 + 1, dx, dy);
    cufftDoubleComplex z11 = phase_corrected(fft, nx, ny, kx0 + 1, ky0 + 1, dx, dy);

    double w00 = (1.0 - tx) * (1.0 - ty);
    double w10 = tx * (1.0 - ty);
    double w01 = (1.0 - tx) * ty;
    double w11 = tx * ty;

    cufftDoubleComplex out;
    out.x = w00*z00.x + w10*z10.x + w01*z01.x + w11*z11.x;
    out.y = w00*z00.y + w10*z10.y + w01*z01.y + w11*z11.y;
    return out;
}

static double lanczos_weight(double x, int radius)
{
    double ax = std::fabs(x);
    if (ax >= radius)
        return 0.0;
    return sinc_pi(x) * sinc_pi(x / radius);
}

static cufftDoubleComplex sample_fft_lanczos(const std::vector<cufftDoubleComplex> &fft,
                                             int nx, int ny,
                                             double dx, double dy,
                                             double fx, double fy,
                                             int radius)
{
    double ux = fx * nx * dx;
    double uy = fy * ny * dy;
    int cx = (int)std::floor(ux);
    int cy = (int)std::floor(uy);
    cufftDoubleComplex out{0.0, 0.0};
    double wsum = 0.0;
    for (int jy = cy - radius + 1; jy <= cy + radius; ++jy)
    {
        double wy = lanczos_weight(uy - jy, radius);
        if (wy == 0.0) continue;
        for (int ix = cx - radius + 1; ix <= cx + radius; ++ix)
        {
            double wx = lanczos_weight(ux - ix, radius);
            double w = wx * wy;
            if (w == 0.0) continue;
            cufftDoubleComplex z = phase_corrected(fft, nx, ny, ix, jy, dx, dy);
            out.x += w * z.x;
            out.y += w * z.y;
            wsum += w;
        }
    }
    if (std::fabs(wsum) > 1e-30)
    {
        out.x /= wsum;
        out.y /= wsum;
    }
    return out;
}

static void usage(const char *argv0)
{
    std::printf(
        "Usage: %s [--nx N] [--ny N] [--width UM] [--height UM] [--pad P] [--samples K]\n"
        "       [--angle-deg A --center-x X --center-y Y --phase-fx F --phase-fy F]\n"
        "       [--beams N --gpu-raster --gpu-sample]\n"
        "       [--angular --nphi N --ntheta N --lambda UM --interp R --hybrid-rel E --no-analytic]\n"
        "\n"
        "Experimental cuFFT aperture probe. Rasterizes a centered rectangular aperture,\n"
        "runs a 2D complex FFT, and compares the fy=0 central cut to the analytic\n"
        "continuous transform of a possibly rotated/shifted rectangle with a phase ramp.\n"
        "\n"
        "Defaults: --nx 2048 --ny 2048 --width 25.12 --height 8.59 --pad 4 --samples 32\n",
        argv0);
}

int main(int argc, char **argv)
{
    int nx = 2048;
    int ny = 2048;
    int samples = 32;
    int nPhi = 600;
    int nTheta = 180;
    int interpRadius = 4;
    int reportWorst = 0;
    double width = 25.12;
    double height = 8.59;
    double pad = 4.0;
    double lambda = 1.55;
    double angleDeg = 0.0;
    double centerX = 0.0;
    double centerY = 0.0;
    double phaseFx = 0.0;
    double phaseFy = 0.0;
    double hybridRel = 0.02;
    int beams = 1;
    bool angular = false;
    bool gpuRaster = false;
    bool gpuSample = false;
    bool noAnalytic = false;

    for (int i = 1; i < argc; ++i)
    {
        if (std::strcmp(argv[i], "--help") == 0 || std::strcmp(argv[i], "-h") == 0)
        {
            usage(argv[0]);
            return 0;
        }
        else if (std::strcmp(argv[i], "--nx") == 0 && i + 1 < argc)
            nx = std::atoi(argv[++i]);
        else if (std::strcmp(argv[i], "--ny") == 0 && i + 1 < argc)
            ny = std::atoi(argv[++i]);
        else if (std::strcmp(argv[i], "--width") == 0 && i + 1 < argc)
            width = std::atof(argv[++i]);
        else if (std::strcmp(argv[i], "--height") == 0 && i + 1 < argc)
            height = std::atof(argv[++i]);
        else if (std::strcmp(argv[i], "--angle-deg") == 0 && i + 1 < argc)
            angleDeg = std::atof(argv[++i]);
        else if (std::strcmp(argv[i], "--center-x") == 0 && i + 1 < argc)
            centerX = std::atof(argv[++i]);
        else if (std::strcmp(argv[i], "--center-y") == 0 && i + 1 < argc)
            centerY = std::atof(argv[++i]);
        else if (std::strcmp(argv[i], "--phase-fx") == 0 && i + 1 < argc)
            phaseFx = std::atof(argv[++i]);
        else if (std::strcmp(argv[i], "--phase-fy") == 0 && i + 1 < argc)
            phaseFy = std::atof(argv[++i]);
        else if (std::strcmp(argv[i], "--beams") == 0 && i + 1 < argc)
            beams = std::atoi(argv[++i]);
        else if (std::strcmp(argv[i], "--pad") == 0 && i + 1 < argc)
            pad = std::atof(argv[++i]);
        else if (std::strcmp(argv[i], "--samples") == 0 && i + 1 < argc)
            samples = std::atoi(argv[++i]);
        else if (std::strcmp(argv[i], "--nphi") == 0 && i + 1 < argc)
            nPhi = std::atoi(argv[++i]);
        else if (std::strcmp(argv[i], "--ntheta") == 0 && i + 1 < argc)
            nTheta = std::atoi(argv[++i]);
        else if (std::strcmp(argv[i], "--lambda") == 0 && i + 1 < argc)
            lambda = std::atof(argv[++i]);
        else if (std::strcmp(argv[i], "--interp") == 0 && i + 1 < argc)
            interpRadius = std::atoi(argv[++i]);
        else if (std::strcmp(argv[i], "--report-worst") == 0 && i + 1 < argc)
            reportWorst = std::atoi(argv[++i]);
        else if (std::strcmp(argv[i], "--hybrid-rel") == 0 && i + 1 < argc)
            hybridRel = std::atof(argv[++i]);
        else if (std::strcmp(argv[i], "--angular") == 0)
            angular = true;
        else if (std::strcmp(argv[i], "--gpu-raster") == 0)
            gpuRaster = true;
        else if (std::strcmp(argv[i], "--gpu-sample") == 0)
            gpuSample = true;
        else if (std::strcmp(argv[i], "--no-analytic") == 0)
            noAnalytic = true;
        else
            die("bad arguments; use --help");
    }

    if (nx <= 0 || ny <= 0 || width <= 0.0 || height <= 0.0 || pad < 1.0 || samples < 0)
        die("invalid numeric arguments");
    if (beams <= 0)
        die("invalid --beams");
    if (hybridRel < 0.0)
        die("invalid --hybrid-rel");

    const double lx = pad * width;
    const double ly = pad * height;
    const double dx = lx / nx;
    const double dy = ly / ny;
    const size_t count = (size_t)nx * (size_t)ny;

    std::vector<RectBeam> beamSpecs;
    beamSpecs.reserve(beams);
    if (beams == 1)
    {
        RectBeam b;
        b.width = width;
        b.height = height;
        b.angleRad = angleDeg * M_PI / 180.0;
        b.centerX = centerX;
        b.centerY = centerY;
        b.phaseFx = phaseFx;
        b.phaseFy = phaseFy;
        b.ampR = 1.0;
        b.ampI = 0.0;
        beamSpecs.push_back(b);
    }
    else
    {
        std::mt19937 rng(12345);
        std::uniform_real_distribution<double> u01(0.0, 1.0);
        double radiusX = 0.35 * width;
        double radiusY = 0.35 * height;
        for (int i = 0; i < beams; ++i)
        {
            double a = 2.0 * M_PI * u01(rng);
            double r = std::sqrt(u01(rng));
            RectBeam b;
            b.width = width * (0.12 + 0.08 * u01(rng));
            b.height = height * (0.18 + 0.10 * u01(rng));
            b.angleRad = (angleDeg + 180.0 * u01(rng)) * M_PI / 180.0;
            b.centerX = centerX + radiusX * r * std::cos(a);
            b.centerY = centerY + radiusY * r * std::sin(a);
            b.phaseFx = phaseFx + 0.04 * (u01(rng) - 0.5);
            b.phaseFy = phaseFy + 0.04 * (u01(rng) - 0.5);
            double ph = 2.0 * M_PI * u01(rng);
            double amp = 1.0 / std::sqrt((double)beams);
            b.ampR = amp * std::cos(ph);
            b.ampI = amp * std::sin(ph);
            beamSpecs.push_back(b);
        }
    }

    bool needHostFft = !gpuSample || samples > 0 || !gpuRaster;
    std::vector<cufftDoubleComplex> h;
    if (needHostFft)
        h.resize(count);
    std::vector<double> gpuAngularAbs;
    cufftDoubleComplex *d = nullptr;
    if (gpuRaster)
    {
        RectBeam *dBeams = nullptr;
        checkCuda(cudaMalloc(&d, count * sizeof(cufftDoubleComplex)), "raster grid malloc");
        checkCuda(cudaMalloc(&dBeams, beamSpecs.size() * sizeof(RectBeam)), "raster beams malloc");
        checkCuda(cudaMemcpy(dBeams, beamSpecs.data(), beamSpecs.size() * sizeof(RectBeam),
                             cudaMemcpyHostToDevice), "raster beams copy");
        int block = 256;
        int grid = (int)((count + block - 1) / block);
        rasterize_rect_beams_kernel<<<grid, block>>>(d, nx, ny, dx, dy,
                                                     dBeams, (int)beamSpecs.size());
        checkCuda(cudaGetLastError(), "raster kernel launch");
        checkCuda(cudaDeviceSynchronize(), "raster kernel sync");
        checkCuda(cudaFree(dBeams), "raster beams free");
    }
    else
    {
        for (int y = 0; y < ny; ++y)
        {
            double yy = (y + 0.5 - ny / 2.0) * dy;
            for (int x = 0; x < nx; ++x)
            {
                double xx = (x + 0.5 - nx / 2.0) * dx;
                size_t idx = (size_t)y * nx + x;
                cufftDoubleComplex value{0.0, 0.0};
                for (const RectBeam &beam : beamSpecs)
                {
                    double ca = std::cos(beam.angleRad), sa = std::sin(beam.angleRad);
                    double rx = xx - beam.centerX;
                    double ry = yy - beam.centerY;
                    double localX = rx * ca + ry * sa;
                    double localY = -rx * sa + ry * ca;
                    bool inside = (std::fabs(localX) < 0.5 * beam.width &&
                                   std::fabs(localY) < 0.5 * beam.height);
                    if (!inside) continue;
                    double phase = 2.0 * M_PI * (beam.phaseFx * xx + beam.phaseFy * yy);
                    cufftDoubleComplex carrier{std::cos(phase), std::sin(phase)};
                    cufftDoubleComplex amp{beam.ampR, beam.ampI};
                    cadd(value, cmul(amp, carrier));
                }
                h[idx] = value;
            }
        }
    }

    if (!d)
    {
        checkCuda(cudaMalloc(&d, count * sizeof(cufftDoubleComplex)), "cudaMalloc");
        checkCuda(cudaMemcpy(d, h.data(), count * sizeof(cufftDoubleComplex),
                             cudaMemcpyHostToDevice), "copy h2d");
    }

    cufftHandle plan;
    checkCufft(cufftPlan2d(&plan, ny, nx, CUFFT_Z2Z), "cufftPlan2d");
    checkCufft(cufftExecZ2Z(plan, d, d, CUFFT_FORWARD), "cufftExecZ2Z");
    checkCuda(cudaDeviceSynchronize(), "fft sync");
    if (angular && gpuSample)
    {
        const int angularCount = (nTheta + 1) * nPhi;
        double *dAngularAbs = nullptr;
        checkCuda(cudaMalloc(&dAngularAbs, angularCount * sizeof(double)),
                  "angular abs malloc");
        int block = 256;
        int grid = (angularCount + block - 1) / block;
        angular_sample_abs_kernel<<<grid, block>>>(d, dAngularAbs, nx, ny, dx, dy,
                                                   lambda, nPhi, nTheta,
                                                   std::max(1, interpRadius));
        checkCuda(cudaGetLastError(), "angular sample kernel launch");
        gpuAngularAbs.resize(angularCount);
        checkCuda(cudaMemcpy(gpuAngularAbs.data(), dAngularAbs,
                             angularCount * sizeof(double),
                             cudaMemcpyDeviceToHost), "angular abs copy d2h");
        checkCuda(cudaFree(dAngularAbs), "angular abs free");
    }
    if (needHostFft)
    {
        checkCuda(cudaMemcpy(h.data(), d, count * sizeof(cufftDoubleComplex),
                             cudaMemcpyDeviceToHost), "copy d2h");
    }
    checkCufft(cufftDestroy(plan), "cufftDestroy");
    checkCuda(cudaFree(d), "cudaFree");

    if (angular && noAnalytic)
    {
        std::printf("FFT_GPU_SAMPLE_SUMMARY nphi=%d ntheta=%d samples=%d gpu_sample=%s analytic=off\n",
                    nPhi, nTheta + 1, (nTheta + 1) * nPhi,
                    gpuSample ? "on" : "off");
        return 0;
    }

    double maxAbs = 0.0;
    double maxRel = 0.0;
    double maxRelNonzero = 0.0;
    double rms = 0.0;
    int n = 0;
    int nNonzero = 0;
    const int ky0 = 0;
    const int yidx = freq_index(ky0, ny);

    std::printf("# nx=%d ny=%d beams=%d raster=%s width=%.10g height=%.10g angle_deg=%.10g center=(%.10g,%.10g) phase=(%.10g,%.10g) dx=%.10g dy=%.10g\n",
                nx, ny, beams, gpuRaster ? "gpu" : "cpu", width, height, angleDeg, centerX, centerY, phaseFx, phaseFy, dx, dy);
    std::printf("# k fx fft_abs analytic_abs abs_err rel_err\n");

    for (int k = -samples; samples > 0 && k <= samples; ++k)
    {
        int xidx = freq_index(k, nx);
        double fx = (double)k / (nx * dx);
        const cufftDoubleComplex &z = h[(size_t)yidx * nx + xidx];
        double fftAbs = std::hypot(z.x, z.y) * dx * dy;
        double analyticAbs = analytic_rect_abs(fx, 0.0, beamSpecs);
        double absErr = std::fabs(fftAbs - analyticAbs);
        double relErr = (analyticAbs > 1e-12) ? absErr / analyticAbs : absErr;
        maxAbs = std::max(maxAbs, absErr);
        maxRel = std::max(maxRel, relErr);
        if (analyticAbs > 1e-6 * width * height)
        {
            maxRelNonzero = std::max(maxRelNonzero, relErr);
            ++nNonzero;
        }
        rms += absErr * absErr;
        ++n;
        if (std::abs(k) <= 8 || k == -samples || k == samples)
            std::printf("%d %.12g %.12g %.12g %.6g %.6g\n",
                        k, fx, fftAbs, analyticAbs, absErr, relErr);
    }

    rms = std::sqrt(rms / std::max(1, n));
    std::printf("FFT_PROBE_SUMMARY max_abs=%.10g max_rel=%.10g max_rel_nonzero=%.10g rms=%.10g samples=%d nonzero=%d\n",
                maxAbs, maxRel, maxRelNonzero, rms, n, nNonzero);

    if (angular)
    {
        double aMaxAbs = 0.0;
        double aMaxRel = 0.0;
        double aRms = 0.0;
        double iMaxAbs = 0.0;
        double iMaxRel = 0.0;
        double iRms = 0.0;
        double iSumFft = 0.0;
        double iSumAnalytic = 0.0;
        double iSumHybrid = 0.0;
        double iSumHybridIrel = 0.0;
        double wSumFft = 0.0;
        double wSumAnalytic = 0.0;
        double wSumHybrid = 0.0;
        double wSumHybridIrel = 0.0;
        struct WorstItem
        {
            double absErr;
            double relErr;
            int ip;
            int it;
            double fftAbs;
            double analyticAbs;
        };
        struct WorstCmp
        {
            bool operator()(const WorstItem &a, const WorstItem &b) const
            {
                return a.absErr > b.absErr;
            }
        };
        std::priority_queue<WorstItem, std::vector<WorstItem>, WorstCmp> worst;
        int aCount = 0;
        int aNonzero = 0;
        int aBad = 0;
        int aSignificant = 0;
        int aBadSignificant1 = 0;
        int aBadSignificant2 = 0;
        int hybridCorrections = 0;
        int hybridIrelCorrections = 0;
        for (int it = 0; it <= nTheta; ++it)
        {
            double theta = M_PI * it / nTheta;
            double st = std::sin(theta);
            double thetaWeight = st;
            for (int ip = 0; ip < nPhi; ++ip)
            {
                double phi = 2.0 * M_PI * ip / nPhi;
                double fx = st * std::cos(phi) / lambda;
                double fy = st * std::sin(phi) / lambda;
                double fftAbs;
                if (gpuSample)
                {
                    fftAbs = gpuAngularAbs[(size_t)it * nPhi + ip];
                }
                else
                {
                    cufftDoubleComplex z = (interpRadius > 1)
                        ? sample_fft_lanczos(h, nx, ny, dx, dy, fx, fy, interpRadius)
                        : sample_fft_bilinear(h, nx, ny, dx, dy, fx, fy);
                    fftAbs = std::hypot(z.x, z.y);
                }
                double analyticAbs = analytic_rect_abs(fx, fy, beamSpecs);
                double absErr = std::fabs(fftAbs - analyticAbs);
                double relErr = (analyticAbs > 1e-10) ? absErr / analyticAbs : absErr;
                double fftI = fftAbs * fftAbs;
                double analyticI = analyticAbs * analyticAbs;
                double iErr = std::fabs(fftI - analyticI);
                double iRel = (analyticI > 1e-10) ? iErr / analyticI : iErr;
                aMaxAbs = std::max(aMaxAbs, absErr);
                aMaxRel = std::max(aMaxRel, relErr);
                aRms += absErr * absErr;
                if (reportWorst > 0)
                {
                    WorstItem item{absErr, relErr, ip, it, fftAbs, analyticAbs};
                    if ((int)worst.size() < reportWorst)
                        worst.push(item);
                    else if (item.absErr > worst.top().absErr)
                    {
                        worst.pop();
                        worst.push(item);
                    }
                }
                iMaxAbs = std::max(iMaxAbs, iErr);
                iMaxRel = std::max(iMaxRel, iRel);
                iRms += iErr * iErr;
                iSumFft += fftI;
                iSumAnalytic += analyticI;
                bool significant = analyticAbs > 1e-3 * width * height;
                bool needsHybridCorrection = significant && relErr > hybridRel;
                bool needsHybridIrelCorrection = significant && iRel > hybridRel;
                iSumHybrid += needsHybridCorrection ? analyticI : fftI;
                iSumHybridIrel += needsHybridIrelCorrection ? analyticI : fftI;
                if (needsHybridCorrection)
                    ++hybridCorrections;
                if (needsHybridIrelCorrection)
                    ++hybridIrelCorrections;
                wSumFft += thetaWeight * fftI;
                wSumAnalytic += thetaWeight * analyticI;
                wSumHybrid += thetaWeight * (needsHybridCorrection ? analyticI : fftI);
                wSumHybridIrel += thetaWeight * (needsHybridIrelCorrection ? analyticI : fftI);
                ++aCount;
                if (analyticAbs > 1e-5 * width * height)
                {
                    ++aNonzero;
                    if (relErr > 0.02)
                        ++aBad;
                }
                if (significant)
                {
                    ++aSignificant;
                    if (relErr > 0.01)
                        ++aBadSignificant1;
                    if (relErr > 0.02)
                        ++aBadSignificant2;
                }
            }
        }
        aRms = std::sqrt(aRms / std::max(1, aCount));
        iRms = std::sqrt(iRms / std::max(1, aCount));
        double iSumRel = (iSumAnalytic > 1e-30)
            ? std::fabs(iSumFft - iSumAnalytic) / iSumAnalytic : 0.0;
        double iHybridRel = (iSumAnalytic > 1e-30)
            ? std::fabs(iSumHybrid - iSumAnalytic) / iSumAnalytic : 0.0;
        double iHybridIrelRel = (iSumAnalytic > 1e-30)
            ? std::fabs(iSumHybridIrel - iSumAnalytic) / iSumAnalytic : 0.0;
        double wSumRel = (wSumAnalytic > 1e-30)
            ? std::fabs(wSumFft - wSumAnalytic) / wSumAnalytic : 0.0;
        double wHybridRel = (wSumAnalytic > 1e-30)
            ? std::fabs(wSumHybrid - wSumAnalytic) / wSumAnalytic : 0.0;
        double wHybridIrelRel = (wSumAnalytic > 1e-30)
            ? std::fabs(wSumHybridIrel - wSumAnalytic) / wSumAnalytic : 0.0;
        std::printf("FFT_ANGULAR_SUMMARY nphi=%d ntheta=%d lambda=%.10g interp=%d max_abs=%.10g max_rel=%.10g rms=%.10g samples=%d nonzero=%d bad_rel_gt_2pct=%d\n",
                    nPhi, nTheta + 1, lambda, interpRadius, aMaxAbs, aMaxRel, aRms,
                    aCount, aNonzero, aBad);
        std::printf("FFT_INTENSITY_SUMMARY max_abs=%.10g max_rel=%.10g rms=%.10g sum_fft=%.10g sum_analytic=%.10g sum_rel=%.10g\n",
                    iMaxAbs, iMaxRel, iRms, iSumFft, iSumAnalytic, iSumRel);
        std::printf("FFT_SIGNIFICANT_SUMMARY threshold=%.10g count=%d bad_rel_gt_1pct=%d bad_rel_gt_2pct=%d\n",
                    1e-3 * width * height, aSignificant,
                    aBadSignificant1, aBadSignificant2);
        std::printf("FFT_HYBRID_SUMMARY rel_threshold=%.10g corrections=%d correction_fraction=%.10g sum_hybrid=%.10g sum_analytic=%.10g sum_rel=%.10g\n",
                    hybridRel, hybridCorrections,
                    (aCount > 0) ? (double)hybridCorrections / (double)aCount : 0.0,
                    iSumHybrid, iSumAnalytic, iHybridRel);
        std::printf("FFT_HYBRID_IREL_SUMMARY rel_threshold=%.10g corrections=%d correction_fraction=%.10g sum_hybrid=%.10g sum_analytic=%.10g sum_rel=%.10g\n",
                    hybridRel, hybridIrelCorrections,
                    (aCount > 0) ? (double)hybridIrelCorrections / (double)aCount : 0.0,
                    iSumHybridIrel, iSumAnalytic, iHybridIrelRel);
        std::printf("FFT_WEIGHTED_SUMMARY weight=sin_theta sum_fft=%.10g sum_analytic=%.10g sum_rel=%.10g hybrid_amp_rel=%.10g hybrid_irel_rel=%.10g\n",
                    wSumFft, wSumAnalytic, wSumRel, wHybridRel, wHybridIrelRel);
        if (reportWorst > 0)
        {
            std::vector<WorstItem> items;
            while (!worst.empty())
            {
                items.push_back(worst.top());
                worst.pop();
            }
            std::sort(items.begin(), items.end(),
                      [](const WorstItem &a, const WorstItem &b) {
                          return a.absErr > b.absErr;
                      });
            std::printf("# worst_abs_errors: rank theta_deg phi_deg fft_abs analytic_abs abs_err rel_err\n");
            for (size_t i = 0; i < items.size(); ++i)
            {
                double thetaDeg = 180.0 * items[i].it / nTheta;
                double phiDeg = 360.0 * items[i].ip / nPhi;
                std::printf("FFT_WORST %zu %.10g %.10g %.12g %.12g %.6g %.6g\n",
                            i + 1, thetaDeg, phiDeg,
                            items[i].fftAbs, items[i].analyticAbs,
                            items[i].absErr, items[i].relErr);
            }
        }
    }
    return 0;
}
