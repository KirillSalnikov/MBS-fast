#include <cuda_runtime.h>

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <vector>

struct Vec4
{
    float x, y, z, d;
};

static inline float dot3(const Vec4 &a, const Vec4 &b)
{
    return a.x * b.x + a.y * b.y + a.z * b.z;
}

static void cpu_project(const Vec4 *points, const Vec4 *dirs,
                        const Vec4 *normals, Vec4 *out, int count)
{
    for (int idx = 0; idx < count; ++idx)
    {
        const Vec4 p = points[idx];
        const Vec4 dir = dirs[idx];
        const Vec4 n = normals[idx];
        float dp0 = dot3(dir, n);
        float t = (dot3(p, n) + n.d) / dp0;
        out[idx] = Vec4{p.x - dir.x * t, p.y - dir.y * t,
                        p.z - dir.z * t, 0.0f};
    }
}

__global__ void project_kernel(const Vec4 *points, const Vec4 *dirs,
                               const Vec4 *normals, Vec4 *out, int count)
{
    int idx = (int)(blockIdx.x * blockDim.x + threadIdx.x);
    if (idx >= count) return;
    Vec4 p = points[idx];
    Vec4 dir = dirs[idx];
    Vec4 n = normals[idx];
    float dp0 = dir.x * n.x + dir.y * n.y + dir.z * n.z;
    float t = (p.x * n.x + p.y * n.y + p.z * n.z + n.d) / dp0;
    out[idx] = Vec4{p.x - dir.x * t, p.y - dir.y * t,
                    p.z - dir.z * t, 0.0f};
}

static bool check(cudaError_t err, const char *where)
{
    if (err == cudaSuccess) return true;
    std::fprintf(stderr, "CUDA error at %s: %s\n", where, cudaGetErrorString(err));
    return false;
}

int main(int argc, char **argv)
{
    int count = argc > 1 ? std::atoi(argv[1]) : 1000000;
    int repeats = argc > 2 ? std::atoi(argv[2]) : 20;
    if (count <= 0 || repeats <= 0)
    {
        std::fprintf(stderr, "usage: %s [count] [repeats]\n", argv[0]);
        return 2;
    }

    std::vector<Vec4> points(count), dirs(count), normals(count);
    std::vector<Vec4> cpuOut(count), gpuOut(count);
    for (int i = 0; i < count; ++i)
    {
        float a = 0.001f * (float)(i % 1009);
        float b = 0.001f * (float)(i % 917);
        points[i] = Vec4{std::sin(a), std::cos(b), std::sin(a + b), 0.0f};
        dirs[i] = Vec4{0.2f + 0.1f * std::sin(b), -0.3f,
                       -0.9f + 0.05f * std::cos(a), 0.0f};
        normals[i] = Vec4{0.1f, 0.2f + 0.01f * std::sin(a),
                          0.97f, -0.05f * std::cos(b)};
    }

    auto cpuStart = std::chrono::high_resolution_clock::now();
    for (int r = 0; r < repeats; ++r)
        cpu_project(points.data(), dirs.data(), normals.data(), cpuOut.data(), count);
    auto cpuEnd = std::chrono::high_resolution_clock::now();

    Vec4 *dPoints = nullptr, *dDirs = nullptr, *dNormals = nullptr, *dOut = nullptr;
    size_t bytes = (size_t)count * sizeof(Vec4);
    if (!check(cudaMalloc(&dPoints, bytes), "cudaMalloc points") ||
        !check(cudaMalloc(&dDirs, bytes), "cudaMalloc dirs") ||
        !check(cudaMalloc(&dNormals, bytes), "cudaMalloc normals") ||
        !check(cudaMalloc(&dOut, bytes), "cudaMalloc out"))
        return 1;

    cudaEvent_t e0, e1, k0, k1;
    cudaEventCreate(&e0);
    cudaEventCreate(&e1);
    cudaEventCreate(&k0);
    cudaEventCreate(&k1);
    cudaEventRecord(e0);
    cudaMemcpy(dPoints, points.data(), bytes, cudaMemcpyHostToDevice);
    cudaMemcpy(dDirs, dirs.data(), bytes, cudaMemcpyHostToDevice);
    cudaMemcpy(dNormals, normals.data(), bytes, cudaMemcpyHostToDevice);
    int block = 256;
    int grid = (count + block - 1) / block;
    cudaEventRecord(k0);
    for (int r = 0; r < repeats; ++r)
        project_kernel<<<grid, block>>>(dPoints, dDirs, dNormals, dOut, count);
    cudaEventRecord(k1);
    cudaMemcpy(gpuOut.data(), dOut, bytes, cudaMemcpyDeviceToHost);
    cudaEventRecord(e1);
    cudaEventSynchronize(e1);
    cudaEventSynchronize(k1);
    float gpuMs = 0.0f;
    float kernelMs = 0.0f;
    cudaEventElapsedTime(&gpuMs, e0, e1);
    cudaEventElapsedTime(&kernelMs, k0, k1);

    double maxErr = 0.0;
    for (int i = 0; i < std::min(count, 1000); ++i)
    {
        maxErr = std::max(maxErr, std::fabs((double)cpuOut[i].x - gpuOut[i].x));
        maxErr = std::max(maxErr, std::fabs((double)cpuOut[i].y - gpuOut[i].y));
        maxErr = std::max(maxErr, std::fabs((double)cpuOut[i].z - gpuOut[i].z));
    }

    double cpuMs = std::chrono::duration<double, std::milli>(cpuEnd - cpuStart).count();
    double ops = (double)count * repeats;
    std::printf("count=%d repeats=%d ops=%.0f\n", count, repeats, ops);
    std::printf("cpu_ms=%.3f gpu_ms_including_copies=%.3f kernel_ms=%.3f "
                "speedup=%.2f kernel_speedup=%.2f max_err=%.3g\n",
                cpuMs, gpuMs, kernelMs, cpuMs / gpuMs, cpuMs / kernelMs, maxErr);

    cudaFree(dPoints);
    cudaFree(dDirs);
    cudaFree(dNormals);
    cudaFree(dOut);
    return 0;
}
