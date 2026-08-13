#include <cuda_runtime.h>

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <vector>

struct Quaternion
{
    double x, y, z, w;
};

struct Matrix3
{
    double v[9];
};

static bool Check(cudaError_t error, const char *where)
{
    if (error == cudaSuccess)
        return true;
    std::fprintf(stderr, "CUDA error at %s: %s\n", where,
                 cudaGetErrorString(error));
    return false;
}

static Matrix3 QuaternionMatrix(const Quaternion &q)
{
    const double xx = q.x*q.x, yy = q.y*q.y, zz = q.z*q.z;
    const double xy = q.x*q.y, xz = q.x*q.z, yz = q.y*q.z;
    const double wx = q.w*q.x, wy = q.w*q.y, wz = q.w*q.z;
    Matrix3 m = {{
        1.0 - 2.0*(yy + zz), 2.0*(xy - wz),       2.0*(xz + wy),
        2.0*(xy + wz),       1.0 - 2.0*(xx + zz), 2.0*(yz - wx),
        2.0*(xz - wy),       2.0*(yz + wx),       1.0 - 2.0*(xx + yy)
    }};
    return m;
}

__global__ void RotateQuaternionKernel(const double3 *input, const Quaternion *q,
                                       double3 *output, int vectorCount,
                                       int vectorsPerOrientation)
{
    const int index = (int)(blockIdx.x * blockDim.x + threadIdx.x);
    if (index >= vectorCount)
        return;
    const Quaternion rotation = q[index / vectorsPerOrientation];
    const double3 value = input[index];
    const double tx = 2.0 * (rotation.y * value.z - rotation.z * value.y);
    const double ty = 2.0 * (rotation.z * value.x - rotation.x * value.z);
    const double tz = 2.0 * (rotation.x * value.y - rotation.y * value.x);
    output[index] = make_double3(
        value.x + rotation.w * tx + rotation.y * tz - rotation.z * ty,
        value.y + rotation.w * ty + rotation.z * tx - rotation.x * tz,
        value.z + rotation.w * tz + rotation.x * ty - rotation.y * tx);
}

__global__ void RotateMatrixKernel(const double3 *input, const Matrix3 *matrix,
                                   double3 *output, int vectorCount,
                                   int vectorsPerOrientation)
{
    const int index = (int)(blockIdx.x * blockDim.x + threadIdx.x);
    if (index >= vectorCount)
        return;
    const Matrix3 m = matrix[index / vectorsPerOrientation];
    const double3 value = input[index];
    output[index] = make_double3(
        m.v[0]*value.x + m.v[1]*value.y + m.v[2]*value.z,
        m.v[3]*value.x + m.v[4]*value.y + m.v[5]*value.z,
        m.v[6]*value.x + m.v[7]*value.y + m.v[8]*value.z);
}

static float TimeKernel(bool quaternion, const double3 *input,
                        const Quaternion *q, const Matrix3 *matrix,
                        double3 *output, int vectorCount,
                        int vectorsPerOrientation, int repeats)
{
    const int block = 256;
    const int grid = (vectorCount + block - 1) / block;
    cudaEvent_t begin, end;
    cudaEventCreate(&begin);
    cudaEventCreate(&end);
    cudaEventRecord(begin);
    for (int repeat = 0; repeat < repeats; ++repeat)
    {
        if (quaternion)
            RotateQuaternionKernel<<<grid, block>>>(input, q, output, vectorCount,
                                                     vectorsPerOrientation);
        else
            RotateMatrixKernel<<<grid, block>>>(input, matrix, output, vectorCount,
                                                 vectorsPerOrientation);
    }
    cudaEventRecord(end);
    cudaEventSynchronize(end);
    float milliseconds = 0.0f;
    cudaEventElapsedTime(&milliseconds, begin, end);
    cudaEventDestroy(begin);
    cudaEventDestroy(end);
    return milliseconds;
}

int main(int argc, char **argv)
{
    const int orientationCount = argc > 1 ? std::atoi(argv[1]) : 65536;
    const int vectorsPerOrientation = argc > 2 ? std::atoi(argv[2]) : 64;
    const int repeats = argc > 3 ? std::atoi(argv[3]) : 50;
    if (orientationCount <= 0 || vectorsPerOrientation <= 0 || repeats <= 0)
    {
        std::fprintf(stderr, "usage: %s [orientations] [vectors-per-orientation] [repeats]\n",
                     argv[0]);
        return 2;
    }

    const int vectorCount = orientationCount * vectorsPerOrientation;
    std::vector<Quaternion> quaternions(orientationCount);
    std::vector<Matrix3> matrices(orientationCount);
    std::vector<double3> input(vectorCount), quaternionOutput(vectorCount),
        matrixOutput(vectorCount);
    for (int i = 0; i < orientationCount; ++i)
    {
        const double u1 = (i + 0.5) / orientationCount;
        const double u2 = std::fmod(i * 0.6180339887498949, 1.0);
        const double u3 = std::fmod(i * 0.7548776662466927, 1.0);
        const double a = 2.0 * M_PI * u2;
        const double b = 2.0 * M_PI * u3;
        quaternions[i] = Quaternion{
            std::sqrt(1.0 - u1) * std::sin(a),
            std::sqrt(1.0 - u1) * std::cos(a),
            std::sqrt(u1) * std::sin(b),
            std::sqrt(u1) * std::cos(b)};
        matrices[i] = QuaternionMatrix(quaternions[i]);
    }
    for (int i = 0; i < vectorCount; ++i)
    {
        const double a = 0.001 * (i % 1009);
        input[i] = make_double3(std::sin(a), std::cos(0.7*a),
                                std::sin(1.3*a));
    }

    double3 *deviceInput = nullptr, *deviceOutput = nullptr;
    Quaternion *deviceQuaternion = nullptr;
    Matrix3 *deviceMatrix = nullptr;
    const size_t vectorBytes = (size_t)vectorCount * sizeof(double3);
    if (!Check(cudaMalloc(&deviceInput, vectorBytes), "cudaMalloc input")
        || !Check(cudaMalloc(&deviceOutput, vectorBytes), "cudaMalloc output")
        || !Check(cudaMalloc(&deviceQuaternion,
                             quaternions.size() * sizeof(Quaternion)), "cudaMalloc quaternion")
        || !Check(cudaMalloc(&deviceMatrix,
                             matrices.size() * sizeof(Matrix3)), "cudaMalloc matrix"))
        return 1;
    Check(cudaMemcpy(deviceInput, input.data(), vectorBytes,
                     cudaMemcpyHostToDevice), "copy input");
    Check(cudaMemcpy(deviceQuaternion, quaternions.data(),
                     quaternions.size() * sizeof(Quaternion),
                     cudaMemcpyHostToDevice), "copy quaternion");
    Check(cudaMemcpy(deviceMatrix, matrices.data(), matrices.size() * sizeof(Matrix3),
                     cudaMemcpyHostToDevice), "copy matrix");

    TimeKernel(true, deviceInput, deviceQuaternion, deviceMatrix, deviceOutput,
               vectorCount, vectorsPerOrientation, 1);
    TimeKernel(false, deviceInput, deviceQuaternion, deviceMatrix, deviceOutput,
               vectorCount, vectorsPerOrientation, 1);
    const float quaternionMs = TimeKernel(
        true, deviceInput, deviceQuaternion, deviceMatrix, deviceOutput,
        vectorCount, vectorsPerOrientation, repeats);
    const auto transferBegin = std::chrono::high_resolution_clock::now();
    Check(cudaMemcpy(quaternionOutput.data(), deviceOutput, vectorBytes,
                     cudaMemcpyDeviceToHost), "copy quaternion output");
    const double deviceToHostMs = std::chrono::duration<double, std::milli>(
        std::chrono::high_resolution_clock::now() - transferBegin).count();
    const float matrixMs = TimeKernel(
        false, deviceInput, deviceQuaternion, deviceMatrix, deviceOutput,
        vectorCount, vectorsPerOrientation, repeats);
    Check(cudaMemcpy(matrixOutput.data(), deviceOutput, vectorBytes,
                     cudaMemcpyDeviceToHost), "copy matrix output");

    double maxError = 0.0;
    for (int i = 0; i < vectorCount; ++i)
    {
        maxError = std::max(maxError, std::fabs(
            quaternionOutput[i].x - matrixOutput[i].x));
        maxError = std::max(maxError, std::fabs(
            quaternionOutput[i].y - matrixOutput[i].y));
        maxError = std::max(maxError, std::fabs(
            quaternionOutput[i].z - matrixOutput[i].z));
    }

    const double transforms = (double)vectorCount * repeats;
    std::printf("precision=double orientations=%d vectors_per_orientation=%d repeats=%d\n",
                orientationCount, vectorsPerOrientation, repeats);
    std::printf("quaternion_ms=%.3f matrix_ms=%.3f quaternion_speedup_vs_matrix=%.3f "
                "quaternion_Gtransform_s=%.3f matrix_Gtransform_s=%.3f max_error=%.3g\n",
                quaternionMs, matrixMs, matrixMs / quaternionMs,
                transforms / (quaternionMs * 1.0e6),
                transforms / (matrixMs * 1.0e6), maxError);
    std::printf("single_quaternion_kernel_ms=%.3f device_to_host_ms=%.3f "
                "transfer_over_kernel=%.1f\n",
                quaternionMs / repeats, deviceToHostMs,
                deviceToHostMs / (quaternionMs / repeats));

    cudaFree(deviceInput);
    cudaFree(deviceOutput);
    cudaFree(deviceQuaternion);
    cudaFree(deviceMatrix);
    return Check(cudaGetLastError(), "final CUDA status") ? 0 : 1;
}
