#include "../src/cuda/GpuPhaseReduction.cuh"
#include <cuda_runtime.h>
#include <algorithm>
#include <cstdio>
#include <random>
#include <vector>

__global__ void evaluate(const double* phases, float2* values, int n)
{
    const int i = blockIdx.x*blockDim.x+threadIdx.x;
    if (i>=n) return;
    double r;
    if (MbsGpuPhase::Reduce(phases[i],r))
        sincosf((float)r,&values[i].x,&values[i].y);
    else {
        double s,c;
        sincos(phases[i],&s,&c);
        values[i]=make_float2((float)s,(float)c);
    }
}

int main()
{
    std::vector<double> x(1<<20);
    std::mt19937_64 rng(30211);
    std::uniform_real_distribution<double> draw(-1048576.,1048576.);
    for (auto& v:x) v=draw(rng);
    for (int i=0;i<4096;++i)
        x[i]=(double)(i-2048)*acos(-1.)+((i%3)-1)*1e-12;
    x[4096]=1048576.;x[4097]=-1048576.;
    x[4098]=1e10;x[4099]=-1e10;x[4100]=1e100;x[4101]=-1e100;
    std::vector<float2> y(x.size());
    double* dx=nullptr;float2* dy=nullptr;
    if (cudaMalloc(&dx,x.size()*sizeof(double))!=cudaSuccess ||
        cudaMalloc(&dy,y.size()*sizeof(float2))!=cudaSuccess) return 1;
    if (cudaMemcpy(dx,x.data(),x.size()*sizeof(double),cudaMemcpyHostToDevice)!=cudaSuccess) return 2;
    evaluate<<<(x.size()+255)/256,256>>>(dx,dy,x.size());
    if (cudaGetLastError()!=cudaSuccess ||
        cudaMemcpy(y.data(),dy,y.size()*sizeof(float2),cudaMemcpyDeviceToHost)!=cudaSuccess) return 3;
    double worst=0;
    for (size_t i=0;i<x.size();++i) {
        const double error=(double)hypotl((long double)y[i].x-sinl((long double)x[i]),
                                         (long double)y[i].y-cosl((long double)x[i]));
        if (!std::isfinite(error)) return 4;
        worst=std::max(worst,error);
    }
    cudaFree(dx);cudaFree(dy);
    std::printf("CUDA reduced phase, including FP64 fallback: max complex error %.9g\n",worst);
    return worst<2e-7 ? 0 : 5;
}
