#include "GpuAnalyticFacet.h"
#include "AnalyticRadialTail.h"
#include <cuda_runtime.h>
#include <map>
#include <memory>
#include <algorithm>

namespace {
constexpr size_t tileSize=65536;
struct Workspace {
    int device;
    cudaStream_t stream=nullptr;
    double *factor=nullptr,*values=nullptr,*derivatives=nullptr;
    explicit Workspace(int d):device(d) {}
    ~Workspace(){cudaSetDevice(device);if(stream)cudaStreamDestroy(stream);
        if(factor)cudaFree(factor);if(values)cudaFree(values);if(derivatives)cudaFree(derivatives);}
    cudaError_t Init(){
        cudaError_t s=cudaStreamCreateWithFlags(&stream,cudaStreamNonBlocking);if(s!=cudaSuccess)return s;
        s=cudaMalloc(reinterpret_cast<void**>(&factor),161*sizeof(double));if(s!=cudaSuccess)return s;
        s=cudaMalloc(reinterpret_cast<void**>(&values),tileSize*sizeof(double));if(s!=cudaSuccess)return s;
        return cudaMalloc(reinterpret_cast<void**>(&derivatives),tileSize*sizeof(double));
    }
};
__global__ void radial_tail(int degree,double step,size_t first,size_t count,
    const double*factor,double z0,double offset,double fluxOffset,double*values,double*derivatives){
    const size_t i=size_t(blockIdx.x)*blockDim.x+threadIdx.x;if(i>=count)return;
    const double z=(first+i)*step;double v,d;
    AnalyticRadialTail::Exact(degree,z,factor,v,d);
    values[i]=v+offset+fluxOffset*log(z/z0);derivatives[i]=d+fluxOffset/z;
}
}

int AnalyticMeanGpuDevice(){int device=-1;return cudaGetDevice(&device)==cudaSuccess?device:-1;}

bool EvaluateAnalyticMeanTailGpu(int device,int degree,double step,size_t first,size_t end,
    const std::vector<double>&factor,double z0,double offset,double fluxOffset,
    double*values,double*derivatives,std::string*error){
    const auto fail=[&](const char*operation,cudaError_t s){if(error)*error=std::string(operation)+": "+cudaGetErrorString(s);return false;};
    if(degree<2||degree>320||degree%2||factor.size()!=size_t(degree/2+1)
        ||!(step>0)||!(z0>0)||first>=end||first*step<=degree+10||!values||!derivatives){
        if(error)*error="invalid radial tail domain";return false;
    }
    cudaError_t s=cudaSetDevice(device);if(s!=cudaSuccess)return fail("cudaSetDevice",s);
    static thread_local std::map<int,std::unique_ptr<Workspace>> workspaces;
    auto& saved=workspaces[device];
    if(!saved){std::unique_ptr<Workspace> candidate(new Workspace(device));s=candidate->Init();
        if(s!=cudaSuccess)return fail("bounded radial workspace",s);saved=std::move(candidate);}
    auto& w=*saved;
    s=cudaMemcpyAsync(w.factor,factor.data(),factor.size()*sizeof(double),cudaMemcpyHostToDevice,w.stream);
    if(s!=cudaSuccess)return fail("radial coefficients upload",s);
    for(size_t begin=first;begin<end;){
        const size_t n=std::min(tileSize,end-begin);
        radial_tail<<<(n+127)/128,128,0,w.stream>>>(degree,step,begin,n,w.factor,z0,offset,fluxOffset,w.values,w.derivatives);
        s=cudaGetLastError();if(s!=cudaSuccess){cudaStreamSynchronize(w.stream);return fail("radial kernel",s);}
        s=cudaMemcpyAsync(values+begin,w.values,n*sizeof(double),cudaMemcpyDeviceToHost,w.stream);
        if(s==cudaSuccess)s=cudaMemcpyAsync(derivatives+begin,w.derivatives,n*sizeof(double),cudaMemcpyDeviceToHost,w.stream);
        if(s!=cudaSuccess){cudaStreamSynchronize(w.stream);return fail("radial download",s);}
        s=cudaStreamSynchronize(w.stream);if(s!=cudaSuccess)return fail("radial stream",s);
        begin+=n;
    }
    return true;
}
