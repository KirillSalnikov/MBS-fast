// Dedicated small regression. Compile GpuAnalyticFacet.cu with
// -DMBS_ANALYTIC_GPU_TEST_TILE_CELLS=31 to exercise partial phi/theta tiles.
// The reference switch only moves the final phi reduction back onto the CPU;
// both paths execute identical physical kernels and orientation reductions.
#include "GpuAnalyticFacet.h"
#include <cuda_runtime.h>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <stdexcept>
static void require(bool ok,const std::string&message){if(!ok)throw std::runtime_error(message);}
int main(){try{
    AnalyticGpuModel model;model.indexReal=1.3116;model.indexImag=0;model.wave=.532;model.circleRadius=3.2;model.shadowMode=1;
    for(int axis=0;axis<3;++axis)for(int sign:{-1,1}){
        AnalyticGpuPatch patch={};patch.normal[axis]=sign;patch.horizontal[(axis+1)%3]=1;patch.vertical[(axis+2)%3]=sign;
        patch.firstVertex=model.verticesXY.size()/2;patch.vertices=4;patch.area=25;
        for(double xy:{-2.5,-2.5,2.5,-2.5,2.5,2.5,-2.5,2.5})model.verticesXY.push_back(xy);
        model.patches.push_back(patch);
    }
    const double pi=std::acos(-1.);
    const std::vector<double>theta={0,1e-8,.1,.4,1.,1.4,2.,2.61,2.7,2.9,pi-1e-8,pi};
    size_t comparisons=0;
    for(int count:{65,9,65})for(int phi:{1,4,7,32})for(int shadowMode:{0,1,2})for(double absorption:{0.,.00042407}){
        model.shadowMode=shadowMode;model.indexImag=absorption;
        std::vector<std::array<AnalyticBackscatter::Vec,3>>frames(count);
        std::vector<double>weights(count);
        for(int i=0;i<count;++i){const double a=i*.039,b=i*.023,ca=cos(a),sa=sin(a),cb=cos(b),sb=sin(b);
            frames[i]={{{ca*cb,sa*cb,-sb},{-sa,ca,0},{ca*sb,sa*sb,cb}}};weights[i]=(i%7+1.)/count;}
        std::vector<double>r,s,oldR,oldS;
        std::vector<AnalyticFacetAverage::Components>samples,oldSamples;
        std::string error;
        require(setenv("MBS_CUDA_ANALYTIC_CPU_SAMPLES","1",1)==0,"setenv failed");
        require(EvaluateAnalyticFacetGpu(model,frames,weights,theta,phi,oldR,oldS,&oldSamples,&error),error);
        require(unsetenv("MBS_CUDA_ANALYTIC_CPU_SAMPLES")==0,"unsetenv failed");
        require(EvaluateAnalyticFacetGpu(model,frames,weights,theta,phi,r,s,&samples,&error),error);
        require(r.size()==oldR.size()&&s.size()==oldS.size()&&samples.size()==oldSamples.size(),"array sizes differ");
        require(!std::memcmp(r.data(),oldR.data(),r.size()*sizeof(double)),"reflection differs bitwise");
        require(!std::memcmp(s.data(),oldS.data(),s.size()*sizeof(double)),"shadow differs bitwise");
        require(!std::memcmp(samples.data(),oldSamples.data(),samples.size()*sizeof(samples[0])),"phi-reduced samples differ bitwise");
        ++comparisons;
    }
    std::printf("PASS analytic GPU sequential phi reduction: %zu bitwise comparisons\n",comparisons);
    return 0;
}catch(const std::exception&error){std::fprintf(stderr,"FAIL: %s\n",error.what());return 1;}}
