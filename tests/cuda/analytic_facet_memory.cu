// GPU regression for angular tiling. Optionally link an unchanged pre-fix
// source with its public function names suffixed Legacy for bitwise checks.
#include "GpuAnalyticFacet.h"
#include <cuda_runtime.h>
#include <cmath>
#include <cstdio>
#include <cstring>
#include <stdexcept>
#include <chrono>
#ifdef MBS_COMPARE_LEGACY
bool EvaluateAnalyticFacetGpuLegacy(const AnalyticGpuModel&,
    const std::vector<std::array<AnalyticBackscatter::Vec,3>>&,
    const std::vector<double>&,const std::vector<double>&,int,
    std::vector<double>&,std::vector<double>&,
    std::vector<AnalyticFacetAverage::Components>*);
#endif
static void require(bool condition,const std::string&message){if(!condition)throw std::runtime_error(message);}
int main(int argc,char**argv){try{
    const std::string mode=argc>1?argv[1]:"compare";
    const bool large=mode=="large";
    const int count=large?4096:1024,rows=large?5470:257,phi=large?64:33;
    const double pi=acos(-1.);
    AnalyticGpuModel model;model.indexReal=1.3116;model.indexImag=0;model.wave=.532;model.circleRadius=73.55;model.shadowMode=1;
    // Six rectangular facets, including grazing directions and both normals.
    for(int axis=0;axis<3;++axis)for(int sign:{-1,1}){
        AnalyticGpuPatch p={};p.normal[axis]=sign;p.horizontal[(axis+1)%3]=1;p.vertical[(axis+2)%3]=sign;
        p.firstVertex=model.verticesXY.size()/2;p.vertices=4;p.area=147.1*446.7;
        for(double v:{-73.55,-223.35,73.55,-223.35,73.55,223.35,-73.55,223.35})model.verticesXY.push_back(v);
        model.patches.push_back(p);
    }
    std::vector<std::array<AnalyticBackscatter::Vec,3>>frames(count);
    std::vector<double>weights(count),theta(rows);
    for(int i=0;i<count;++i){const double a=i*.039,b=i*.023,ca=cos(a),sa=sin(a),cb=cos(b),sb=sin(b);
        frames[i]={{{ca*cb,sa*cb,-sb},{-sa,ca,0},{ca*sb,sa*sb,cb}}};weights[i]=(i%7+1.)/count;}
    for(int j=0;j<rows;++j)theta[j]=pi*j/(rows-1);
    std::vector<double>r,s,oldR,oldS;
    std::vector<AnalyticFacetAverage::Components>samples,oldSamples;
    std::string error;
    if(large){
        // Reproduce the failed fullauto grid without spending time on facet
        // trigonometry. Nonzero refractive-index cases are checked below.
        model.indexReal=1;model.indexImag=0;
#ifdef MBS_COMPARE_LEGACY
        require(!EvaluateAnalyticFacetGpuLegacy(model,frames,weights,theta,phi,oldR,oldS,nullptr),"expected legacy allocation failure on 16 GB GPU");
        cudaGetLastError();
        std::puts("legacy failed as expected at 4096 x 5470 x 64");
#endif
        require(EvaluateAnalyticFacetGpu(model,frames,weights,theta,phi,r,s,nullptr,&error),error);
        for(size_t i=0;i<r.size();++i)require(r[i]==0&&s[i]==0,"vacuum result changed");
    }else{
        for(int shadowMode=0;shadowMode<=2;++shadowMode)for(double absorption:{0.,.00042407}){
            model.shadowMode=shadowMode;model.indexImag=absorption;
#ifdef MBS_COMPARE_LEGACY
            require(EvaluateAnalyticFacetGpuLegacy(model,frames,weights,theta,phi,oldR,oldS,&oldSamples),"legacy evaluation failed");
#endif
            require(EvaluateAnalyticFacetGpu(model,frames,weights,theta,phi,r,s,&samples,&error),error);
#ifdef MBS_COMPARE_LEGACY
            require(r.size()==oldR.size()&&s.size()==oldS.size()&&samples.size()==oldSamples.size(),"output dimensions changed");
            require(!memcmp(r.data(),oldR.data(),r.size()*sizeof(double)),"reflection differs bitwise");
            require(!memcmp(s.data(),oldS.data(),s.size()*sizeof(double)),"shadow differs bitwise");
            require(!memcmp(samples.data(),oldSamples.data(),samples.size()*sizeof(samples[0])),"samples differ bitwise");
#endif
            std::printf("matched shadow=%d absorption=%.8g cells=%zu samples=%zu\n",shadowMode,absorption,r.size(),samples.size());
        }
        if(mode=="pressure"){
            oldR=r;oldS=s;oldSamples=samples;
            size_t free=0,total=0;require(cudaMemGetInfo(&free,&total)==cudaSuccess,"meminfo failed");
            void*pressure=nullptr;require(free>64*1024*1024,"need 64 MiB for pressure fixture");
            require(cudaMalloc(&pressure,free-32*1024*1024)==cudaSuccess,"pressure allocation failed");
            const bool ok=EvaluateAnalyticFacetGpu(model,frames,weights,theta,phi,r,s,&samples,&error);
            cudaFree(pressure);require(ok,error);
            require(!memcmp(r.data(),oldR.data(),r.size()*sizeof(double))&&!memcmp(s.data(),oldS.data(),s.size()*sizeof(double))
                &&!memcmp(samples.data(),oldSamples.data(),samples.size()*sizeof(samples[0])),"memory pressure changed results");
            std::puts("32 MiB free: results unchanged");
        }
        frames[0][0].x=2;
        require(!EvaluateAnalyticFacetGpu(model,frames,weights,theta,phi,r,s,nullptr,&error)&&!error.empty(),"invalid frame lacks diagnostic");
    }
    std::puts("PASS analytic CUDA bounded workspace");return 0;
}catch(const std::exception&e){std::fprintf(stderr,"FAIL: %s\n",e.what());return 1;}}
