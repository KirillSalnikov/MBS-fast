// Build with GpuAnalyticFacet.cu and MBS_ANALYTIC_GPU_TEST_TILE_CELLS=31.
// Tests identical optical kernels/reductions and asynchronous workspace ownership.
#include "GpuAnalyticFacet.h"
#include <cuda_runtime.h>
#include <cmath>
#include <cstdio>
#include <cstring>
#include <stdexcept>
#include <thread>

static void require(bool ok,const std::string&message)
{if(!ok)throw std::runtime_error(message);}
static void same(const std::vector<double>&a,const std::vector<double>&b)
{require(a.size()==b.size()&&!std::memcmp(a.data(),b.data(),a.size()*sizeof(double)),"controls differ bitwise");}

int main(){try{
    AnalyticGpuModel model;model.indexReal=1.3116;model.indexImag=0;
    model.wave=.532;model.circleRadius=3.2;model.shadowMode=1;
    for(int axis=0;axis<3;++axis)for(int sign:{-1,1}){
        AnalyticGpuPatch patch={};patch.normal[axis]=sign;
        patch.horizontal[(axis+1)%3]=1;patch.vertical[(axis+2)%3]=sign;
        patch.firstVertex=model.verticesXY.size()/2;patch.vertices=4;patch.area=25;
        for(double xy:{-2.5,-2.5,2.5,-2.5,2.5,2.5,-2.5,2.5})model.verticesXY.push_back(xy);
        model.patches.push_back(patch);
    }
    const std::vector<double>theta={0.,1e-8,.4,2.7,std::acos(-1.)-1e-8,std::acos(-1.)};
    std::vector<std::array<AnalyticBackscatter::Vec,3>>frames(32);
    std::vector<double>weights(frames.size());
    for(size_t i=0;i<frames.size();++i){
        const double a=i*.039,b=i*.023,ca=cos(a),sa=sin(a),cb=cos(b),sb=sin(b);
        frames[i]={{{ca*cb,sa*cb,-sb},{-sa,ca,0},{ca*sb,sa*sb,cb}}};weights[i]=(i%7+1.)/frames.size();
    }
    std::string error;size_t comparisons=0;
    for(int shadowMode:{0,1,2})for(double absorption:{0.,.00042407}){
        model.shadowMode=shadowMode;model.indexImag=absorption;
        std::vector<double>referenceR,referenceS,r(1,17.),s(1,19.);
        require(EvaluateAnalyticFacetGpu(model,frames,weights,theta,4,referenceR,referenceS,nullptr,&error),error);
        AnalyticFacetGpuPending pending;bool supported=false;
        require(pending.Prepare(model,frames,weights,theta,4,supported,&error)&&supported,error);
        require(!pending.Collect(r,s,&error),"collect before submit accepted");
        require(r==std::vector<double>(1,17.)&&s==std::vector<double>(1,19.),"failed collect published output");
        std::vector<double>blockedR(1,23.),blockedS(1,29.);
        require(!EvaluateAnalyticFacetGpu(model,frames,weights,theta,4,blockedR,blockedS,nullptr,&error),"pending workspace reused");
        require(blockedR==std::vector<double>(1,23.)&&blockedS==std::vector<double>(1,29.),"busy workspace changed output");
        bool wrongThreadAccepted=true;
        std::thread wrongThread([&](){std::string why;wrongThreadAccepted=pending.Submit(&why);});wrongThread.join();
        require(!wrongThreadAccepted,"submit on another host thread accepted");
        require(pending.Submit(&error),error);
        require(!pending.Submit(&error),"duplicate submit accepted");
        require(pending.Collect(r,s,&error),error);same(r,referenceR);same(s,referenceS);
        require(!pending.Collect(r,s,&error),"duplicate collect accepted");
        same(r,referenceR);same(s,referenceS);++comparisons;
        // Exceptional scope exit must drain and unlock the same workspace.
        try{
            AnalyticFacetGpuPending abandoned;
            require(abandoned.Prepare(model,frames,weights,theta,4,supported,&error)&&supported,error);
            require(abandoned.Submit(&error),error);
            throw 7;
        }catch(int value){require(value==7,"unexpected exception");}
        require(EvaluateAnalyticFacetGpu(model,frames,weights,theta,4,r,s,nullptr,&error),error);
        same(r,referenceR);same(s,referenceS);
        // Preparation without submission also releases the reservation.
        {AnalyticFacetGpuPending unused;require(unused.Prepare(model,frames,weights,theta,4,supported,&error)&&supported,error);}
    }
#ifdef MBS_ANALYTIC_GPU_TEST_TILE_CELLS
    std::vector<double>largeTheta(20);for(size_t i=0;i<largeTheta.size();++i)largeTheta[i]=std::acos(-1.)*i/(largeTheta.size()-1);
    std::vector<double>r,s,referenceR,referenceS;
    require(EvaluateAnalyticFacetGpu(model,frames,weights,largeTheta,4,referenceR,referenceS,nullptr,&error),error);
    {AnalyticFacetGpuPending tiled;bool supported=true;
        require(tiled.Prepare(model,frames,weights,largeTheta,4,supported,&error)&&!supported,"multi-tile case was submitted");
        require(!tiled.Submit(&error),"unsupported operation accepted submit");
        // Unsupported preflight holds no workspace and permits the old path.
        require(EvaluateAnalyticFacetGpu(model,frames,weights,largeTheta,4,r,s,nullptr,&error),error);
        same(r,referenceR);same(s,referenceS);}
#endif
    std::printf("PASS pending analytic controls: %zu bitwise optical comparisons; ownership, abort and tiled fallback\n",comparisons);
    return 0;
}catch(const std::exception&error){std::fprintf(stderr,"FAIL: %s\n",error.what());return 1;}}
