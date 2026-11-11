// Link the prior source with EvaluateAnalyticFacetGpu renamed to ...Baseline.
// This real hexagonal-column aperture fixture compares identical scientific
// work while orientation chunk sizes grow and shrink. No artificial GPU work.
#include "GpuAnalyticFacet.h"
#include <cuda_runtime.h>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <chrono>
#include <stdexcept>
bool EvaluateAnalyticFacetGpuBaseline(const AnalyticGpuModel&,
    const std::vector<std::array<AnalyticBackscatter::Vec,3>>&,
    const std::vector<double>&,const std::vector<double>&,int,
    std::vector<double>&,std::vector<double>&,
    std::vector<AnalyticFacetAverage::Components>*,std::string*);
static void require(bool ok,const std::string&message){if(!ok)throw std::runtime_error(message);}
int main(){try{
    const double pi=std::acos(-1.),radius=69.6/2,length=100;
    AnalyticGpuModel model;model.indexReal=1.3116;model.indexImag=0;model.wave=.532;model.circleRadius=radius;model.shadowMode=1;
    for(int sign:{-1,1}){
        AnalyticGpuPatch patch={};patch.normal[2]=sign;patch.horizontal[0]=1;patch.vertical[1]=sign;
        patch.firstVertex=model.verticesXY.size()/2;patch.vertices=6;patch.area=3*std::sqrt(3.)*radius*radius/2;
        for(int j=0;j<6;++j){model.verticesXY.push_back(radius*std::cos(j*pi/3));model.verticesXY.push_back(radius*std::sin(j*pi/3));}
        model.patches.push_back(patch);
    }
    for(int j=0;j<6;++j){
        AnalyticGpuPatch patch={};const double angle=(j+.5)*pi/3,c=std::cos(angle),s=std::sin(angle);
        patch.normal[0]=c;patch.normal[1]=s;patch.horizontal[0]=-s;patch.horizontal[1]=c;patch.vertical[2]=1;
        patch.firstVertex=model.verticesXY.size()/2;patch.vertices=4;patch.area=length*radius;
        for(double xy:{-radius/2,-length/2,radius/2,-length/2,radius/2,length/2,-radius/2,length/2})model.verticesXY.push_back(xy);
        model.patches.push_back(patch);
    }
    std::vector<double>theta(65);for(size_t j=0;j<theta.size();++j)theta[j]=(150+20.*j/(theta.size()-1))*pi/180;
    double total[2]={};size_t checks=0;
    // Repeat with orientation-row samples, also checking accumulation across
    // changed tile boundaries. The sample case is a regression, not a claim
    // that fullauto currently requests these diagnostics.
    for(bool withSamples:{false,true})for(int count:{32,256,512,1024,2048,2947,4096,1205,3055,8192,1024,2947}){
        std::vector<std::array<AnalyticBackscatter::Vec,3>>frames(count);
        std::vector<double>weights(count);
        for(int i=0;i<count;++i){const double a=i*.039,b=i*.023,ca=cos(a),sa=sin(a),cb=cos(b),sb=sin(b);
            frames[i]={{{ca*cb,sa*cb,-sb},{-sa,ca,0},{ca*sb,sa*sb,cb}}};weights[i]=(i%7+1.)/count;}
        std::vector<double>r[2],s[2];std::vector<AnalyticFacetAverage::Components>samples[2];
        // Alternate order to expose rather than conceal a warmup/order effect.
        for(int pass=0;pass<2;++pass){const int method=(checks%2) ? 1-pass:pass;
            std::printf("CASE method=%s count=%d rows=65 phi=32 samples=%d\n",method ? "reuse":"baseline",count,int(withSamples));std::fflush(stdout);
            std::string error;const auto start=std::chrono::steady_clock::now();
            const bool ok=method ? EvaluateAnalyticFacetGpu(model,frames,weights,theta,32,r[method],s[method],withSamples ? &samples[method]:nullptr,&error)
                : EvaluateAnalyticFacetGpuBaseline(model,frames,weights,theta,32,r[method],s[method],withSamples ? &samples[method]:nullptr,&error);
            require(ok,error);const double seconds=std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count();total[method]+=seconds;
            std::printf("CASE_RESULT method=%s count=%d samples=%d wall_seconds=%.9f\n",method ? "reuse":"baseline",count,int(withSamples),seconds);
        }
        require(r[0].size()==r[1].size()&&s[0].size()==s[1].size()&&samples[0].size()==samples[1].size(),"dimensions differ");
        require(!std::memcmp(r[0].data(),r[1].data(),r[0].size()*sizeof(double)),"reflection differs bitwise");
        require(!std::memcmp(s[0].data(),s[1].data(),s[0].size()*sizeof(double)),"shadow differs bitwise");
        if(withSamples)require(!std::memcmp(samples[0].data(),samples[1].data(),samples[0].size()*sizeof(samples[0][0])),"sample rows differ bitwise");
        ++checks;
    }
    std::printf("PASS analytic retained tiles: %zu bitwise comparisons baseline_seconds=%.9f reuse_seconds=%.9f\n",checks,total[0],total[1]);
    return 0;
}catch(const std::exception&error){std::fprintf(stderr,"FAIL: %s\n",error.what());return 1;}}
