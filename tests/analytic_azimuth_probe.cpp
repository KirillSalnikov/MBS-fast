#include "AnalyticAzimuthGaussian.h"
#include <cmath>
#include <iomanip>
#include <iostream>
#include <string>

int main(int argc,char**argv)
{
    using namespace AnalyticAzimuthGaussian;
    std::cout<<std::setprecision(17);
    if(argc==2 && std::string(argv[1])=="i0")
    {double x;while(std::cin>>x)std::cout<<ScaledI0(x)<<'\n';return 0;}
    double dx,dy,dz,area,norm,wave,theta,phi;
    while(std::cin>>dx>>dy>>dz>>area>>norm>>wave>>theta>>phi)
    {
        const Beam b=Aperture(AnalyticBackscatter::Vec(dx,dy,dz),area,norm,wave);
        const Values v=Evaluate(b,theta,phi);
        std::cout<<b.peak<<' '<<b.kappa<<' '<<v.point<<' '<<v.mean;
#ifdef MBS_GPU_AZIMUTH_PROBE
        const std::vector<std::vector<Beam>> poses{{b}};std::vector<double>p,m;
        if(!EvaluateAnalyticAzimuthGaussianGpu(poses,{1},{theta},1,p,m,nullptr))return 3;
        // The CUDA grid observes phi=0.
        const Values zero=Evaluate(b,theta,0);
        std::cout<<' '<<zero.point<<' '<<p[0]<<' '<<m[0];
#endif
        std::cout<<'\n';
    }
    return 0;
}
