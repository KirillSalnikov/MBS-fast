#include "AnalyticAzimuthGaussian.h"
#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace AnalyticAzimuthGaussian
{
namespace { const double pi=3.1415926535897932384626433832795; }
double ScaledI0(double x)
{
    x=std::fabs(x);
    if(x<50)
    {
        double term=1,sum=1;
        for(int j=1;j<140;++j)
        {term*=x*x/(4*j*j);sum+=term;if(term<sum*1e-16)break;}
        return std::exp(-x)*sum;
    }
    double term=1,sum=1;
    for(int j=1;j<=12;++j){term*=(2*j-1.)*(2*j-1.)/(8*j*x);sum+=term;}
    return sum/std::sqrt(2*pi*x);
}
Beam Aperture(const AnalyticBackscatter::Vec &direction,double area,double norm,double wave)
{
    const double length=std::hypot(std::hypot(direction.x,direction.y),direction.z);
    if(!(length>0) || !std::isfinite(length) || !(area>=0) || !std::isfinite(area)
        || !(norm>=0) || !std::isfinite(norm) || !(wave>0) || !std::isfinite(wave))
        throw std::invalid_argument("invalid physical Gaussian beam aperture");
    Beam b;b.dx=direction.x/length;b.dy=direction.y/length;b.dz=direction.z/length;
    b.peak=.5*norm*area*area/(wave*wave);b.kappa=2*pi*area/(wave*wave);
    b.sine=std::hypot(b.dx,b.dy);b.theta=std::atan2(b.sine,-b.dz);return b;
}
Values Evaluate(const Beam &b,double theta,double azimuth)
{
    const double st=(theta==0 || theta==pi) ? 0:std::sin(theta),ct=std::cos(theta);
    const double x=b.dx-st*std::cos(azimuth),y=b.dy-st*std::sin(azimuth),z=b.dz+ct;
    const double half=std::sin((theta-b.theta)/2),B=b.kappa*st*b.sine;
    if(2*b.kappa*half*half>PolarExponentLimit)return Values{0,0};
    return Values{b.peak*std::exp(-.5*b.kappa*(x*x+y*y+z*z)),
        b.peak*std::exp(-2*b.kappa*half*half)*ScaledI0(B)};
}
bool EvaluateCells(const std::vector<std::vector<Beam>>&poses,
    const std::vector<double>&weights,const std::vector<double>&theta,int phi,
    std::vector<double>&point,std::vector<double>&mean,std::vector<Values>*samples)
{
    if(poses.empty() || weights.size()!=poses.size() || theta.empty() || phi<1)return false;
    const size_t rows=theta.size(),cells=rows*phi;
    point.assign(cells,0);mean.assign(cells,0);
    if(samples)samples->assign(poses.size()*rows,Values{0,0});
    // Each cell has its own sum; diagnostic samples are partitioned by theta.
    #pragma omp parallel for schedule(static)
    for(size_t t=0;t<rows;++t)
    {
        for(size_t i=0;i<poses.size();++i)for(int p=0;p<phi;++p)
        {
            Values value{0,0};
            for(const Beam &b:poses[i]){const auto v=Evaluate(b,theta[t],2*pi*p/phi);value.point+=v.point;value.mean+=v.mean;}
            point[p*rows+t]+=weights[i]*value.point;mean[p*rows+t]+=weights[i]*value.mean;
            if(samples){(*samples)[i*rows+t].point+=value.point/phi;(*samples)[i*rows+t].mean+=value.mean/phi;}
        }
    }
    return true;
}
}
