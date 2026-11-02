#include "../src/AnalyticFacetAverage.cpp"
#include "AnalyticRadialTail.h"
#include <iomanip>

int main(){
    double worst=0;
    for(int degree:{128,192,320}){
        AnalyticFacetAverage::Basis basis(degree,65,.06);
        for(double angle:{.01,21.,90.,150.,159.37,170.,179.999,180.}){
            const auto a=AnalyticFacetAverage::Control::ReflectionMoments({1.3116,0},degree,angle*std::acos(-1.)/180);
            std::vector<double> factor(degree/2+1);double p=1;
            for(size_t l=0;l<factor.size();++l){if(l)p*=-(2*l-1.)/(2*l);factor[l]=a[2*l]*(l%2?-1:1)*p;}
            for(double z:{degree+10.000001,400.,1000.,20000.,100000.}){
                const auto expected=basis.Exact(z,factor);double v,d;
                AnalyticRadialTail::Exact(degree,z,factor.data(),v,d);
                const double error=std::max(std::fabs(v-expected[0])/std::max(1e-30,std::fabs(expected[0])),
                    std::fabs(d-expected[1])/std::max(1e-30,std::fabs(expected[1])));
                worst=std::max(worst,error);if(error>1e-12)throw std::runtime_error("streamed radial primitive mismatch");
            }
        }
    }
    std::cout<<"PASS: streamed primitive vs existing Bessel recurrence, max relative error="<<std::setprecision(17)<<worst<<'\n';
}
