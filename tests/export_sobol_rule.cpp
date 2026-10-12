#include "math/Sobol.h"
#include <cstdlib>
#include <iomanip>
#include <iostream>
int main(int argc,char** argv) {
    if(argc!=3) return 2;
    const int n=std::atoi(argv[1]);
    if(n<=0) return 2;
    Sobol3D points(static_cast<unsigned>(std::strtoul(argv[2],nullptr,10)));
    std::cout<<"# BETA_DEG GAMMA_DEG ALPHA_DEG WEIGHT\n"<<std::setprecision(17);
    for(int i=0;i<n;++i) {
        double u,v,w;points.next(u,v,w);
        std::cout<<std::acos(1-2*u)*180/M_PI<<' '<<360*v<<' '<<360*w<<' '<<1./n<<'\n';
    }
}
