#pragma once
#include <cmath>

#ifdef __CUDACC__
#define MBS_RADIAL_HD __host__ __device__
#else
#define MBS_RADIAL_HD
#endif

// Same exact primitive as Basis::Exact, in the stable upward-Bessel domain
// z > degree+10. Stream two Bessel values instead of allocating degree+1
// values for each radial sample. No truncation of the Legendre expansion.
namespace AnalyticRadialTail {
MBS_RADIAL_HD inline void Exact(int degree,double z,const double* factor,
                                double& value,double& derivative) {
    double a=1,b=1,ta=1,tb=1;
    for(int k=1;k<40;++k){
        const double na=-ta*(2*k)*(2*k-1)/(z*z),nb=-tb*(2*k)*(2*k+1)/(z*z);
        if(fabs(na)>fabs(ta)||fabs(nb)>fabs(tb))break;
        a+=na;b+=nb;ta=na;tb=nb;
        if(fmax(fabs(ta),fabs(tb))<1e-17)break;
    }
    const double sine=sin(z),cosine=cos(z);
    const double g0=.5772156649015328606065120900824024+log(z)
        -(sine/z*a-cosine/(z*z)*b);
    double even=sine/z,odd=even/z-cosine/z;
    value=factor[0]*g0;derivative=factor[0]*(1-cosine)/z;
    double H=1-even,A=1,T=0;const double K=g0-1+even;
    for(int l=2;l<=degree;l+=2){
        if(l>2){
            even=(2*l-5)*odd/z-even;
            odd=(2*l-3)*even/z-odd;
            const int q=l-1;const double next=double(q-1)/q*H-double(2*q-1)/q*even;
            T=double(q-1)/q*T-(H+next)/q;A*=double(q-1)/q;H=next;
        }
        value+=factor[l/2]*(l*A*K+l*T-H);
        derivative+=factor[l/2]*(l*H/z-odd);
    }
}
}
#undef MBS_RADIAL_HD
