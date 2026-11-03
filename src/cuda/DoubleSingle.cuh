#pragma once
#include <cmath>
#include <cfloat>
#ifdef __CUDACC__
#include <cuda_runtime.h>
#define DS_INLINE __host__ __device__ __forceinline__
#else
#define DS_INLINE inline
#endif

// Normalized two-word arithmetic for bounded diffraction intermediates. No fast-math/reassociation:
// each primitive explicitly retains the rounding required by TwoSum/TwoProd.
namespace MbsDS {
DS_INLINE float add_rn(float a,float b) {
#ifdef __CUDA_ARCH__
    return __fadd_rn(a,b);
#else
    volatile float r=a+b; return r;
#endif
}
DS_INLINE float mul_rn(float a,float b) {
#ifdef __CUDA_ARCH__
    return __fmul_rn(a,b);
#else
    volatile float r=a*b; return r;
#endif
}
DS_INLINE float fma_rn(float a,float b,float c) {
#ifdef __CUDA_ARCH__
    return __fmaf_rn(a,b,c);
#else
    return ::fmaf(a,b,c);
#endif
}
DS_INLINE float div_rn(float a,float b) {
#ifdef __CUDA_ARCH__
    return __fdiv_rn(a,b);
#else
    volatile float r=a/b; return r;
#endif
}
struct Number {
    float hi,lo;
    DS_INLINE Number():hi(0),lo(0){}
    DS_INLINE Number(float x):hi(x),lo(0){}
    DS_INLINE Number(int x):hi((float)x),lo(0){}
    DS_INLINE Number(double x):hi((float)x),lo((float)(x-(double)hi)){}
    DS_INLINE Number(float h,float l):hi(h),lo(l){}
    DS_INLINE explicit operator float()const {return add_rn(hi,lo);}
    DS_INLINE explicit operator double()const {return (double)hi+(double)lo;}
};
DS_INLINE Number two_sum(float a,float b) {
    float s=add_rn(a,b),bb=add_rn(s,-a);
    return Number(s,add_rn(add_rn(a,-add_rn(s,-bb)),add_rn(b,-bb)));
}
DS_INLINE Number quick_sum(float a,float b) {
    float s=add_rn(a,b); return Number(s,add_rn(b,-add_rn(s,-a)));
}
DS_INLINE Number operator-(Number a){return Number(-a.hi,-a.lo);}
DS_INLINE Number operator+(Number a,Number b){
    Number s=two_sum(a.hi,b.hi),t=two_sum(a.lo,b.lo);
    s.lo=add_rn(s.lo,t.hi); s=quick_sum(s.hi,s.lo);
    s.lo=add_rn(s.lo,t.lo); return quick_sum(s.hi,s.lo);
}
DS_INLINE Number operator-(Number a,Number b){return a+(-b);}
DS_INLINE Number operator*(Number a,Number b){
    float p=mul_rn(a.hi,b.hi),e=fma_rn(a.hi,b.hi,-p);
    e=fma_rn(a.hi,b.lo,e);e=fma_rn(a.lo,b.hi,e);
    return quick_sum(p,e);
}
DS_INLINE Number operator/(Number a,Number b){
    float q=div_rn(a.hi,b.hi);Number r=a-b*Number(q);
    return quick_sum(q,div_rn(r.hi,b.hi));
}
DS_INLINE Number& operator+=(Number&a,Number b){a=a+b;return a;}
DS_INLINE Number& operator-=(Number&a,Number b){a=a-b;return a;}
DS_INLINE Number& operator*=(Number&a,Number b){a=a*b;return a;}
DS_INLINE bool operator<(Number a,Number b){return a.hi<b.hi||(a.hi==b.hi&&a.lo<b.lo);}
DS_INLINE bool operator>(Number a,Number b){return b<a;}
DS_INLINE bool operator<=(Number a,Number b){return a<b || (a.hi==b.hi && a.lo==b.lo);}
DS_INLINE bool operator>=(Number a,Number b){return b<=a;}
DS_INLINE Number fabs(Number a){return a.hi<0||(a.hi==0&&a.lo<0)?-a:a;}
DS_INLINE Number fmax(Number a,Number b){return a<b?b:a;}
DS_INLINE Number fmin(Number a,Number b){return a<b?a:b;}

// Range reduction uses a pair for pi/2; normal use is |phase| <= 2^20.
// Larger/nonfinite arguments use native FP64 trigonometry of the completed
// DS argument. This does not recover rounding in earlier DS operations.
// This type has FP32 exponent range; it is not a general IEEE FP64 replacement.
DS_INLINE Number reduce_quadrant(Number x,int &quadrant) {
    float q=::nearbyintf(mul_rn(x.hi,0.63661977236758134308f));
    quadrant=((int)q)&3;
    return x-Number(q)*Number(0x1.921fb60000000p+0f,-0x1.777a5c0000000p-25f);
}
DS_INLINE void sincos_pair(Number x,Number &s,Number &c) {
    if (!(::fabsf(x.hi)<=1048576.0f)) {
        double arg=(double)x;
        s=Number(::sin(arg));c=Number(::cos(arg));return;
    }
    int quadrant;Number r=reduce_quadrant(x,quadrant),z=r*r;
    Number ps=Number(0x1.952c780000000p-49f,-0x1.f9ea560000000p-74f);
    ps=Number(-0x1.ae7f3e0000000p-41f,-0x1.ccee080000000p-67f)+z*ps;
    ps=Number(0x1.6124620000000p-33f,-0x1.8af25e0000000p-58f)+z*ps;
    ps=Number(-0x1.ae64560000000p-26f,-0x1.fd51380000000p-52f)+z*ps;
    ps=Number(0x1.71de3a0000000p-19f,0x1.55b1cc0000000p-45f)+z*ps;
    ps=Number(-0x1.a01a020000000p-13f,0x1.7f97fa0000000p-39f)+z*ps;
    ps=Number(0x1.1111120000000p-7f,-0x1.ddddde0000000p-32f)+z*ps;
    ps=Number(-0x1.5555560000000p-3f,0x1.5555560000000p-28f)+z*ps;
    Number pc=Number(-0x1.6827860000000p-53f,-0x1.dcbecc0000000p-80f);
    pc=Number(0x1.ae7f3e0000000p-45f,0x1.ccee080000000p-71f)+z*pc;
    pc=Number(-0x1.93974a0000000p-37f,-0x1.180f940000000p-62f)+z*pc;
    pc=Number(0x1.1eed8e0000000p-29f,0x1.ff1b140000000p-54f)+z*pc;
    pc=Number(-0x1.27e4fc0000000p-22f,0x1.10ec140000000p-47f)+z*pc;
    pc=Number(0x1.a01a020000000p-16f,-0x1.7f97fa0000000p-42f)+z*pc;
    pc=Number(-0x1.6c16c20000000p-10f,0x1.27d27e0000000p-35f)+z*pc;
    pc=Number(0x1.5555560000000p-5f,-0x1.5555560000000p-30f)+z*pc;
    pc=Number(-0x1.0000000000000p-1f,0x0.0p+0f)+z*pc;
    Number sr=r+r*z*ps,cr=Number(1.0f)+z*pc;
    switch(quadrant){case 0:s=sr;c=cr;break;case 1:s=cr;c=-sr;break;
        case 2:s=-sr;c=-cr;break;default:s=-cr;c=sr;break;}
}
DS_INLINE Number sin(Number x){Number s,c;sincos_pair(x,s,c);return s;}
DS_INLINE void sincos_float(Number x,float *s,float *c,bool precise) {
    if (precise || !(::fabsf(x.hi)<=1048576.0f)) {
        Number sd,cd;sincos_pair(x,sd,cd);*s=(float)sd;*c=(float)cd;return;
    }
    int quadrant;Number r=reduce_quadrant(x,quadrant);
    // Bounded polynomial instead of approximate hardware sin/cos: the latter
    // exceeded the predeclared 3e-7 complex-error limit in the CUDA unit test.
    // |r.hi| <= pi/4 plus a small quotient-rounding margin.
    float z=mul_rn(r.hi,r.hi);
    float ps=fma_rn(z,1.0f/362880.0f,-1.0f/5040.0f);
    ps=fma_rn(z,ps,1.0f/120.0f);ps=fma_rn(z,ps,-1.0f/6.0f);
    float sr=fma_rn(mul_rn(r.hi,z),ps,r.hi);
    float pc=fma_rn(z,-1.0f/3628800.0f,1.0f/40320.0f);
    pc=fma_rn(z,pc,-1.0f/720.0f);pc=fma_rn(z,pc,1.0f/24.0f);
    pc=fma_rn(z,pc,-0.5f);float cr=fma_rn(z,pc,1.0f);
    float sc=fma_rn(cr,r.lo,sr),cc=fma_rn(-sr,r.lo,cr);
    switch(quadrant){case 0:*s=sc;*c=cc;break;case 1:*s=cc;*c=-sc;break;
        case 2:*s=-sc;*c=-cc;break;default:*s=-cc;*c=sc;break;}
}
}
using DsReal=MbsDS::Number;
#undef DS_INLINE
