#include "DoubleSingle.cuh"
#include <cuda_runtime.h>
#include <vector>
#include <random>
#include <cstdio>
#include <algorithm>
using MbsDS::Number;
struct Input{Number a,b,p;};
struct Output{Number sum,product,quotient,s,c;float sf,cf,sp,cp;};
__global__ void evaluate(const Input*in,Output*out,int n){
 int i=blockIdx.x*blockDim.x+threadIdx.x;if(i>=n)return;
 Input a=in[i];Output r;r.sum=a.a+a.b;r.product=a.a*a.b;r.quotient=a.a/a.b;
 MbsDS::sincos_pair(a.p,r.s,r.c);MbsDS::sincos_float(a.p,&r.sf,&r.cf,false);
 MbsDS::sincos_float(a.p,&r.sp,&r.cp,true);out[i]=r;
}
long double val(Number x){return (long double)x.hi+(long double)x.lo;}
int main(){
 int n=200000;std::vector<Input>in(n);std::vector<Output>out(n);std::mt19937_64 rng(73913);std::uniform_real_distribution<double>d(-1,1);
 for(int i=0;i<n;i++){
  in[i].a=Number(d(rng)*std::pow(10.,(i%19)-9));in[i].b=Number(d(rng)*std::pow(10.,((i/19)%19)-9));
  if(i%7==0)in[i].b=-in[i].a+Number(float(std::ldexp(d(rng),-35)));
  in[i].p=Number(d(rng)*1048576.0);
 }
 Input*di;Output*doo;
 if(cudaMalloc(&di,n*sizeof(Input))!=cudaSuccess||cudaMalloc(&doo,n*sizeof(Output))!=cudaSuccess)return 2;
 if(cudaMemcpy(di,in.data(),n*sizeof(Input),cudaMemcpyHostToDevice)!=cudaSuccess)return 3;
 evaluate<<<(n+255)/256,256>>>(di,doo,n);
 if(cudaDeviceSynchronize()!=cudaSuccess||cudaMemcpy(out.data(),doo,n*sizeof(Output),cudaMemcpyDeviceToHost)!=cudaSuccess)return 4;
 long double ae=0,me=0,de=0,te=0,fe=0,pe=0;
 for(int i=0;i<n;i++){
  long double a=val(in[i].a),b=val(in[i].b),p=val(in[i].p);Output r=out[i];
  ae=std::max(ae,fabsl(val(r.sum)-(a+b))/std::max(std::max(fabsl(a),fabsl(b)),1e-100L));
  me=std::max(me,fabsl(val(r.product)-a*b)/std::max(fabsl(a*b),1e-100L));
  de=std::max(de,fabsl(val(r.quotient)-a/b)/std::max(fabsl(a/b),1e-100L));
  long double s=sinl(p),c=cosl(p);te=std::max(te,hypotl(val(r.s)-s,val(r.c)-c));
  fe=std::max(fe,hypotl((long double)r.sf-s,(long double)r.cf-c));
  pe=std::max(pe,hypotl((long double)r.sp-s,(long double)r.cp-c));
 }
 printf("{\"samples\":%d,\"add_scaled_error\":%.12Lg,\"mul_relative_error\":%.12Lg,\"div_relative_error\":%.12Lg,\"pair_trig_absolute_error\":%.12Lg,\"fast_trig_absolute_error\":%.12Lg,\"rounded_pair_trig_absolute_error\":%.12Lg}\n",n,ae,me,de,te,fe,pe);
 cudaFree(di);cudaFree(doo);
 return !(ae<3e-14L&&me<3e-14L&&de<3e-14L&&te<1e-8L&&fe<3e-7L&&pe<1e-7L);
}
