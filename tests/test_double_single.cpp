#include "DoubleSingle.cuh"
#include <random>
#include <cstdio>
#include <algorithm>
using MbsDS::Number;
long double value(Number x){return (long double)x.hi+(long double)x.lo;}
int main(){
 std::mt19937_64 rng(73913);std::uniform_real_distribution<double>d(-1,1);
 long double add_error=0,mul_error=0,div_error=0,trig_error=0,fast_error=0,rounded_error=0;
 for(int i=0;i<200000;i++){
  double x=d(rng)*std::pow(10.0,(i%19)-9),y=d(rng)*std::pow(10.0,((i/19)%19)-9);
  Number a(x),b(y);if(i%7==0)b=-a+Number(float(std::ldexp(d(rng),-35)));
  long double av=value(a),bv=value(b),scale=std::max(fabsl(av),fabsl(bv));
  add_error=std::max(add_error,fabsl(value(a+b)-(av+bv))/std::max(scale,1e-100L));
  mul_error=std::max(mul_error,fabsl(value(a*b)-av*bv)/std::max(fabsl(av*bv),1e-100L));
  if(bv!=0)div_error=std::max(div_error,fabsl(value(a/b)-av/bv)/std::max(fabsl(av/bv),1e-100L));
  Number phase(d(rng)*1048576.0),s,c;MbsDS::sincos_pair(phase,s,c);
  long double ref_s=sinl(value(phase)),ref_c=cosl(value(phase));
  trig_error=std::max(trig_error,hypotl(value(s)-ref_s,value(c)-ref_c));
  float sf,cf;MbsDS::sincos_float(phase,&sf,&cf,false);
  fast_error=std::max(fast_error,hypotl((long double)sf-ref_s,(long double)cf-ref_c));
  MbsDS::sincos_float(phase,&sf,&cf,true);
  rounded_error=std::max(rounded_error,hypotl((long double)sf-ref_s,(long double)cf-ref_c));
 }
 printf("{\"samples\":200000,\"add_scaled_error\":%.12Lg,\"mul_relative_error\":%.12Lg,\"div_relative_error\":%.12Lg,\"pair_trig_absolute_error\":%.12Lg,\"fast_trig_absolute_error\":%.12Lg,\"rounded_pair_trig_absolute_error\":%.12Lg}\n",add_error,mul_error,div_error,trig_error,fast_error,rounded_error);
 return !(add_error<3e-14L&&mul_error<3e-14L&&div_error<3e-14L&&trig_error<1e-8L&&fast_error<3e-7L&&rounded_error<1e-7L);
}
