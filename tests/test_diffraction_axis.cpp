#include "../src/handler/HandlerPO_fast.h"
#include <random>
#include <cassert>
#include <iostream>
int main(){std::mt19937 rng(12);auto random=[&](){return (int(rng()%20001)-10000)*.001;};
for(int n=0;n<10000;++n){double vc[32],vs[32],a[4],b[4],sx[8],sy[8],r[4],i[4],rr[4],ii[4];bool vx[8],vy[8];for(int j=0;j<32;++j){vc[j]=random();vs[j]=random();}for(int j=0;j<4;++j){a[j]=random();b[j]=random();}for(int j=0;j<8;++j){sx[j]=random();sy[j]=random();vx[j]=rng()%2;vy[j]=rng()%2;}int nv=3+rng()%6;bool expected=diffract_theta_4_regular_impl<-1>(vc,vs,0,4,nv,a,b,sx,vx,sy,vy,.001,.3,.4,r,i);for(bool specialize:{false,true}){bool actual=diffract_theta_4_regular(vc,vs,0,4,nv,a,b,sx,vx,sy,vy,.001,.3,.4,rr,ii,specialize);assert(actual==expected);if(actual)for(int j=0;j<4;++j){assert(r[j]==rr[j]);assert(i[j]==ii[j]);}}}
std::cout<<"PASS 10000 uniform/mixed branch vectors, validity and singular fallback; exact outputs\n";}
