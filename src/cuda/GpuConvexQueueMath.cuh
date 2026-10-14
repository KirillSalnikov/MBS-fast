#pragma once
#include "GpuConvexQueue.h"
#include <cmath>
#include <cfloat>

#ifdef __CUDACC__
#define CQ_INLINE __host__ __device__ inline
#else
#define CQ_INLINE inline
#endif
namespace GpuConvexQueue {
constexpr double Epsilon = 0x1p-40;
CQ_INLINE double maximum(double a,double b) { return a>b?a:b; }
CQ_INLINE double minimum(double a,double b) { return a<b?a:b; }
CQ_INLINE Vec add(Vec a,Vec b) { return {a.x+b.x,a.y+b.y,a.z+b.z}; }
CQ_INLINE Vec sub(Vec a,Vec b) { return {a.x-b.x,a.y-b.y,a.z-b.z}; }
CQ_INLINE Vec mul(Vec a,double s) { return {a.x*s,a.y*s,a.z*s}; }
CQ_INLINE double dot(Vec a,Vec b) { return a.x*b.x+a.y*b.y+a.z*b.z; }
CQ_INLINE double norm(Vec a) { return dot(a,a); }
CQ_INLINE Vec cross(Vec a,Vec b) { return {a.y*b.z-a.z*b.y,a.z*b.x-a.x*b.z,a.x*b.y-a.y*b.x}; }
CQ_INLINE Vec unit(Vec a) { double l=sqrt(norm(a));if(!(l>DBL_MIN)||!isfinite(l))return {0,0,0};return {a.x/l,a.y/l,a.z/l}; }
CQ_INLINE Cx add(Cx a,Cx b) { return {a.re+b.re,a.im+b.im}; }
CQ_INLINE Cx sub(Cx a,Cx b) { return {a.re-b.re,a.im-b.im}; }
CQ_INLINE Cx mul(Cx a,double s) { return {a.re*s,a.im*s}; }
CQ_INLINE Cx mul(Cx a,Cx b) { return {a.re*b.re-a.im*b.im,a.im*b.re+a.re*b.im}; }
CQ_INLINE Cx divide(Cx a,Cx b) {
    if(fabs(b.re)>=fabs(b.im)) {
        double r=b.im/b.re,d=b.re+r*b.im;
        return {(a.re+r*a.im)/d,(a.im-r*a.re)/d};
    }
    double r=b.re/b.im,d=b.im+r*b.re;
    return {(r*a.re+a.im)/d,(r*a.im-a.re)/d};
}
CQ_INLINE Vec center(const Beam& b) {
    Vec result={0,0,0};for(int i=0;i<b.count;++i)result=add(result,b.vertices[i]);
    return mul(result,1.0/b.count);
}
CQ_INLINE double area(const Beam& b) {
    if(b.count<3)return 0;
    Vec p1=sub(b.vertices[1],b.vertices[0]),sum={0,0,0};
    for(int i=2;i<b.count;++i){Vec p2=sub(b.vertices[i],b.vertices[0]);sum=add(sum,cross(p1,p2));p1=p2;}
    return .5*sqrt(norm(sum));
}
CQ_INLINE double jnorm(const Beam& b) {
    double sum=0;for(int i=0;i<4;++i)sum+=b.j[i].re*b.j[i].re+b.j[i].im*b.j[i].im;
    return sum;
}
CQ_INLINE bool pruned(const Input& in,const Beam& b,Status& st) {
    if(!b.depth)return false;
    bool active[3];for(int k=0;k<3;++k)active[k]=in.cutoff[k]>0&&in.reference[k]>0;
    if(!active[0]&&!active[1]&&!active[2])return false;
    ++st.evaluated;
    double j=jnorm(b),a=area(b);
    bool small[3]={active[0]&&j<in.cutoff[0]*in.reference[0],
        active[1]&&a<in.cutoff[1]*in.reference[1],active[2]&&j*a<in.cutoff[2]*in.reference[2]};
    if(!small[0]&&!small[1]&&!small[2])return false;
    ++st.rejected;st.rejectedJones+=small[0];st.rejectedArea+=small[1];st.rejectedImportance+=small[2];return true;
}
CQ_INLINE bool inside(Vec x,Vec p1,Vec p2,Vec n) {
    Vec e=sub(p2,p1),q=sub(x,p1);double side=dot(cross(e,q),n);
    if(side>=0)return true;
    double e2=norm(e);return side*side<=Epsilon*Epsilon*(e2*maximum(e2,norm(q))*norm(n));
}
CQ_INLINE bool intersection(Vec a1,Vec a2,Vec b1,Vec b2,Vec n,Vec& result) {
    Vec va=sub(a2,a1),tr=cross(sub(b2,b1),n);
    double den=dot(va,tr),scale=sqrt(norm(va)*norm(tr));
    if(scale<=DBL_MIN||fabs(den)<=Epsilon*maximum(fabs(scale),DBL_MIN))return false;
    result=sub(a1,mul(va,dot(sub(a1,b1),tr)/den));
    return isfinite(result.x)&&isfinite(result.y)&&isfinite(result.z);
}
CQ_INLINE bool append(Vec* p,int& count,Vec v,Status& st) {
    if(count>=MaxVertices){st.failure=1;return false;}
    p[count++]=v;return true;
}
CQ_INLINE void clean(const Vec* p,int n,Beam& b,const Input& in) {
    b.count=0;if(n<=0)return;
    double tolerance2=in.lengthTolerance*in.lengthTolerance;
    if(n==3) {
        Vec a=sub(p[1],p[0]),c=sub(p[0],p[2]);
        if(norm(a)<=tolerance2||norm(sub(p[2],p[1]))<=tolerance2||norm(c)<=tolerance2
           ||norm(cross(a,mul(c,-1)))<=4*in.areaTolerance*in.areaTolerance)return;
        for(int i=0;i<3;++i)b.vertices[i]=p[i];b.count=3;return;
    }
    Vec prev=p[n-1],areaVector={0,0,0},previous={0,0,0};
    for(int i=0;i<n;++i){
        if(norm(sub(p[i],prev))>tolerance2){
            if(b.count){Vec v=sub(p[i],b.vertices[0]);if(b.count>1)areaVector=add(areaVector,cross(previous,v));previous=v;}
            b.vertices[b.count++]=p[i];
        }
        prev=p[i];
    }
    if(b.count<3||norm(areaVector)<=4*in.areaTolerance*in.areaTolerance)b.count=0;
}
CQ_INLINE void intersect(const Input& in,const Beam& b,int face,Beam& result,Status& st) {
    result.count=0;const Face& f=in.faces[face];
    double dp=dot(b.direction,f.inward),scale=sqrt(norm(b.direction)*norm(f.inward));
    if(scale<=DBL_MIN||fabs(dp)<=Epsilon*maximum(fabs(scale),DBL_MIN))return;
    Vec storage[2][MaxVertices];Vec* out=storage[0];Vec* buffer=storage[1];int n=b.count;
    double inv=1.0/dp;
    for(int i=0;i<n;++i)out[i]=sub(b.vertices[i],mul(b.direction,(dot(b.vertices[i],f.inward)+f.plane)*inv));
    Vec p2=f.vertices[f.count-1],normal=mul(f.inward,-1);
    for(int edge=0;edge<f.count;++edge){
        Vec p1=p2;p2=f.vertices[edge];int size=n;n=0;Vec* temp=out;out=buffer;buffer=temp;
        Vec start=buffer[size-1];bool s=inside(start,p1,p2,normal);
        for(int i=0;i<size;++i){Vec end=buffer[i];bool e=inside(end,p1,p2,normal);Vec x;
            if(e!=s&&intersection(start,end,p1,p2,normal,x)&&!append(out,n,x,st))return;
            if(e&&!append(out,n,end,st))return;
            start=end;s=e;
        }
        if(n<3)return;
    }
    clean(out,n,result,in);if(result.count<3)return;
    const Face& source=in.faces[b.facet];
    double den=dot(b.direction,source.inward),denScale=sqrt(norm(b.direction)*norm(source.inward));
    if(denScale>DBL_MIN&&fabs(den)>Epsilon*maximum(fabs(denScale),DBL_MIN)){
        double sign=den>=0?1.0:-1.0;
        Vec start=result.vertices[result.count-1];
        double sv=(dot(start,source.inward)+source.plane)*sign;
        if(fabs(sv)<=in.depthTolerance)sv=0;
        bool s=sv>=0;n=0;
        for(int i=0;i<result.count;++i){
            Vec end=result.vertices[i];double ev=(dot(end,source.inward)+source.plane)*sign;
            if(fabs(ev)<=in.depthTolerance)ev=0;
            bool e=ev>=0;
            if(s!=e&&fabs(sv-ev)>DBL_MIN){double t=maximum(0.0,minimum(1.0,sv/(sv-ev)));
                if(!append(out,n,add(start,mul(sub(end,start),t)),st))return;}
            if(e&&!append(out,n,end,st))return;
            start=end;sv=ev;s=e;
        }
        clean(out,n,result,in);
    }
}
CQ_INLINE void rows(Beam& b,Cx v,Cx h) {
    b.j[0]=mul(b.j[0],v);b.j[1]=mul(b.j[1],v);b.j[2]=mul(b.j[2],h);b.j[3]=mul(b.j[3],h);
}
CQ_INLINE bool finiteBeam(const Beam& b) {
    if(!isfinite(b.path)||!isfinite(b.front)||!isfinite(norm(b.direction))||!isfinite(norm(b.basis)))return false;
    for(int i=0;i<4;++i)if(!isfinite(b.j[i].re)||!isfinite(b.j[i].im))return false;
    for(int i=0;i<b.count;++i)if(!isfinite(norm(b.vertices[i])))return false;
    return true;
}
CQ_INLINE void rotate(Beam& b,Vec normal) {
    Vec basis=unit(cross(normal,b.direction));
    double cs=maximum(-1.0,minimum(1.0,dot(basis,b.basis)));
    if(fabs(1-cs)>=64*DBL_EPSILON){
        if(fabs(1+cs)<64*DBL_EPSILON){for(int i=0;i<4;++i)b.j[i]=mul(b.j[i],-1);}
        else {Vec k=unit(cross(b.basis,basis));double angle=acos(cs);
            if(norm(add(k,b.direction))<=.5)angle=-angle;
            double sn=sin(angle);
            Cx a=add(mul(b.j[0],cs),mul(b.j[2],sn)),c=add(mul(b.j[1],cs),mul(b.j[3],sn));
            b.j[2]=sub(mul(b.j[2],cs),mul(b.j[0],sn));b.j[3]=sub(mul(b.j[3],cs),mul(b.j[1],sn));b.j[0]=a;b.j[1]=c;
        }
    }
    b.basis=basis;
}
// Returns whether a transmitted beam exists. The reflected beam always exists.
CQ_INLINE bool split(const Input& input,const Beam& original,Vec normal,double cosA,Beam& reflected,Beam& transmitted) {
    const double ri=input.ri,eff=ri*ri;
    Vec r=sub(original.direction,mul(normal,cosA));double s=1.0/eff-norm(r);
    bool normalIncidence=cosA>=0.99999999998254670756866631966593;
    Beam incident=original;if(!normalIncidence)rotate(incident,normal);
    transmitted=reflected;
    reflected.direction=normalIncidence?mul(incident.direction,-1):unit(sub(incident.direction,mul(normal,2*cosA)));
    reflected.basis=transmitted.basis=incident.basis;
    for(int j=0;j<4;++j)reflected.j[j]=transmitted.j[j]=incident.j[j];
    bool outgoing=normalIncidence||s>0;
    if(normalIncidence){
        transmitted.direction=incident.direction;
        Cx t=divide(Cx{1-ri,0},Cx{1+ri,0});rows(reflected,t,mul(t,-1));
        Cx u=divide(Cx{2*ri,0},Cx{1+ri,0});rows(transmitted,u,u);
    }else if(outgoing){
        transmitted.direction=unit(add(r,mul(normal,sqrt(maximum(0,s)))));
        double cosG=dot(normal,transmitted.direction),t0=ri*cosA,t1=ri*cosG;
        rows(reflected,divide(Cx{cosA-t1,0},Cx{t1+cosA,0}),divide(Cx{t0-cosG,0},Cx{t0+cosG,0}));
        rows(transmitted,divide(Cx{2*t0,0},Cx{t1+cosA,0}),divide(Cx{2*t0,0},Cx{t0+cosG,0}));
    }else{
        Cx sq={0,sqrt(maximum(eff*(1-cosA*cosA)-1,0))},t0={ri*cosA,0},t1=mul(sq,ri),c={cosA,0};
        rows(reflected,divide(sub(c,t1),add(t1,c)),divide(sub(t0,sq),add(t0,sq)));
    }
    double path=(dot(original.direction,center(reflected))+original.front)*sqrt(eff)+original.path;
    reflected.path=path;reflected.front=dot(mul(reflected.direction,-1),center(reflected));
    transmitted.front=outgoing?dot(mul(transmitted.direction,-1),center(transmitted)):0;
    transmitted.path=path+transmitted.front;
    return outgoing;
}
}
#undef CQ_INLINE
