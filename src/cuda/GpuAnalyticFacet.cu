#include "GpuAnalyticFacet.h"
#include <cuda_runtime.h>
#include <cmath>
#include <cstdio>
#include <cstring>
#include <map>
#include <memory>

namespace
{
const double pi=3.1415926535897932384626433832795;
__device__ double sinc_value(double x)
{
    if(fabs(x)<1e-3){double q=x*x;return 1-q/6+q*q/120-q*q*q/5040;}
    return sin(x)/x;
}
__device__ double reflectance(double mu,double re,double im)
{
    mu=fmin(1.,fmax(0.,mu));double nr=re*re-im*im,ni=2*re*im;
    double a=nr-1+mu*mu,b=ni,length=hypot(a,b);
    double gr,gi;
    if(a>=0){gr=sqrt(.5*(length+a));gi=gr>0 ? b/(2*gr):0;}
    else{gi=copysign(sqrt(.5*(length-a)),b);gr=gi!=0 ? b/(2*gi):0;}
    double rs=((mu-gr)*(mu-gr)+gi*gi)/((mu+gr)*(mu+gr)+gi*gi);
    double pr=nr*mu-gr,ps=ni*mu-gi,qr=nr*mu+gr,qs=ni*mu+gi;
    return .5*(rs+(pr*pr+ps*ps)/(qr*qr+qs*qs));
}
__device__ double dot3(const double*a,const double*b){return a[0]*b[0]+a[1]*b[1]+a[2]*b[2];}
__device__ double fourier_norm(const AnalyticGpuPatch &p,const double*xy,double qx,double qy)
{
    xy+=2*p.firstVertex;double maximum=0;
    for(int i=0;i<p.vertices;++i)maximum=fmax(maximum,fabs(qx*xy[2*i]+qy*xy[2*i+1]));
    double sr=0,si=0;
    if(maximum<1e-3)
    {
        double p0=qx*xy[0]+qy*xy[1];
        for(int i=1;i+1<p.vertices;++i)
        {
            double area=.5*((xy[2*i]-xy[0])*(xy[2*(i+1)+1]-xy[1])-(xy[2*i+1]-xy[1])*(xy[2*(i+1)]-xy[0]));
            double p1=qx*xy[2*i]+qy*xy[2*i+1],p2=qx*xy[2*(i+1)]+qy*xy[2*(i+1)+1];
            sr+=area*(1-(p0*p0+p1*p1+p2*p2+p0*p1+p0*p2+p1*p2)/12);
            si+=area*(p0+p1+p2)/3;
        }
        return sr*sr+si*si;
    }
    for(int i=0;i<p.vertices;++i)
    {
        int j=(i+1)%p.vertices;
        double dx=xy[2*j]-xy[2*i],dy=xy[2*j+1]-xy[2*i+1];
        double phase=.5*(qx*(xy[2*i]+xy[2*j])+qy*(xy[2*i+1]+xy[2*j+1]));
        double value=(qx*dy-qy*dx)*sinc_value(.5*(qx*dx+qy*dy));
        double s,c;sincos(phase,&s,&c);sr+=value*c;si+=value*s;
    }
    double q2=qx*qx+qy*qy;return (sr*sr+si*si)/(q2*q2);
}
__global__ void controls_kernel(const AnalyticGpuPatch*patches,int count,const double*xy,
    const double*frames,const double*theta,int rows,int phi,int orientations,
    double nr,double ni,double wave,double radius,int shadowMode,double*reflection,double*shadow)
{
    size_t index=static_cast<size_t>(blockIdx.x)*blockDim.x+threadIdx.x;
    size_t cells=static_cast<size_t>(rows)*phi,total=cells*orientations;if(index>=total)return;
    int pose=index/cells,cell=index%cells,row=cell%rows,p=cell/rows;
    const double*f=frames+pose*9;const double*u=f+6;double angle=theta[row];
    if(nr==1 && ni==0){reflection[index]=0;shadow[index]=0;return;}
    double st,ct,sp,cp;sincos(angle,&st,&ct);sincos(2*pi*p/phi,&sp,&cp);
    if(angle==0 || angle==pi)st=0;
    double v[3],rq[3],sq[3],k=2*pi/wave,h2=0,c2=0;
    for(int i=0;i<3;++i)
    {
        v[i]=st*(cp*f[i]+sp*f[3+i])-ct*u[i];
        double sum=u[i]+v[i],difference=u[i]-v[i];
        h2+=.25*sum*sum;c2+=.25*difference*difference;rq[i]=k*sum;
    }
    h2=fmin(1.,h2);c2=fmin(1.,c2);
    const double t=dot3(u,v);
    for(int i=0;i<3;++i)sq[i]=k*(v[i]-t*u[i]);
    double r=0,s=0,area=0;
    for(int i=0;i<count;++i)
    {
        const AnalyticGpuPatch &patch=patches[i];double mu=dot3(patch.normal,u);if(mu<=AnalyticFacetAverage::GrazingTolerance)continue;
        area+=mu*patch.area;
        if(h2>0)
        {
            double nu=dot3(patch.normal,v),weight=h2*fmax(0.,c2+mu*nu)*reflectance(mu,nr,ni);
            r+=weight*fourier_norm(patch,xy,dot3(rq,patch.horizontal),dot3(rq,patch.vertical))/(wave*wave);
        }
        if(shadowMode==0 && c2>0 && h2>0)
            s+=c2*c2*mu*mu*fourier_norm(patch,xy,dot3(sq,patch.horizontal),dot3(sq,patch.vertical))/(wave*wave);
    }
    if(shadowMode!=2 && (h2==0 || shadowMode==1))
    {
        double x=k*radius*st,env;
        if(fabs(x)<1e-3){double q=x*x;env=1-q/8+q*q/192-q*q*q/9216;}
        else env=2*j1(x)/x;
        s=c2*c2*env*env*area*area/(wave*wave);
    }
    reflection[index]=r;shadow[index]=s;
}
__global__ void sum_kernel(const double*r,const double*s,const double*weights,int orientations,int cells,double*sumR,double*sumS)
{
    int cell=blockIdx.x*blockDim.x+threadIdx.x;if(cell>=cells)return;double a=0,b=0;
    for(int i=0;i<orientations;++i){size_t j=static_cast<size_t>(i)*cells+cell;a+=weights[i]*r[j];b+=weights[i]*s[j];}
    sumR[cell]=a;sumS[cell]=b;
}
__device__ double scaled_i0(double x)
{
    x=fabs(x);double term=1,sum=1;
    if(x<50)
    {
        for(int j=1;j<140;++j){term*=x*x/(4*j*j);sum+=term;if(term<sum*1e-16)break;}
        return exp(-x)*sum;
    }
    for(int j=1;j<=12;++j){term*=(2*j-1.)*(2*j-1.)/(8*j*x);sum+=term;}
    return sum/sqrt(2*pi*x);
}
__global__ void gaussian_kernel(const AnalyticAzimuthGaussian::Beam*beams,const int*offset,
    const double*theta,int rows,int phi,int orientations,double*point,double*mean)
{
    const size_t index=static_cast<size_t>(blockIdx.x)*blockDim.x+threadIdx.x;
    const size_t cells=static_cast<size_t>(rows)*phi;if(index>=cells*orientations)return;
    const int pose=index/cells,cell=index%cells,t=cell%rows,p=cell/rows;
    double st,ct,sp,cp;sincos(theta[t],&st,&ct);sincos(2*pi*p/phi,&sp,&cp);
    if(theta[t]==0 || theta[t]==pi)st=0;
    double a=0,m=0;
    for(int i=offset[pose];i<offset[pose+1];++i)
    {
        const auto &b=beams[i];const double x=b.dx-st*cp,y=b.dy-st*sp,z=b.dz+ct;
        const double half=sin((theta[t]-b.theta)/2);
        if(2*b.kappa*half*half>AnalyticAzimuthGaussian::PolarExponentLimit)continue;
        a+=b.peak*exp(-.5*b.kappa*(x*x+y*y+z*z));
        m+=b.peak*exp(-2*b.kappa*half*half)*scaled_i0(b.kappa*st*b.sine);
    }
    point[index]=a;mean[index]=m;
}
template<class T>struct Memory
{
    T*pointer=nullptr;
    size_t capacity=0;
    ~Memory(){if(pointer)cudaFree(pointer);}
    bool Allocate(size_t count)
    {
        if(count<=capacity)return true;
        if(pointer){cudaFree(pointer);pointer=nullptr;capacity=0;}
        if(cudaMalloc(reinterpret_cast<void**>(&pointer),count*sizeof(T))!=cudaSuccess)return false;
        capacity=count;return true;
    }
    bool Load(const T*data,size_t count)
    {return Allocate(count) && (!count || cudaMemcpy(pointer,data,count*sizeof(T),cudaMemcpyHostToDevice)==cudaSuccess);}
};
struct FacetWorkspace
{
    Memory<AnalyticGpuPatch>patch;
    Memory<double>xy,angles,frames,weights,reflection,shadow,sumReflection,sumShadow;
    std::vector<AnalyticGpuPatch>hostPatches;
    std::vector<double>hostXY,hostTheta;
    bool Bind(const AnalyticGpuModel &model,const std::vector<double> &theta)
    {
        const bool same=hostPatches.size()==model.patches.size()
            && (hostPatches.empty() || std::memcmp(hostPatches.data(),model.patches.data(),hostPatches.size()*sizeof(AnalyticGpuPatch))==0)
            && hostXY==model.verticesXY;
        if(!same)
        {
            if(!patch.Load(model.patches.data(),model.patches.size()) || !xy.Load(model.verticesXY.data(),model.verticesXY.size()))return false;
            hostPatches=model.patches;hostXY=model.verticesXY;
        }
        if(hostTheta!=theta)
        {
            if(!angles.Load(theta.data(),theta.size()))return false;
            hostTheta=theta;
        }
        return true;
    }
};
}

bool EvaluateAnalyticFacetGpu(const AnalyticGpuModel &model,
    const std::vector<std::array<AnalyticBackscatter::Vec,3>>&poses,
    const std::vector<double>&weights,const std::vector<double>&theta,int phi,
    std::vector<double>&reflected,std::vector<double>&shadow,
    std::vector<AnalyticFacetAverage::Components>*samples)
{
    if(poses.size()!=weights.size() || poses.empty() || theta.empty() || phi<1)return false;
    for(double angle:theta)if(!std::isfinite(angle) || angle<0 || angle>pi)return false;
    for(const auto &pose:poses)
    {
        for(int a=0;a<3;++a)for(int b=a;b<3;++b)
        {
            const auto &u=pose[a],&v=pose[b];const double dot=u.x*v.x+u.y*v.y+u.z*v.z;
            if(!std::isfinite(dot) || std::fabs(dot-(a==b ? 1.:0.))>1e-10)return false;
        }
        const auto &a=pose[0],&b=pose[1],&u=pose[2];
        const double determinant=(a.y*b.z-a.z*b.y)*u.x+(a.z*b.x-a.x*b.z)*u.y+(a.x*b.y-a.y*b.x)*u.z;
        if(!(determinant>0) || !std::isfinite(determinant))return false;
    }
    const size_t count=poses.size(),rows=theta.size(),cells=rows*phi,total=count*cells;
    std::vector<double>packed(count*9);
    for(size_t i=0;i<count;++i)for(int a=0;a<3;++a)
    {packed[i*9+a*3]=poses[i][a].x;packed[i*9+a*3+1]=poses[i][a].y;packed[i*9+a*3+2]=poses[i][a].z;}
    // Geometry/theta are constant across streamed chunks; retain buffers per
    // thread/device, also avoiding malloc/free synchronization on every chunk.
    int device=0;if(cudaGetDevice(&device)!=cudaSuccess)return false;
    static thread_local std::map<int,std::unique_ptr<FacetWorkspace>> workspaces;
    auto &entry=workspaces[device];if(!entry)entry.reset(new FacetWorkspace);
    auto &cache=*entry;
    if(!cache.Bind(model,theta))return false;
    auto &patch=cache.patch;auto &xy=cache.xy;auto &angles=cache.angles;
    auto &frames=cache.frames;auto &w=cache.weights;auto &r=cache.reflection;auto &s=cache.shadow;
    auto &sr=cache.sumReflection;auto &ss=cache.sumShadow;
    if(!frames.Load(packed.data(),packed.size()) || !w.Load(weights.data(),count)
        || !r.Allocate(total) || !s.Allocate(total) || !sr.Allocate(cells) || !ss.Allocate(cells))return false;
    controls_kernel<<<(total+127)/128,128>>>(patch.pointer,model.patches.size(),xy.pointer,frames.pointer,angles.pointer,
        rows,phi,count,model.indexReal,model.indexImag,model.wave,model.circleRadius,model.shadowMode,r.pointer,s.pointer);
    if(cudaGetLastError()!=cudaSuccess)return false;
    sum_kernel<<<(cells+127)/128,128>>>(r.pointer,s.pointer,w.pointer,count,cells,sr.pointer,ss.pointer);
    if(cudaGetLastError()!=cudaSuccess)return false;
    reflected.resize(cells);shadow.resize(cells);
    if(cudaMemcpy(reflected.data(),sr.pointer,cells*sizeof(double),cudaMemcpyDeviceToHost)!=cudaSuccess
        || cudaMemcpy(shadow.data(),ss.pointer,cells*sizeof(double),cudaMemcpyDeviceToHost)!=cudaSuccess)return false;
    if(samples)
    {
        std::vector<double>hr(total),hs(total);
        if(cudaMemcpy(hr.data(),r.pointer,total*sizeof(double),cudaMemcpyDeviceToHost)!=cudaSuccess
            || cudaMemcpy(hs.data(),s.pointer,total*sizeof(double),cudaMemcpyDeviceToHost)!=cudaSuccess)return false;
        samples->assign(count*rows,AnalyticFacetAverage::Components{0,0});
        for(size_t i=0;i<count;++i)for(int p=0;p<phi;++p)for(size_t t=0;t<rows;++t)
        {size_t j=i*cells+p*rows+t;(*samples)[i*rows+t].reflection+=hr[j]/phi;(*samples)[i*rows+t].shadow+=hs[j]/phi;}
    }
    return true;
}

bool EvaluateAnalyticAzimuthGaussianGpu(
    const std::vector<std::vector<AnalyticAzimuthGaussian::Beam>>&poses,
    const std::vector<double>&weights,const std::vector<double>&theta,int phi,
    std::vector<double>&point,std::vector<double>&mean,
    std::vector<AnalyticAzimuthGaussian::Values>*samples)
{
    if(poses.empty() || poses.size()!=weights.size() || theta.empty() || phi<1)return false;
    const size_t count=poses.size(),rows=theta.size(),cells=rows*phi,total=count*cells;
    std::vector<AnalyticAzimuthGaussian::Beam> beams;std::vector<int>offset(1,0);
    for(const auto &pose:poses){beams.insert(beams.end(),pose.begin(),pose.end());offset.push_back(beams.size());}
    if(beams.empty())
    {point.assign(cells,0);mean.assign(cells,0);if(samples)samples->assign(count*rows,AnalyticAzimuthGaussian::Values{0,0});return true;}
    Memory<AnalyticAzimuthGaussian::Beam>packed;Memory<int>indices;Memory<double>angles,w,r,s,sr,ss;
    if(!packed.Load(beams.data(),beams.size()) || !indices.Load(offset.data(),offset.size())
        || !angles.Load(theta.data(),rows) || !w.Load(weights.data(),count)
        || !r.Allocate(total) || !s.Allocate(total) || !sr.Allocate(cells) || !ss.Allocate(cells))return false;
    gaussian_kernel<<<(total+127)/128,128>>>(packed.pointer,indices.pointer,angles.pointer,rows,phi,count,r.pointer,s.pointer);
    if(cudaGetLastError()!=cudaSuccess)return false;
    sum_kernel<<<(cells+127)/128,128>>>(r.pointer,s.pointer,w.pointer,count,cells,sr.pointer,ss.pointer);
    if(cudaGetLastError()!=cudaSuccess)return false;
    point.resize(cells);mean.resize(cells);
    if(cudaMemcpy(point.data(),sr.pointer,cells*sizeof(double),cudaMemcpyDeviceToHost)!=cudaSuccess
        || cudaMemcpy(mean.data(),ss.pointer,cells*sizeof(double),cudaMemcpyDeviceToHost)!=cudaSuccess)return false;
    if(samples)
    {
        std::vector<double>hr(total),hs(total);
        if(cudaMemcpy(hr.data(),r.pointer,total*sizeof(double),cudaMemcpyDeviceToHost)!=cudaSuccess
            || cudaMemcpy(hs.data(),s.pointer,total*sizeof(double),cudaMemcpyDeviceToHost)!=cudaSuccess)return false;
        samples->assign(count*rows,AnalyticAzimuthGaussian::Values{0,0});
        for(size_t i=0;i<count;++i)for(int p=0;p<phi;++p)for(size_t t=0;t<rows;++t)
        {const size_t j=i*cells+p*rows+t;(*samples)[i*rows+t].point+=hr[j]/phi;(*samples)[i*rows+t].mean+=hs[j]/phi;}
    }
    return true;
}
