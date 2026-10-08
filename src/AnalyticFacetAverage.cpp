#include "AnalyticFacetAverage.h"
#include "cuda/GpuAnalyticFacet.h"
#include <algorithm>
#include <cmath>
#include <functional>
#include <stdexcept>
#include <cstdint>
#include <cstring>
#include <fstream>
#include <iomanip>
#include <sstream>
#include <cstdio>
#include <fcntl.h>
#include <sys/file.h>
#include <unistd.h>

namespace AnalyticFacetAverage
{
std::vector<std::array<double,2>> ReadControlWeights(
    const std::string &path,const std::vector<double> &thetaRadians)
{
    std::vector<std::array<double,2>> result(thetaRadians.size(),{{1.,1.}});
    if(path.empty())return result;
    std::ifstream input(path.c_str());
    if(!input)throw std::runtime_error("cannot read analytic control weights: "+path);
    std::string line;size_t row=0;bool header=false;
    while(std::getline(input,line))
    {
        if(line.empty() || line[0]=='#')continue;
        if(!header)
        {
            std::istringstream fields(line);std::string a,b,c,extra;
            if(!(fields>>a>>b>>c) || fields>>extra || a!="theta_deg"
                || b!="reflection_weight" || c!="shadow_weight")
                throw std::runtime_error("invalid analytic control weight header");
            header=true;continue;
        }
        std::istringstream fields(line);double theta,reflection,shadow;std::string extra;
        if(!(fields>>theta>>reflection>>shadow) || fields>>extra || row>=result.size()
            || !std::isfinite(theta) || !std::isfinite(reflection) || !std::isfinite(shadow)
            || std::fabs(reflection)>4 || std::fabs(shadow)>4
            || std::fabs(theta-thetaRadians[row]*180/3.14159265358979323846)>1e-8)
            throw std::runtime_error("analytic control weights must match every theta row and lie in [-4,4]");
        result[row++]={{reflection,shadow}};
    }
    if(!header || row!=result.size())throw std::runtime_error("incomplete analytic control weight file");
    return result;
}
namespace
{
const double pi=3.1415926535897932384626433832795;
using C=std::complex<double>;
Vec Sub(Vec a,Vec b){return Vec(a.x-b.x,a.y-b.y,a.z-b.z);}
double Dot(Vec a,Vec b){return a.x*b.x+a.y*b.y+a.z*b.z;}
double Length(Vec a){return std::sqrt(Dot(a,a));}
Vec Unit(Vec a)
{
    if(!std::isfinite(a.x) || !std::isfinite(a.y) || !std::isfinite(a.z))
        throw std::invalid_argument("analytic facet: direction must be finite and nonzero");
    const double scale=std::max(std::fabs(a.x),std::max(std::fabs(a.y),std::fabs(a.z)));
    if(!(scale>0))throw std::invalid_argument("analytic facet: direction must be finite and nonzero");
    a=Vec(a.x/scale,a.y/scale,a.z/scale);
    const double norm=std::hypot(std::hypot(a.x,a.y),a.z);
    return Vec(a.x/norm,a.y/norm,a.z/norm);
}
Vec Cross(Vec a,Vec b){return Vec(a.y*b.z-a.z*b.y,a.z*b.x-a.x*b.z,a.x*b.y-a.y*b.x);}
void Gauss(int n,std::vector<double>&x,std::vector<double>&w)
{
    x.resize(n);w.resize(n);
    for(int i=0;i<(n+1)/2;++i)
    {
        double z=std::cos(pi*(i+.75)/(n+.5)),d=0;
        for(int it=0;it<100;++it)
        {
            double a=1,b=z;
            for(int l=2;l<=n;++l){double p=((2*l-1)*z*b-(l-1)*a)/l;a=b;b=p;}
            d=n*(z*b-a)/(z*z-1);double next=z-b/d;
            if(std::fabs(next-z)<4e-16){z=next;break;} z=next;
        }
        double a=1,b=z;
        for(int l=2;l<=n;++l){double p=((2*l-1)*z*b-(l-1)*a)/l;a=b;b=p;}
        d=n*(z*b-a)/(z*z-1);
        x[i]=-z;x[n-1-i]=z;w[i]=w[n-1-i]=2/((1-z*z)*d*d);
    }
}
double Reflectance(double mu,C index)
{
    mu=std::max(0.,std::min(1.,mu));
    C g=std::sqrt(index*index-1.+mu*mu);
    C rs=(mu-g)/(mu+g),rp=(index*index*mu-g)/(index*index*mu+g);
    return .5*(std::norm(rs)+std::norm(rp));
}
double Sinc(double x)
{if(std::fabs(x)<1e-3){double q=x*x;return 1-q/6+q*q/120-q*q*q/5040;}return std::sin(x)/x;}
double Jinc(double x)
{if(std::fabs(x)<1e-3){double q=x*x;return 1-q/8+q*q/192-q*q*q/9216;}return 2*::j1(x)/x;}
void Angles(double theta,double &h2,double &c2)
{
    if(std::fabs(theta)<1e-13){h2=0;c2=1;}
    else if(std::fabs(theta-pi)<1e-13){h2=1;c2=0;}
    else {h2=std::pow(std::sin(theta/2),2);c2=1-h2;}
}
struct Moments
{
    std::vector<double>A,B,Cprev;
    Moments(C index,int degree)
    {
        A.assign(degree+1,0);B=A;Cprev=A;
        std::vector<double>x,w;Gauss(std::max(512,2*degree+64),x,w);
        for(size_t q=0;q<x.size();++q)
        {
            double mu=(x[q]+1)/2,r=.5*w[q]*Reflectance(mu,index),prev=1,p=1;
            for(int l=0;l<=degree;++l)
            {
                A[l]+=r*p;B[l]+=r*mu*mu*p;
                if(l)Cprev[l]+=r*mu*prev;
                double next=l==0 ? mu : ((2*l+1)*mu*p-l*prev)/(l+1);
                prev=p;p=next;
            }
        }
    }
    std::vector<double> At(double theta) const
    {
        double h2,c2;Angles(theta,h2,c2);double h=std::sqrt(h2),t=2*h2-1;
        std::vector<double>a(A.size(),0);double prev=1,p=1;
        for(size_t l=0;l<a.size();++l)
        {
            double spin=l ? 2*h*l/(l+1)*(prev-h*p)*(Cprev[l]-B[l]) : 0;
            a[l]=(2*l+1)/2.*h2*(p*(c2*A[l]+t*B[l])+spin);
            double next=l==0 ? h : ((2*l+1)*h*p-l*prev)/(l+1);
            prev=p;p=next;
        }
        return a;
    }
};
std::vector<double> Cumulative(const std::vector<double>&f,double step)
{
    std::vector<double>result(f.size(),0);
    for(size_t i=2;i<f.size();i+=2)
    {
        result[i]=result[i-2]+step/3*(f[i-2]+4*f[i-1]+f[i]);
        result[i-1]=result[i-2]+step/12*(5*f[i-2]+8*f[i-1]-f[i]);
    }
    return result;
}
// Cubic Hermite radial Green function, with exact polynomial antiderivatives
// for coincident-edge integrals. Phase evaluation never samples orientations.
struct Green
{
    double step;
    std::vector<double>g,d,integral,moment;
    Green(double s,std::vector<double>values,std::vector<double>derivatives)
        :step(s),g(std::move(values)),d(std::move(derivatives)),integral(g.size(),0),moment(g.size(),0)
    {
        for(size_t i=0;i+1<g.size();++i)
        {
            double a,b,c,e;Polynomial(i,a,b,c,e);
            double h=a+b/2+c/3+e/4;
            integral[i+1]=integral[i]+step*h;
            moment[i+1]=moment[i]+step*(i*step)*h+step*step*(a/2+b/3+c/4+e/5);
        }
    }
    void Polynomial(size_t i,double&a,double&b,double&c,double&e)const
    {a=g[i];b=step*d[i];c=3*(g[i+1]-g[i])-step*(2*d[i]+d[i+1]);e=2*(g[i]-g[i+1])+step*(d[i]+d[i+1]);}
    double Value(double z)const
    {
        size_t i=std::min(size_t(std::max(0.,z)/step),g.size()-2);double t=(z-i*step)/step,a,b,c,e;
        Polynomial(i,a,b,c,e);return a+t*(b+t*(c+t*e));
    }
    double Same(double z)const
    {
        if(z==0)return 0;
        size_t i=std::min(size_t(z/step),g.size()-2);double t=(z-i*step)/step,a,b,c,e;
        Polynomial(i,a,b,c,e);
        double h=t*(a+t*(b/2+t*(c/3+t*e/4)));
        double j=t*t*(a/2+t*(b/3+t*(c/4+t*e/5)));
        double H=integral[i]+step*h,J=moment[i]+step*(i*step)*h+step*step*j;
        return 2/z*(H-J/z);
    }
};
struct Basis
{
    int degree;size_t count,terms;double step;
    std::vector<double>g,d;
    Basis(int degree_,double maxz,double step_):degree(degree_),step(step_)
    {
        count=std::max(size_t(4),size_t(std::ceil(maxz/step/2))*2+1);terms=degree/2+1;
        // Two FP64 tables: at most ~1 GiB. The previous 32M-entry limit
        // excluded the 891.3 um test prism at the required refinement.
        if(count>1000001 || count*terms>64000000)
            throw std::runtime_error("analytic facet kernel exceeds radial-table budget.\n  Fix: use a smaller particle or remove --analytic-facet-average.");
        g.resize(count*terms);d.resize(count*terms);
        for(size_t i=0;i<count;++i)
        {
            auto j=AnalyticBackscatter::SphericalBessel(degree,i*step);
            for(size_t l=0;l<terms;++l)g[i*terms+l]=(i*step)*j[2*l];
        }
        std::vector<double>f(count);
        for(size_t l=0;l<terms;++l)
        {
            for(size_t i=0;i<count;++i)f[i]=g[i*terms+l];
            auto flux=Cumulative(f,step);
            for(size_t i=0;i<count;++i)f[i]=i ? flux[i]/(i*step) : 0;
            auto values=Cumulative(f,step);
            for(size_t i=0;i<count;++i){g[i*terms+l]=values[i];d[i*terms+l]=f[i];}
        }
    }
    Green Combine(const std::vector<double>&a,double maxz)const
    {
        size_t n=std::min(count,std::max(size_t(4),size_t(std::ceil(maxz/step))+2));
        std::vector<double>values(n,0),derivatives(n,0),factor(terms);
        double p=1;
        for(size_t l=0;l<terms;++l)
        {
            if(l)p*=-(2*l-1.)/(2*l);
            factor[l]=a[2*l]*(l%2 ? -1:1)*p;
        }
        for(size_t i=0;i<n;++i)
            for(size_t l=0;l<terms;++l)
            {values[i]+=factor[l]*g[i*terms+l];derivatives[i]+=factor[l]*d[i*terms+l];}
        return Green(step,std::move(values),std::move(derivatives));
    }
};
double EdgeMean(const std::vector<double>&px,const std::vector<double>&py,double Q,int nodes,
                const std::function<double(double)>&value,const std::function<double(double)>&same)
{
    std::vector<double>x,w;Gauss(nodes,x,w);
    for(int i=0;i<nodes;++i){x[i]=(x[i]+1)/2;w[i]/=2;}
    double result=0;size_t m=px.size();
    for(size_t i=0;i<m;++i)
    {
        size_t in=(i+1)%m;double ex=px[in]-px[i],ey=py[in]-py[i];
        for(size_t j=i;j<m;++j)
        {
            size_t jn=(j+1)%m;double fx=px[jn]-px[j],fy=py[jn]-py[j];
            double dot=ex*fx+ey*fy;if(dot==0)continue;
            double integral=0;
            if(i==j)integral=same(Q*std::hypot(ex,ey));
            else for(int a=0;a<nodes;++a)for(int b=0;b<nodes;++b)
            {
                double dx=px[i]+x[a]*ex-px[j]-x[b]*fx,dy=py[i]+x[a]*ey-py[j]-x[b]*fy;
                integral+=w[a]*w[b]*value(Q*std::hypot(dx,dy));
            }
            result-=(i==j ? 1:2)*dot*integral/(Q*Q);
        }
    }
    return result;
}
C Fourier(const std::vector<double>&x,const std::vector<double>&y,double qx,double qy)
{
    double maxphase=0;for(size_t i=0;i<x.size();++i)maxphase=std::max(maxphase,std::fabs(qx*x[i]+qy*y[i]));
    if(maxphase<1e-3)
    {
        C result=0;double p0=qx*x[0]+qy*y[0];
        for(size_t i=1;i+1<x.size();++i)
        {
            double area=.5*((x[i]-x[0])*(y[i+1]-y[0])-(y[i]-y[0])*(x[i+1]-x[0]));
            double p1=qx*x[i]+qy*y[i],p2=qx*x[i+1]+qy*y[i+1];
            result+=area*C(1-(p0*p0+p1*p1+p2*p2+p0*p1+p0*p2+p1*p2)/12,(p0+p1+p2)/3);
        }
        return result;
    }
    C result=0;
    for(size_t i=0;i<x.size();++i)
    {
        size_t j=(i+1)%x.size();double dx=x[j]-x[i],dy=y[j]-y[i];
        double dp=qx*dx+qy*dy,phase=.5*(qx*(x[i]+x[j])+qy*(y[i]+y[j]));
        result+=(qx*dy-qy*dx)*std::exp(C(0,phase))*Sinc(dp/2);
    }
    return result/C(0,qx*qx+qy*qy);
}
double ShadowGreen(double z)
{
    if(std::fabs(z)<1e-3){double q=z*z;return .25*(q/6-q*q/120+q*q*q/5040);}
    return .25*(1-Sinc(z));
}
double ShadowSame(double z)
{
    if(std::fabs(z)<1e-3){double q=z*z;return q/144-q*q/7200+q*q*q/564480;}
    return .25-AnalyticBackscatter::FiniteSincIntegral(z/2)/8;
}
}

std::vector<double> Control::ReflectionMoments(C index,int degree,double theta)
{return Moments(index,degree).At(theta);}

AnalyticGpuModel Control::GpuModel() const
{
    AnalyticGpuModel result;result.indexReal=index.real();result.indexImag=index.imag();
    result.wave=wave;result.circleRadius=circleRadius;
    result.shadowMode=shadowMode=="facets" ? 0:(shadowMode=="circular" ? 1:2);
    for(const Patch &p:patches)
    {
        AnalyticGpuPatch target;
        target.normal[0]=p.n.x;target.normal[1]=p.n.y;target.normal[2]=p.n.z;
        target.horizontal[0]=p.h.x;target.horizontal[1]=p.h.y;target.horizontal[2]=p.h.z;
        target.vertical[0]=p.v.x;target.vertical[1]=p.v.y;target.vertical[2]=p.v.z;
        target.area=p.area;target.firstVertex=static_cast<int>(result.verticesXY.size()/2);target.vertices=static_cast<int>(p.x.size());
        for(size_t i=0;i<p.x.size();++i){result.verticesXY.push_back(p.x[i]);result.verticesXY.push_back(p.y[i]);}
        result.patches.push_back(target);
    }
    return result;
}

Control::Control(const std::vector<Face>&input,C refractive,double wavelength,
                 const std::vector<double>&theta,const std::string &mode,const std::string &meanCache)
    :faces(input),index(refractive),wave(wavelength),refinementError(0),forwardShadowMean(0),circleRadius(0),shadowMode(mode)
{
    if(faces.empty() || !(wave>0) || index.real()<1 || !std::isfinite(index.real()+index.imag()+wave))
        throw std::invalid_argument("analytic facet average requires finite index REAL>=1 and positive wavelength");
    if(mode!="facets" && mode!="circular" && mode!="off")throw std::invalid_argument("invalid analytic shadow model");
    // The full transparent-particle field vanishes by shadow/transmission
    // cancellation. Avoid adding control noise to this exact null solution.
    if(index==C(1,0))shadowMode="off";
    double totalArea=0;
    for(const Face&f:faces)
    {
        if(f.vertices.size()<3)throw std::runtime_error("analytic facet average: degenerate facet");
        Patch p;p.n=Unit(f.normal);p.h=Unit(Sub(f.vertices[1],f.vertices[0]));p.v=Unit(Cross(p.n,p.h));p.area=0;p.diameter=0;
        Vec center;for(Vec r:f.vertices){center.x+=r.x/f.vertices.size();center.y+=r.y/f.vertices.size();center.z+=r.z/f.vertices.size();}
        for(Vec r:f.vertices){Vec d=Sub(r,center);p.x.push_back(Dot(d,p.h));p.y.push_back(Dot(d,p.v));}
        for(size_t i=0;i<p.x.size();++i)
        {
            size_t j=(i+1)%p.x.size();p.area+=.5*(p.x[i]*p.y[j]-p.y[i]*p.x[j]);
            for(size_t t=0;t<p.x.size();++t)p.diameter=std::max(p.diameter,std::hypot(p.x[i]-p.x[t],p.y[i]-p.y[t]));
        }
        p.area=std::fabs(p.area);if(!(p.area>0))throw std::runtime_error("analytic facet average: zero polygon area");
        totalArea+=p.area;patches.push_back(p);
    }
    circleRadius=std::sqrt(totalArea/(4*pi));
    for(const Patch&a:patches)for(const Patch&b:patches)
    {
        double dot=std::max(-1.,std::min(1.,Dot(a.n,b.n))),angle=std::acos(dot);
        forwardShadowMean+=a.area*b.area*(std::sin(angle)+(pi-angle)*dot)/(6*pi*wave*wave);
    }
    PrepareCachedMeans(theta,meanCache);
}

void Control::PrepareCachedMeans(const std::vector<double>&theta,const std::string &path)
{
    if(path.empty()){PrepareMeans(theta);return;}
    // Cache EXACT physical inputs and algorithm revision. Never use a fitted
    // mean from another particle, wavelength, grid or Fresnel/shadow model.
    std::uint64_t fingerprint=1469598103934665603ULL;
    const auto mix=[&](const void *data,size_t bytes)
    {
        const unsigned char *p=static_cast<const unsigned char*>(data);
        for(size_t i=0;i<bytes;++i){fingerprint^=p[i];fingerprint*=1099511628211ULL;}
    };
    const char version[]="MBS_PHYSICAL_FACET_MEAN_V2";mix(version,sizeof(version));
    double re=index.real(),im=index.imag();mix(&re,sizeof(re));mix(&im,sizeof(im));mix(&wave,sizeof(wave));
    mix(shadowMode.data(),shadowMode.size());
    const size_t count=faces.size();mix(&count,sizeof(count));
    for(const Face &face:faces)
    {
        mix(&face.normal.x,sizeof(double));mix(&face.normal.y,sizeof(double));mix(&face.normal.z,sizeof(double));
        const size_t vertices=face.vertices.size();mix(&vertices,sizeof(vertices));
        for(const Vec &v:face.vertices){mix(&v.x,sizeof(double));mix(&v.y,sizeof(double));mix(&v.z,sizeof(double));}
    }
    const size_t rows=theta.size();mix(&rows,sizeof(rows));
    for(double angle:theta)mix(&angle,sizeof(angle));
    struct Lock
    {
        int fd;
        explicit Lock(const std::string &file):fd(::open((file+".lock").c_str(),O_CREAT|O_RDWR,0600))
        {
            if(fd<0 || ::flock(fd,LOCK_EX)!=0)
            {if(fd>=0)::close(fd);throw std::runtime_error("cannot lock analytic mean cache.\n  Fix: use a writable cache path with an existing parent directory.");}
        }
        ~Lock(){::flock(fd,LOCK_UN);::close(fd);}
    } lock(path);
    {
        std::ifstream file(path.c_str());std::string magic;std::uint64_t saved=0;size_t n=0;double error=0;
        if(file>>magic>>saved>>n>>error && magic=="MBS_PHYSICAL_FACET_MEAN_V2" && saved==fingerprint
            && n==theta.size() && std::isfinite(error) && error>=0 && error<1e-4)
        {
            std::vector<Components> loaded(n);bool valid=true;
            for(size_t i=0;i<n;++i)
            {
                double angle=0;
                if(!(file>>angle>>loaded[i].reflection>>loaded[i].shadow) || angle!=theta[i]
                    || !std::isfinite(loaded[i].reflection+loaded[i].shadow) || loaded[i].reflection<0 || loaded[i].shadow<0)
                {valid=false;break;}
            }
            if(valid){means.swap(loaded);refinementError=error;return;}
        }
    }
    PrepareMeans(theta);
    const std::string temporary=path+".tmp."+std::to_string(::getpid());
    std::ofstream file(temporary.c_str());file<<std::setprecision(17)
        <<"MBS_PHYSICAL_FACET_MEAN_V2 "<<fingerprint<<' '<<theta.size()<<' '<<refinementError<<'\n';
    for(size_t i=0;i<theta.size();++i)file<<theta[i]<<' '<<means[i].reflection<<' '<<means[i].shadow<<'\n';
    file.flush();if(!file)throw std::runtime_error("cannot save analytic mean cache.\n  Fix: check disk space and cache permissions.");
    file.close();if(::rename(temporary.c_str(),path.c_str())!=0)
        throw std::runtime_error("cannot publish analytic mean cache.\n  Fix: check cache path permissions.");
}

std::array<Vec,3> SourceFrame(const Vec &source,const Vec &azimuthReference)
{
    const Vec u=Unit(source),reference=Unit(azimuthReference);
    const double projection=Dot(reference,u);
    const Vec transverse(reference.x-projection*u.x,reference.y-projection*u.y,reference.z-projection*u.z);
    if(Length(transverse)<1e-12)
        throw std::invalid_argument("analytic facet: azimuth reference must not be parallel to source");
    const Vec horizontal=Unit(transverse);
    return {{horizontal,Unit(Cross(u,horizontal)),u}};
}

Vec Observer(const std::array<Vec,3> &frame,double theta,double azimuth)
{
    if(!std::isfinite(theta) || theta<0 || theta>pi || !std::isfinite(azimuth))
        throw std::invalid_argument("analytic facet: theta must be in [0, pi] and azimuth finite");
    const double st=(theta==0 || theta==pi) ? 0:std::sin(theta),ct=std::cos(theta),cp=std::cos(azimuth),sp=std::sin(azimuth);
    return Vec(st*(cp*frame[0].x+sp*frame[1].x)-ct*frame[2].x,
               st*(cp*frame[0].y+sp*frame[1].y)-ct*frame[2].y,
               st*(cp*frame[0].z+sp*frame[1].z)-ct*frame[2].z);
}

Components Control::Evaluate(const Vec&source,const Vec&observer)const
{return EvaluateUnit(Unit(source),Unit(observer));}

Components Control::Evaluate(const Vec&source,const Vec&observer,double theta)const
{
    const Vec u=Unit(source),v=Unit(observer);
    if(!std::isfinite(theta) || theta<0 || theta>pi || std::fabs(Dot(u,v)+std::cos(theta))>1e-10)
        throw std::invalid_argument("analytic facet: theta is inconsistent with source/observer directions");
    return EvaluateUnit(u,v);
}

Components Control::EvaluateUnit(const Vec&u,const Vec&v)const
{
    if(index==C(1,0))return Components{0,0};
    // Half-angle lengths remain accurate at both poles, unlike 1 +/- u.v.
    const Vec sum(u.x+v.x,u.y+v.y,u.z+v.z),difference=Sub(u,v);
    const double h2=std::min(1.,Dot(sum,sum)/4),c2=std::min(1.,Dot(difference,difference)/4);
    const double k=2*pi/wave,t=Dot(u,v),sine=Length(Cross(u,v));
    Vec reflectedQ(k*(u.x+v.x),k*(u.y+v.y),k*(u.z+v.z));
    Vec shadowQ(k*(v.x-t*u.x),k*(v.y-t*u.y),k*(v.z-t*u.z));
    Components result{0,0};double projectedArea=0;
    for(const Patch&p:patches)
    {
        double mu=Dot(p.n,u);if(mu<=GrazingTolerance)continue;
        projectedArea+=mu*p.area;
        if(h2>0)
        {
            double nu=Dot(p.n,v),weight=h2*std::max(0.,c2+mu*nu)*Reflectance(mu,index);
            result.reflection+=weight*std::norm(Fourier(p.x,p.y,Dot(reflectedQ,p.h),Dot(reflectedQ,p.v)))/(wave*wave);
        }
        if(shadowMode=="facets" && c2>0 && h2>0)
            result.shadow+=c2*c2*mu*mu*std::norm(Fourier(p.x,p.y,Dot(shadowQ,p.h),Dot(shadowQ,p.v)))/(wave*wave);
    }
    if(shadowMode!="off" && (h2==0 || shadowMode=="circular"))
    {
        double envelope=Jinc(k*circleRadius*sine);
        result.shadow=c2*c2*envelope*envelope*projectedArea*projectedArea/(wave*wave);
    }
    return result;
}

void Control::PrepareMeans(const std::vector<double>&theta)
{
    if(index==C(1,0)){means.assign(theta.size(),Components{0,0});return;}
    double maxDiameter=0,maxQ=0;
    for(const Patch&p:patches)maxDiameter=std::max(maxDiameter,p.diameter);
    for(double angle:theta)maxQ=std::max(maxQ,4*pi/wave*std::sin(angle/2));
    // Haar self means depend only on the ordered polygon's pair distances.
    // Congruent facets share one boundary integral; point controls and the
    // coherent projected-area mean still retain every physical facet normal.
    std::vector<size_t> representatives,multiplicity;
    const auto congruent=[](const Patch&a,const Patch&b)
    {
        const size_t m=a.x.size();if(m!=b.x.size())return false;
        for(size_t shift=0;shift<m;++shift)for(int sign:{-1,1})
        {
            bool same=true;
            for(size_t i=0;i<m && same;++i)for(size_t j=i+1;j<m;++j)
            {
                const size_t ii=(shift+m+(sign<0 ? m-i:i))%m,jj=(shift+m+(sign<0 ? m-j:j))%m;
                const double da=std::hypot(a.x[i]-a.x[j],a.y[i]-a.y[j]),db=std::hypot(b.x[ii]-b.x[jj],b.y[ii]-b.y[jj]);
                if(std::fabs(da-db)>1e-13*std::max(da,db)){same=false;break;}
            }
            if(same)return true;
        }
        return false;
    };
    for(size_t i=0;i<patches.size();++i)
    {
        size_t g=0;for(;g<representatives.size();++g)if(congruent(patches[i],patches[representatives[g]]))break;
        if(g==representatives.size()){representatives.push_back(i);multiplicity.push_back(1);}else ++multiplicity[g];
    }
    means.assign(theta.size(),Components{0,0});
    std::vector<Components>previous(theta.size(),Components{0,0});
    double lastReflection=0,lastShadow=0,lastAngle=0;
    for(int pass=0;pass<5;++pass)
    {
        int degree=pass==0 ? 128:(pass==1 ? 192:320),nodes=pass==0 ? 64:(pass==1 ? 128:(pass==2 ? 192:(pass==3 ? 384:768)));
        double step=pass==0 ? .06:(pass==1 ? .03:.015);
        Moments moments(index,degree);Basis basis(degree,maxQ*maxDiameter+step*4,step);
        double worst=0;
        for(size_t row=0;row<theta.size();++row)
        {
            double h2,c2;Angles(theta[row],h2,c2);double reflection=0,shadow=0,Q=4*pi/wave*std::sqrt(h2);
            if(Q>0)
            {
                Green green=basis.Combine(moments.At(theta[row]),Q*maxDiameter+step);
                for(size_t g=0;g<representatives.size();++g)
                {
                    const Patch &p=patches[representatives[g]];
                    reflection+=multiplicity[g]*EdgeMean(p.x,p.y,Q,nodes,[&green](double z){return green.Value(z);},[&green](double z){return green.Same(z);})/(wave*wave);
                }
            }
            if(shadowMode!="off")
            {
                if(h2==0 || shadowMode=="circular")
                {double j=Jinc(2*pi/wave*circleRadius*std::sin(theta[row]));shadow=c2*c2*j*j*forwardShadowMean;}
                else if(c2>0)
                {
                    double qs=2*pi/wave*std::sin(theta[row]);
                    for(size_t g=0;g<representatives.size();++g)
                    {
                        const Patch &p=patches[representatives[g]];
                        shadow+=multiplicity[g]*c2*c2*EdgeMean(p.x,p.y,qs,nodes,ShadowGreen,ShadowSame)/(wave*wave);
                    }
                }
            }
            if(pass)
            {
                double re=std::fabs(reflection-previous[row].reflection)/std::max(std::fabs(reflection),1e-100);
                // Shadow self terms can be tiny far from forward scattering.
                // Control their absolute error relative to the FULL control;
                // reflection remains independently checked for no-shadow use.
                double se=std::fabs(shadow-previous[row].shadow)
                    /std::max(std::fabs(reflection)+std::fabs(shadow),1e-100);
                if(std::max(re,se)>worst){lastReflection=re;lastShadow=se;lastAngle=theta[row]*180/pi;}
                worst=std::max(worst,std::max(re,se));
            }
            if(!std::isfinite(reflection+shadow) || reflection<0 || shadow<0)
                throw std::runtime_error("analytic facet mean is nonfinite or negative.\n  Fix: reduce size or remove --analytic-facet-average.");
            means[row]=Components{reflection,shadow};previous[row]=means[row];
        }
        if(pass && worst<1e-4){refinementError=worst;return;}
        if(pass==4)
            throw std::runtime_error("analytic facet mean did not converge within1e-4 at theta="
                +std::to_string(lastAngle)+", reflection change="+std::to_string(lastReflection)
                +", shadow change="+std::to_string(lastShadow)
                +".\n  Fix: reduce particle size/grid or remove --analytic-facet-average.");
    }
}

bool Control::SupportsDomain(double betaSym,double gammaSym)const
{
    bool half=std::fabs(betaSym-pi/2)<1e-10;
    if(!half && std::fabs(betaSym-pi)>1e-10)return false;
    if(!(gammaSym>0) || gammaSym>2*pi+1e-10)return false;
    double periods=2*pi/gammaSym;
    if(std::fabs(periods-std::round(periods))>1e-9 || periods>64)return false;
    Vec center;size_t count=0;double extent=0;
    for(const Face&f:faces)for(Vec v:f.vertices){center.x+=v.x;center.y+=v.y;center.z+=v.z;++count;}
    center.x/=count;center.y/=count;center.z/=count;
    for(const Face&f:faces)for(Vec v:f.vertices)extent=std::max(extent,Length(Sub(v,center)));
    auto invariant=[&](double rotation,bool flipXZ)
    {
        auto transform=[&](Vec v){if(flipXZ){v.x=-v.x;v.z=-v.z;}return Vec(std::cos(rotation)*v.x-std::sin(rotation)*v.y,std::sin(rotation)*v.x+std::cos(rotation)*v.y,v.z);};
        std::vector<bool>used(faces.size(),false);
        for(const Face&a:faces)
        {
            bool found=false;
            for(size_t j=0;j<faces.size();++j)
            {
                const Face&b=faces[j];if(used[j] || a.vertices.size()!=b.vertices.size() || Length(Sub(transform(Unit(a.normal)),Unit(b.normal)))>1e-8)continue;
                bool match=true;
                for(Vec v:a.vertices)
                {
                    Vec target=transform(Sub(v,center));bool vertex=false;
                    for(Vec w:b.vertices)if(Length(Sub(target,Sub(w,center)))<1e-8*extent){vertex=true;break;}
                    if(!vertex){match=false;break;}
                }
                if(match){used[j]=true;found=true;break;}
            }
            if(!found)return false;
        }
        return true;
    };
    return (periods<1.5 || invariant(gammaSym,false)) && (!half || invariant(0,true));
}
}
