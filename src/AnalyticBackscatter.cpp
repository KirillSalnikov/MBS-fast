#include "AnalyticBackscatter.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <utility>

namespace AnalyticBackscatter
{
namespace
{
const double pi = 3.1415926535897932384626433832795;
typedef std::complex<double> C;
Vec Sub(const Vec &a, const Vec &b) { return Vec(a.x-b.x, a.y-b.y, a.z-b.z); }
double Dot(const Vec &a, const Vec &b) { return a.x*b.x+a.y*b.y+a.z*b.z; }
double Length(const Vec &a) { return std::sqrt(Dot(a,a)); }
Vec Cross(const Vec &a, const Vec &b)
{ return Vec(a.y*b.z-a.z*b.y, a.z*b.x-a.x*b.z, a.x*b.y-a.y*b.x); }
Vec Unit(const Vec &a)
{
    const double r = Length(a);
    if (!(r > 0)) throw std::runtime_error("analytic backscatter: zero normal");
    return Vec(a.x/r,a.y/r,a.z/r);
}

void Gauss(int count, std::vector<double> &nodes, std::vector<double> &weights)
{
    nodes.resize(count); weights.resize(count);
    for (int i=0; i<(count+1)/2; ++i)
    {
        double z=std::cos(pi*(i+.75)/(count+.5)), derivative=0;
        for (int iteration=0; iteration<100; ++iteration)
        {
            double p0=1, p1=z;
            for (int l=2; l<=count; ++l)
            { const double next=((2*l-1)*z*p1-(l-1)*p0)/l; p0=p1; p1=next; }
            derivative=count*(z*p1-p0)/(z*z-1);
            const double next=z-p1/derivative;
            if (std::fabs(next-z)<4e-16) { z=next; break; }
            z=next;
        }
        // Recompute the derivative at the final root, including n=1.
        double p0=1, p1=z;
        for (int l=2; l<=count; ++l)
        { const double next=((2*l-1)*z*p1-(l-1)*p0)/l; p0=p1; p1=next; }
        derivative=count*(z*p1-p0)/(z*z-1);
        const double w=2/((1-z*z)*derivative*derivative);
        nodes[i]=-z; nodes[count-1-i]=z;
        weights[i]=weights[count-1-i]=w;
    }
}

std::pair<C,C> Fresnel(double mu, double relativeIndex)
{
    const C g=std::sqrt(C(relativeIndex*relativeIndex-1+mu*mu,0));
    return std::make_pair((mu-g)/(mu+g),
                         (relativeIndex*relativeIndex*mu-g)/(relativeIndex*relativeIndex*mu+g));
}
double Window(double displacement, double width)
{ return std::max(0.,std::min(displacement,width)-std::max(0.,displacement-width)); }
struct Amplitude { double self, path, jacobian; C cross; };
Amplitude PhysicalAmplitude(double beta, double n, double B, double L)
{
    const double s=std::sin(beta), c=std::cos(beta);
    const double gz=std::sqrt(n*n-s*s), gx=std::sqrt(n*n-c*c);
    const auto eb=Fresnel(c,n), ex=Fresnel(s,n);
    const auto rb=Fresnel(gz/n,1/n), rs=Fresnel(s/n,1/n);
    const auto rl=Fresnel(gx/n,1/n), rz=Fresnel(c/n,1/n);
    const C abS=(1.-eb.first*eb.first)*rb.first*rs.first;
    const C abP=(1.-eb.second*eb.second)*rb.second*rs.second;
    const C axS=(1.-ex.first*ex.first)*rl.first*rz.first;
    const C axP=(1.-ex.second*ex.second)*rl.second*rz.second;
    const double pb=c*Window(2*L*s/gz,B), px=s*Window(2*B*c/gx,L);
    Amplitude a;
    a.self=.5*(pb*pb*(std::norm(abS)+std::norm(abP))+
               px*px*(std::norm(axS)+std::norm(axP)));
    a.cross=.5*pb*px*(abS*std::conj(axS)+abP*std::conj(axP));
    a.path=L*(gz-c)-B*(gx-s);
    a.jacobian=L*s*(1-c/gz)+B*c*(1-s/gx);
    return a;
}

std::vector<double> Breakpoints(double n, double B, double L)
{
    std::vector<double> p={0,pi/2};
    for (double v : {B/(2*L), B/L})
    { const double s=n*v/std::sqrt(1+v*v); if (s<1) p.push_back(std::asin(s)); }
    for (double v : {L/(2*B), L/B})
    { const double c=n*v/std::sqrt(1+v*v); if (c<1) p.push_back(std::acos(c)); }
    if (n<std::sqrt(2.))
    { const double t=std::asin(std::sqrt(n*n-1)); p.push_back(t); p.push_back(pi/2-t); }
    std::sort(p.begin(),p.end());
    p.erase(std::unique(p.begin(),p.end()),p.end());
    return p;
}

double TopSelf(double beta, double n, double B, double L, int p, int q)
{
    const double s=std::sin(beta), c=std::cos(beta), g=std::sqrt(n*n-s*s);
    const double width=Window(2*p*L*s/g-2*(q-1)*B,B);
    if (width==0) return 0;
    const auto entry=Fresnel(c,n), longitudinal=Fresnel(g/n,1/n), lateral=Fresnel(s/n,1/n);
    const double ss=std::norm(1.-entry.first*entry.first)
        *std::pow(std::norm(longitudinal.first),2*p-1)*std::pow(std::norm(lateral.first),2*q-1);
    const double pp=std::norm(1.-entry.second*entry.second)
        *std::pow(std::norm(longitudinal.second),2*p-1)*std::pow(std::norm(lateral.second),2*q-1);
    return .5*c*c*width*width*(ss+pp);
}

double HigherSelf(double beta, double n, double B, double L, int order)
{
    double value=0;
    for (int p=1; p<=order; ++p)
        for (int q=1; p+q<=order+1; ++q)
        {
            if (p==1 && q==1) continue;
            value+=TopSelf(beta,n,B,L,p,q)+TopSelf(pi/2-beta,n,L,B,p,q);
        }
    return value;
}

std::vector<double> HigherBreakpoints(double n, double B, double L, int order)
{
    auto points=Breakpoints(n,B,L);
    for (int p=1; p<=order; ++p)
        for (int q=1; p+q<=order+1; ++q)
            for (int t=2*q-2; t<=2*q; ++t)
            {
                const double xb=t*B/(2*p*L), xl=t*L/(2*p*B);
                const double s=n*xb/std::sqrt(1+xb*xb), c=n*xl/std::sqrt(1+xl*xl);
                if (s<1) points.push_back(std::asin(s));
                if (c<1) points.push_back(std::acos(c));
            }
    std::sort(points.begin(),points.end());
    points.erase(std::unique(points.begin(),points.end()),points.end());
    return points;
}

double SineIntegral(double x)
{
    if (x<=2)
    {
        double factorialTerm=x, sum=x;
        for (int j=1; j<100; ++j)
        {
            factorialTerm *= -x*x/((2*j)*(2*j+1.));
            const double term=factorialTerm/(2*j+1.);
            sum+=term;
            if (std::fabs(term)<2e-16*std::max(1.,std::fabs(sum))) break;
        }
        return sum;
    }
    // Continued fraction of the exponential integral; Si(x)=pi/2+Im(...).
    // Avoid quadrature of the rapidly oscillating sinc profile.
    C b(1,x), c(1e100,0), d=1./b, h=d;
    for (int j=1; j<10000; ++j)
    {
        const double a=-double(j)*j;
        b+=2.; d=1./(a*d+b); c=b+a/c;
        const C delta=c*d; h*=delta;
        if (std::abs(delta-1.)<3e-15) return pi/2+(h*std::exp(C(0,-x))).imag();
    }
    throw std::runtime_error("analytic backscatter: Si continued fraction failed");
}
}

double FiniteSincIntegral(double a)
{
    if (!std::isfinite(a) || a<0) throw std::invalid_argument("sinc argument must be finite and nonnegative");
    if (a<1e-3)
    { const double q=a*a; return 2*(1-q/9+2*q*q/225-q*q*q/2205); }
    return 2/a*(SineIntegral(2*a)+(std::cos(2*a)-1)/(2*a));
}

std::vector<double> SphericalBessel(int degree, double x)
{
    if (degree<0 || !std::isfinite(x) || x<0) throw std::invalid_argument("invalid spherical Bessel arguments");
    std::vector<double> j(degree+1,0);
    if (x<1e-8)
    {
        double term=1;
        for (int l=0; l<=degree; ++l)
        { j[l]=term*(1-x*x/(2*(2*l+3.))); term*=x/(2*l+3.); }
        return j;
    }
    const double j0=std::sin(x)/x, j1=j0/x-std::cos(x)/x;
    if (x>degree+10)
    {
        j[0]=j0; if (degree) j[1]=j1;
        for (int l=1; l<degree; ++l) j[l+1]=(2*l+1)*j[l]/x-j[l-1];
        return j;
    }
    const int top=std::max(degree+60,static_cast<int>(x)+60);
    std::vector<double> work(top+2,0); work[top]=1;
    for (int l=top; l>0; --l)
    {
        work[l-1]=(2*l+1)*work[l]/x-work[l+1];
        if (std::fabs(work[l-1])>1e100)
            for (int t=l-1; t<=top+1; ++t) work[t]*=1e-100;
    }
    const double scale=std::fabs(j0)>std::fabs(j1) ? j0/work[0] : j1/work[1];
    for (int l=0; l<=degree; ++l) j[l]=scale*work[l];
    return j;
}

ReturnStrip::ReturnStrip(double index, double b, double l, double h, double wavelength, int degree, int returnOrder)
    : n(index), B(b), L(l), H(h), wave(wavelength), mean(0), leadingMean(0), order(returnOrder)
{
    if (!(n>1) || !std::isfinite(n) || !(B>0) || !(L>0) || !(H>0) || !(wave>0)
        || !std::isfinite(B+L+H+wave) || degree<8 || degree>512 || order<1 || order>8)
        throw std::invalid_argument("analytic backscatter requires real n>1 and positive finite dimensions/wavelength");
    const double size=std::max(B,L), bn=B/size, ln=L/size, k=2*pi/wave;
    const auto breaks=Breakpoints(n,bn,ln);
    std::vector<double> nodes, weights;
    Gauss(std::max(256,2*degree+64),nodes,weights);
    double self=0; C cross=0;
    for (size_t segment=1; segment<breaks.size(); ++segment)
    {
        const double a=breaks[segment-1], b1=breaks[segment];
        if (b1-a<1e-14) continue;
        // A sin^2 map regularizes square-root Fresnel endpoints for self terms.
        for (size_t q=0; q<nodes.size(); ++q)
        {
            const double t=pi/4*(nodes[q]+1);
            const double beta=a+(b1-a)*std::sin(t)*std::sin(t);
            self+=weights[q]*pi/4*(b1-a)*std::sin(2*t)*PhysicalAmplitude(beta,n,bn,ln).self;
        }
        const double va=PhysicalAmplitude(a,n,bn,ln).path, vb=PhysicalAmplitude(b1,n,bn,ln).path;
        if (vb-va<1e-14) continue;
        std::vector<C> coefficients(degree+1,0);
        for (size_t q=0; q<nodes.size(); ++q)
        {
            const double v=.5*(va+vb)+.5*(vb-va)*nodes[q];
            double left=a, right=b1;
            for (int iteration=0; iteration<52; ++iteration)
            {
                const double middle=.5*(left+right);
                if (PhysicalAmplitude(middle,n,bn,ln).path<v) left=middle; else right=middle;
            }
            const Amplitude amp=PhysicalAmplitude(.5*(left+right),n,bn,ln);
            const C f=amp.cross/amp.jacobian;
            double p0=1, p1=nodes[q];
            coefficients[0]+=.5*weights[q]*f;
            coefficients[1]+=1.5*weights[q]*f*p1;
            for (int ell=2; ell<=degree; ++ell)
            {
                const double next=((2*ell-1)*nodes[q]*p1-(ell-1)*p0)/ell;
                coefficients[ell]+=.5*(2*ell+1)*weights[q]*f*next;
                p0=p1; p1=next;
            }
        }
        const auto j=SphericalBessel(degree,k*size*(vb-va));
        C sum=0, power=1;
        for (int ell=0; ell<=degree; ++ell)
        { sum+=coefficients[ell]*power*j[ell]; power*=C(0,1); }
        cross+=(vb-va)*std::exp(C(0,k*size*(va+vb)))*sum;
    }
    if (order>1)
    {
        const auto points=HigherBreakpoints(n,bn,ln,order);
        for (size_t segment=1; segment<points.size(); ++segment)
        {
            const double a=points[segment-1], b1=points[segment];
            if (b1-a<1e-14) continue;
            for (size_t q=0; q<nodes.size(); ++q)
            {
                const double t=pi/4*(nodes[q]+1);
                const double beta=a+(b1-a)*std::sin(t)*std::sin(t);
                self+=weights[q]*pi/4*(b1-a)*std::sin(2*t)*HigherSelf(beta,n,bn,ln,order);
            }
        }
    }
    leadingMean=H*size*size/(8*pi*wave)*(self+2*cross.real());
    mean=leadingMean*FiniteSincIntegral(k*H)/(pi/(k*H));
    if (!std::isfinite(mean) || mean<0)
        throw std::runtime_error("analytic backscatter: invalid integrated physical return profile");
}

double ReturnStrip::Evaluate(double sx, double sz, double se) const
{
    if (sx<=0 || sz<=0) return 0;
    const Amplitude a=PhysicalAmplitude(std::atan2(sx,sz),n,B,L);
    const double k=2*pi/wave, x=k*H*se;
    const double sinc=std::fabs(x)<1e-8 ? 1-x*x/6 : std::sin(x)/x;
    const double value=a.self+2*(a.cross*std::exp(C(0,2*k*a.path))).real()
        +HigherSelf(std::atan2(sx,sz),n,B,L,order);
    return H*H/(wave*wave)*std::max(0.,value)*sinc*sinc;
}

std::vector<Strip> FindReturnStrips(const std::vector<Face> &faces)
{
    std::vector<Vec> normals; double extent=0;
    for (const auto &f:faces)
    {
        if (f.vertices.size()<3) throw std::runtime_error("analytic backscatter: degenerate facet");
        normals.push_back(Unit(f.normal));
        for (const Vec &v:f.vertices)
            extent=std::max(extent,Length(Sub(v,faces.front().vertices.front())));
    }
    if (!(extent>0)) throw std::runtime_error("analytic backscatter: zero particle extent");
    const double tolerance=1e-9*extent;
    std::vector<double> widths(faces.size(),0);
    for (size_t i=0; i<faces.size(); ++i)
    {
        const double plane=Dot(normals[i],faces[i].vertices[0]);
        for (size_t j=0; j<faces.size(); ++j)
        {
            for (const Vec &v:faces[j].vertices)
                if (Dot(normals[i],v)>plane+tolerance)
                    throw std::runtime_error("--analytic-backscatter requires convex outward-oriented facets.\n  Fix: use a convex prism or remove the flag.");
            if (Dot(normals[i],normals[j]) < -1+1e-10)
                widths[i]=std::max(widths[i],plane-Dot(normals[i],faces[j].vertices[0]));
        }
    }
    std::vector<Strip> result;
    for (size_t i=0; i<faces.size(); ++i)
        for (size_t j=i+1; j<faces.size(); ++j)
        {
            if (std::fabs(Dot(normals[i],normals[j]))>1e-10 || widths[i]<=tolerance || widths[j]<=tolerance) continue;
            for (size_t ei=0; ei<faces[i].vertices.size(); ++ei)
            {
                const Vec &a=faces[i].vertices[ei], &b=faces[i].vertices[(ei+1)%faces[i].vertices.size()];
                bool shared=false;
                for (size_t ej=0; ej<faces[j].vertices.size(); ++ej)
                {
                    const Vec &c=faces[j].vertices[ej], &d=faces[j].vertices[(ej+1)%faces[j].vertices.size()];
                    if ((Length(Sub(a,c))<tolerance && Length(Sub(b,d))<tolerance) ||
                        (Length(Sub(a,d))<tolerance && Length(Sub(b,c))<tolerance)) shared=true;
                }
                if (!shared || Length(Sub(a,b))<=tolerance) continue;
                Strip strip;
                const size_t x=widths[i]<=widths[j] ? i : j, z=x==i ? j : i;
                strip.normalX=normals[x]; strip.normalZ=normals[z];
                strip.edge=Unit(Cross(strip.normalX,strip.normalZ));
                strip.B=widths[x]; strip.L=widths[z]; strip.H=Length(Sub(a,b));
                result.push_back(strip);
            }
        }
    if (result.empty())
        throw std::runtime_error("--analytic-backscatter found no right-angle edges with parallel opposing faces.\n  Fix: use a rectangular/hexagonal prism or remove the flag.");
    return result;
}

double ReturnSelfIntensity(double beta, double index, double B, double L,
                           double H, double wavelength, int p, int q)
{
    if (p<1 || q<1) throw std::invalid_argument("return image orders must be positive");
    return H*H/(wavelength*wavelength)*TopSelf(beta,index,B,L,p,q);
}

ReturnControl::ReturnControl(const std::vector<Face> &faces, double index, double wavelength, int returnOrder)
    : strips(FindReturnStrips(faces)), mean(0), refinementError(0)
{
    std::vector<Strip> unique;
    for (const Strip &s:strips)
    {
        size_t match=unique.size();
        for (size_t i=0; i<unique.size(); ++i)
            if (std::fabs(s.B/unique[i].B-1)<1e-10 && std::fabs(s.L/unique[i].L-1)<1e-10 && std::fabs(s.H/unique[i].H-1)<1e-10)
            { match=i; break; }
        if (match==unique.size())
        {
            unique.push_back(s);
            const ReturnStrip coarse(index,s.B,s.L,s.H,wavelength,96,returnOrder);
            ReturnStrip fine(index,s.B,s.L,s.H,wavelength,160,returnOrder);
            double error=std::fabs(fine.Mean()-coarse.Mean())/std::max(fine.Mean(),1e-300);
            if (error>1e-6)
            {
                const ReturnStrip finer(index,s.B,s.L,s.H,wavelength,320,returnOrder);
                error=std::fabs(finer.Mean()-fine.Mean())/std::max(finer.Mean(),1e-300);
                fine=finer;
            }
            if (error>1e-6)
                throw std::runtime_error("analytic backscatter moments did not converge within 1e-6.\n  Fix: remove the flag for this size/index; the ordinary solver remains available.");
            refinementError=std::max(refinementError,error);
            models.push_back(fine);
        }
        modelIndex.push_back(static_cast<int>(match));
        mean+=models[match].Mean();
    }
}

double ReturnControl::Evaluate(const Vec &source) const
{
    double value=0;
    for (size_t i=0; i<strips.size(); ++i)
        value+=models[modelIndex[i]].Evaluate(Dot(source,strips[i].normalX),
                                             Dot(source,strips[i].normalZ),Dot(source,strips[i].edge));
    return value;
}

bool ReturnControl::SupportsDomain(double betaSym, double gammaSym) const
{
    const bool halfBeta=std::fabs(betaSym-pi/2)<1e-10;
    if (!halfBeta && std::fabs(betaSym-pi)>1e-10) return false;
    if (!(gammaSym>0) || gammaSym>2*pi+1e-10) return false;
    const double periods=2*pi/gammaSym;
    if (std::fabs(periods-std::round(periods))>1e-9 || periods>64) return false;
    const auto invariant = [this](double rotation, bool mirrorZ)
    {
        const auto transform = [rotation,mirrorZ](const Vec &v)
        {
            return Vec(std::cos(rotation)*v.x-std::sin(rotation)*v.y,
                       std::sin(rotation)*v.x+std::cos(rotation)*v.y,
                       mirrorZ ? -v.z : v.z);
        };
        std::vector<bool> used(strips.size(),false);
        for (const Strip &a:strips)
        {
            const Vec x=transform(a.normalX), z=transform(a.normalZ);
            bool found=false;
            for (size_t i=0; i<strips.size(); ++i)
            {
                if (used[i]) continue;
                const Strip &b=strips[i];
                const bool direct=Length(Sub(x,b.normalX))<1e-9 && Length(Sub(z,b.normalZ))<1e-9
                    && std::fabs(a.B/b.B-1)<1e-9 && std::fabs(a.L/b.L-1)<1e-9;
                const bool swapped=Length(Sub(x,b.normalZ))<1e-9 && Length(Sub(z,b.normalX))<1e-9
                    && std::fabs(a.B/b.L-1)<1e-9 && std::fabs(a.L/b.B-1)<1e-9;
                if ((direct || swapped) && std::fabs(a.H/b.H-1)<1e-9)
                { used[i]=true; found=true; break; }
            }
            if (!found) return false;
        }
        return true;
    };
    return (periods<1.5 || invariant(gammaSym,false)) && (!halfBeta || invariant(0,true));
}
}
