#include "AnalyticFacetAverage.h"
#include "Sobol.h"
#include "HandlerPO_fast.h"
#include <cmath>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <string>
#include <algorithm>
#ifdef MBS_GPU_DIRECTION_PROBE
#include "cuda/GpuAnalyticFacet.h"
#endif

using AnalyticBackscatter::Vec;
Vec cross(Vec a,Vec b){return Vec(a.y*b.z-a.z*b.y,a.z*b.x-a.x*b.z,a.x*b.y-a.y*b.x);}
AnalyticBackscatter::Face Polygon(std::string shape,double scale)
{
    AnalyticBackscatter::Face f;f.normal=Vec(0,0,1);
    if(shape=="rectangle")f.vertices={Vec(0,0,0),Vec(scale,0,0),Vec(scale,.6*scale,0),Vec(0,.6*scale,0)};
    else if(shape=="triangle")f.vertices={Vec(0,0,0),Vec(scale,0,0),Vec(.2*scale,.7*scale,0)};
    else for(int i=0;i<6;++i){double phi=i*std::acos(-1.)/3;f.vertices.push_back(Vec(5.0510479784*scale*std::cos(phi),5.0510479784*scale*std::sin(phi),0));}
    return f;
}
int main(int argc,char**argv)
{
    const double pi=std::acos(-1.);std::cout<<std::setprecision(17);
    if(argc==2 && std::string(argv[1])=="frame")
    {
        try
        {
        Vec u,x;double theta,phi;
        while(std::cin>>u.x>>u.y>>u.z>>x.x>>x.y>>x.z>>theta>>phi)
        {
            auto frame=AnalyticFacetAverage::SourceFrame(u,x);
            const Vec v=AnalyticFacetAverage::Observer(frame,theta,phi);
            for(Vec axis:frame)std::cout<<axis.x<<' '<<axis.y<<' '<<axis.z<<' ';
            std::cout<<v.x<<' '<<v.y<<' '<<v.z<<'\n';
        }
        return 0;
        }
        catch(const std::exception &error){std::cerr<<error.what()<<'\n';return 2;}
    }
    if(argc==12 && (std::string(argv[1])=="directions" || std::string(argv[1])=="gpu-directions"))
    {
        const double x=std::atof(argv[8]),y=std::atof(argv[9]),z=std::atof(argv[10]),w=std::atof(argv[11]);
        const double R[3][3]={{1-2*(y*y+z*z),2*(x*y-w*z),2*(x*z+w*y)},
            {2*(x*y+w*z),1-2*(x*x+z*z),2*(y*z-w*x)},
            {2*(x*z-w*y),2*(y*z+w*x),1-2*(x*x+y*y)}};
        const auto rotate=[&](Vec p){return Vec(R[0][0]*p.x+R[0][1]*p.y+R[0][2]*p.z,
            R[1][0]*p.x+R[1][1]*p.y+R[1][2]*p.z,R[2][0]*p.x+R[2][1]*p.y+R[2][2]*p.z);};
        auto face=Polygon(argv[2],std::atof(argv[3])),rotated=face;
        rotated.normal=rotate(face.normal);for(Vec &p:rotated.vertices)p=rotate(p);
        const double theta=std::atof(argv[7])*pi/180;
#ifndef MBS_GPU_DIRECTION_PROBE
        (void)theta;
#endif
        const std::complex<double> index(std::atof(argv[4]),std::atof(argv[5]));
        AnalyticFacetAverage::Control a({face},index,.532,{0},argv[6]),b({rotated},index,.532,{0},argv[6]);
        Vec u,v;
        if(std::string(argv[1])=="directions")
        {
            while(std::cin>>u.x>>u.y>>u.z>>v.x>>v.y>>v.z)
            {
                auto p=a.Evaluate(u,v),q=b.Evaluate(rotate(u),rotate(v));
                std::cout<<p.reflection<<' '<<p.shadow<<' '<<q.reflection<<' '<<q.shadow<<'\n';
            }
            return 0;
        }
#ifdef MBS_GPU_DIRECTION_PROBE
        std::vector<std::array<Vec,3>> frames,rotatedFrames;
        while(std::cin>>u.x>>u.y>>u.z>>v.x>>v.y>>v.z)
        {
            auto frame=AnalyticFacetAverage::SourceFrame(u,v);frames.push_back(frame);
            for(Vec &axis:frame)axis=rotate(axis);rotatedFrames.push_back(frame);
        }
        std::vector<double> weights(frames.size(),1),r,s,rr,ss;
        std::vector<AnalyticFacetAverage::Components> samples,rotatedSamples;
        if(!EvaluateAnalyticFacetGpu(a.GpuModel(),frames,weights,{theta},1,r,s,&samples)
            || !EvaluateAnalyticFacetGpu(b.GpuModel(),rotatedFrames,weights,{theta},1,rr,ss,&rotatedSamples))return 3;
        for(size_t i=0;i<frames.size();++i)
        {
            auto p=a.Evaluate(frames[i][2],AnalyticFacetAverage::Observer(frames[i],theta,0));
            std::cout<<p.reflection<<' '<<p.shadow<<' '<<samples[i].reflection<<' '<<samples[i].shadow
                <<' '<<rotatedSamples[i].reflection<<' '<<rotatedSamples[i].shadow<<'\n';
        }
        return 0;
#else
        return 3;
#endif
    }
    if(argc==5 && std::string(argv[1])=="sobol")
    {
        uint32_t seed=std::strtoul(argv[2],nullptr,10);Sobol3D s(seed,std::atoi(argv[3]));Sobol2D old(seed);
        for(int i=0;i<std::atoi(argv[4]);++i){double x,y,z,a,b;s.next(x,y,z);old.next(a,b);std::cout<<x<<' '<<y<<' '<<z<<' '<<a<<' '<<b<<'\n';}
        return 0;
    }
    if(argc==6 && std::string(argv[1])=="moments")
    {
        auto coefficients=AnalyticFacetAverage::Control::ReflectionMoments(std::complex<double>(std::atof(argv[2]),std::atof(argv[3])),std::atoi(argv[5]),std::atof(argv[4])*pi/180);
        for(double a:coefficients)std::cout<<a<<'\n';
        return 0;
    }
    if(argc==5 && std::string(argv[1])=="groups")
    {
        const double theta=std::atof(argv[4])*pi/180;
        auto rectangle=Polygon("rectangle",1),translated=rectangle,tilted=rectangle,triangle=Polygon("triangle",1.3);
        std::rotate(translated.vertices.begin(),translated.vertices.begin()+1,translated.vertices.end());
        for(Vec &p:translated.vertices)p=Vec(p.x+3,p.y-2,p.z+1);
        const double c=std::cos(.37),s=std::sin(.37);
        tilted.normal=Vec(s,0,c);for(Vec &p:tilted.vertices)p=Vec(c*p.x+s*p.z,p.y,-s*p.x+c*p.z);
        std::reverse(tilted.vertices.begin(),tilted.vertices.end());
        const std::complex<double> ri(std::atof(argv[2]),std::atof(argv[3]));
        AnalyticFacetAverage::Control grouped({rectangle,translated,tilted,triangle},ri,.532,{theta});
        AnalyticFacetAverage::Control a({rectangle},ri,.532,{theta}),b({triangle},ri,.532,{theta});
        auto g=grouped.Mean(0),x=a.Mean(0),y=b.Mean(0);
        std::cout<<g.reflection<<' '<<g.shadow<<' '<<3*x.reflection+y.reflection<<' '<<3*x.shadow+y.shadow<<'\n';
        return 0;
    }
    if(argc>=8 && std::string(argv[1])=="mean")
    {
        std::vector<double>theta;for(int i=7;i<argc;++i)theta.push_back(std::atof(argv[i])*pi/180);
        auto f=Polygon(argv[2],std::atof(argv[3]));
        AnalyticFacetAverage::Control control({f},std::complex<double>(std::atof(argv[4]),std::atof(argv[5])),.532,theta,argv[6]);
        for(size_t i=0;i<theta.size();++i){auto value=control.Mean(i);std::cout<<theta[i]*180/pi<<' '<<value.reflection<<' '<<value.shadow<<' '<<control.RefinementError()<<'\n';}
        return 0;
    }
    if(argc==8 && std::string(argv[1])=="point")
    {
        double theta=std::atof(argv[6])*pi/180;
        AnalyticFacetAverage::Control control({Polygon(argv[2],std::atof(argv[3]))},std::complex<double>(std::atof(argv[4]),std::atof(argv[5])),.532,{theta},argv[7]);
        Vec u,v;while(std::cin>>u.x>>u.y>>u.z>>v.x>>v.y>>v.z){auto value=control.Evaluate(u,v,theta);std::cout<<value.reflection<<' '<<value.shadow<<'\n';}
        return 0;
    }
    if(argc==2 && std::string(argv[1])=="projection")
    {
        double mu,phi,theta,az;
        while(std::cin>>mu>>phi>>theta>>az)
        {
            double root=std::sqrt(1-mu*mu);Vec n(root*std::cos(phi),root*std::sin(phi),mu),d(2*mu*n.x,2*mu*n.y,2*mu*mu-1),s(-std::sin(phi),std::cos(phi),0),p=cross(d,s);
            Vec nt=cross(n,p),np=cross(n,s),a=cross(n,cross(d,p)),b=cross(n,cross(d,s));
            Vec obs(std::sin(theta)*std::cos(az),std::sin(theta)*std::sin(az),-std::cos(theta)),vf(-std::sin(az),std::cos(az),0);
            double r00,r01,r10,r11;rotate_jones_inline(nt.x,nt.y,nt.z,np.x,np.y,np.z,a.x,a.y,a.z,b.x,b.y,b.z,vf.x,vf.y,vf.z,obs.x,obs.y,obs.z,r00,r01,r10,r11);
            std::cout<<r00*r00+r10*r10<<' '<<r01*r01+r11*r11<<'\n';
        }
        return 0;
    }
    std::cerr<<"Unknown probe mode\n";return 2;
}
