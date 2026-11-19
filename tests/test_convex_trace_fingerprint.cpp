#include "scattering/ScatteringConvex.h"
#include "particle/Hexagonal.h"
#include "particle/Bullet.h"
#include "particle/Droxtal.h"
#ifdef NDEBUG
#error "Trace validation/benchmark requires assertions: compile without -DNDEBUG"
#endif
#include <memory>
#include <chrono>
#include <cstdint>
#include <cstring>
static uint64_t checksum=1469598103934665603ULL;
static void Hash(double v){uint64_t bits;std::memcpy(&bits,&v,sizeof(bits));checksum=(checksum^bits)*1099511628211ULL;}
#include <cassert>
#include <cstdlib>
#include <cmath>
#include <iostream>
#include <random>
#include <vector>
static std::vector<Beam> Trace(bool fast,double scale,double beta,double gamma,double alpha,double absorption,int shape=0){
 (void)fast; // Trace fingerprint can be compared across release binaries.
 std::unique_ptr<Particle> particle;
 if(shape==0)particle.reset(new Hexagonal(complex(1.3116,absorption),34.89*scale,100*scale));
 else if(shape==1)particle.reset(new Bullet(complex(1.3116,absorption),34.89*scale,100*scale,20*scale));
 else particle.reset(new Droxtal(complex(1.3116,absorption),.56,1.25,50*scale));
 particle->Rotate(beta,gamma,alpha);
 Light light;light.direction=Point3f(0,0,-1);light.polarizationBasis=Point3f(0,1,0);
 light.direction.d_param=particle->LongRadius();
 ScatteringConvex scatter(particle.get(),&light,true,18);std::vector<Beam> out;
 if(const char* cutoff=std::getenv("MBS_TEST_TRACE_CUTOFF")){
 const double value=std::strtod(cutoff,nullptr);scatter.m_traceCutoffJRel=value;
 scatter.m_traceCutoffAreaRel=value;scatter.m_traceCutoffImportanceRel=value;
 }
 assert(scatter.ScatterLight(beta,gamma,out));return out;
}
// Match the production preparation loop: reuse one rotated particle, tracer
// and output capacity across orientations. The fresh-allocation benchmark is
// retained separately because its page faults can dominate sampling profiles.
struct ReusedTrace {
 std::unique_ptr<Particle> particle;
 Light light;
 std::unique_ptr<ScatteringConvex> scatter;
 std::vector<Beam> out;
 ReusedTrace(int shape,double scale,double absorption){
  if(shape==0)particle.reset(new Hexagonal(complex(1.3116,absorption),34.89*scale,100*scale));
  else if(shape==1)particle.reset(new Bullet(complex(1.3116,absorption),34.89*scale,100*scale,20*scale));
  else particle.reset(new Droxtal(complex(1.3116,absorption),.56,1.25,50*scale));
  light.direction=Point3f(0,0,-1);light.polarizationBasis=Point3f(0,1,0);light.direction.d_param=particle->LongRadius();
  scatter.reset(new ScatteringConvex(particle.get(),&light,true,18));
  if(const char* cutoff=std::getenv("MBS_TEST_TRACE_CUTOFF")){
   const double value=std::strtod(cutoff,nullptr);scatter->m_traceCutoffJRel=value;
   scatter->m_traceCutoffAreaRel=value;scatter->m_traceCutoffImportanceRel=value;
  }
 }
 const std::vector<Beam>& Run(double beta,double gamma,double alpha){
  particle->Rotate(beta,gamma,alpha);out.clear();assert(scatter->ScatterLight(0,0,out));return out;
 }
};
static void SamePoint(const Point3f& a,const Point3f& b){assert(a.cx==b.cx&&a.cy==b.cy&&a.cz==b.cz);Hash(a.cx);Hash(a.cy);Hash(a.cz);}
int main(int argc,char** argv){
 if(argc>1){
  const bool reuse=std::strcmp(argv[1],"bench-reuse")==0;
  if(!reuse && std::strcmp(argv[1],"bench")!=0){std::cerr<<"usage: test_convex_trace_fingerprint [bench|bench-reuse]\n";return 2;}
  for(int shape=0;shape<3;++shape)for(int pass=0;pass<3;++pass){
   std::unique_ptr<ReusedTrace> steady;if(reuse)steady.reset(new ReusedTrace(shape,1.,.001));
   const auto start=std::chrono::steady_clock::now();size_t count=0;
   for(int i=0;i<192;++i)count+=reuse?steady->Run(.01+i*.016,.02+i*.061,.1+i*.043).size():Trace(pass==1,1.,.01+i*.016,.02+i*.061,.1+i*.043,.001,shape).size();
   std::cout<<"shape="<<shape<<" pass="<<pass<<" reuse="<<reuse<<" seconds="<<std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count()<<" beams="<<count<<std::endl;
  }return 0;
 }
 std::mt19937_64 rng(49081);std::uniform_real_distribution<double> u(0,1);size_t beams=0;
 for(int shape=0;shape<3;++shape)for(double scale:{.001,1.,1000.})for(double absorption:{0.,.001}){
  ReusedTrace steady(shape,scale,absorption);
  for(int pose=0;pose<24;++pose){
  const double beta=pose==0?0:pose==1?std::acos(-1.)/2:u(rng)*std::acos(-1.);
  const double gamma=pose<2?0:u(rng)*2*std::acos(-1.),alpha=pose<2?0:u(rng)*2*std::acos(-1.);
  setenv("MBS_CONVEX_FACET_BOUNDS","0",1);
  const auto a=Trace(false,scale,beta,gamma,alpha,absorption,shape);
  setenv("MBS_CONVEX_FACET_BOUNDS","1",1);
  const auto& b=steady.Run(beta,gamma,alpha);
  assert(a.size()==b.size());beams+=a.size();
  for(size_t k=0;k<a.size();++k){const auto& x=a[k];const auto& y=b[k];
   assert(x.nVertices==y.nVertices&&x.nActs==y.nActs&&x.lastFacetId==y.lastFacetId&&x.id==y.id&&x.opticalPath==y.opticalPath);
   Hash(x.opticalPath);Hash(x.nActs);Hash(x.lastFacetId);SamePoint(x.direction,y.direction);SamePoint(x.polarizationBasis,y.polarizationBasis);
   for(int v=0;v<x.nVertices;++v)SamePoint(x.arr[v],y.arr[v]);
   const complex xx[]={x.J.m11,x.J.m12,x.J.m21,x.J.m22},yy[]={y.J.m11,y.J.m12,y.J.m21,y.J.m22};
   for(int j=0;j<4;++j){assert(real(xx[j])==real(yy[j])&&imag(xx[j])==imag(yy[j]));Hash(real(xx[j]));Hash(imag(xx[j]));}
  }
 }
 }
 std::cout<<"PASS 432 fresh-versus-reused traces (column, bullet, droxtal): "<<beams<<" beams, exact vertices/Jones/path/order across scales, absorption and boundary orientations; checksum="<<checksum<<"\n";
}
