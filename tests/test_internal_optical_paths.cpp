#include "scattering/ScatteringConvex.h"
#include "scattering/ScatteringNonConvex.h"
#include "particle/Hexagonal.h"
#include "particle/Bullet.h"
#include "particle/Droxtal.h"
#include "Tracks.h"
#include <array>
#include <cassert>
#include <memory>
#include <iostream>
#include <limits>
#include <cmath>
class Convex:public ScatteringConvex {public:using ScatteringConvex::ScatteringConvex;void General(){m_unitFacetNormals=false;}};
class NonConvex:public ScatteringNonConvex {public:using ScatteringNonConvex::ScatteringNonConvex;void General(){m_unitFacetNormals=false;}};
static size_t checked=0,different=0;static double worst=0;
void Check(Scattering&s,const Beam&b,const std::vector<int>&track,std::array<Point3f,3> points){
 double actual[3];s.ComputeInternalOpticalPaths3(b,points.data(),track,actual);
 for(int j=0;j<3;++j){double ref=s.ComputeInternalOpticalPath(b,points[j],track);assert(std::isfinite(ref)==std::isfinite(actual[j]));if(std::isfinite(ref)){++checked;if(ref!=actual[j])++different;worst=std::max(worst,std::fabs(ref-actual[j])/std::max(1.,std::fabs(ref)));}}
}
int main(int argc,char**argv){assert(argc==2);
 for(int shape=0;shape<4;++shape)for(double absorption:{0.,.001})for(int pose=0;pose<3;++pose){
  std::unique_ptr<Particle> p;
  if(shape==0)p.reset(new Hexagonal(complex(1.31,absorption),35,100));
  else if(shape==1)p.reset(new Bullet(complex(1.31,absorption),35,100,20));
  else if(shape==2)p.reset(new Droxtal(complex(1.31,absorption),.56,1.25,50));
  else{p.reset(new Particle);p->SetFromFile(argv[1]);p->Resize(15.85);p->SetRefractiveIndex(complex(1.31,absorption));}
  const double beta=pose==0?0:pose==1?std::acos(-1.)/2:.73,gamma=pose<2?0:.42;p->Rotate(beta,gamma,.17);
  Light light;light.direction={0,0,-1};light.polarizationBasis={0,1,0};light.direction.d_param=p->LongRadius();
  std::unique_ptr<Scattering>s;int depth=shape==3?6:18;
  if(shape==3)s.reset(new NonConvex(p.get(),&light,true,depth));else s.reset(new Convex(p.get(),&light,true,depth));
  s->SetTrackIdsRequired(true);std::vector<Beam>beams;assert(s->ScatterLight(beta,gamma,beams));
  for(int general=0;general<2;++general){
   if(general){if(shape==3)static_cast<NonConvex*>(s.get())->General();else static_cast<Convex*>(s.get())->General();}
   bool malformed=false;
   for(const auto&b:beams){if(b.nVertices<3||b.lastFacetId==__INT_MAX__||b.nActs==0)continue;std::vector<int>track;Tracks::RecoverTrack(b,p->nFacets,track);std::array<Point3f,3>pts{{b.arr[0],b.arr[1],b.AreaCentroid()}};Check(*s,b,track,pts);
    if(!malformed){malformed=true;for(double x:{std::numeric_limits<double>::infinity(),std::numeric_limits<double>::quiet_NaN()}){auto bad=pts;bad[1].cx=x;Check(*s,b,track,bad);}}
   }
  }
  std::cout<<"shape="<<shape<<" absorption="<<absorption<<" pose="<<pose<<" beams="<<beams.size()<<std::endl;
 }
 std::cout<<"checked="<<checked<<" differing="<<different<<" max_relative="<<worst<<std::endl;return worst<1e-12?0:1;
}
