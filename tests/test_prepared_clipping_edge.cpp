#include "geometry/Intersection.h"
#include <cassert>
#include <iostream>
#include <random>

int main(){
    std::mt19937_64 rng(318751);std::uniform_real_distribution<double> random(-1,1);
    size_t checked=0;
    for(double scale:{1e-9,1e-3,1.,1e3,1e9})for(int trial=0;trial<20000;++trial){
        auto point=[&](){return Point3f(scale*random(rng),scale*random(rng),scale*random(rng));};
        const Point3f p1=point(),p2=trial%17?p1+point():p1;
        const Point3f normal=trial%29?Point3f(random(rng),random(rng),random(rng)):Point3f(0,0,0);
        const PreparedClippingEdge prepared(p1,p2,normal);
        for(const Point3f& p:{point(),p1,p2,p1+(p2-p1)*.5,p1+point()*1e-14}){
            assert(prepared.Contains(p)==is_inside_i(p,p1,p2,normal));++checked;
        }
        const Point3f a=point(),b=trial%13?point():a+(p2-p1);
        bool oldOk=false,newOk=false;const Point3f old=intersect_i(a,b,p1,p2,normal,oldOk),now=prepared.Intersect(a,b,newOk);
        assert(oldOk==newOk);
        if(oldOk)assert(old.cx==now.cx&&old.cy==now.cy&&old.cz==now.cz);
    }
    std::cout<<"PASS "<<checked<<" side tests and 100000 exact line intersections across scales, degeneracy, near-edge and parallel cases\n";
}
