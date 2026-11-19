#pragma once
#include <algorithm>
#include <cmath>
#include <limits>
#include "../Beam.h"
#include "../geometry/Facet.h"
#include "../geometry/geometry_lib.h"

// A projected ray keeps u = x[j] - d[j]/d[k] * x[k] constant.
// Compare the incident polygon's u interval with a box containing the facet.
// Exact clipping still handles near contacts, grazing rays and acute facets.
namespace FacetRayBounds {
struct Box {
    Point3f center{0,0,0}, half{0,0,0};
    double margin=0;
    bool usable=false;
};
inline Box Limits(const Facet& facet) {
    Box box;
    if(facet.nVertices<3 || facet.nVertices>MAX_VERTEX_NUM)return box;
    Point3f lo=facet.arr[0],hi=lo;
    double minSin=1.,scale=0.;
    for(int i=0;i<facet.nVertices;++i){
        const auto& v=facet.arr[i];
        for(int j=0;j<3;++j){
            if(!std::isfinite(v.coordinates[j]))return box;
            lo.coordinates[j]=std::min(lo.coordinates[j],v.coordinates[j]);
            hi.coordinates[j]=std::max(hi.coordinates[j],v.coordinates[j]);
            scale=std::max(scale,std::abs(v.coordinates[j]));
        }
        const auto a=v-facet.arr[i?i-1:facet.nVertices-1];
        const auto b=facet.arr[(i+1)%facet.nVertices]-v;
        double denom=Norm(a)*Norm(b);if(!(denom>0) || !std::isfinite(denom))return box;
        minSin=std::min(minSin,std::sqrt(Norm(CrossProduct(a,b))/denom));
    }
    if(!(minSin>.01) || !(scale>0))return box;
    box.center=(hi+lo)*.5;box.half=(hi-lo)*.5;box.margin=1e-5*scale;box.usable=true;
    return box;
}
struct Footprint {
    int k,u,v;
    double du,dv,loU,hiU,loV,hiV,scale;
    explicit Footprint(const Beam& beam) {
        const auto& d=beam.direction;
        k=std::fabs(d.cx)>std::fabs(d.cy)?0:1;
        if(std::fabs(d.cz)>std::fabs(d.coordinates[k]))k=2;
        u=(k+1)%3;v=(k+2)%3;
        du=d.coordinates[u]/d.coordinates[k];dv=d.coordinates[v]/d.coordinates[k];
        loU=loV=std::numeric_limits<double>::infinity();hiU=hiV=-loU;scale=0.;
        for(int i=0;i<beam.nVertices;++i){
            const auto& x=beam.arr[i];const double a=x.coordinates[u]-du*x.coordinates[k];
            const double b=x.coordinates[v]-dv*x.coordinates[k];
            loU=std::min(loU,a);hiU=std::max(hiU,a);
            loV=std::min(loV,b);hiV=std::max(hiV,b);
            scale=std::max(scale,std::fabs(x.coordinates[k]));
        }
        scale=1e-5*std::max(scale,std::max(std::max(std::fabs(loU),std::fabs(hiU)),std::max(std::fabs(loV),std::fabs(hiV))));
    }
};
}
