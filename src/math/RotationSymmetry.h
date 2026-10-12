#pragma once
#include <array>
#include <vector>
#include <cmath>
#include <algorithm>

namespace RotationSymmetry {
using Point=std::array<double,3>;
using Face=std::vector<Point>;
// Check a bijection of whole oriented polygons about the body origin.
// A translated or differently aligned symmetric particle is deliberately
// left unreduced. Comparing cyclic vertex order also preserves normals.
inline bool Preserves(const std::vector<Face>& faces,double angle,bool aroundX){
    double extent=0;for(const auto& face:faces)for(const auto& p:face)
        extent=std::max(extent,std::hypot(std::hypot(p[0],p[1]),p[2]));
    if(!(extent>0)||faces.empty())return false;
    const double tolerance=extent*1e-11,c=std::cos(angle),s=std::sin(angle);
    std::vector<bool> used(faces.size(),false);
    for(const auto& face:faces){bool found=false;if(face.size()<3)return false;
        Face transformed;for(const auto& p:face)transformed.push_back(aroundX?
            Point{{p[0],c*p[1]-s*p[2],s*p[1]+c*p[2]}}:
            Point{{c*p[0]-s*p[1],s*p[0]+c*p[1],p[2]}});
        for(size_t f=0;f<faces.size()&&!found;++f){const auto& target=faces[f];if(used[f]||target.size()!=face.size())continue;
            for(size_t shift=0;shift<face.size()&&!found;++shift){bool match=true;
                for(size_t v=0;v<face.size();++v){const auto& a=transformed[v];const auto& b=target[(v+shift)%face.size()];
                    if(std::hypot(std::hypot(a[0]-b[0],a[1]-b[1]),a[2]-b[2])>tolerance){match=false;break;}}
                if(match){used[f]=true;found=true;}}
        }if(!found)return false;
    }return true;
}
inline std::array<int,2> Detect(const std::vector<Face>& faces){
    const double pi=std::acos(-1.);std::array<int,2> factors{{1,1}};
    // Bound the conservative detector's cost for very large arbitrary meshes.
    if(faces.size()>4096)return factors;
    for(int q:{12,8,6,4,3,2})if(Preserves(faces,2*pi/q,false)){factors[1]=q;break;}
    if(Preserves(faces,pi,true))factors[0]=2;
    return factors;
}
}
