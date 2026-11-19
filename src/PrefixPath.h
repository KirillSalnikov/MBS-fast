#pragma once
#include "Beam.h"
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <vector>

// The trace and preparation must use the same per-orientation output vector.
// Copies, shadow rays, resident GPU traces and nonconvex traces use normal recovery.
namespace PrefixPath {
struct Step {
    Point3f ray, normal;
    double reciprocal;
    int parent;
};
inline const std::vector<int>& deferredTrack(){static const std::vector<int> empty;return empty;}
struct Scope;
inline Scope*& current() { static thread_local Scope* active=nullptr; return active; }
struct Scope {
    std::vector<Beam>& output;
    Scope* previous;
    std::vector<Step> steps;
    std::vector<int> emitted;
    bool recording=false;
    unsigned long long hits=0;
    explicit Scope(std::vector<Beam>& beams) :output(beams),previous(current()) {
        static const bool enabled=[](){ const char* v=std::getenv("MBS_CONVEX_PATH_PREFIX");return v && std::strcmp(v,"1")==0; }();
        current()=enabled?this:nullptr;
    }
    ~Scope(){current()=previous;}
    int extend(int parent,const Point3f& ray,const Point3f& normal) {
        const double dot=DotProduct(ray,normal);
        if(!recording || steps.size()>=65536 || parent==-2 || !std::isfinite(dot) || std::fabs(dot)<1e-5)return -2;
        steps.push_back({ray,normal,1./dot,parent});
        return static_cast<int>(steps.size()-1);
    }
    void bindExit(int link,const Point3f& reverse) {
        if(link>=0){
            const Point3f mismatch=reverse+steps[link].ray;
            if(!(Norm(mismatch)<1e-24))link=-2;
        }
        if(output.size()>262144){recording=false;return;}
        emitted.resize(output.size(),-2);
        emitted.back()=link;
    }
    int lookup(const Beam& beam) const {
        if(!recording || emitted.size()!=output.size() || output.empty())return -2;
        const uintptr_t addr=reinterpret_cast<uintptr_t>(&beam),start=reinterpret_cast<uintptr_t>(output.data());
        if(addr<start || addr-start>=output.size()*sizeof(Beam) || (addr-start)%sizeof(Beam)!=0)return -2;
        return emitted[(addr-start)/sizeof(Beam)];
    }
    template<unsigned N> bool evaluate(int node,const Point3f* points,double* lengths) {
        if(node<0)return false;
        Point3f xyz[N];
        for(unsigned j=0;j<N;++j){xyz[j]=points[j];lengths[j]=0;}
        for(;node>=0;node=steps[node].parent){
            const Step& s=steps[node];
            for(unsigned j=0;j<N;++j){
                const double t=(DotProduct(xyz[j],s.normal)+s.normal.d_param)*s.reciprocal;
                if(!std::isfinite(t))return false;
                lengths[j]+=std::fabs(t);
                xyz[j]=Point3f(xyz[j].cx-s.ray.cx*t,xyz[j].cy-s.ray.cy*t,xyz[j].cz-s.ray.cz*t);
            }
        }
        hits+=N;
        return true;
    }
};
}
