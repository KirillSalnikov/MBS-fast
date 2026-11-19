#pragma once
#include <vector>
#include <algorithm>
#include "math/compl.hpp"
namespace DepthExtinction {
struct Capture;
inline Capture*& current(){static thread_local Capture* ptr=nullptr;return ptr;}
struct Capture {
    const std::vector<int>& depths;
    std::vector<complex> diag0,diag1;
    std::vector<double>& result;
    Capture* previous;
    explicit Capture(const std::vector<int>& levels,std::vector<double>& values):
        depths(levels),diag0(levels.size()),diag1(levels.size()),result(values),previous(current()){
        current()=levels.empty()?nullptr:this;
    }
    ~Capture(){current()=previous;}
    void add(int depth,const complex& a,const complex& b){
        for(size_t i=0;i<depths.size();++i)if(depth<=depths[i]){diag0[i]+=a;diag1[i]+=b;}
    }
};
}
