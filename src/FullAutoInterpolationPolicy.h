#pragma once

#include "FullAuto.h"
#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace FullAuto {

// A fixed linear operator is shared by every replica, matrix component and
// paired residual. Choosing slopes from replica values (e.g. PCHIP) would not
// commute with the sample mean and would change the control-variate estimator.
inline bool UseLinearPOInterpolation(bool go,bool coherent) {
    return !go&&!coherent;
}

inline const char* POOutputInterpolationPolicy(bool go,bool coherent) {
    return go?"go-angular-bins-v1":coherent?"po-local-cubic-v1":"po-incoherent-piecewise-linear-v1";
}

inline void ValidatePiecewiseLinearGrid(const std::vector<double>& x) {
    if(x.size()<2)throw std::runtime_error("fullauto: insufficient piecewise-linear support");
    for(size_t i=0;i<x.size();++i)
        if(!std::isfinite(x[i])||(i&&!(x[i]>x[i-1])))
            throw std::runtime_error("fullauto: invalid piecewise-linear coordinates");
}

struct PiecewiseLinearLocation {size_t left,right;double fraction;};

inline PiecewiseLinearLocation LocatePiecewiseLinear(const std::vector<double>& x,double q) {
    if(!std::isfinite(q)||q<x.front()-1e-12||q>x.back()+1e-12)
        throw std::runtime_error("fullauto: piecewise-linear extrapolation forbidden");
    // Retain the existing endpoint-roundoff allowance, but never extrapolate
    // with negative weights. Exact knots are returned without arithmetic.
    if(q<=x.front())return {0,0,0};
    if(q>=x.back())return {x.size()-1,x.size()-1,0};
    const size_t right=std::lower_bound(x.begin(),x.end(),q)-x.begin();
    if(x[right]==q)return {right,right,0};
    const size_t left=right-1;
    const double u=(q-x[left])/(x[right]-x[left]);
    return {left,right,u};
}

inline std::vector<std::pair<size_t,double>> PiecewiseLinearWeights(const std::vector<double>& x,double q) {
    ValidatePiecewiseLinearGrid(x);
    const auto at=LocatePiecewiseLinear(x,q);
    if(at.left==at.right)return {{at.left,1.}};
    return {{at.left,1-at.fraction},{at.right,at.fraction}};
}

inline Matrix PiecewiseLinearInterpolate(const std::vector<double>& x,const Matrix& y,const std::vector<double>& query) {
    ValidatePiecewiseLinearGrid(x);
    if(y.size()!=x.size())throw std::runtime_error("fullauto: invalid piecewise-linear matrix dimensions");
    Matrix output(query.size());
    for(size_t r=0;r<query.size();++r){
        const auto at=LocatePiecewiseLinear(x,query[r]);
        if(at.left==at.right){output[r]=y[at.left];continue;}
        for(int j=0;j<16;++j)
            output[r][j]=(1-at.fraction)*y[at.left][j]+at.fraction*y[at.right][j];
    }
    return output;
}

inline Matrix POOutputInterpolate(bool go,bool coherent,const std::vector<double>& x,const Matrix& y,const std::vector<double>& query) {
    return UseLinearPOInterpolation(go,coherent)?PiecewiseLinearInterpolate(x,y,query):LocalCubicInterpolate(x,y,query);
}

inline std::vector<std::pair<size_t,double>> POOutputWeights(bool go,bool coherent,const std::vector<double>& x,double q) {
    return UseLinearPOInterpolation(go,coherent)?PiecewiseLinearWeights(x,q):LocalCubicWeights(x,q);
}

}
