#pragma once
#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <vector>

// Shared component budgets and scheduling-only noise forecasts. These helpers
// neither modify a Mueller matrix nor accept accuracy. Production acceptance
// remains in CheckAngularSampling.
namespace FullAutoAngularPilotCost {
inline double AngularFraction(double epsilon,double interpolation) {
    return std::min(interpolation,epsilon/4);
}
inline double AbsoluteCap(double epsilon,double interpolation) {
    return std::min(interpolation,epsilon/32);
}
template<class Matrix>
std::array<double,16> ComponentEnvelope(const Matrix& mean,double epsilon) {
    std::array<double,16> envelope;envelope.fill(epsilon);
    for(const auto& row:mean)if(std::isfinite(row[0])&&row[0]>0)
        for(size_t j=0;j<16;++j)if(std::isfinite(row[j]))
            envelope[j]=std::max(envelope[j],std::fabs(row[j]/row[0]));
    return envelope;
}
inline std::array<double,16> ComponentAllowances(const std::array<double,16>& envelope,
                                               double epsilon,double interpolation) {
    std::array<double,16> allowance;
    for(size_t j=0;j<16;++j)
        allowance[j]=std::min(AbsoluteCap(epsilon,interpolation),
                              AngularFraction(epsilon,interpolation)*envelope[j]);
    return allowance;
}
struct Result {
    bool resolved=true,finite=true;
    double maximumRatio=0,maximumNoise95=0;
    size_t pairs=0,shapes=0;
    int worstElement=-1,worstSpan=0;
    bool worstNormalized=false;
    size_t worstLeft=0,worstCenter=0,worstRight=0;
    const char* reason="resolved";
};
inline Result Invalid(const char* reason,bool finite) {
    Result result;result.resolved=false;result.finite=finite;
    result.maximumRatio=1e50;result.maximumNoise95=1e50;result.reason=reason;return result;
}
// x is the already calculated pilot grid. All of its rows are supports for
// this cost forecast, so no interpolation or extra optical jobs are required.
// Samples/Statistics deliberately use structural templates to avoid coupling
// the production controller's anonymous-namespace storage to this pure helper.
template<class Samples,class Statistics>
Result Evaluate(const Samples& samples,const Statistics& stats,const std::vector<double>& x,
                double epsilon,double interpolation,bool ratios,double t95) {
    if(!std::isfinite(epsilon)||!std::isfinite(interpolation)||!std::isfinite(t95))
        return Invalid("invalid angular pilot tolerance",false);
    if(!(epsilon>0&&epsilon<1&&interpolation>0&&t95>0))
        return Invalid("invalid angular pilot tolerance",true);
    if(samples.size()!=8||x.size()<2||stats.mean.size()!=x.size()
            ||(ratios&&stats.ratioInterval.size()!=x.size()))
        return Invalid("invalid angular pilot dimensions",true);
    for(size_t r=0;r<x.size();++r){
        if(!std::isfinite(x[r]))return Invalid("invalid angular pilot grid",false);
        if(r&&!(x[r]>x[r-1]))return Invalid("invalid angular pilot grid",true);
        if(!std::isfinite(stats.mean[r][0]))return Invalid("nonfinite angular pilot intensity",false);
        if(!(stats.mean[r][0]>0))return Invalid("nonpositive angular pilot intensity",true);
        for(size_t j=0;j<16;++j)if(!std::isfinite(stats.mean[r][j]))return Invalid("nonfinite angular pilot mean",false);
    }
    for(const auto& sample:samples){
        if(sample.size()!=x.size())return Invalid("invalid angular pilot sample rows",true);
        for(const auto& row:sample)for(size_t j=0;j<16;++j)
            if(!std::isfinite(row[j]))return Invalid("nonfinite angular pilot sample",false);
    }
    const auto allowance=ComponentAllowances(ComponentEnvelope(stats.mean,epsilon),epsilon,interpolation);
    for(double a:allowance)if(!(a>0)||!std::isfinite(a))return Invalid("invalid angular pilot allowance",false);
    Result result;
    auto check=[&](const std::array<size_t,3>& nodes,const std::array<double,3>& weights,size_t size,int span){
        if(span)++result.shapes;else ++result.pairs;
        double scale=stats.mean[nodes[0]][0];
        for(size_t k=1;k<size;++k)scale=std::min(scale,stats.mean[nodes[k]][0]);
        for(size_t j=0;j<16;++j)for(int normalized=0;normalized<=(ratios&&j?1:0);++normalized){
            std::array<double,8> values{};
            for(size_t seed=0;seed<8;++seed)for(size_t k=0;k<size;++k){
                const size_t row=nodes[k];const auto& v=samples[seed][row];const auto& m=stats.mean[row];
                values[seed]+=weights[k]*(normalized?(v[j]-(m[j]/m[0])*v[0])/m[0]:v[j]/scale);
            }
            bool finite=true,unresolved=false;double mean=0,variance=0;
            for(double v:values){finite=finite&&std::isfinite(v);mean+=v/8;}
            if(finite)for(double v:values)variance+=(v-mean)*(v-mean)/7;
            double noise=finite?t95*std::sqrt(variance/8):1e50;
            finite=finite&&std::isfinite(noise);
            if(normalized)for(size_t k=0;k<size;++k){
                const double interval=stats.ratioInterval[nodes[k]][j];
                finite=finite&&std::isfinite(interval);
                unresolved=unresolved||!std::isfinite(interval)||interval>=1e50||interval<0;
            }
            double score=noise/allowance[j];
            if(!finite||!std::isfinite(score)||unresolved||score>=1e50){
                result.resolved=false;result.finite=result.finite&&finite&&std::isfinite(score);
                result.reason=unresolved?"unresolved angular pilot ratio":"unresolved angular pilot contrast";
                noise=std::min(1e50,std::isfinite(noise)?noise:1e50);score=1e50;
            }
            result.maximumNoise95=std::max(result.maximumNoise95,noise);
            if(result.worstElement<0||score>result.maximumRatio){
                result.maximumRatio=score;result.worstElement=int(j);result.worstNormalized=normalized;
                result.worstSpan=span;result.worstLeft=nodes[0];result.worstCenter=size==3?nodes[1]:nodes[0];
                result.worstRight=nodes[size-1];
            }
        }
    };
    for(size_t k=1;k<x.size();++k)check({{k-1,k,k}},{{-1.,1.,0.}},2,0);
    for(size_t span:{size_t(1),size_t(2)})for(size_t k=span;k+span<x.size();++k){
        const size_t a=k-span,b=k,c=k+span;const double w=(x[b]-x[a])/(x[c]-x[a]);
        check({{a,b,c}},{{-(1-w),1.,-w}},3,int(span));
    }
    return result;
}
}
