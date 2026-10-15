// Scheduling forecast only. Scale the eight centered *joint* level samples,
// preserving all observed cross-level and numerator/denominator covariances.
// Future accuracy is always measured again from native matrices.
namespace {
struct AllocationForecast {
    std::vector<Samples> centered;std::vector<int> counts;Matrix reference;bool ratios;
    AllocationForecast(const std::vector<Samples>& data,const std::vector<int>& ns,const Matrix& mean,bool useRatios)
        :centered(data),counts(ns),reference(mean),ratios(useRatios){
        Require(data.size()==ns.size()&&!data.empty()&&mean.size()==1&&mean[0][0]>0,"invalid covariance allocation model");
        for(size_t l=0;l<data.size();++l){const auto m=Mean(data[l]);Require(m.size()==1&&ns[l]>0,"allocation expects one angular row");
            for(auto& sample:centered[l])for(int j=0;j<16;++j)sample[0][j]-=m[0][j];}
        // An unbounded Fieller interval cannot rank small improvements.
        // First resolve the common denominator using the raw matrix metrics;
        // the next native decision still tests every requested ratio.
        if(ratios&&Predict(ns).interval[0][0]>=1)ratios=false;
    }
    Statistics Predict(const std::vector<int>& next)const{
        Require(next.size()==counts.size(),"allocation forecast dimensions differ");Samples samples(8,reference);
        for(size_t l=0;l<counts.size();++l){Require(next[l]>=counts[l],"allocation cannot shrink a live prefix");const double scale=std::sqrt(double(counts[l])/next[l]);
            for(int s=0;s<8;++s)for(int j=0;j<16;++j)samples[s][0][j]+=scale*centered[l][s][0][j];}
        return Estimate(samples);
    }
    double Loss(const std::vector<int>& next,double eps)const{
        const auto stats=Predict(next);double loss=0;
        auto add=[&](double v){const double r=std::min(1e90,v/eps);loss+=std::max(0.,r*r-1);};
        for(double v:stats.interval[0])add(v);if(ratios)for(int j=1;j<16;++j)add(stats.ratioInterval[0][j]);return loss;
    }
};
}
std::vector<int> AllocateCovariant(const std::vector<Samples>& levels,const std::vector<double>& costs,
        const std::vector<int>& counts,const Matrix& reference,double eps,bool ratios,int maximum,double decisionTolerance){
    Require(levels.size()==costs.size()&&levels.size()<=8&&eps>0,"invalid covariance allocation dimensions");
    AllocationForecast model(levels,counts,reference,ratios);auto result=counts;std::vector<int> cap;
    for(size_t l=0;l<counts.size();++l){Require(costs[l]>0&&std::isfinite(costs[l])&&counts[l]<=maximum,"invalid covariance allocation cost/limit");cap.push_back(counts[l]>maximum/4?maximum:counts[l]*4);}
    // At most two dyadic increments per level in a production round. Joint
    // candidates retain cancellations that a one-coordinate greedy step loses.
    for(size_t step=0;step<2*counts.size();++step){const double before=model.Loss(result,eps);if(before==0)break;
        auto best=result;double score=0;
        for(unsigned mask=1;mask<(1u<<counts.size());++mask){auto next=result;double seconds=0;
            for(size_t l=0;l<counts.size();++l)if(mask&(1u<<l)){next[l]=std::min(cap[l],result[l]>maximum/2?maximum:2*result[l]);seconds+=costs[l]*(next[l]-result[l]);}
            if(seconds<=0)continue;const double gain=before-model.Loss(next,eps),efficiency=gain/seconds;
            if(efficiency>score){score=efficiency;best=next;}}
        if(best==result){ // Preserve progress when the small pilot cannot rank it.
            if(result!=counts)break;
            for(size_t l=0;l<counts.size();++l)best[l]=std::min(cap[l],result[l]>maximum/2?maximum:2*result[l]);
            if(best==result)break;
        }
        result=best;
    }
    // Keep the committed independent-level policy as a cost ceiling. A new
    // covariance model is used only when its *next* step is cheaper and its
    // measured joint-noise reduction per unit work is at least as good.
    auto legacy=Allocate(levels,costs,counts,reference,eps);
    if(ratios){auto transformed=levels;for(auto& level:transformed)for(auto& sample:level)for(int j=1;j<16;++j)sample[0][j]-=reference[0][j]/reference[0][0]*sample[0][0];
        const auto normalized=Allocate(transformed,costs,counts,reference,eps);for(size_t l=0;l<legacy.size();++l)legacy[l]=std::max(legacy[l],normalized[l]);}
    for(size_t l=0;l<legacy.size();++l)legacy[l]=std::max(counts[l],std::min(cap[l],legacy[l]));
    if(!(decisionTolerance>0))decisionTolerance=eps;
    legacy=PrioritizeAllocation(levels,costs,counts,legacy,decisionTolerance,ratios,false);
    result=PrioritizeAllocation(levels,costs,counts,result,decisionTolerance,ratios,true);
    double oldWork=0,newWork=0;for(size_t l=0;l<counts.size();++l){oldWork+=costs[l]*(legacy[l]-counts[l]);newWork+=costs[l]*(result[l]-counts[l]);}
    const double before=model.Loss(counts,eps),oldGain=before-model.Loss(legacy,eps),newGain=before-model.Loss(result,eps);
    const bool accepted=newGain>0&&(oldWork==0||newWork<=oldWork*(1+1e-12))
        &&newGain/std::max(1e-30,newWork)>=oldGain/std::max(1e-30,oldWork);
    return accepted?result:legacy;
}
