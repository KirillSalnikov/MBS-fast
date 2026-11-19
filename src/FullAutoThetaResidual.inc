// Included inside Engine. This estimator proposes interpolation knots only.
// The full-resolution prefix and sparse tail use disjoint Sobol indices of
// the same scramble. The angular map is fixed before reading any samples.
using ResidualStencil=std::vector<std::pair<size_t,double>>;
ResidualStencil ResidualWeights(const std::vector<double>& x,double q)const{
    const auto upper=std::lower_bound(x.begin(),x.end(),q);
    if(upper!=x.end()&&*upper==q)return {{size_t(upper-x.begin()),1.}};
    if(x.size()<8||upper==x.begin()||upper==x.end())return LocalCubicWeights(x,q);
    const size_t right=size_t(upper-x.begin());
    const size_t first=std::min(right>4?right-4:0,x.size()-8);
    ResidualStencil weights;
    double amplification=0;
    for(size_t j=first;j<first+8;++j){
        double w=1.;
        for(size_t k=first;k<first+8;++k)if(k!=j)w*=(q-x[k])/(x[j]-x[k]);
        weights.emplace_back(j,w);amplification+=std::fabs(w);
    }
    // A cluster of anchors can greatly amplify sample noise. Select this
    // fixed map by geometry alone; the subsequent independent guard remains.
    return amplification<=3.?weights:LocalCubicWeights(x,q);
}
#include "FullAutoPrefixAnglePlan.inc"
bool PairedThetaPilot(int count,Samples& data,Matrix& reference){
    if(!pairedThetaPilot||theta.size()<128||count<32||count>256||(count&(count-1))||
       phiCutoffPilot||(controls&&!sharedMeans))return false;
    const int prefix=count/4,tail=count-prefix;
    const double fraction=double(tail)/count;
    std::vector<size_t> anchors,probes,missing,measured;
    std::vector<double> anchorAngles,measuredAngles;
    for(size_t r=0;r<theta.size();++r){
        const bool anchor=r<4||r+4>=theta.size()||r%4==0;
        const size_t probeOffset=4*((r/16)%4)+1+(r/16)%3;
        const bool probe=!anchor&&r%16==probeOffset;
        if(anchor){anchors.push_back(r);anchorAngles.push_back(theta[r]);}
        if(probe)probes.push_back(r);
        if(anchor||probe){measured.push_back(r);measuredAngles.push_back(theta[r]);}
        else missing.push_back(r);
    }
    if(probes.empty()||missing.empty())return false;
    const auto head=DirectStage(prefix,angularValidation,8,target,theta);
    const Samples low=Data(head);
    std::vector<ResidualStencil> stencils;
    if(prefixAnglePlan){
        std::vector<double> scale(theta.size(),0);
        for(size_t i=0;i<theta.size();++i)for(size_t s=0;s<8;++s)
            scale[i]+=(controls?head[s].controls[i][0]:low[s][i][0])/8;
        const size_t old=measured.size();
        if(PrefixAnglePlan(low,scale,anchors,probes,stencils)){
            measured=anchors;measured.insert(measured.end(),probes.begin(),probes.end());std::sort(measured.begin(),measured.end());
            anchorAngles.clear();for(auto i:anchors)anchorAngles.push_back(theta[i]);
            missing.clear();measuredAngles.clear();
            for(size_t i=0;i<theta.size();++i){
                if(std::binary_search(measured.begin(),measured.end(),i))measuredAngles.push_back(theta[i]);
                else missing.push_back(i);
            }
        }
        Log("theta_tail_plan","{\"original_rows\":"+Number(old)+",\"planned_rows\":"+Number(measured.size())+",\"probes\":"+Number(probes.size())+"}");
    }
    // Native randomized Sobol blocks have power-of-two lengths. The tail
    // [N/4,N) is exactly two aligned blocks of lengths N/4 and N/2.
    const auto runTail=[&](const std::vector<double>& angles,std::vector<double>& raw){
        Samples combined(8,Matrix(angles.size()));raw.assign(angles.size(),0.);
        for(int offset=prefix;offset<count;offset*=2){
            const auto jobs=DirectStage(offset,angularValidation,8,target,angles,true,false,offset);
            const auto samples=Data(jobs);const double weight=double(offset)/tail;
            for(size_t s=0;s<8;++s)for(size_t r=0;r<angles.size();++r){
                if(controls)raw[r]+=weight*jobs[s].controls[r][0]/8;
                for(size_t c=0;c<16;++c)combined[s][r][c]+=weight*samples[s][r][c];
            }
        }
        return combined;
    };
    std::vector<double> measuredRaw;
    const Samples high=runTail(measuredAngles,measuredRaw);
    Matrix rawScale=Mean(low);
    if(controls)for(size_t r=0;r<theta.size();++r){
        rawScale[r][0]=0;
        for(const auto& j:head)rawScale[r][0]+=j.controls[r][0]/8;
    }
    data=low;
    Samples held(8,Matrix(probes.size()));
    if(highOrderThetaResidual&&stencils.empty()){
        stencils.reserve(theta.size());
        for(double q:theta)stencils.push_back(ResidualWeights(anchorAngles,q));
    }
    for(size_t s=0;s<8;++s){
        Matrix change;size_t i=0;
        for(size_t row:anchors){
            while(measured.at(i)<row)++i;
            Mueller difference{};
            for(size_t c=0;c<16;++c)difference[c]=fraction*(high[s][i][c]-low[s][row][c]);
            change.push_back(difference);
        }
        if(highOrderThetaResidual){
            for(size_t r=0;r<theta.size();++r){
                Mueller delta{};
                for(const auto& term:stencils[r])for(size_t c=0;c<16;++c)delta[c]+=term.second*change[term.first][c];
                for(size_t c=0;c<16;++c)data[s][r][c]+=delta[c];
            }
        }else{
            const auto predicted=LocalCubicInterpolate(anchorAngles,change,theta);
            for(size_t r=0;r<theta.size();++r)for(size_t c=0;c<16;++c)data[s][r][c]+=predicted[r][c];
        }
        size_t probe=0;
        for(size_t k=0;k<measured.size();++k){
            const size_t r=measured[k];const bool heldOut=probe<probes.size()&&r==probes[probe];
            for(size_t c=0;c<16;++c){
                const double full=low[s][r][c]+fraction*(high[s][k][c]-low[s][r][c]);
                if(heldOut)held[s][probe][c]=data[s][r][c]-full;
                data[s][r][c]=full;
            }
            if(heldOut)++probe;
        }
    }
    const auto bias=Mean(held),variance=Variance(held);double maximum=0;
    for(size_t r=0;r<probes.size();++r)for(size_t c=0;c<16;++c){
        const double scale=rawScale[probes[r]][0];
        const double value=scale>0?(std::fabs(bias[r][c])+t95*std::sqrt(variance[r][c]/8))/scale:1e100;
        maximum=std::max(maximum,value);
    }
    reference=Mean(data);
    // Use the directly evaluated physical intensity for pilot normalization.
    // A signed analytical control can produce a negative finite mean.
    if(controls)for(size_t r=0;r<theta.size();++r)reference[r][0]=rawScale[r][0];
    const bool complete=maximum>interpolation/4||std::any_of(reference.begin(),reference.end(),[](const Mueller& a){return !(a[0]>0);});
    if(complete){
        std::vector<double> angles;for(size_t r:missing)angles.push_back(theta[r]);
        std::vector<double> extraRaw;
        const auto extra=runTail(angles,extraRaw);
        for(size_t s=0;s<8;++s)for(size_t k=0;k<missing.size();++k){
            const size_t r=missing[k];
            for(size_t c=0;c<16;++c)data[s][r][c]=low[s][r][c]+fraction*(extra[s][k][c]-low[s][r][c]);
        }
        reference=Mean(data);
        if(controls)for(size_t k=0,m=0,r=0;r<theta.size();++r){
            const bool known=k<measured.size()&&r==measured[k];
            const double raw=known?measuredRaw[k]:extraRaw[m];
            reference[r][0]=rawScale[r][0]+fraction*(raw-rawScale[r][0]);
            if(known)++k;else ++m;
        }
    }
    Log("paired_theta_pilot","{\"prefix_count\":"+Number(prefix)+",\"tail_count\":"+Number(tail)
        +",\"full_grid_rows\":"+Number(theta.size())+",\"anchor_rows\":"+Number(anchors.size())
        +",\"check_rows\":"+Number(probes.size())+",\"residual_order\":"+Number(highOrderThetaResidual?7:3)+",\"check_bound_scaled_M11\":"+Number(maximum)
        +",\"check_budget\":"+Number(interpolation/4)+",\"completed_all_tail_rows\":"+(complete?"true":"false")
        +",\"production_accuracy_accepted\":false}");
    return true;
}
