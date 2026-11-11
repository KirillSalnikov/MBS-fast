// Shared automatic physical-depth selection. Sampling/geometry and exact
// worker commands retain the same cache identity as production.
namespace {
void SaveSelectedDepth(Engine& e,int selected){
    const std::string payload=e.signature+" "+Number(selected)+"\n";
    Atomic(e.root+"/automatic_depth.tsv",Hash(payload)+"\n"+payload);
    e.target=selected;e.SetMirror(e.mirror);
    e.Log("physical_depth_selected","{\"target_depth\":"+Number(selected)+",\"start\":"+Number(e.depthStart)+",\"maximum\":"+Number(e.depthMaximum)+",\"budget\":"+Number(e.depthBudget)+"}");
}
struct DepthResidual {double upper=0,lower=0,mean=0,interval=0;};
DepthResidual DepthDifference(const Samples& first,const Samples& second,bool go,size_t bins,bool ratios=true,const Matrix* requestedScale=nullptr,std::vector<DepthResidual>* perRow=nullptr){
    const auto a=Mean(first),b=Mean(second);const auto& scale=requestedScale?*requestedScale:b;
    const auto delta=Subtract(second,first);const auto mean=Mean(delta),variance=Variance(delta);
    DepthResidual total;if(perRow)perRow->clear();
    for(size_t row=0;row<scale.size();++row){
        DepthResidual result;
        auto accumulate=[&](double bias,double ci){
            if(!std::isfinite(bias)||!std::isfinite(ci)){result.upper=1e100;return;}
            result.upper=std::max(result.upper,bias+ci);result.lower=std::max(result.lower,std::max(0.,bias-ci));
            result.mean=std::max(result.mean,bias);result.interval=std::max(result.interval,ci);};
        auto finish=[&](){total.upper=std::max(total.upper,result.upper);total.lower=std::max(total.lower,result.lower);
            total.mean=std::max(total.mean,result.mean);total.interval=std::max(total.interval,result.interval);if(perRow)perRow->push_back(result);};
        if(!std::isfinite(scale[row][0])){result.upper=1e100;finish();continue;}
        if(!(scale[row][0]>0)){
            if(!go||row<bins||mean[row][0]!=0||variance[row][0]!=0)result.upper=1e100;
            finish();continue;
        }
        if(ratios){
            double va=0,vb=0;for(int seed=0;seed<8;++seed){va+=std::pow(first[seed][row][0]-a[row][0],2)/7;vb+=std::pow(second[seed][row][0]-b[row][0],2)/7;}
            if(t95*std::sqrt(va/8)>=a[row][0]||t95*std::sqrt(vb/8)>=b[row][0])result.upper=1e100;
        }
        for(int j=0;j<16;++j){
            accumulate(std::fabs(mean[row][j])/scale[row][0],t95*std::sqrt(variance[row][j]/8)/scale[row][0]);
            if(!ratios||j==0)continue;
            if(!(a[row][0]>0&&b[row][0]>0)){result.upper=1e100;continue;}
            const double ra=a[row][j]/a[row][0],rb=b[row][j]/b[row][0];double v=0;
            // Paired influence function for a difference of ratios of means.
            // Retains covariance both between depths and with M11.
            for(int seed=0;seed<8;++seed){
                const double influence=(second[seed][row][j]-rb*second[seed][row][0])/b[row][0]
                    -(first[seed][row][j]-ra*first[seed][row][0])/a[row][0];v+=influence*influence/7;
            }
            accumulate(std::fabs(rb-ra),t95*std::sqrt(v/8));
        }
        finish();
    }return total;
}
double PhysicalDepthBound(const Samples& first,const Samples& second,const Matrix& reference,bool go,size_t bins,bool ratios){
    return DepthDifference(first,second,go,bins,ratios,&reference).upper;
}
void PrepareDepthPhi(Engine& e,const std::vector<double>& x){
    if(e.fixedPhi)return;
    const unsigned seeds[8]={1800001001u,1800001013u,1800001031u,1800001043u,1800001057u,1800001069u,1800001081u,1800001099u};
    for(unsigned seed:seeds)Require(std::find(e.seeds.begin(),e.seeds.end(),seed)==e.seeds.end(),"startup phi/production seed overlap");
    std::set<int> candidates;for(int candidate:{4,32})candidates.insert(std::max(e.minPhi,std::min(e.maxPhi,candidate)));
    const int n=std::min(e.pilot,64);int winner=e.phi;double best=std::numeric_limits<double>::max();bool resolvedWinner=false;
    for(int candidate:candidates){
        e.phi=candidate;const auto jobs=e.Stage(n,seeds,8,e.target,x,false);const auto stats=Estimate(e.Data(jobs,false));
        const double uncertainty=e.ratios?std::max(stats.maxInterval,stats.maxRatioInterval):stats.maxInterval;
        const bool resolved=FullAutoScheduling::Resolved(uncertainty);
        const double error=resolved?uncertainty:(FullAutoScheduling::Resolved(stats.maxInterval)?stats.maxInterval:1e25);
        double seconds=0;for(const auto& job:jobs)seconds+=std::max(1e-6,job.seconds-job.controlSetup);
        const double score=seconds*std::max(double(e.depthPilot)/n,std::pow(error/e.depthBudget,2));
        if((resolved&&!resolvedWinner)||(resolved==resolvedWinner&&score<best)){winner=candidate;best=score;resolvedWinner=resolved;}
        e.Log("startup_phi_candidate","{\"phi\":"+Number(candidate)+",\"count_per_seed\":"+Number(n)
            +",\"theta_rows\":"+Number(x.size())+",\"seconds\":"+Number(seconds)+",\"forecast_score\":"+Number(score)
            +",\"uncertainty95\":"+Number(uncertainty)+",\"resolved\":"+(resolved?"true":"false")
            +",\"accuracy_accepted\":false,\"adaptive_pilot_growth\":false,\"full_production_selection_follows\":true}");
    }
    e.phi=winner;
}
void PrepareDepthStartup(Engine& e){
    if(!e.autoDepth||e.go||e.depthOrder=="depth-first"||(e.depthOrder=="auto"&&!e.config.useGpu))return;
    const std::string checkpoint=e.root+"/depth_startup.tsv";
    if(Exists(checkpoint)){
        const auto text=Read(checkpoint);const auto split=text.find('\n');const auto payload=text.substr(split+1);Require(Hash(payload)==text.substr(0,split),"depth startup changed");
        std::istringstream input(payload);std::string signature;Require(bool(input>>signature>>e.phi>>e.kernel.block>>e.kernel.pipeline>>e.kernel.warp>>e.controls)&&signature==e.signature,"depth startup inputs changed");
        e.startupKernel=true;e.SetMirror(false);e.WarmCache();return;
    }
    std::vector<double> x;const size_t size=std::min(size_t(9),std::min(e.theta.size(),size_t(e.maxTheta)));
    for(size_t row=0;row<size;++row)x.push_back(e.theta[row*(e.theta.size()-1)/(size-1)]);
    const int oldPilot=e.pilot,oldKernel=e.kernelCount;const bool oldControls=e.controls;
    e.pilot=std::min(e.pilot,64);e.kernelCount=std::min(e.kernelCount,64);e.controls=false;
    const double begin=Seconds();e.State("startup_sampling_tuning",0,{e.pilot},x.size(),0,0,0);
    // The depth pilot uses paired, unbiased Haar estimates. It does not need
    // a full search over azimuth rules or control coefficients at a depth
    // which has not yet been selected. Calibrate the execution kernel only;
    // production still performs the complete phi/control validation below.
    try{e.Calibrate(x,false,true);PrepareDepthPhi(e,x);}catch(...){e.pilot=oldPilot;e.kernelCount=oldKernel;e.controls=oldControls;throw;}
    e.pilot=oldPilot;e.kernelCount=oldKernel;e.controls=oldControls;e.startupKernel=true;
    std::ostringstream payload;payload<<e.signature<<' '<<e.phi<<' '<<e.kernel.block<<' '<<e.kernel.pipeline<<' '<<e.kernel.warp<<' '<<e.controls<<'\n';
    Atomic(checkpoint,Hash(payload.str())+"\n"+payload.str());
    e.Record(e.root+"/depth_startup.json","{\"order\":\"kernel-first\",\"phi\":"+Number(e.phi)+",\"kernel_pipeline\":"+(e.kernel.pipeline?std::string("true"):"false")+",\"pilot_count\":"+Number(std::min(oldPilot,64))+",\"theta_rows\":"+Number(x.size())+",\"seconds\":"+Number(Seconds()-begin)+",\"full_phi_search_deferred\":true,\"scope\":\"bounded kernel pilot before paired depth; complete phi/control selection and resource tuning follow at the selected physical depth\"}");
}
bool DepthPilotFailureResolved(double lower,double interval,double budget,bool terminal){
    return lower>budget&&(!terminal||interval<=budget/4);
}
#include "FullAutoDepthPromotion.inc"
#include "FullAutoDepthSearch.inc"
bool PromoteAutomaticDepth(Engine& e,int code){
    if(code!=4||!e.autoDepth||e.target>=e.depthMaximum)return false;
    const auto validation=CompactJson(Read(e.root+"/physics_validation.json"));
    if(validation.find("\"depth_pass\":false")==std::string::npos||validation.find("\"cutoff_pass\":false")!=std::string::npos)return false;
    const int old=e.target,next=std::min(e.depthMaximum,old+2);
    BeginDepthPromotion(e,old,next);return ResumeDepthPromotion(e);
}
}
