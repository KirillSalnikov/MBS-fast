// A facet control changes only M11 in a native unweighted calibration job.
// Keep bounded raw validation samples from the last controlled grid. A fresh
// process or different physics automatically uses the normal calibration.
bool rawPhiPilots=false,rawPhiFallback=false;
size_t rawPhiBytes=0;
std::string rawPhiDomain;
std::map<std::pair<int,int>,Samples> rawPhiSamples;
std::string RawPhiDomain(const std::vector<double>& x,int depth)const{
    if(!rawPhiPilots)return "";
    std::ostringstream recipe;
    recipe<<signature<<' '<<ExecutionHash()<<' '<<depth<<' '<<mirror<<' '<<kernel.block
          <<' '<<kernel.warp<<' '<<kernel.pipeline<<' '<<dispatchPipeline<<'\n';
    for(double angle:x)recipe<<Number(angle)<<'\n';
    return Hash(recipe.str());
}
void RememberRawPhi(int n,int rule,const std::vector<Job>& valid){
    const auto key=std::make_pair(n,rule);
    if(rawPhiSamples.count(key))return;
    const size_t bytes=valid.size()*valid.front().matrix.size()*sizeof(Mueller);
    const size_t limit=24u*1024*1024;
    if(bytes>limit-rawPhiBytes)return;
    Samples raw;
    for(size_t s=0;s<valid.size();++s){
        const auto& job=valid[s];
        Require(job.seed==validation[s]&&job.hasControls&&job.matrix.size()==job.controls.size(),"invalid raw phi source");
        raw.push_back(job.matrix);
        for(size_t r=0;r<job.controls.size();++r)raw.back()[r][0]=job.controls[r][0];
    }
    rawPhiSamples.emplace(key,std::move(raw));rawPhiBytes+=bytes;
}
bool RestoreRawPhi(const std::string& domain,int n,int rule,int depth,const std::vector<double>& x,std::vector<Job>& valid){
    if(!rawPhiPilots||controls||!rawPhiFallback||domain!=rawPhiDomain||phiCutoffPilot)return false;
    const auto it=rawPhiSamples.find({n,rule});
    if(it==rawPhiSamples.end())return false;
    const auto& raw=it->second;
    const size_t timed=PhiTimingSeeds();
    auto directlyTimed=DirectStage(n,validation,timed,depth,x,false);
    double difference=0.;
    for(size_t s=0;s<timed;++s){
        const auto& ref=directlyTimed[s].matrix;
        for(size_t r=0;r<x.size();++r){
            const auto& a=ref[r];const auto& b=raw[s][r];
            if(!(a[0]>0&&b[0]>0))return false;
            for(size_t j=0;j<16;++j){
                const double absolute=std::fabs(a[j]-b[j])/a[0];
                const double normalized=std::fabs(a[j]/a[0]-b[j]/b[0]);
                if(!std::isfinite(absolute)||!std::isfinite(normalized))return false;
                difference=std::max(difference,std::max(absolute,normalized));
            }
        }
    }
    const bool passed=difference<=eps/1024;
    Log("raw_phi_fallback_check","{\"phi\":"+Number(rule)+",\"depth\":"+Number(depth)+",\"count\":"+Number(n)
        +",\"theta_rows\":"+Number(x.size())+",\"timed_all16_and_ratio_difference\":"+Number(difference)
        +",\"budget\":"+Number(eps/1024)+",\"reused_scrambles\":"+Number(passed?8-timed:0)+",\"directly_timed_scrambles\":"+Number(timed)+"}");
    if(!passed)return false;
    valid=std::move(directlyTimed);
    for(size_t s=timed;s<raw.size();++s){
        Job view;SetSampling(view,n,validation[s],depth,x,false);
        view.matrix=raw[s];view.angles=x;valid.push_back(std::move(view));
    }
    return true;
}
