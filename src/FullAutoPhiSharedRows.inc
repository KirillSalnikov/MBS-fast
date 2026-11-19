// Retain only observed per-row uncertainties and timing measurements. A record
// may guide a working-grid fit only for matching native settings and angles.
struct StartupPhiRows {
    std::vector<double> angles, error;
    std::vector<std::pair<double,double>> times;
    int count=0;
    double setup=0;
    std::string recipe;
};
std::map<int,StartupPhiRows> startupPhiRows;
std::string PhiRowRecipe(int depth) const {
    std::ostringstream key;
    key<<ExecutionHash()<<' '<<depth<<' '<<controls<<' '<<sharedMeans<<' '<<polarControls
       <<' '<<ratios<<' '<<mirror<<' '<<kernel.block<<' '<<kernel.warp<<' '<<kernel.pipeline;
    return key.str();
}
void RememberStartupPhiRows(int rule,int depth,int n,const std::vector<double>& x,
        const Statistics& stats,const std::vector<Job>& valid,size_t timed,double setup) {
    if(depth!=target||valid.size()!=8||!timed||timed>valid.size())return;
    StartupPhiRows rows;
    rows.angles=x;rows.count=n;rows.setup=setup;rows.recipe=PhiRowRecipe(depth);
    for(size_t i=0;i<x.size();++i){
        double error=*std::max_element(stats.interval[i].begin(),stats.interval[i].end());
        if(ratios)error=std::max(error,*std::max_element(stats.ratioInterval[i].begin(),stats.ratioInterval[i].end()));
        rows.error.push_back(error);
    }
    for(size_t i=0;i<timed;++i){
        const double wall=valid[i].seconds-valid[i].controlSetup,dif=valid[i].diffractionSeconds;
        if(!(wall>0&&dif>0&&std::isfinite(wall)&&std::isfinite(dif)))return;
        rows.times.emplace_back(wall,std::min(wall,dif));
    }
    startupPhiRows[rule]=std::move(rows);
}
struct PhiSharedProposal {int winner=0;double runner=0;size_t overlap=0;};
PhiSharedProposal SharedPhiProposal(const std::vector<double>& work) const {
    PhiSharedProposal proposal;
    if(fixedPhi||refitPhiWinner||phiCutoffPilot||work.empty())return proposal;
    struct Forecast {int rule;double lower,upper;};
    std::vector<Forecast> observed;
    for(int rule:{1,2,4,8,16,32,64}){
        if(rule<minPhi||rule>maxPhi||(rule!=phi&&rule!=phi/2&&rule!=2*phi&&rule!=minPhi))continue;
        const auto it=startupPhiRows.find(rule);
        if(it==startupPhiRows.end()||it->second.recipe!=PhiRowRecipe(target))return {};
        const auto& r=it->second;
        if(r.angles.empty()||r.angles.front()!=work.front()||r.angles.back()!=work.back())return {};
        double peak=0;size_t overlap=0;
        for(size_t i=0;i<r.angles.size();++i)if(std::binary_search(work.begin(),work.end(),r.angles[i])){
            peak=std::max(peak,r.error[i]);++overlap;
        }
        if(overlap<8||!FullAutoScheduling::Resolved(peak))return {};
        proposal.overlap=overlap;
        double lo=std::numeric_limits<double>::max(),hi=0;
        const double extra=double(work.size())/r.angles.size()-1.;
        if(extra<0)return {};
        // A shared angle supplies a lower estimate of the full-grid worst
        // uncertainty. For a possible winner include every measured startup
        // angle, even when the new grid has no row at that exact coordinate.
        const double allPeak=*std::max_element(r.error.begin(),r.error.end());
        if(!FullAutoScheduling::Resolved(allPeak))return {};
        const auto orientations=[&](double uncertainty){
            return 8*std::max(double(initial)/r.count,std::pow(uncertainty/statistical,2));
        };
        const double lowScale=orientations(peak),highScale=orientations(allPeak);
        for(const auto& sample:r.times){
            const double time=sample.first+extra*sample.second;
            const double lower=time*lowScale+r.setup,upper=time*highScale+r.setup;
            if(!(lower>0&&std::isfinite(lower)&&std::isfinite(upper)))return {};
            lo=std::min(lo,lower);hi=std::max(hi,upper);
        }
        observed.push_back({rule,lo,hi});
    }
    if(observed.size()<2)return {};
    for(const auto& v:observed){
        double runner=std::numeric_limits<double>::max();
        for(const auto& c:observed)if(c.rule!=v.rule)runner=std::min(runner,c.lower);
        if(v.upper<runner){proposal.winner=v.rule;proposal.runner=runner;return proposal;}
    }
    return {};
}
