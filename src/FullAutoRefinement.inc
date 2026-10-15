// Scheduling decisions use observed errors; acceptance always uses Estimate.
struct ErrorMetric {int element;bool ratio;double error;ErrorMetric(int j=0,bool r=false,double v=0):element(j),ratio(r),error(v){}};
ErrorMetric WorstMetric(const Statistics& stats,size_t row,bool ratios){
    ErrorMetric m;for(int j=0;j<16;++j){if(stats.interval[row][j]>m.error)m={j,false,stats.interval[row][j]};
        if(ratios&&stats.ratioInterval[row][j]>m.error)m={j,true,stats.ratioInterval[row][j]};}
    // An unbounded Fieller interval first needs a resolved denominator.
    if(m.ratio&&stats.interval[row][0]>=1)m={0,false,stats.interval[row][0]};return m;
}
std::vector<Samples> WeightedLevels(const AdaptiveData& data,const std::vector<size_t>& ids,const std::vector<std::pair<size_t,double>>& weights){
    std::vector<Samples> result(data.levels.size(),Samples(8,Matrix(1)));
    for(const auto& w:weights){const auto pos=std::lower_bound(ids.begin(),ids.end(),w.first);Require(pos!=ids.end()&&*pos==w.first,"missing refinement support");const size_t row=pos-ids.begin();
        for(size_t l=0;l<result.size();++l)for(int s=0;s<8;++s)for(int j=0;j<16;++j)result[l][s][0][j]+=w.second*data.levels[l][s][row][j];}return result;
}
std::vector<std::pair<size_t,double>> QueryWeights(const Engine& e,const std::vector<size_t>& knots,double q){
    auto weights=LocalCubicWeights(Angles(e.theta,knots),q);for(auto& w:weights)w.first=knots[w.first];return weights;
}
int Doubled(const Engine& e,int n){return n>e.maxCount/2?e.maxCount:2*n;}
// All support nodes of one selected level grow together. This preserves the
// observed cancellation between their seed-wise contributions. There is no
// assumption that variances of neighbouring angular rows are independent.
bool RefineError(Engine& e,AdditionalWork& work,const std::vector<int>& depths,
        const Budgets& current,Budgets& proposed,const std::vector<size_t>& evaluated,
        const AdaptiveData& data,const std::vector<std::pair<size_t,double>>& weights,
        const ErrorMetric& metric,const Mueller& reference,double query,const std::string& reason,
        std::ostringstream& log,const std::vector<double>& discrepancy={},const Budgets& observedCounts={}){
    auto values=WeightedLevels(data,evaluated,weights);
    if(!observedCounts.empty()){
        Require(observedCounts.size()==depths.size(),"audit refinement budget level mismatch");
        for(auto& level:values)for(auto& seed:level)seed[0].fill(0);
        for(const auto& w:weights){const size_t row=std::lower_bound(evaluated.begin(),evaluated.end(),w.first)-evaluated.begin();
            for(size_t l=0;l<depths.size();++l){
                const double projection=std::sqrt(double(observedCounts[l][w.first])/current[l][w.first]);
                for(int s=0;s<8;++s)for(int j=0;j<16;++j)values[l][s][0][j]+=w.second*projection*data.levels[l][s][row][j];
            }
        }
    }
    std::vector<size_t> nodes;for(const auto& w:weights)nodes.push_back(w.first);
    if(!e.local)nodes=evaluated;
    std::vector<std::array<double,8>> centered(depths.size());std::array<double,8> total{};
    const double ratio=metric.ratio?reference[metric.element]/reference[0]:0;
    for(size_t l=0;l<depths.size();++l){double mean=0;
        for(int s=0;s<8;++s){centered[l][s]=(values[l][s][0][metric.element]-ratio*values[l][s][0][0]);mean+=centered[l][s]/8;}
        for(int s=0;s<8;++s){centered[l][s]-=mean;total[s]+=centered[l][s];}}
    struct Choice {double variance=0,reduction=0,seconds=0;bool profiled=false,eligible=false;};std::vector<Choice> choices(depths.size());
    bool positive=false;
    for(size_t l=0;l<depths.size();++l){auto& choice=choices[l];double shrink=0;size_t growing=0;
        for(size_t id:nodes){const int n=current[l][id],next=Doubled(e,n);shrink=std::max(shrink,std::sqrt(double(n)/next));growing+=next>n;}
        if(!growing)continue;choice.eligible=true;
        double covariance=0;for(int s=0;s<8;++s){choice.variance+=centered[l][s]*centered[l][s]/7;covariance+=centered[l][s]*total[s]/7;}
        choice.reduction=2*(1-shrink)*covariance-(1-shrink)*(1-shrink)*choice.variance;
        if(!discrepancy.empty())choice.reduction+=(1-shrink*shrink)*discrepancy[l]*discrepancy[l];
        positive=positive||choice.reduction>0;
        auto requests=[&](bool grow){std::map<int,std::vector<double>> groups;for(size_t id:nodes)groups[grow?Doubled(e,current[l][id]):current[l][id]].push_back(e.theta[id]);
            std::vector<Engine::WorkRequest> result;for(auto& item:groups){auto& x=item.second;if(x.size()==1){x.push_back(e.theta[evaluated.front()]==x[0]?e.theta[evaluated.back()]:e.theta[evaluated.front()]);std::sort(x.begin(),x.end());}
                result.push_back(Engine::WorkRequest{item.first,0,depths[l],l&&e.depthSharing?depths[l-1]:0,x});
                if(l&&!e.depthSharing)result.push_back(Engine::WorkRequest{item.first,0,depths[l-1],0,x});}return result;};
        const auto before=work.Forecast(requests(false)),after=work.Forecast(requests(true));choice.profiled=before.profiled&&after.profiled;
        if(choice.profiled)choice.seconds=std::max(0.,after.seconds-before.seconds);
        else{double cost=0;for(const auto& w:weights){const size_t row=std::lower_bound(evaluated.begin(),evaluated.end(),w.first)-evaluated.begin();cost+=data.costs[l][row]/weights.size();}
            int delta=0;for(size_t id:nodes)delta=std::max(delta,Doubled(e,current[l][id])-current[l][id]);choice.seconds=8*cost*delta;}
    }
    size_t best=depths.size();double score=-1;
    for(size_t l=0;l<choices.size();++l){const auto& c=choices[l];if(!c.eligible)continue;
        const double gain=positive?std::max(0.,c.reduction):c.variance;
        const double efficiency=gain/std::max(1e-12,c.seconds);if(efficiency>score){score=efficiency;best=l;}}
    for(size_t l=0;l<choices.size();++l){const auto& c=choices[l];log<<reason<<','<<Number(query)<<",M"<<metric.element/4+1<<metric.element%4+1<<(metric.ratio?"/M11":"/scale(M11)")<<','<<l
        <<','<<Number(metric.error)<<','<<Number(c.variance)<<','<<Number(c.reduction)<<','<<Number(c.seconds)<<','<<c.profiled<<','<<(l==best)<<','<<weights.size()<<'\n';}
    if(best==depths.size())return false;
    bool grew=false;for(size_t id:nodes){const int next=std::max(proposed[best][id],Doubled(e,current[best][id]));grew=grew||next>proposed[best][id];proposed[best][id]=next;}return grew;
}
std::string RefinementHeader(){return "reason,theta_deg,metric,level,error95,variance,predicted_variance_reduction,additional_worker_seconds,cost_profiled,selected,support_nodes\n";}
void RefinementLog(Engine& e,int round,const std::string& label,const std::ostringstream& rows){
    const auto filename=label+"_round"+Number(round)+"_"+Hash(rows.str())+".csv";Atomic(e.root+"/"+filename,rows.str());
    e.Log(label,"{\"details\":"+Quoted(filename)+",\"scope\":\"error and additional-cost scheduling only; observed accuracy gates unchanged\"}");
}
