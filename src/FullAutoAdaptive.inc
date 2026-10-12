// Included only by FullAuto.cpp: native adaptive scheduling and durable state.
namespace {
using Budgets=std::vector<std::vector<int>>;
struct AdaptiveData {std::vector<Samples> levels;std::vector<std::vector<double>> costs;std::vector<double> analytic;std::string sources;};
std::vector<int> PeakBudgets(const Budgets& n){std::vector<int> result;for(const auto& level:n)result.push_back(*std::max_element(level.begin(),level.end()));return result;}
double StatisticalError(const Engine& e,const Statistics& stats){return e.ratios?std::max(stats.maxInterval,stats.maxRatioInterval):stats.maxInterval;}
bool VaryingBudgets(const Budgets& ns,const std::vector<size_t>& ids){for(const auto& level:ns)for(size_t id:ids)if(level[id]!=level[ids.front()])return true;return false;}
void GateLocalCost(Engine& e,Budgets& proposed,const std::vector<int>& depths,const std::vector<size_t>& ids){
    if(e.localPolicy!="auto")return;
    std::ostringstream report;report<<"{\"policy\":\"auto\",\"levels\":[";
    for(size_t l=0;l<proposed.size();++l){double trace=0,diffraction=0;size_t observations=0;
        for(const auto& entry:e.cache){const auto& job=entry.second;if(job.assembled||job.offset||job.depth!=depths[l]||job.phi!=e.phi||job.traceSeconds<=0||job.diffractionSeconds<=0)continue;
            trace+=job.traceSeconds/job.count;diffraction+=job.diffractionSeconds/job.count/job.inputAngles.size();++observations;}
        bool unified=false;double grouped=0,global=0;
        if(observations){trace/=observations;diffraction/=observations;std::map<int,size_t> groups;int maximum=0;for(size_t id:ids){++groups[proposed[l][id]];maximum=std::max(maximum,proposed[l][id]);}
            for(const auto& group:groups)grouped+=group.first*(trace+diffraction*group.second);
            global=maximum*(trace+diffraction*ids.size());
            if(grouped>global*1.10){for(size_t id:ids)proposed[l][id]=maximum;unified=true;}}
        report<<(l?",":"")<<"{\"level\":"<<l<<",\"observations\":"<<observations<<",\"projected_local_seconds_per_seed\":"<<Number(grouped)<<",\"projected_global_seconds_per_seed\":"<<Number(global)<<",\"unified\":"<<(unified?"true":"false")<<"}";
    }report<<"]}";Atomic(e.root+"/local_cost_gate.json",report.str());
}
std::vector<size_t> AdaptiveGuards(const std::vector<size_t>& knots,size_t size,bool final=false){
    if(final)return Guards(knots,size);
    std::set<size_t> all(knots.begin(),knots.end());
    for(size_t i=1;i<knots.size();++i){size_t q=(knots[i-1]+knots[i])/2;if(q>knots[i-1]&&q<knots[i])all.insert(q);}return {all.begin(),all.end()};
}
AdaptiveData EvaluateBudget(Engine& e,const std::vector<int>& depths,const Budgets& counts,const std::vector<size_t>& ids,const unsigned* seeds){
    const auto x=Angles(e.theta,ids);AdaptiveData result;result.analytic.assign(ids.size(),0);std::ostringstream sources;
    for(size_t level=0;level<depths.size();++level){Samples data(8,Matrix(ids.size()));std::vector<double> cost(ids.size());std::map<int,std::vector<size_t>> groups;
        for(size_t r=0;r<ids.size();++r)groups[counts[level][ids[r]]].push_back(r);
        for(const auto& group:groups){std::vector<double> subset;for(size_t r:group.second)subset.push_back(x[r]);
            if(subset.size()==1){subset.push_back(x[0]==subset[0]?x.back():x[0]);std::sort(subset.begin(),subset.end());}
            std::vector<Job> upper,lower;if(level&&e.depthSharing){auto pair=e.SharedStage(group.first,seeds,8,depths[level],depths[level-1],subset);upper=pair.first;lower=pair.second;}else upper=e.Stage(group.first,seeds,8,depths[level],subset);
            Samples value=e.Data(upper);double seconds=e.Cost(upper,group.first);
            for(size_t s=0;s<8;++s)sources<<level<<" 1 "<<group.first<<' '<<seeds[s]<<' '<<upper[s].folder<<'\n';
            if(level){if(lower.empty())lower=e.Stage(group.first,seeds,8,depths[level-1],subset);value=Subtract(value,e.Data(lower));seconds+=e.Cost(lower,group.first);
                for(size_t s=0;s<8;++s)sources<<level<<" -1 "<<group.first<<' '<<seeds[s]<<' '<<lower[s].folder<<'\n';}
            for(size_t r:group.second){const size_t row=std::lower_bound(subset.begin(),subset.end(),x[r])-subset.begin();cost[r]=seconds;for(int s=0;s<8;++s)data[s][r]=value[s][row];
                if(!level&&e.controls){const auto polar=e.Polar(x[r]);const auto scalar=e.ScalarWeights(x[r]);for(const auto& job:upper)for(int j=0;j<16;++j){const double wr=j?polar[0][j]:scalar[0],ws=j?polar[1][j]:scalar[1];const auto& c=job.controls[row];result.analytic[r]=std::max(result.analytic[r],job.controlRefinement*(std::fabs(wr*c[1])+std::fabs(ws*c[2])));}}
            }
        }
        result.levels.push_back(data);result.costs.push_back(cost);
    }
    result.sources=sources.str();return result;
}
std::vector<int> AllocateRow(const Engine& e,const AdaptiveData& data,const Budgets& counts,const std::vector<size_t>& ids,size_t row,const Matrix& mean){
    std::vector<Samples> values;std::vector<double> costs;std::vector<int> ns;
    for(size_t l=0;l<data.levels.size();++l){Samples v(8,Matrix(1));for(int s=0;s<8;++s)v[s][0]=data.levels[l][s][row];values.push_back(v);costs.push_back(data.costs[l][row]);ns.push_back(counts[l][ids[row]]);}
    const Matrix reference={mean[row]};auto proposed=Allocate(values,costs,ns,reference,e.statistical*.85);
    if(e.ratios){for(auto& level:values)for(auto& seed:level)for(int j=1;j<16;++j)seed[0][j]-=mean[row][j]/mean[row][0]*seed[0][0];
        auto ratio=Allocate(values,costs,ns,reference,e.statistical*.85);for(size_t l=0;l<proposed.size();++l)proposed[l]=std::max(proposed[l],ratio[l]);}
    return proposed;
}
struct AdaptiveState {int round=0,streak=0;std::vector<size_t> knots;Budgets counts;Matrix previous;std::vector<bool> valid,coarse;};
void SaveAdaptive(Engine& e,const AdaptiveState& state,const std::string& phase="active"){
    std::ostringstream payload;payload<<std::setprecision(17)<<e.signature<<' '<<phase<<' '<<state.round<<' '<<state.streak<<' '<<e.phi<<' '<<e.controls<<' '<<e.sharedMeans<<' '<<e.mirror<<' '<<e.mirrorError<<' '<<e.kernel.block<<' '<<e.kernel.pipeline<<' '<<e.kernel.warp<<'\n';
    payload<<e.depthSharing<<' '<<e.depthVerified<<' '<<e.local<<'\n';payload<<e.weightTheta.size()<<'\n';for(size_t r=0;r<e.weightTheta.size();++r){payload<<e.weightTheta[r];
        const auto b=e.Polar(e.weightTheta[r]);for(int j=0;j<16;++j)for(int c=0;c<2;++c)payload<<' '<<(j?b[c][j]:(e.weights.empty()?1:e.weights[r][c]));payload<<'\n';}
    payload<<state.knots.size();for(size_t k:state.knots)payload<<' '<<k;payload<<'\n';
    const auto evaluated=AdaptiveGuards(state.knots,e.theta.size());std::set<size_t> saved(evaluated.begin(),evaluated.end());
    for(size_t i=0;i<e.theta.size();++i)if(state.valid[i])saved.insert(i);
    payload<<state.counts.size()<<' '<<saved.size()<<'\n';for(size_t i:saved){payload<<i<<' '<<state.valid[i]<<' '<<state.coarse[i];for(const auto& level:state.counts)payload<<' '<<level[i];for(double v:state.previous[i])payload<<' '<<v;payload<<'\n';}
    Atomic(e.root+"/adaptive_checkpoint.tsv","2 "+Hash(payload.str())+"\n"+payload.str());
}
bool RestoreAdaptive(Engine& e,AdaptiveState& state){
    const std::string path=e.root+"/adaptive_checkpoint.tsv";if(!Exists(path))return false;
    const std::string text=Read(path);const size_t newline=text.find('\n');std::istringstream header(text.substr(0,newline));int version;std::string hash;Require(bool(header>>version>>hash)&&version==2,"unsupported checkpoint schema");
    const std::string payload=text.substr(newline+1);Require(Hash(payload)==hash,"controller checkpoint changed");std::istringstream file(payload);std::string signature,phase;
    Require(bool(file>>signature>>phase>>state.round>>state.streak>>e.phi>>e.controls>>e.sharedMeans>>e.mirror>>e.mirrorError>>e.kernel.block>>e.kernel.pipeline>>e.kernel.warp)&&signature==e.signature,"checkpoint inputs changed");
    file>>e.depthSharing>>e.depthVerified>>e.local;size_t rows;file>>rows;e.weightTheta.resize(rows);e.weights.resize(rows);e.polarWeights.resize(rows);
    for(size_t r=0;r<rows;++r){file>>e.weightTheta[r];for(int j=0;j<16;++j)for(int c=0;c<2;++c){double v;Require(bool(file>>v)&&std::isfinite(v),"invalid checkpoint coefficient");if(j)e.polarWeights[r][c][j]=v;else e.weights[r][c]=v;}}
    if(!e.controls){e.weights.clear();e.polarWeights.clear();}const double savedMirror=e.mirrorError;e.SetMirror(e.mirror);e.mirrorError=savedMirror;size_t knots;file>>knots;state.knots.resize(knots);for(auto& k:state.knots)file>>k;
    Require(knots>=2&&std::is_sorted(state.knots.begin(),state.knots.end())&&state.knots.back()<e.theta.size(),"invalid checkpoint support");
    size_t levels,saved;file>>levels>>saved;Require(levels==state.counts.size()&&saved<=e.theta.size(),"checkpoint budget dimensions changed");
    for(size_t r=0;r<saved;++r){size_t i;bool valid,coarse;Require(bool(file>>i>>valid>>coarse)&&i<e.theta.size(),"invalid checkpoint row");state.valid[i]=valid;state.coarse[i]=coarse;
        for(auto& level:state.counts){int n;Require(bool(file>>n)&&n>0&&n<=e.maxCount&&(n&(n-1))==0,"invalid checkpoint count");level[i]=n;}
        for(double& v:state.previous[i])Require(bool(file>>v)&&std::isfinite(v),"invalid checkpoint matrix");}
    std::string extra;Require(!(file>>extra),"extra checkpoint fields");
    if(phase=="converged"&&Exists(e.root+"/final_outputs.fnv")){std::istringstream hashes(Read(e.root+"/final_outputs.fnv"));std::string path,digest;
        while(hashes>>path>>digest)Require(Hash(Read(e.root+"/"+path))==digest,"completed output changed");
        std::set<std::string> visited;std::istringstream sources(Read(e.root+"/final_stage_sources.tsv"));size_t level;int sign,n;unsigned seed;
        while(sources>>level>>sign>>n>>seed>>path)ValidateStoredJob(path,visited);
        e.completedCheckpoint=true;
    }else e.WarmCache();return true;
}
void BudgetFile(Engine& e,const AdaptiveState& state,const std::vector<size_t>& evaluated){std::ostringstream file;file<<"theta_deg";for(size_t l=0;l<state.counts.size();++l)file<<",N_level"<<l;file<<"\n";
    for(size_t i:evaluated){file<<Number(e.theta[i]);for(const auto& level:state.counts)file<<','<<level[i];file<<'\n';}Atomic(e.root+"/orientation_budget.csv",file.str());}
void HistoryRecord(Engine& e,int round,const std::string& record){
    std::ostringstream history;if(Exists(e.root+"/refinements.tsv")){std::istringstream old(Read(e.root+"/refinements.tsv"));std::string line;while(std::getline(old,line)){std::istringstream row(line);int id=-1;row>>id;if(id!=round)history<<line<<'\n';}}
    history<<record<<'\n';Atomic(e.root+"/refinements.tsv",history.str());
}
std::vector<double> WorkingGrid(const Engine& e,const std::vector<size_t>& evaluated){std::vector<double> grid;const size_t n=std::min(size_t(257),evaluated.size());
    for(size_t i=0;i<n;++i)grid.push_back(e.theta[evaluated[i*(evaluated.size()-1)/(n-1)]]);
    return grid;}
struct InterpolationResult {Samples output;std::vector<double> errors;std::vector<size_t> bad;double maximum=0;};
InterpolationResult InterpolateBudget(const Engine& e,const std::vector<size_t>& knots,const std::vector<size_t>& evaluated,const Samples& samples){
    InterpolationResult result;result.output.resize(8);result.errors.resize(evaluated.size());Samples predicted(8);const auto x=Angles(e.theta,evaluated),sx=Angles(e.theta,knots);const auto support=Support(samples,evaluated,knots);
    for(int s=0;s<8;++s){predicted[s]=LocalCubicInterpolate(sx,support[s],x);result.output[s]=LocalCubicInterpolate(sx,support[s],e.theta);}
    const auto delta=Subtract(predicted,samples);const auto mean=Mean(delta),variance=Variance(delta),reference=Mean(samples),prediction=Mean(predicted);
    for(size_t r=0;r<evaluated.size();++r){Require(reference[r][0]>0,"nonpositive interpolation reference");for(int j=0;j<16;++j)result.errors[r]=std::max(result.errors[r],(std::fabs(mean[r][j])+t95*std::sqrt(variance[r][j]/8))/reference[r][0]);
        if(e.ratios&&prediction[r][0]>0)for(int j=1;j<16;++j){const double a=prediction[r][j]/prediction[r][0],b=reference[r][j]/reference[r][0];double variance=0;
            for(int s=0;s<8;++s){const double influence=(predicted[s][r][j]-a*predicted[s][r][0])/prediction[r][0]-(samples[s][r][j]-b*samples[s][r][0])/reference[r][0];variance+=influence*influence/7;}
            result.errors[r]=std::max(result.errors[r],std::fabs(a-b)+t95*std::sqrt(variance/8));}
        result.maximum=std::max(result.maximum,result.errors[r]);if(result.errors[r]>e.interpolation)result.bad.push_back(evaluated[r]);}
    return result;
}
InterpolationResult CommonPrefixInterpolation(Engine& e,const AdaptiveState& state,const std::vector<int>& depths,const std::vector<size_t>& evaluated,const unsigned* seeds){
    InterpolationResult result;result.errors.assign(evaluated.size(),0);const size_t width=std::min(size_t(4),state.knots.size());
    using Key=std::vector<int>;std::map<Key,std::set<size_t>> groups;
    for(size_t interval=0;interval+1<state.knots.size();++interval){const size_t left=std::min(interval?interval-1:0,state.knots.size()-width),right=std::min(interval,state.knots.size()-width);
        std::set<size_t> nodes;for(size_t k=left;k<right+width;++k)nodes.insert(state.knots[k]);for(size_t id:evaluated)if(id>=state.knots[interval]&&id<=state.knots[interval+1])nodes.insert(id);
        Key counts(depths.size(),0);for(size_t l=0;l<depths.size();++l)for(size_t id:nodes)counts[l]=std::max(counts[l],state.counts[l][id]);groups[counts].insert(nodes.begin(),nodes.end());}
    std::ostringstream sources,budgets;budgets<<"group,theta_deg";for(size_t l=0;l<depths.size();++l)budgets<<",N_level"<<l;budgets<<'\n';size_t groupNumber=0;
    for(const auto& group:groups){const std::vector<size_t> ids(group.second.begin(),group.second.end());Budgets ns(depths.size());for(size_t l=0;l<depths.size();++l)ns[l].assign(e.theta.size(),group.first[l]);
        auto stage=EvaluateBudget(e,depths,ns,ids,seeds);const Samples data=Combine(stage.levels);const Matrix reference=Mean(data);const auto x=Angles(e.theta,ids);
        for(size_t interval=0;interval+1<state.knots.size();++interval){const size_t left=std::min(interval?interval-1:0,state.knots.size()-width),right=std::min(interval,state.knots.size()-width);
            Key counts(depths.size(),0);std::set<size_t> nodes;for(size_t k=left;k<right+width;++k)nodes.insert(state.knots[k]);for(size_t id:evaluated)if(id>=state.knots[interval]&&id<=state.knots[interval+1])nodes.insert(id);
            for(size_t l=0;l<depths.size();++l)for(size_t id:nodes)counts[l]=std::max(counts[l],state.counts[l][id]);
            if(counts!=group.first)continue;
            std::vector<double> knots;Samples supports(8);for(size_t k=left;k<right+width;++k){const size_t row=std::lower_bound(ids.begin(),ids.end(),state.knots[k])-ids.begin();knots.push_back(e.theta[state.knots[k]]);for(int s=0;s<8;++s)supports[s].push_back(data[s][row]);}
            // Use the globally selected endpoint derivatives: their local
            // four-node windows are entirely contained in this stencil.
            std::vector<double> query;std::vector<size_t> rows,positions;for(size_t r=0;r<evaluated.size();++r)if(evaluated[r]>=state.knots[interval]&&evaluated[r]<=state.knots[interval+1]){query.push_back(e.theta[evaluated[r]]);rows.push_back(std::lower_bound(ids.begin(),ids.end(),evaluated[r])-ids.begin());positions.push_back(r);}
            Samples predicted(8);for(int s=0;s<8;++s)predicted[s]=LocalCubicInterpolate(knots,supports[s],query);const Matrix pm=Mean(predicted);
            for(size_t r=0;r<rows.size();++r){double error=0;const auto& ref=reference[rows[r]];Require(ref[0]>0,"nonpositive common-prefix interpolation reference");
                for(int j=0;j<16;++j){double mean=0,var=0;for(int s=0;s<8;++s)mean+=(predicted[s][r][j]-data[s][rows[r]][j])/8;for(int s=0;s<8;++s)var+=std::pow(predicted[s][r][j]-data[s][rows[r]][j]-mean,2)/7;error=std::max(error,(std::fabs(mean)+t95*std::sqrt(var/8))/ref[0]);
                    if(e.ratios&&j&&pm[r][0]>0){const double a=pm[r][j]/pm[r][0],b=ref[j]/ref[0];double v=0;for(int s=0;s<8;++s){const double influence=(predicted[s][r][j]-a*predicted[s][r][0])/pm[r][0]-(data[s][rows[r]][j]-b*data[s][rows[r]][0])/ref[0];v+=influence*influence/7;}error=std::max(error,std::fabs(a-b)+t95*std::sqrt(v/8));}}
                result.errors[positions[r]]=std::max(result.errors[positions[r]],error);}
        }
        std::istringstream records(stage.sources);std::string line;while(std::getline(records,line))sources<<groupNumber<<' '<<line<<'\n';
        for(size_t id:ids){budgets<<groupNumber<<','<<Number(e.theta[id]);for(int n:group.first)budgets<<','<<n;budgets<<'\n';}++groupNumber;
    }
    for(size_t r=0;r<evaluated.size();++r){result.maximum=std::max(result.maximum,result.errors[r]);if(result.errors[r]>e.interpolation)result.bad.push_back(evaluated[r]);}
    std::ostringstream seedNames;for(int s=0;s<8;++s)seedNames<<seeds[s]<<' ';const std::string label=seeds[0]==e.seeds[0]?"production":"validation";
    Atomic(e.root+"/interpolation_"+label+"_sources.tsv",sources.str());Atomic(e.root+"/interpolation_"+label+"_budgets.csv",budgets.str());return result;
}
void SaveReplicaMatrices(Engine& e,const std::string& label,const Samples& data,const std::vector<double>& x,const unsigned* seeds){Directory(e.root+"/replicas");for(int s=0;s<8;++s)WriteMatrix(e.root+"/replicas/"+label+"_seed"+Number(seeds[s])+".dat",x,data[s]);}
bool IndependentAudit(Engine& e,AdaptiveState& state,const std::vector<int>& depths,const Samples& productionSamples,const std::vector<size_t>& productionIds){
    if(!e.finalAudit){Atomic(e.root+"/independent_final_validation.json","{\"enabled\":false,\"passed\":false}");return true;}
    std::array<unsigned,8> seeds;for(int s=0;s<8;++s){seeds[s]=1000000000u+unsigned(state.round)*1024u+unsigned(s)*17u;Require(std::find(e.seeds.begin(),e.seeds.end(),seeds[s])==e.seeds.end(),"final/production seed overlap");}
    auto ids=AdaptiveGuards(state.knots,e.theta.size(),true);
    if(e.theta.size()<=4097){ids.resize(e.theta.size());for(size_t i=0;i<ids.size();++i)ids[i]=i;}
    for(size_t i:ids)if(!std::binary_search(productionIds.begin(),productionIds.end(),i)){
        auto upper=std::upper_bound(productionIds.begin(),productionIds.end(),i);const size_t right=upper==productionIds.end()?productionIds.back():*upper,left=upper==productionIds.begin()?productionIds.front():*(upper-1);
        for(auto& level:state.counts)level[i]=std::max(level[left],level[right]);}
    auto data=EvaluateBudget(e,depths,state.counts,ids,seeds.data());const Samples combined=Combine(data.levels);const auto stats=Estimate(combined);auto interp=InterpolateBudget(e,state.knots,ids,combined);
    const bool common=e.local&&VaryingBudgets(state.counts,ids);
    Atomic(e.root+"/interpolation_validation_mode.json","{\"common_prefix\":"+std::string(common?"true":"false")+"}");
    if(common){auto checked=CommonPrefixInterpolation(e,state,depths,ids,seeds.data());interp.errors=checked.errors;interp.maximum=checked.maximum;interp.bad=checked.bad;}
    const auto outputStats=Estimate(interp.output);const auto production=Mean(productionSamples);
    // Compare frozen estimators at common direct nodes; no seed-wise ratios.
    Matrix matched;for(size_t id:productionIds){size_t r=std::lower_bound(ids.begin(),ids.end(),id)-ids.begin();matched.push_back(stats.mean[r]);}
    const double change=Difference(production,matched,production);
    const bool pass=StatisticalError(e,stats)<=e.statistical&&StatisticalError(e,outputStats)<=e.statistical&&interp.maximum<=e.interpolation&&change<=e.eps;
    std::ostringstream report;report<<"{\"enabled\":true,\"passed\":"<<(pass?"true":"false")<<",\"fresh_seeds_per_attempt\":true,\"seeds\":[";
    for(int s=0;s<8;++s)report<<(s?",":"")<<seeds[s];
    report<<"],\"max_pointwise95_scaled_M11\":"<<Number(std::max(stats.maxInterval,outputStats.maxInterval))<<",\"max_ratio95\":"<<Number(std::max(stats.maxRatioInterval,outputStats.maxRatioInterval))<<",\"max_interpolation_residual95\":"<<Number(interp.maximum)<<",\"difference_to_production_scaled_M11\":"<<Number(change)<<",\"checked_theta\":"<<ids.size()<<",\"scope\":\"fresh-seed frozen-estimator pointwise validation; all requested nodes up to 4097, otherwise quarter interval guards; no simultaneous-confidence claim\"}";
    Atomic(e.root+"/independent_final_validation.json",report.str());Atomic(e.root+"/independent_final_sources.tsv",data.sources);SaveReplicaMatrices(e,"validation_evaluated",combined,Angles(e.theta,ids),seeds.data());
    WriteGrid(e.root+"/independent_theta.csv",Angles(e.theta,ids));
    if(!pass){std::set<size_t> next(state.knots.begin(),state.knots.end());next.insert(interp.bad.begin(),interp.bad.end());state.knots.assign(next.begin(),next.end());}
    e.independentPassed=pass;return pass;
}
bool PhysicsAudit(Engine& e,const AdaptiveState& state,const std::vector<size_t>& evaluated,const Matrix& reference){
    if(!e.physicsAudit){Atomic(e.root+"/physics_validation.json","{\"enabled\":false,\"passed\":false,\"requested_physics_retained\":true}");return true;}
    std::vector<size_t> selected;const size_t size=std::min(size_t(65),evaluated.size());for(size_t r=0;r<size;++r)selected.push_back(evaluated[r*(evaluated.size()-1)/(size-1)]);
    const auto x=Angles(e.theta,selected);Matrix scale;for(size_t id:selected)scale.push_back(reference[std::lower_bound(evaluated.begin(),evaluated.end(),id)-evaluated.begin()]);
    std::array<unsigned,8> seeds;for(int s=0;s<8;++s){seeds[s]=1500000000u+unsigned(state.round)*1024u+unsigned(s)*17u;Require(std::find(e.seeds.begin(),e.seeds.end(),seeds[s])==e.seeds.end(),"physics/production seed overlap");}
    int count=e.minimum;for(size_t l=1;l<state.counts.size();++l)for(size_t id:selected)count=std::max(count,state.counts[l][id]);
    const auto original=e.base;const auto savedCutoff=e.selectedCutoff;bool pass=false;double depthError=0,cutoffError=0;
    auto bound=[&](const Samples& difference){const auto mean=Mean(difference),variance=Variance(difference);double error=0;for(size_t r=0;r<x.size();++r)for(int j=0;j<16;++j)error=std::max(error,(std::fabs(mean[r][j])+t95*std::sqrt(variance[r][j]/8))/scale[r][0]);return error;};
    for(int attempt=0;attempt<5;++attempt){
        const auto target=e.DirectStage(count,seeds.data(),8,e.target,x),deeper=e.DirectStage(count,seeds.data(),8,e.target+2,x);
        depthError=bound(Subtract(e.Data(deeper),e.Data(target)));
        e.SetCutoff("off");
        try{const auto uncut=e.DirectStage(count,seeds.data(),8,e.target,x);e.base=original;e.selectedCutoff=savedCutoff;e.SetMirror(e.mirror);cutoffError=bound(Subtract(e.Data(uncut),e.Data(target)));}
        catch(...){e.base=original;e.selectedCutoff=savedCutoff;e.SetMirror(e.mirror);throw;}
        pass=depthError<=e.eps&&cutoffError<=e.eps;if(pass||count>=e.maxCount||attempt==4)break;count=std::min(e.maxCount,count>e.maxCount/2?e.maxCount:2*count);
    }
    Atomic(e.root+"/physics_validation.json","{\"enabled\":true,\"passed\":"+std::string(pass?"true":"false")+",\"depth_plus2_residual95_scaled_M11\":"+Number(depthError)+",\"cutoffs_off_residual95_scaled_M11\":"+Number(cutoffError)+",\"count_per_seed\":"+Number(count)+",\"checked_theta\":"+Number(x.size())+",\"requested_physics_retained\":true,\"scope\":\"representative angles; depth+2 and no-cutoff sensitivity estimates; no ADDA or infinite-depth guarantee\"}");return pass;
}
int RunAdaptive(Engine& e){
    std::vector<int> depths;for(int d:{std::max(1,e.target-4),std::max(1,e.target-2),e.target})if(depths.empty()||depths.back()!=d)depths.push_back(d);
    AdaptiveState state;state.counts=Budgets(depths.size(),std::vector<int>(e.theta.size(),e.minimum));state.counts[0].assign(e.theta.size(),e.initial);
    state.previous=Matrix(e.theta.size());state.valid.assign(e.theta.size(),false);state.coarse.assign(e.theta.size(),false);
    std::vector<double> tuningAngles;const size_t tuningRows=std::min(e.theta.size(),size_t(33));
    for(size_t i=0;i<tuningRows;++i)tuningAngles.push_back(e.theta[i*(e.theta.size()-1)/(tuningRows-1)]);
    e.PrepareTuning(tuningAngles);
    const bool resumed=RestoreAdaptive(e,state);
    if(e.completedCheckpoint){
        const std::string path=e.root+"/fullauto_status.json";if(Exists(path)){std::string text=Read(path);const size_t key=text.find("\"status\":\"");if(key!=std::string::npos){const size_t begin=key+10,end=text.find('"',begin);text.replace(begin,end-begin,"converged");Atomic(path,text);}}
        std::cout<<"fullauto: completed outputs and dependencies verified; no calculation repeated"<<std::endl;return 0;
    }
    if(!resumed){const size_t initial=std::min(e.theta.size(),size_t(e.startTheta));for(size_t i=0;i<initial;++i)state.knots.push_back(i*(e.theta.size()-1)/(initial-1));
        auto ids=AdaptiveGuards(state.knots,e.theta.size());auto x=Angles(e.theta,ids);e.State("calibrating",0,PeakBudgets(state.counts),ids.size(),0,0,0);
        auto calibrate=[&](){e.PrepareMirror(x);e.Calibrate(x);if(!e.VerifyMirror(x)){e.SetMirror(false);e.Calibrate(x);}};
        try{calibrate();if(e.thetaAuditCount)state.knots=e.AuditThetaSupport(state.knots);}catch(const std::exception& error){if(!e.DisableUnsupportedControls(error))throw;calibrate();if(e.thetaAuditCount)state.knots=e.AuditThetaSupport(state.knots);}
        ids=AdaptiveGuards(state.knots,e.theta.size());if(e.phiRecheck&&!e.args.IsCatched("nphi")&&ids.size()>2*x.size()){
            Atomic(e.root+"/initial_calibration.json",Read(e.root+"/calibration.json"));e.Calibrate(WorkingGrid(e,ids),true);if(!e.VerifyMirror(WorkingGrid(e,ids))){e.SetMirror(false);e.Calibrate(WorkingGrid(e,ids),true);}}
        // On small grids startup, matrix assembly and tracing dominate. The
        // measured regressions do not justify grouping there; retain a
        // common budget. Explicit 'on' remains available for experiments.
        if(e.localPolicy=="auto"&&e.theta.size()*size_t(e.phi)<=256)e.local=false;
        SaveAdaptive(e,state);
    }
    double lastCI=0,lastInterpolation=0;size_t lastNodes=0;
    for(;state.round<e.roundLimit;++state.round){
        const auto ids=AdaptiveGuards(state.knots,e.theta.size());Require(ids.size()<=size_t(e.maxTheta),"theta budget exhausted");const auto x=Angles(e.theta,ids);lastNodes=ids.size();SaveAdaptive(e,state);
        e.State("active",state.round,PeakBudgets(state.counts),ids.size(),lastCI,lastInterpolation,state.streak);
        if(!e.VerifyMirror(x)){e.SetMirror(false);e.Calibrate(WorkingGrid(e,ids),true);std::fill(state.valid.begin(),state.valid.end(),false);std::fill(state.coarse.begin(),state.coarse.end(),false);state.streak=0;SaveAdaptive(e,state);--state.round;continue;}
        e.VerifyDepthSharing(WorkingGrid(e,ids));
        AdaptiveData stage;Samples data;Statistics stats;
        try{stage=EvaluateBudget(e,depths,state.counts,ids,e.seeds.data());data=Combine(stage.levels);stats=Estimate(data);
            e.lastAnalyticError=0;for(size_t r=0;r<ids.size();++r)e.lastAnalyticError=std::max(e.lastAnalyticError,stage.analytic[r]/stats.mean[r][0]);
            Require(!e.controls||e.lastAnalyticError<=e.eps/64,"analytic mean numerical budget exceeded; use complete sampled PO or refine the analytic provider");}
        catch(const std::exception& error){if(!e.DisableUnsupportedControls(error))throw;e.Calibrate(WorkingGrid(e,ids),true);std::fill(state.coarse.begin(),state.coarse.end(),false);std::fill(state.valid.begin(),state.valid.end(),false);state.streak=0;SaveAdaptive(e,state);--state.round;continue;}
        Budgets half=state.counts;std::vector<size_t> check;for(size_t id:ids)if(!state.coarse[id]&&state.counts[0][id]>=2*e.initial){half[0][id]/=2;check.push_back(id);}
        if(!check.empty()){auto previous=EvaluateBudget(e,{depths[0]},{half[0]},check,e.seeds.data());const auto pm=Mean(previous.levels[0]),cm=Mean(stage.levels[0]);
            for(size_t k=0;k<check.size();++k){const size_t r=std::lower_bound(ids.begin(),ids.end(),check[k])-ids.begin();double change=0;for(int j=0;j<16;++j)change=std::max(change,std::fabs(cm[r][j]-pm[k][j])/stats.mean[r][0]);state.coarse[check[k]]=change<=e.statistical;}}
        auto interp=InterpolateBudget(e,state.knots,ids,data);
        const bool common=e.local&&VaryingBudgets(state.counts,ids);
        Atomic(e.root+"/interpolation_production_mode.json","{\"common_prefix\":"+std::string(common?"true":"false")+"}");
        if(common){auto checked=CommonPrefixInterpolation(e,state,depths,ids,e.seeds.data());interp.errors=checked.errors;interp.maximum=checked.maximum;interp.bad=checked.bad;}
        bool positive=true;for(const auto& m:Mean(interp.output))positive=positive&&m[0]>0;
        std::vector<Matrix> levelVariance;for(const auto& level:stage.levels)levelVariance.push_back(Variance(level));
        Statistics output;if(positive)output=Estimate(interp.output);double outputCI=positive?StatisticalError(e,output):1e100;
        double change=0;bool previous=true,coarse=true;for(size_t r=0;r<ids.size();++r){const size_t id=ids[r];previous=previous&&state.valid[id];coarse=coarse&&state.coarse[id];if(state.valid[id])for(int j=0;j<16;++j)change=std::max(change,std::fabs(stats.mean[r][j]-state.previous[id][j])/stats.mean[r][0]);}
        if(!previous)change=1e100;
        const bool pass=StatisticalError(e,stats)<=e.statistical&&outputCI<=e.statistical&&interp.maximum<=e.interpolation&&change<=e.statistical&&coarse;
        state.streak=pass?state.streak+1:0;lastCI=std::max(stats.maxInterval,positive?output.maxInterval:1e100);lastInterpolation=interp.maximum;e.lastRatio=std::max(stats.maxRatioInterval,positive?output.maxRatioInterval:1e100);
        BudgetFile(e,state,ids);Atomic(e.root+"/final_stage_sources.tsv",stage.sources);
        std::ostringstream record;record<<state.round<<' '<<PeakBudgets(state.counts)[0]<<' '<<ids.size()<<' '<<state.knots.size()<<' '<<Number(stats.maxInterval)<<' '<<Number(outputCI)<<' '<<Number(interp.maximum)<<' '<<Number(change)<<' '<<state.streak;HistoryRecord(e,state.round,record.str());
        if(state.streak>=e.stablePasses){e.State("validating",state.round,PeakBudgets(state.counts),ids.size(),lastCI,lastInterpolation,state.streak);
            if(IndependentAudit(e,state,depths,data,ids)){
                if(!e.AuditAutomaticCutoffs(x,PeakBudgets(state.counts)[0],state.round)){
                    e.Calibrate(x,true);state.streak=0;std::fill(state.valid.begin(),state.valid.end(),false);
                    std::fill(state.coarse.begin(),state.coarse.end(),false);SaveAdaptive(e,state);continue;
                }
                WriteMatrix(e.root+"/mueller_fullauto.dat",e.theta,output.mean);WriteMatrix(e.root+"/mueller_evaluated.dat",x,stats.mean);WriteGrid(e.root+"/theta_support.csv",Angles(e.theta,state.knots));WriteGrid(e.root+"/theta_evaluated.csv",x);
                SaveReplicaMatrices(e,"production_evaluated",data,x,e.seeds.data());SaveReplicaMatrices(e,"production_interpolated",interp.output,e.theta,e.seeds.data());
                std::ostringstream accuracy,ratio,guards;accuracy<<"theta_deg";ratio<<"theta_deg";for(int j=0;j<16;++j){accuracy<<",estimated95_M"<<j/4+1<<j%4+1<<"_scaled_M11";ratio<<",estimated95_M"<<j/4+1<<j%4+1<<"_over_M11";}accuracy<<'\n';ratio<<'\n';
                for(size_t i=0;i<e.theta.size();++i){accuracy<<Number(e.theta[i]);ratio<<Number(e.theta[i]);for(int j=0;j<16;++j){accuracy<<','<<Number(output.interval[i][j]);ratio<<','<<Number(output.ratioInterval[i][j]);}accuracy<<'\n';ratio<<'\n';}
                guards<<"theta_deg,is_support,M11,estimated_guard_error_plus95_scaled_M11\n";for(size_t r=0;r<ids.size();++r)guards<<Number(x[r])<<','<<std::binary_search(state.knots.begin(),state.knots.end(),ids[r])<<','<<Number(stats.mean[r][0])<<','<<Number(interp.errors[r])<<'\n';
                Atomic(e.root+"/all_mueller_confidence.csv",accuracy.str());Atomic(e.root+"/normalized_mueller_confidence.csv",ratio.str());Atomic(e.root+"/theta_validation.csv",guards.str());
                const bool physics=PhysicsAudit(e,state,ids,stats.mean);
                Atomic(e.root+"/precision_contract.json","{\"statistics\":\"eight-scramble pointwise Student-t estimates\",\"normalized_ratios\":\"paired-covariance Fieller intervals for ratios of means\",\"interpolation_budget\":"+Number(e.interpolation)+",\"execution_budget\":"+Number((e.resourcePolicy=="auto"||e.config.useGpu)?e.eps/1024:0)+",\"rotation_budget\":"+Number((e.rotationBeta>1||e.rotationGamma>1)?e.eps/256:0)+",\"cutoff_budget\":"+Number(e.selectedCutoff!="requested"&&e.selectedCutoff!="off"?e.eps/8:0)+",\"selected_cutoff\":"+Quoted(e.selectedCutoff)+",\"automatic_cutoff_audit\":"+std::string(e.tuneCutoffs&&e.selectedCutoff!="off"?"true":"false")+",\"physical_depth\":"+Number(e.target)+",\"physical_accuracy_audited\":"+std::string(e.physicsAudit?"true":"false")+",\"analytic_mean_error\":\"reported internal refinement, not an absolute bound\",\"global_simultaneous_confidence\":false}");
                if(!physics){SaveAdaptive(e,state,"physics_validation_failed");e.State("physics_validation_failed",state.round,PeakBudgets(state.counts),ids.size(),lastCI,lastInterpolation,state.streak);return 4;}
                std::ostringstream finalHashes;for(const char* path:{"mueller_fullauto.dat","mueller_evaluated.dat","all_mueller_confidence.csv","normalized_mueller_confidence.csv","theta_support.csv","theta_evaluated.csv","orientation_budget.csv","final_stage_sources.tsv","independent_final_validation.json","precision_contract.json","physics_validation.json","execution_tuning.tsv","resource_calibration.json","rotation_verification.json","cutoff_calibration.json","automatic_cutoff_validation.json"})finalHashes<<path<<' '<<Hash(Read(e.root+"/"+path))<<'\n';Atomic(e.root+"/final_outputs.fnv",finalHashes.str());
                SaveAdaptive(e,state,"converged");e.State("converged",state.round,PeakBudgets(state.counts),ids.size(),lastCI,lastInterpolation,state.streak);return 0;
            }
            state.streak=0;std::fill(state.valid.begin(),state.valid.end(),false);for(auto& level:state.counts)for(size_t id:ids)level[id]=std::min(e.maxCount,level[id]>e.maxCount/2?e.maxCount:2*level[id]);SaveAdaptive(e,state);continue;
        }
        for(size_t r=0;r<ids.size();++r){state.previous[ids[r]]=stats.mean[r];state.valid[ids[r]]=true;}
        if(!interp.bad.empty()||!positive){std::set<size_t> next(state.knots.begin(),state.knots.end());next.insert(interp.bad.begin(),interp.bad.end());if(!positive)next.insert(ids.begin(),ids.end());
            state.knots.assign(next.begin(),next.end());for(size_t id:AdaptiveGuards(state.knots,e.theta.size()))if(!std::binary_search(ids.begin(),ids.end(),id)){
                auto upper=std::upper_bound(ids.begin(),ids.end(),id);const size_t right=upper==ids.end()?ids.back():*upper,left=upper==ids.begin()?ids.front():*(upper-1);for(auto& level:state.counts)level[id]=std::max(level[left],level[right]);}state.streak=0;SaveAdaptive(e,state);continue;}
        Budgets proposed=state.counts;bool grew=false;
        for(size_t r=0;r<ids.size();++r){const size_t id=ids[r];double error=0;for(int j=0;j<16;++j){error=std::max(error,stats.interval[r][j]);if(e.ratios)error=std::max(error,stats.ratioInterval[r][j]);}
            if(error>e.statistical){auto ns=AllocateRow(e,stage,state.counts,ids,r,stats.mean);for(size_t l=0;l<ns.size();++l)proposed[l][id]=std::max(state.counts[l][id],std::min(e.maxCount,std::min(state.counts[l][id]>e.maxCount/4?e.maxCount:4*state.counts[l][id],ns[l])));}
            if(!state.coarse[id])proposed[0][id]=std::min(e.maxCount,std::max(proposed[0][id],state.counts[0][id]>e.maxCount/2?e.maxCount:2*state.counts[0][id]));
        }
        if(!e.local){for(size_t l=0;l<proposed.size();++l){int maximum=0;for(size_t id:ids)maximum=std::max(maximum,proposed[l][id]);for(size_t id:ids)proposed[l][id]=maximum;}}
        for(size_t l=0;l<proposed.size();++l)for(size_t id:ids)grew=grew||proposed[l][id]>state.counts[l][id];
        if(outputCI>e.statistical&&positive){
            for(size_t q=0;q<e.theta.size();++q){double error=0;for(int j=0;j<16;++j){error=std::max(error,output.interval[q][j]);if(e.ratios)error=std::max(error,output.ratioInterval[q][j]);}if(error<=e.statistical)continue;
                size_t near=std::upper_bound(state.knots.begin(),state.knots.end(),q)-state.knots.begin();const size_t begin=near>2?near-2:0,end=std::min(state.knots.size(),near+2);
                for(size_t k=begin;k<end;++k){const size_t id=state.knots[k],row=std::lower_bound(ids.begin(),ids.end(),id)-ids.begin();size_t chosen=0;double largest=-1;
                    for(size_t l=0;l<stage.levels.size();++l){double score=0;for(int j=0;j<16;++j)score=std::max(score,levelVariance[l][row][j]);if(score>largest){largest=score;chosen=l;}}
                    proposed[chosen][id]=std::max(proposed[chosen][id],std::min(e.maxCount,state.counts[chosen][id]>e.maxCount/2?e.maxCount:2*state.counts[chosen][id]));}
            }
            if(!e.local)for(size_t l=0;l<proposed.size();++l){int maximum=0;for(size_t id:ids)maximum=std::max(maximum,proposed[l][id]);for(size_t id:ids)proposed[l][id]=maximum;}
            for(size_t l=0;l<proposed.size();++l)for(size_t id:ids)grew=grew||proposed[l][id]>state.counts[l][id];
        }
        if(!grew){ // A fresh prefix is required for a stability pass or spline uncertainty.
            size_t worst=0;double maximum=-1;for(size_t r=0;r<ids.size();++r)for(int j=0;j<16;++j)if(stats.interval[r][j]>maximum){maximum=stats.interval[r][j];worst=r;}
            const size_t id=ids[worst];size_t chosen=0;double largest=-1;for(size_t l=0;l<stage.levels.size();++l){double v=0;for(int j=0;j<16;++j)v=std::max(v,levelVariance[l][worst][j]);if(v>largest){largest=v;chosen=l;}}
            auto grow=[&](size_t i){proposed[chosen][i]=std::min(e.maxCount,state.counts[chosen][i]>e.maxCount/2?e.maxCount:2*state.counts[chosen][i]);};if(e.local)grow(id);else for(size_t i:ids)grow(i);
        }
        GateLocalCost(e,proposed,depths,ids);
        if(proposed==state.counts)break;
        state.counts=proposed;SaveAdaptive(e,state);
    }
    e.State("requires_further_sampling",state.round,PeakBudgets(state.counts),lastNodes,lastCI,lastInterpolation,state.streak);return 3;
}
} // anonymous namespace
