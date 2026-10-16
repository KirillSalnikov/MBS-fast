// Included only by FullAuto.cpp: native adaptive scheduling and durable state.
namespace {
using Budgets=std::vector<std::vector<int>>;
struct AdaptiveData {std::vector<Samples> levels;std::vector<std::vector<double>> costs;std::vector<double> analytic;std::string sources;};
std::vector<int> PeakBudgets(const Budgets& n){std::vector<int> result;for(const auto& level:n)result.push_back(*std::max_element(level.begin(),level.end()));return result;}
double StatisticalError(const Engine& e,const Statistics& stats){return e.ratios?std::max(stats.maxInterval,stats.maxRatioInterval):stats.maxInterval;}
bool VaryingBudgets(const Budgets& ns,const std::vector<size_t>& ids){for(const auto& level:ns)for(size_t id:ids)if(level[id]!=level[ids.front()])return true;return false;}
#include "FullAutoWork.inc"
void GateLocalCost(Engine& e,AdditionalWork& work,Budgets& proposed,const std::vector<int>& depths,const std::vector<size_t>& ids){
    if(e.localPolicy!="auto"||!e.local)return;
    std::ostringstream report;report<<"{\"policy\":\"additional native work\",\"scope\":\"sum of worker seconds; excludes controller assembly; not wall time\",\"levels\":[";
    for(size_t l=0;l<proposed.size();++l){Budgets uniform=proposed;int maximum=0;for(size_t id:ids)maximum=std::max(maximum,proposed[l][id]);for(size_t id:ids)uniform[l][id]=maximum;
        const auto local=work.Forecast(BudgetRequests(e,depths,proposed,ids)),global=work.Forecast(BudgetRequests(e,depths,uniform,ids));
        const bool unified=local.profiled&&global.profiled&&local.seconds>global.seconds*1.10;
        if(unified)proposed=uniform;
        report<<(l?",":"")<<"{\"level\":"<<l<<",\"profiled\":"<<(local.profiled&&global.profiled?"true":"false")
            <<",\"additional_local_worker_seconds\":"<<Number(local.seconds)<<",\"additional_global_worker_seconds\":"<<Number(global.seconds)
            <<",\"local_native_jobs\":"<<local.jobs<<",\"global_native_jobs\":"<<global.jobs<<",\"unified\":"<<(unified?"true":"false")<<"}";
    }report<<"]}";e.Record(e.root+"/local_cost_gate.json",report.str());
}
std::vector<size_t> AdaptiveGuards(const std::vector<size_t>& knots,size_t size,bool final=false){
    if(final)return Guards(knots,size);
    std::set<size_t> all(knots.begin(),knots.end());
    for(size_t i=1;i<knots.size();++i){size_t q=(knots[i-1]+knots[i])/2;if(q>knots[i-1]&&q<knots[i])all.insert(q);}return {all.begin(),all.end()};
}
AdaptiveData EvaluateBudget(Engine& e,const std::vector<int>& depths,const Budgets& counts,const std::vector<size_t>& ids,const unsigned* seeds){
    const auto x=Angles(e.theta,ids);AdaptiveData result;result.analytic.assign(ids.size(),0);std::ostringstream sources;
    e.Prefetch(BudgetRequests(e,depths,counts,ids),seeds);
    for(size_t level=0;level<depths.size();++level){Samples data(8,Matrix(ids.size()));std::vector<double> cost(ids.size());std::map<int,std::vector<size_t>> groups;
        for(size_t r=0;r<ids.size();++r)groups[counts[level][ids[r]]].push_back(r);
        for(const auto& group:groups){std::vector<double> subset;for(size_t r:group.second)subset.push_back(x[r]);
            if(subset.size()==1){subset.push_back(x.size()>1?(x[0]==subset[0]?x.back():x[0]):(e.theta.front()==subset[0]?e.theta.back():e.theta.front()));std::sort(subset.begin(),subset.end());}
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
std::vector<int> AllocateRow(const Engine& e,const AdaptiveData& data,const Budgets& counts,const std::vector<size_t>& ids,size_t row,const Matrix& mean,double tolerance){
    std::vector<Samples> values;std::vector<double> costs;std::vector<int> ns;
    for(size_t l=0;l<data.levels.size();++l){Samples v(8,Matrix(1));for(int s=0;s<8;++s)v[s][0]=data.levels[l][s][row];values.push_back(v);costs.push_back(data.costs[l][row]);ns.push_back(counts[l][ids[row]]);}
    return AllocateCovariant(values,costs,ns,{mean[row]},tolerance*.85,e.ratios,e.maxCount,tolerance);
}
#include "FullAutoRefinement.inc"
struct AdaptiveState {int round=0,streak=0;unsigned auditSerial=0;std::vector<size_t> knots;Budgets counts;Matrix previous;std::vector<bool> valid,coarse;};
std::string SaveRoundDiagnostics(Engine& e,const AdaptiveState& state,const AdaptiveData& stage,const Statistics& stats,const std::vector<size_t>& ids,const std::vector<double>& errors,const std::vector<double>& x){
    const std::string folder=e.root+"/round_snapshots/round"+Number(state.round)+"_"+Hash(stage.sources);
    Directory(folder);WriteMatrix(folder+"/mueller_evaluated.dat",x,stats.mean);WriteGrid(folder+"/theta_support.csv",Angles(e.theta,state.knots));
    Atomic(folder+"/stage_sources.tsv",stage.sources);
    std::ostringstream table;table<<"theta_deg,is_support,interpolation_residual95,previous_available,half_prefix_stable";
    for(size_t l=0;l<state.counts.size();++l)table<<",N_level"<<l;
    for(int j=0;j<16;++j)table<<",estimated95_M"<<j/4+1<<j%4+1<<"_scaled_M11";
    for(int j=1;j<16;++j)table<<",estimated95_M"<<j/4+1<<j%4+1<<"_over_M11";table<<'\n';
    for(size_t r=0;r<ids.size();++r){const size_t id=ids[r];table<<Number(x[r])<<','<<std::binary_search(state.knots.begin(),state.knots.end(),id)<<','<<Number(errors[r])<<','<<state.valid[id]<<','<<state.coarse[id];
        for(const auto& level:state.counts)table<<','<<level[id];for(double v:stats.interval[r])table<<','<<Number(v);for(int j=1;j<16;++j)table<<','<<Number(stats.ratioInterval[r][j]);table<<'\n';}
    Atomic(folder+"/parameters_accuracy.csv",table.str());
    e.Log("round_diagnostics","{\"round\":"+Number(state.round)+",\"folder\":"+Quoted(folder)+",\"sources_fnv1a64\":"+Quoted(Hash(stage.sources))+"}");
    return folder;
}
const char* BudgetPolicy();
void UpdatePassStreak(Engine& e,AdaptiveState& state,bool passed,const AdaptiveData& data,const std::string& snapshot){
    const std::string key=Hash(e.ExecutionHash()+data.sources+Number(e.statistical)+" "+Number(e.interpolation)+" "+BudgetPolicy()+" "+Number(passed));
    const std::string path=e.root+"/production_validation_receipt.tsv";bool reused=false;int result=std::min(state.streak,e.stablePasses);
    if(Exists(path)){const auto text=Read(path);const auto end=text.find('\n');Require(end!=std::string::npos&&Hash(text.substr(end+1))==text.substr(0,end),"production validation receipt changed");
        std::istringstream in(text.substr(end+1));std::string signature,savedKey,extra;int count;bool savedPass;
        Require(bool(in>>signature>>savedKey>>count>>savedPass)&&!(in>>extra)&&signature==e.signature&&count>=0,"invalid production validation receipt");
        if(savedKey==key&&savedPass==passed){result=std::min(count,e.stablePasses);reused=true;}
    }else if(passed&&state.streak>0&&Exists(e.root+"/intermediate_current.json")){
        // Legacy checkpoints can already contain this decision's streak. A
        // matching saved estimator must not count again after an upgrade.
        const auto old=CompactJson(Read(e.root+"/intermediate_current.json"));
        reused=old.find("\"snapshot_folder\":"+Quoted(snapshot))!=std::string::npos
            &&old.find("\"numerical_round_passed\":true")!=std::string::npos;
    }
    if(!passed)result=0;else if(!reused)result=std::min(e.stablePasses,result+1);
    state.streak=result;std::ostringstream payload;payload<<e.signature<<' '<<key<<' '<<result<<' '<<passed<<'\n';Atomic(path,Hash(payload.str())+"\n"+payload.str());
    e.Log("production_validation_decision","{\"round\":"+Number(state.round)+",\"reused_saved_estimator\":"+(reused?"true":"false")+",\"passed\":"+(passed?"true":"false")+",\"distinct_passes\":"+Number(result)+",\"required_passes\":"+Number(e.stablePasses)+"}");
}
void SaveAdaptive(Engine& e,const AdaptiveState& state,const std::string& phase="active"){
    std::ostringstream payload;payload<<std::setprecision(17)<<e.signature<<' '<<phase<<' '<<state.round<<' '<<state.streak<<' '<<e.phi<<' '<<e.controls<<' '<<e.sharedMeans<<' '<<e.mirror<<' '<<e.mirrorError<<' '<<e.kernel.block<<' '<<e.kernel.pipeline<<' '<<e.kernel.warp<<'\n';
    payload<<state.auditSerial<<'\n';payload<<e.depthSharing<<' '<<e.depthVerified<<' '<<e.local<<'\n';payload<<e.weightTheta.size()<<'\n';for(size_t r=0;r<e.weightTheta.size();++r){payload<<e.weightTheta[r];
        const auto b=e.Polar(e.weightTheta[r]);for(int j=0;j<16;++j)for(int c=0;c<2;++c)payload<<' '<<(j?b[c][j]:(e.weights.empty()?1:e.weights[r][c]));payload<<'\n';}
    payload<<state.knots.size();for(size_t k:state.knots)payload<<' '<<k;payload<<'\n';
    const auto evaluated=AdaptiveGuards(state.knots,e.theta.size());std::set<size_t> saved(evaluated.begin(),evaluated.end());
    for(size_t i=0;i<e.theta.size();++i)if(state.valid[i])saved.insert(i);
    payload<<state.counts.size()<<' '<<saved.size()<<'\n';for(size_t i:saved){payload<<i<<' '<<state.valid[i]<<' '<<state.coarse[i];for(const auto& level:state.counts)payload<<' '<<level[i];for(double v:state.previous[i])payload<<' '<<v;payload<<'\n';}
    Atomic(e.root+"/adaptive_checkpoint.tsv","3 "+Hash(payload.str())+"\n"+payload.str());
    std::ostringstream history;history<<"{\"phase\":"<<Quoted(phase)<<",\"round\":"<<state.round<<",\"pass_streak\":"<<state.streak<<",\"phi\":"<<e.phi<<",\"theta_support_indices\":[";
    for(size_t i=0;i<state.knots.size();++i)history<<(i?",":"")<<state.knots[i];history<<"],\"counts_by_level\":[";
    for(size_t level=0;level<state.counts.size();++level){if(level)history<<',';std::map<int,std::vector<size_t>> groups;for(size_t id:saved)groups[state.counts[level][id]].push_back(id);
        history<<'[';bool first=true;for(const auto& group:groups){history<<(first?"":",")<<"{\"count_per_seed\":"<<group.first<<",\"theta_indices\":[";
            for(size_t i=0;i<group.second.size();++i)history<<(i?",":"")<<group.second[i];history<<"]}";first=false;}history<<']';}
    history<<"]}";e.Log("allocation_checkpoint",history.str());
}
bool RestoreAdaptive(Engine& e,AdaptiveState& state){
    const std::string path=e.root+"/adaptive_checkpoint.tsv";if(!Exists(path))return false;
    const std::string text=Read(path);const size_t newline=text.find('\n');std::istringstream header(text.substr(0,newline));int version;std::string hash;Require(bool(header>>version>>hash)&&(version==2||version==3),"unsupported checkpoint schema");
    const std::string payload=text.substr(newline+1);Require(Hash(payload)==hash,"controller checkpoint changed");std::istringstream file(payload);std::string signature,phase;
    Require(bool(file>>signature>>phase>>state.round>>state.streak>>e.phi>>e.controls>>e.sharedMeans>>e.mirror>>e.mirrorError>>e.kernel.block>>e.kernel.pipeline>>e.kernel.warp)&&signature==e.signature,"checkpoint inputs changed");
    if(version>=3)file>>state.auditSerial;else state.auditSerial=unsigned(state.round)*8u+8u;
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
        ValidateIntegralSources(e.root);
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
struct InterpolationResult {Samples output;Matrix matrixErrors,ratioErrors;std::vector<double> errors;std::vector<size_t> bad;double maximum=0;};
#include "FullAutoCorrectionProbe.inc"
#include "FullAutoWorkingPhi.inc"
InterpolationResult InterpolateBudget(const Engine& e,const std::vector<size_t>& knots,const std::vector<size_t>& evaluated,const Samples& samples){
    InterpolationResult result;result.output.resize(8);result.errors.resize(evaluated.size());result.matrixErrors.resize(evaluated.size());result.ratioErrors.resize(evaluated.size());Samples predicted(8);const auto x=Angles(e.theta,evaluated),sx=Angles(e.theta,knots);const auto support=Support(samples,evaluated,knots);
    for(int s=0;s<8;++s){predicted[s]=LocalCubicInterpolate(sx,support[s],x);result.output[s]=LocalCubicInterpolate(sx,support[s],e.theta);}
    const auto delta=Subtract(predicted,samples);const auto mean=Mean(delta),variance=Variance(delta),reference=Mean(samples),prediction=Mean(predicted);
    for(size_t r=0;r<evaluated.size();++r){Require(reference[r][0]>0,"nonpositive interpolation reference");for(int j=0;j<16;++j){result.matrixErrors[r][j]=(std::fabs(mean[r][j])+t95*std::sqrt(variance[r][j]/8))/reference[r][0];result.errors[r]=std::max(result.errors[r],result.matrixErrors[r][j]);}
        if(e.ratios&&prediction[r][0]>0)for(int j=1;j<16;++j){const double a=prediction[r][j]/prediction[r][0],b=reference[r][j]/reference[r][0];double variance=0;
            for(int s=0;s<8;++s){const double influence=(predicted[s][r][j]-a*predicted[s][r][0])/prediction[r][0]-(samples[s][r][j]-b*samples[s][r][0])/reference[r][0];variance+=influence*influence/7;}
            result.ratioErrors[r][j]=std::fabs(a-b)+t95*std::sqrt(variance/8);result.errors[r]=std::max(result.errors[r],result.ratioErrors[r][j]);}
        result.maximum=std::max(result.maximum,result.errors[r]);if(result.errors[r]>e.interpolation)result.bad.push_back(evaluated[r]);}
    return result;
}
InterpolationResult CommonPrefixInterpolation(Engine& e,const AdaptiveState& state,const std::vector<int>& depths,const std::vector<size_t>& evaluated,const unsigned* seeds){
    InterpolationResult result;result.errors.assign(evaluated.size(),0);result.matrixErrors.resize(evaluated.size());result.ratioErrors.resize(evaluated.size());const size_t width=std::min(size_t(4),state.knots.size());
    const auto groups=FullAutoAuditPlan::Stencils(state.knots,state.counts,evaluated);
    std::ostringstream sources,budgets;budgets<<"group,theta_deg";for(size_t l=0;l<depths.size();++l)budgets<<",N_level"<<l;budgets<<'\n';size_t groupNumber=0;
    for(const auto& group:groups){const std::vector<size_t> ids(group.second.rows.begin(),group.second.rows.end());Budgets ns(depths.size());for(size_t l=0;l<depths.size();++l)ns[l].assign(e.theta.size(),group.first[l]);
        auto stage=EvaluateBudget(e,depths,ns,ids,seeds);const Samples data=Combine(stage.levels);const Matrix reference=Mean(data);const auto x=Angles(e.theta,ids);
        for(size_t interval:group.second.intervals){const size_t left=std::min(interval?interval-1:0,state.knots.size()-width),right=std::min(interval,state.knots.size()-width);
            std::vector<double> knots;Samples supports(8);for(size_t k=left;k<right+width;++k){const size_t row=std::lower_bound(ids.begin(),ids.end(),state.knots[k])-ids.begin();knots.push_back(e.theta[state.knots[k]]);for(int s=0;s<8;++s)supports[s].push_back(data[s][row]);}
            // Use the globally selected endpoint derivatives: their local
            // four-node windows are entirely contained in this stencil.
            std::vector<double> query;std::vector<size_t> rows,positions;for(size_t r=0;r<evaluated.size();++r)if(evaluated[r]>=state.knots[interval]&&evaluated[r]<=state.knots[interval+1]){query.push_back(e.theta[evaluated[r]]);rows.push_back(std::lower_bound(ids.begin(),ids.end(),evaluated[r])-ids.begin());positions.push_back(r);}
            Samples predicted(8);for(int s=0;s<8;++s)predicted[s]=LocalCubicInterpolate(knots,supports[s],query);const Matrix pm=Mean(predicted);
            for(size_t r=0;r<rows.size();++r){double error=0;const auto& ref=reference[rows[r]];Require(ref[0]>0,"nonpositive common-prefix interpolation reference");
                for(int j=0;j<16;++j){double mean=0,var=0;for(int s=0;s<8;++s)mean+=(predicted[s][r][j]-data[s][rows[r]][j])/8;for(int s=0;s<8;++s)var+=std::pow(predicted[s][r][j]-data[s][rows[r]][j]-mean,2)/7;const double raw=(std::fabs(mean)+t95*std::sqrt(var/8))/ref[0];result.matrixErrors[positions[r]][j]=std::max(result.matrixErrors[positions[r]][j],raw);error=std::max(error,raw);
                    if(e.ratios&&j&&pm[r][0]>0){const double a=pm[r][j]/pm[r][0],b=ref[j]/ref[0];double v=0;for(int s=0;s<8;++s){const double influence=(predicted[s][r][j]-a*predicted[s][r][0])/pm[r][0]-(data[s][rows[r]][j]-b*data[s][rows[r]][0])/ref[0];v+=influence*influence/7;}const double ratio=std::fabs(a-b)+t95*std::sqrt(v/8);result.ratioErrors[positions[r]][j]=std::max(result.ratioErrors[positions[r]][j],ratio);error=std::max(error,ratio);}}
                result.errors[positions[r]]=std::max(result.errors[positions[r]],error);}
        }
        std::istringstream records(stage.sources);std::string line;while(std::getline(records,line))sources<<groupNumber<<' '<<line<<'\n';
        for(size_t id:ids){budgets<<groupNumber<<','<<Number(e.theta[id]);for(int n:group.first)budgets<<','<<n;budgets<<'\n';}++groupNumber;
    }
    for(size_t r=0;r<evaluated.size();++r){result.maximum=std::max(result.maximum,result.errors[r]);if(result.errors[r]>e.interpolation)result.bad.push_back(evaluated[r]);}
    std::ostringstream seedNames;for(int s=0;s<8;++s)seedNames<<seeds[s]<<' ';const std::string label=seeds[0]==e.seeds[0]?"production":"validation";
    Atomic(e.root+"/interpolation_"+label+"_sources.tsv",sources.str());Atomic(e.root+"/interpolation_"+label+"_budgets.csv",budgets.str());return result;
}
#include "FullAutoBudget.inc"
void SaveReplicaMatrices(Engine& e,const std::string& label,const Samples& data,const std::vector<double>& x,const unsigned* seeds){Directory(e.root+"/replicas");for(int s=0;s<8;++s)WriteMatrix(e.root+"/replicas/"+label+"_seed"+Number(seeds[s])+".dat",x,data[s]);}
void RefineFailedAudit(Engine& e,AdaptiveState& state,const std::vector<int>& depths,
        const AdaptiveData& productionData,const std::vector<size_t>& productionIds,const AdaptiveData& audit,
        const std::vector<size_t>& ids,const Budgets& observedCounts,const Statistics& stats,
        const Statistics& output,const InterpolationResult& interp){
    const Budgets before=state.counts;Budgets proposed=before;AdditionalWork work(e,e.seeds.data());
    std::set<size_t> knots(state.knots.begin(),state.knots.end());knots.insert(interp.bad.begin(),interp.bad.end());
    std::ostringstream log;log<<RefinementHeader();
    for(size_t r=0;r<ids.size();++r){const auto residual=ResidualAt(e,interp,ids,ids[r]);const auto metric=BudgetMetric(e,stats,r,residual);if(MetricWithinBudget(e,metric,residual))continue;
        RefineError(e,work,depths,before,proposed,ids,audit,{{ids[r],1}},metric,stats.mean[r],e.theta[ids[r]],"audit_sampled",log,{},observedCounts);knots.insert(ids[r]);}
    for(size_t q=0;q<e.theta.size();++q){const auto residual=ResidualAt(e,interp,ids,q);const auto metric=BudgetMetric(e,output,q,residual);if(MetricWithinBudget(e,metric,residual))continue;
        RefineError(e,work,depths,before,proposed,ids,audit,QueryWeights(e,state.knots,e.theta[q]),metric,output.mean[q],e.theta[q],"audit_interpolated",log,{},observedCounts);}
    const auto production=Mean(Combine(productionData.levels));std::vector<Matrix> pm,am;
    for(const auto& level:productionData.levels)pm.push_back(Mean(level));for(const auto& level:audit.levels)am.push_back(Mean(level));
    for(size_t r=0;r<productionIds.size();++r){const size_t id=productionIds[r],row=std::lower_bound(ids.begin(),ids.end(),id)-ids.begin();ErrorMetric metric;
        for(int j=0;j<16;++j){const double delta=std::fabs(production[r][j]-stats.mean[row][j])/production[r][0];if(delta>metric.error)metric=ErrorMetric(j,false,delta);}
        if(metric.error<=e.eps)continue;
        std::vector<double> discrepancy(depths.size());for(size_t l=0;l<depths.size();++l)discrepancy[l]=pm[l][r][metric.element]-am[l][row][metric.element];
        RefineError(e,work,depths,before,proposed,ids,audit,{{id,1}},metric,stats.mean[row],e.theta[id],"audit_mean_difference",log,discrepancy,observedCounts);knots.insert(id);
    }
    state.knots.assign(knots.begin(),knots.end());
    // New midpoint guards inherit the surrounding requested budget. Explicit
    // refinements at failed audit nodes must not be overwritten by inheritance.
    const auto nextIds=AdaptiveGuards(state.knots,e.theta.size());
    for(size_t id:nextIds)if(!std::binary_search(ids.begin(),ids.end(),id)){
        const auto upper=std::upper_bound(ids.begin(),ids.end(),id);const size_t right=upper==ids.end()?ids.back():*upper,left=upper==ids.begin()?ids.front():*(upper-1);
        for(auto& level:proposed)level[id]=std::max(level[id],std::max(level[left],level[right]));}
    if(!e.local)for(auto& level:proposed){int maximum=0;for(size_t id:nextIds)maximum=std::max(maximum,level[id]);for(size_t id:nextIds)level[id]=maximum;}
    GateLocalCost(e,work,proposed,depths,nextIds);
    std::ostringstream budgets;budgets<<"theta_deg,level,previous_count,next_count\n";
    for(size_t id:nextIds)for(size_t l=0;l<depths.size();++l)budgets<<Number(e.theta[id])<<','<<l<<','<<before[l][id]<<','<<proposed[l][id]<<'\n';
    const auto filename="audit_refinement_budget_"+Number(state.auditSerial)+".csv";Atomic(e.root+"/"+filename,budgets.str());
    RefinementLog(e,state.round,"audit_error_growth",log);e.Log("audit_refinement","{\"counts\":"+Quoted(filename)+",\"next_audit_requires_fresh_seeds\":true}");state.counts=proposed;
}
std::array<unsigned,8> ReserveAuditSeeds(Engine& e,AdaptiveState& state){
    Require(state.auditSerial<976562u,"independent audit seed namespace exhausted");
    const unsigned serial=state.auditSerial++;std::array<unsigned,8> seeds;
    for(int s=0;s<8;++s){seeds[s]=1000000000u+serial*512u+unsigned(s)*17u;Require(std::find(e.seeds.begin(),e.seeds.end(),seeds[s])==e.seeds.end(),"final/production seed overlap");}
    // Reserve before reading any audit result. A crash can skip an unused
    // holdout, but cannot reuse one whose data influenced an allocation.
    SaveAdaptive(e,state);return seeds;
}
std::string AuditInputFingerprint(const Engine& e,const AdaptiveState& state,const AdaptiveData& production,const std::vector<size_t>& ids){
    std::ostringstream out;out<<e.signature<<'\n'<<e.ExecutionHash()<<'\n'<<state.round<<' '<<e.local<<' '<<e.initial<<' '<<e.minimum<<' '<<e.orientationBlock
        <<' '<<Number(e.statistical)<<' '<<Number(e.interpolation)<<' '<<BudgetPolicy()<<'\n'<<production.sources;
    for(size_t k:state.knots)out<<k<<' ';out<<'\n';
    for(size_t id:ids){out<<id;for(const auto& level:state.counts)out<<' '<<level[id];out<<'\n';}
    return Hash(out.str());
}
struct AuditReservation {std::array<unsigned,8> seeds;unsigned serial=0;int legacyFloor=0;bool resumed=false;};
AuditReservation ReserveOrResumeAudit(Engine& e,AdaptiveState& state,const AdaptiveData& production,const std::vector<size_t>& ids,int& attempt,double& fraction,const std::string& allocation=""){
    const std::string fingerprint=Hash(AuditInputFingerprint(e,state,production,ids)+allocation),path=e.root+"/independent_audit_reservation.tsv";
    AuditReservation result;bool restored=false;
    if(Exists(path)){
        const auto text=Read(path);const auto end=text.find('\n');Require(end!=std::string::npos,"invalid audit reservation");
        const auto payload=text.substr(end+1);Require(Hash(payload)==text.substr(0,end),"audit reservation changed");
        std::istringstream input(payload);std::string signature,identity,extra;unsigned serial;int savedAttempt,floor;double savedFraction;
        Require(bool(input>>signature>>identity>>serial>>savedAttempt>>savedFraction>>floor)&&!(input>>extra),"invalid audit reservation fields");
        Require(signature==e.signature&&savedAttempt>=0&&savedAttempt<2&&savedFraction>0&&savedFraction<=1&&floor>=0,"invalid frozen audit reservation");
        if(identity==fingerprint&&serial+1==state.auditSerial&&!Exists(e.root+"/independent_final_validation_serial"+Number(serial)+".json")){
            result.serial=serial;result.legacyFloor=floor;attempt=savedAttempt;fraction=savedFraction;restored=true;
        }
    }
    if(!restored&&attempt==0&&e.pendingLegacyAuditSerial>=0&&state.round==e.pendingLegacyAuditRound
        &&unsigned(e.pendingLegacyAuditSerial)+1==state.auditSerial&&Hash(production.sources)==e.pendingLegacyAuditSources
        &&!Exists(e.root+"/independent_final_validation_serial"+Number(e.pendingLegacyAuditSerial)+".json")){
        result.serial=unsigned(e.pendingLegacyAuditSerial);result.legacyFloor=e.orientationBlock;restored=true;
        // A deterministic one-block floor retains a legacy correction block
        // without choosing its budget from any holdout values or timings.
        e.Log("legacy_pending_audit_import","{\"audit_serial\":"+Number(result.serial)+",\"native_block_floor\":"+Number(result.legacyFloor)+",\"unexposed_holdout_only\":true}");
    }
    e.pendingLegacyAuditSerial=-1;
    if(restored){for(int s=0;s<8;++s)result.seeds[s]=1000000000u+result.serial*512u+unsigned(s)*17u;
        result.resumed=true;e.Log("pending_audit_resumed","{\"audit_serial\":"+Number(result.serial)+",\"attempt\":"+Number(attempt)+",\"allocation_frozen_from_production\":true}");
    }else{result.seeds=ReserveAuditSeeds(e,state);result.serial=state.auditSerial-1;}
    std::ostringstream payload;payload<<e.signature<<' '<<fingerprint<<' '<<result.serial<<' '<<attempt<<' '<<Number(fraction)<<' '<<result.legacyFloor<<'\n';
    Atomic(path,Hash(payload.str())+"\n"+payload.str());return result;
}
void PublishAuditLowerViews(Engine& e,const std::array<unsigned,8>& seeds){
    // Completed legacy paired blocks contain both physical matrices. Make the
    // lower one visible before planning base work, without another ray trace.
    std::vector<Job> paired;for(const auto& entry:e.cache){const auto& job=entry.second;
        if(job.pairedDepth&&!job.assembled&&std::find(seeds.begin(),seeds.end(),job.seed)!=seeds.end()
            &&e.SameSampling(job,job.count,job.seed,job.depth,job.offset))paired.push_back(job);
    }
    for(const auto& job:paired)e.LowerView(job,job.pairedDepth);
    if(!paired.empty())e.Log("audit_paired_cache_reuse","{\"paired_blocks\":"+Number(paired.size())+",\"lower_matrices_reused_without_tracing\":true}");
}
std::vector<Engine::WorkRequest> AuditRequests(const Engine& e,const std::vector<int>& depths,const Budgets& counts,const std::vector<size_t>& ids,const std::vector<size_t>& knots){
    auto requests=BudgetRequests(e,depths,counts,ids);
    if(VaryingBudgets(counts,ids))for(const auto& group:FullAutoAuditPlan::Stencils(knots,counts,ids)){
        const std::vector<size_t> nodes(group.second.rows.begin(),group.second.rows.end());Budgets prefix;
        for(int n:group.first)prefix.push_back(std::vector<int>(e.theta.size(),n));
        auto extra=BudgetRequests(e,depths,prefix,nodes);requests.insert(requests.end(),extra.begin(),extra.end());
    }return requests;
}
bool UseLocalAudit(Engine& e,const AdaptiveState& state,const std::vector<int>& depths,const std::vector<size_t>& ids,double fraction){
    if(!e.local||e.localPolicy!="auto"||!e.angularCache)return e.local;
    // Freeze this scheduling choice using production profiles only. Neither
    // partial holdout values, their timing, nor cache coverage can select it.
    const std::string saved=e.workQueuePolicy;if(saved=="auto")e.workQueuePolicy="on";
    WorkForecast local,uniform;double extraAliases=0;
    auto aliases=[&](const std::vector<Engine::WorkRequest>& requests){double n=0;for(const auto& r:requests)n+=8.*(r.pairedDepth?2:1)*(1+(e.incremental?(r.count+e.orientationBlock-1)/e.orientationBlock:1));return n;};
    try{AdditionalWork work(e,e.seeds.data(),true,true);
        const auto a=FullAutoAuditPlan::Scale(state.counts,ids,fraction,e.initial,e.minimum,true);
        const auto b=FullAutoAuditPlan::Scale(state.counts,ids,fraction,e.initial,e.minimum,false);
        const auto ar=AuditRequests(e,depths,a,ids,state.knots),br=AuditRequests(e,depths,b,ids,state.knots);
        local=work.Forecast(ar);uniform=work.Forecast(br);extraAliases=std::max(0.,aliases(ar)-aliases(br));
    }catch(...){e.workQueuePolicy=saved;throw;}e.workQueuePolicy=saved;
    // Account for serial assembly, measured on production-derived matrices.
    // The alias count is a conservative upper bound; shared views can reduce
    // it. No accuracy budget depends on this execution-cost heuristic.
    double assembly=0;size_t observations=0;
    for(const auto& item:e.cache){const auto& j=item.second;if(j.assembled&&j.assemblySeconds>0
        &&std::find(e.seeds.begin(),e.seeds.end(),j.seed)!=e.seeds.end()){assembly+=j.assemblySeconds;++observations;}}
    if(observations)assembly/=observations;
    const double extraAssembly=extraAliases*assembly,concurrency=std::max(size_t(1),std::min(size_t(8),e.WorkerGpus().size()));
    const bool selected=!(local.profiled&&uniform.profiled&&local.seconds/concurrency+extraAssembly>=uniform.seconds/concurrency*.90);
    e.Record(e.root+"/independent_audit_cost_gate.json","{\"policy\":\"production profiles only; ignore holdout coverage\",\"profiled\":"+std::string(local.profiled&&uniform.profiled?"true":"false")
        +",\"predicted_local_worker_seconds\":"+Number(local.seconds)+",\"predicted_uniform_worker_seconds\":"+Number(uniform.seconds)
        +",\"production_assembly_seconds_per_alias\":"+Number(assembly)+",\"extra_alias_upper_bound\":"+Number(extraAliases)+",\"predicted_extra_assembly_seconds\":"+Number(extraAssembly)+",\"minimum_saving_fraction\":0.1,\"selected_local\":"+(selected?"true":"false")+",\"scope\":\"scheduling forecast, not measured speed or accuracy criterion\"}");
    return selected;
}
#include "FullAutoAngularAudit.inc"
bool IndependentAudit(Engine& e,AdaptiveState& state,const std::vector<int>& depths,const AdaptiveData& productionData,const std::vector<size_t>& productionIds,const InterpolationResult* productionInterpolation=nullptr){
    if(!e.finalAudit){e.Record(e.root+"/independent_final_validation.json","{\"enabled\":false,\"passed\":false}");return true;}
    auto ids=AdaptiveGuards(state.knots,e.theta.size(),true);
    if(e.theta.size()<=size_t(std::min(4097,e.maxTheta))){ids.resize(e.theta.size());for(size_t i=0;i<ids.size();++i)ids[i]=i;}
    FullAutoAuditPlan::Extend(state.counts,productionIds,ids);
    // Freeze local budgets from production observations, before any holdout
    // is read. New guard nodes inherit only neighbouring production budgets.
    const Samples productionSamples=Combine(productionData.levels);const auto productionStats=Estimate(productionSamples);
    double relative=StatisticalError(e,productionStats)/e.statistical;
    if(productionInterpolation){relative=0;for(size_t r=0;r<productionIds.size();++r){
        const auto residual=ResidualAt(e,*productionInterpolation,productionIds,productionIds[r]);
        for(int j=0;j<16;++j){relative=std::max(relative,productionStats.interval[r][j]/std::max(e.statistical,NumericalBudget(e)-residual.matrix[j]));
            if(e.ratios&&j)relative=std::max(relative,productionStats.ratioInterval[r][j]/std::max(e.statistical,NumericalBudget(e)-residual.ratio[j]));}}}
    const double projected=std::pow(relative/.5,2);
    double fraction=std::min(1.,std::max(1./16.,projected));
    if(e.angularAuditPolicy=="on")fraction=1; // explicit forced-control validation/benchmark
    const bool auditLocal=UseLocalAudit(e,state,depths,ids,fraction);
    const auto angularPlan=PlanAngularAudit(e,state,depths,productionData,productionIds,ids,fraction,auditLocal);
    const std::string allocationIdentity=(auditLocal?"local":"uniform")+(angularPlan.enabled?std::string("\nangular-residual-v1 ")+angularPlan.identity:std::string());
    const auto production=Mean(productionSamples);bool pass=false;
    for(int attempt=0;attempt<2;++attempt){
        const auto reservation=ReserveOrResumeAudit(e,state,productionData,ids,attempt,fraction,allocationIdentity);const auto& seeds=reservation.seeds;
        PublishAuditLowerViews(e,seeds);
        Budgets validationCounts=FullAutoAuditPlan::Scale(state.counts,ids,fraction,e.initial,e.minimum,auditLocal);
        if(reservation.legacyFloor)for(size_t l=0;l<validationCounts.size();++l)for(size_t id:ids)
            validationCounts[l][id]=std::max(validationCounts[l][id],std::min(reservation.legacyFloor,state.counts[l][id]));
        std::vector<int> chosen;for(const auto& level:validationCounts){int n=0;for(size_t id:ids)n=std::max(n,level[id]);chosen.push_back(n);}
        const bool angular=angularPlan.enabled&&attempt==0;
        const bool common=auditLocal&&VaryingBudgets(validationCounts,ids);
        std::ostringstream plan;plan<<"{\"schema\":1,\"audit_serial\":"<<state.auditSerial-1<<",\"allocation_from_production\":true,\"local_orientation_budgets\":"<<(auditLocal?"true":"false")
            <<",\"resumed_pending_audit\":"<<(reservation.resumed?"true":"false")<<",\"legacy_block_floor\":"<<reservation.legacyFloor<<",\"common_prefix_interpolation\":"<<(common?"true":"false")<<",\"fraction\":"<<Number(fraction)<<",\"checked_theta\":"<<ids.size()<<",\"levels\":[";
        const auto stencils=FullAutoAuditPlan::Stencils(state.knots,validationCounts,ids);
        for(size_t l=0;l<depths.size();++l){uint64_t work=0,stencilWork=0,nativeWork=0;int lo=chosen[l];std::vector<int> envelope=validationCounts[l];
            for(size_t id:ids){work+=validationCounts[l][id];lo=std::min(lo,validationCounts[l][id]);}
            for(const auto& group:stencils){stencilWork+=uint64_t(group.first[l])*group.second.rows.size();
                if(common)for(size_t id:group.second.rows)envelope[id]=std::max(envelope[id],group.first[l]);}
            if(angular){std::fill(envelope.begin(),envelope.end(),0);
                for(size_t id:state.knots)envelope[id]=angularPlan.anchors[l][id];
                for(const auto& group:FullAutoAuditPlan::Stencils(state.knots,angularPlan.residuals,ids)){
                    bool interior=false;for(size_t id:group.second.rows)interior=interior||!std::binary_search(state.knots.begin(),state.knots.end(),id);
                    if(interior)for(size_t id:group.second.rows)envelope[id]=std::max(envelope[id],group.first[l]);}
            }
            for(size_t id:ids)nativeWork+=envelope[id];
            plan<<(l?",":"")<<"{\"depth\":"<<depths[l]<<",\"minimum_per_seed\":"<<lo<<",\"maximum_per_seed\":"<<chosen[l]
                <<",\"orientation_angle_work_per_seed\":"<<work<<",\"stencil_orientation_angle_work_per_seed\":"<<stencilWork
                <<",\"native_envelope_orientation_angle_work_per_seed\":"<<nativeWork<<",\"uniform_orientation_angle_work_per_seed\":"<<uint64_t(chosen[l])*ids.size()<<"}";}
        plan<<"],\"estimator\":"<<Quoted(angular?"angular_residual_v1":"direct")<<",\"angular_allocation_fnv1a64\":"<<Quoted(angular?angularPlan.identity:"")<<"}";e.Record(e.root+"/independent_audit_plan.json",plan.str());
        std::ostringstream frozen;frozen<<"theta_deg";for(size_t l=0;l<depths.size();++l)frozen<<",N_level"<<l;frozen<<'\n';
        for(size_t id:ids){frozen<<Number(e.theta[id]);for(const auto& level:validationCounts)frozen<<','<<level[id];frozen<<'\n';}
        Atomic(e.root+"/independent_orientation_budget.csv",frozen.str());
        std::cout<<"fullauto independent audit: attempt="<<attempt<<" theta="<<ids.size()<<" allocation="<<(auditLocal?"local":"uniform")<<" maximum_counts=";for(int n:chosen)std::cout<<n<<' ';std::cout<<std::endl;
        // Union angular groups within each orientation window even on one GPU.
        // Explicit --fullauto-work-queue off is preserved. This avoids tracing
        // the same orientations separately for every local count group.
        const std::string savedQueue=e.workQueuePolicy;
        const bool queueAuto=auditLocal&&e.workQueuePolicy=="auto"&&e.angularCache;
        if(queueAuto)e.workQueuePolicy="on";
        AdaptiveData data;InterpolationResult interp;Samples combined;Statistics stats;
        try{
            // Plan point estimates and common-prefix interpolation together.
            // One traced orientation window supplies all of their angular
            // subsets, avoiding a second trace for each stencil group.
            if(angular){data=EvaluateAngularAudit(e,depths,angularPlan,ids,state.knots,seeds.data());}
            else if(common&&e.QueueEnabled()){
                e.Prefetch(AuditRequests(e,depths,validationCounts,ids,state.knots),seeds.data());
            }
            if(!angular)data=EvaluateBudget(e,depths,validationCounts,ids,seeds.data());combined=Combine(data.levels);
            if(angular){bool positive=true;for(const auto& row:Mean(combined))positive=positive&&row[0]>0;
                if(!positive){e.workQueuePolicy=savedQueue;const std::string report="{\"enabled\":true,\"passed\":false,\"estimator\":\"angular_residual_v1\",\"reason\":\"nonpositive independent M11; fresh direct fallback required\"}";
                    e.Record(e.root+"/independent_final_validation.json",report);e.Record(e.root+"/independent_final_validation_serial"+Number(state.auditSerial-1)+".json",report);
                    e.Log("independent_audit_fallback",report);fraction=1;continue;}}
            stats=Estimate(combined);
            interp=InterpolateBudget(e,state.knots,ids,combined);
            if(common&&!angular){AdaptiveState checked=state;checked.counts=validationCounts;
                const auto guards=CommonPrefixInterpolation(e,checked,depths,ids,seeds.data());
                interp.errors=guards.errors;interp.matrixErrors=guards.matrixErrors;interp.ratioErrors=guards.ratioErrors;interp.maximum=guards.maximum;interp.bad=guards.bad;
            }
        }catch(...){e.workQueuePolicy=savedQueue;throw;}
        e.workQueuePolicy=savedQueue;
        e.Record(e.root+"/interpolation_validation_mode.json","{\"angular_residual_estimator\":"+std::string(angular?"true":"false")+",\"common_prefix\":"+std::string(common?"true":"false")+",\"uniform_level_budgets\":"+(auditLocal?"false":"true")+",\"work_queue_auto_for_local_audit\":"+(queueAuto?"true":"false")+"}");
        if(angular){bool positive=true;for(const auto& row:Mean(interp.output))positive=positive&&row[0]>0;
            if(!positive){const std::string report="{\"enabled\":true,\"passed\":false,\"estimator\":\"angular_residual_v1\",\"reason\":\"nonpositive interpolated independent M11; fresh direct fallback required\"}";
                e.Record(e.root+"/independent_final_validation.json",report);e.Record(e.root+"/independent_final_validation_serial"+Number(state.auditSerial-1)+".json",report);
                e.Log("independent_audit_fallback",report);fraction=1;continue;}}
        const auto outputStats=Estimate(interp.output);
        Matrix matched;for(size_t id:productionIds){size_t r=std::lower_bound(ids.begin(),ids.end(),id)-ids.begin();matched.push_back(stats.mean[r]);}
        const double change=Difference(production,matched,production);
        const auto budget=CheckPointwiseBudget(e,stats,outputStats,interp,ids,"independent_serial"+Number(state.auditSerial-1));
        e.independentNeedsSampling=!budget.directPass||!budget.outputPass||change>e.eps;
        pass=!e.independentNeedsSampling&&interp.maximum<=e.interpolation;
        std::ostringstream report;report<<"{\"enabled\":true,\"passed\":"<<(pass?"true":"false")<<",\"fresh_seeds_per_attempt\":true,\"seeds\":[";
        for(int s=0;s<8;++s)report<<(s?",":"")<<seeds[s];
        report<<"],\"counts_per_seed\":[";for(size_t l=0;l<chosen.size();++l)report<<(l?",":"")<<chosen[l];
        report<<"],\"error_budget_policy\":"<<Quoted(BudgetPolicy())<<",\"max_combined_direct95\":"<<Number(budget.direct)<<",\"max_combined_output95\":"<<Number(budget.output)<<",\"fixed_error_reserve\":"<<Number(FixedBudgetReserve(e))<<",\"estimator\":"<<Quoted(angular?"angular_residual_v1":"direct")<<",\"angular_allocation_fnv1a64\":"<<Quoted(angular?angularPlan.identity:"")<<",\"local_orientation_budgets\":"<<(auditLocal?"true":"false")<<",\"allocation_from_production\":true,\"attempt\":"<<attempt<<",\"max_pointwise95_scaled_M11\":"<<Number(std::max(stats.maxInterval,outputStats.maxInterval))<<",\"max_ratio95\":"<<Number(std::max(stats.maxRatioInterval,outputStats.maxRatioInterval))<<",\"max_interpolation_residual95\":"<<Number(interp.maximum)<<",\"difference_to_production_scaled_M11\":"<<Number(change)<<",\"checked_theta\":"<<ids.size()<<",\"scope\":\"fresh-seed frozen-estimator pointwise validation; all requested nodes up to 4097, otherwise quarter interval guards; no simultaneous-confidence claim\"}";
        e.Record(e.root+"/independent_final_validation.json",report.str());e.Record(e.root+"/independent_final_validation_attempt"+Number(attempt)+".json",report.str());
        e.Record(e.root+"/independent_final_validation_serial"+Number(state.auditSerial-1)+".json",report.str());
        Atomic(e.root+"/independent_final_sources.tsv",data.sources);SaveReplicaMatrices(e,"validation_evaluated",combined,Angles(e.theta,ids),seeds.data());
        WriteGrid(e.root+"/independent_theta.csv",Angles(e.theta,ids));
        // A low-cost audit that lacks statistics gets a fresh-seed full
        // allocation fallback. Interpolation-only failure refines theta;
        // it must not blindly double all production orientation counts.
        if(!pass&&attempt==0&&(angular||(e.independentNeedsSampling&&fraction<1))){
            e.Log("independent_audit_fallback","{\"from\":"+Quoted(angular?"angular_residual_v1":"reduced_direct")+",\"fresh_seeds_required\":true,\"next\":\"full direct audit\"}");fraction=1;continue;}
        if(!pass)RefineFailedAudit(e,state,depths,productionData,productionIds,data,ids,validationCounts,stats,outputStats,interp);
        break;
    }
    e.independentPassed=pass;return pass;
}

bool PhysicsAudit(Engine& e,const AdaptiveState& state,const std::vector<size_t>& evaluated,const Matrix& reference){
    if(!e.physicsAudit&&!e.autoDepth){e.Record(e.root+"/physics_validation.json","{\"enabled\":false,\"passed\":false,\"requested_physics_retained\":true}");return true;}
    std::vector<size_t> selected;const size_t size=std::min(size_t(e.autoDepth?4097:65),evaluated.size());for(size_t r=0;r<size;++r)selected.push_back(evaluated[r*(evaluated.size()-1)/(size-1)]);
    const auto x=Angles(e.theta,selected);Matrix scale;for(size_t id:selected)scale.push_back(reference[std::lower_bound(evaluated.begin(),evaluated.end(),id)-evaluated.begin()]);
    std::array<unsigned,8> seeds;for(int s=0;s<8;++s){seeds[s]=1500000000u+unsigned(state.round)*1024u+unsigned(s)*17u;Require(std::find(e.seeds.begin(),e.seeds.end(),seeds[s])==e.seeds.end(),"physics/production seed overlap");}
    int count=e.minimum;for(size_t l=1;l<state.counts.size();++l)for(size_t id:selected)count=std::max(count,state.counts[l][id]);
    const auto original=e.base;const auto savedCutoff=e.selectedCutoff;bool pass=false;double depthError=0,cutoffError=0;
    for(int attempt=0;attempt<5;++attempt){
        std::vector<Job> target,deeper;
        if(e.depthSharing&&e.depthVerified){const auto pair=e.SharedStage(count,seeds.data(),8,e.target+2,e.target,x);deeper=pair.first;target=pair.second;}
        else{target=e.PrefixStage(count,seeds.data(),8,e.target,x);deeper=e.PrefixStage(count,seeds.data(),8,e.target+2,x);}
        depthError=PhysicalDepthBound(e.Data(target),e.Data(deeper),scale,false,e.theta.size(),e.ratios);
        if(e.physicsAudit){e.SetCutoff("off");
        try{const auto uncut=e.DirectStage(count,seeds.data(),8,e.target,x);e.base=original;e.selectedCutoff=savedCutoff;e.SetMirror(e.mirror);cutoffError=PhysicalDepthBound(e.Data(target),e.Data(uncut),scale,false,e.theta.size(),e.ratios);}
        catch(...){e.base=original;e.selectedCutoff=savedCutoff;e.SetMirror(e.mirror);throw;}}
        pass=depthError<=(e.autoDepth?e.depthBudget:e.eps)&&cutoffError<=e.eps;if(pass||count>=e.maxCount||attempt==4)break;count=std::min(e.maxCount,count>e.maxCount/2?e.maxCount:2*count);
    }
    e.Record(e.root+"/physics_validation.json","{\"enabled\":true,\"passed\":"+std::string(pass?"true":"false")+",\"depth_plus2_residual95_scaled_M11\":"+Number(depthError)+",\"cutoffs_off_residual95_scaled_M11\":"+Number(cutoffError)+",\"depth_pass\":"+std::string(depthError<=(e.autoDepth?e.depthBudget:e.eps)?"true":"false")+",\"cutoff_pass\":"+std::string(cutoffError<=e.eps?"true":"false")+",\"cutoffs_checked\":"+std::string(e.physicsAudit?"true":"false")+",\"depth_budget\":"+Number(e.autoDepth?e.depthBudget:e.eps)+",\"count_per_seed\":"+Number(count)+",\"checked_theta\":"+Number(x.size())+",\"requested_physics_retained\":true,\"scope\":\"requested evaluated angles (up to 4097 in automatic depth); depth+2 and optional cutoff sensitivity estimates; no ADDA or infinite-depth guarantee\"}");return pass;
}
int RunAdaptive(Engine& e){
    std::vector<int> depths;for(int d:{std::max(1,e.target-4),std::max(1,e.target-2),e.target})if(depths.empty()||depths.back()!=d)depths.push_back(d);
    AdaptiveState state;state.counts=Budgets(depths.size(),std::vector<int>(e.theta.size(),e.minimum));state.counts[0].assign(e.theta.size(),e.initial);
    state.previous=Matrix(e.theta.size());state.valid.assign(e.theta.size(),false);state.coarse.assign(e.theta.size(),false);
    std::vector<double> tuningAngles;const size_t tuningRows=std::min(e.theta.size(),size_t(std::min(33,e.maxTheta)));
    for(size_t i=0;i<tuningRows;++i)tuningAngles.push_back(e.theta[i*(e.theta.size()-1)/(tuningRows-1)]);
    std::ostringstream levels;levels<<"{\"depths\":[";for(size_t i=0;i<depths.size();++i)levels<<(i?",":"")<<depths[i];levels<<"],\"seeds\":[";
    for(size_t i=0;i<e.seeds.size();++i)levels<<(i?",":"")<<e.seeds[i];levels<<"]}";e.Log("multilevel_definition",levels.str());
    e.PrepareTuning(tuningAngles);
    const bool resumed=RestoreAdaptive(e,state);
    if(e.completedCheckpoint){
        CompleteSavedPOIntegrals(e);
        const std::string path=e.root+"/fullauto_status.json";if(Exists(path)){std::string text=Read(path);const size_t key=text.find("\"status\":\"");if(key!=std::string::npos){const size_t begin=key+10,end=text.find('"',begin);text.replace(begin,end-begin,"converged");Atomic(path,text);}}
        std::cout<<"fullauto: completed outputs and dependencies verified; no calculation repeated"<<std::endl;return 0;
    }
    if(resumed){RestoreWorkingPhi(e,state);RestoreCorrectionProbe(e,state);}
    PreparePointwisePolicy(e,state,resumed);
    e.Record(e.root+"/allocation_policy.json","{\"model\":\"centered joint eight-scramble covariance and Fieller forecast\",\"max_growth_per_round\":4,\"correction_subrounds_per_round\":2,\"correction_subround_max_theta\":33,\"correction_subround_max_base_cost_fraction\":0.1,\"working_phi_recheck\":"+std::string(e.phiRecheck&&!e.fixedPhi&&!e.args.IsCatched("fullauto_pilot")?"true":"false")+",\"forecast_accepts_accuracy\":false,\"pointwise_error_budget_and_independent_audit\":true}");
    if(!resumed){const size_t initial=std::min(e.theta.size(),size_t(e.startTheta));for(size_t i=0;i<initial;++i)state.knots.push_back(i*(e.theta.size()-1)/(initial-1));
        auto ids=AdaptiveGuards(state.knots,e.theta.size());auto x=Angles(e.theta,ids);e.State("calibrating",0,PeakBudgets(state.counts),ids.size(),0,0,0);
        auto calibrate=[&](){e.PrepareMirror(x);e.Calibrate(x);if(!e.VerifyMirror(x)){e.SetMirror(false);e.Calibrate(x);}};
        try{calibrate();if(e.thetaAuditCount)state.knots=e.AuditThetaSupport(state.knots);}catch(const std::exception& error){if(!e.DisableUnsupportedControls(error))throw;calibrate();if(e.thetaAuditCount)state.knots=e.AuditThetaSupport(state.knots);}
        ids=AdaptiveGuards(state.knots,e.theta.size());if(e.phiRecheck&&!e.fixedPhi&&ids.size()>2*x.size()){
            Atomic(e.root+"/initial_calibration.json",Read(e.root+"/calibration.json"));
            try{e.Calibrate(WorkingGrid(e,ids),true);if(!e.VerifyMirror(WorkingGrid(e,ids))){e.SetMirror(false);e.Calibrate(WorkingGrid(e,ids),true);}}
            catch(const std::exception& error){if(!e.DisableUnsupportedControls(error))throw;e.Calibrate(WorkingGrid(e,ids),true);}}
        // On small grids startup, matrix assembly and tracing dominate. The
        // measured regressions do not justify grouping there; retain a
        // common budget. Explicit 'on' remains available for experiments.
        if(e.localPolicy=="auto"&&e.theta.size()*size_t(e.phi)<=256)e.local=false;
        e.TuneTracing(WorkingGrid(e,ids));
        SaveAdaptive(e,state);
    }
    double lastCI=0,lastInterpolation=0;size_t lastNodes=0;
    for(;state.round<e.roundLimit;++state.round){
        const auto ids=AdaptiveGuards(state.knots,e.theta.size());Require(ids.size()<=size_t(e.maxTheta),"theta budget exhausted");const auto x=Angles(e.theta,ids);lastNodes=ids.size();SaveAdaptive(e,state);
        e.State("evaluating",state.round,PeakBudgets(state.counts),ids.size(),lastCI,lastInterpolation,state.streak);
        if(!e.VerifyMirror(x)){e.SetMirror(false);e.Calibrate(WorkingGrid(e,ids),true);std::fill(state.valid.begin(),state.valid.end(),false);std::fill(state.coarse.begin(),state.coarse.end(),false);state.streak=0;SaveAdaptive(e,state);--state.round;continue;}
        e.VerifyDepthSharing(WorkingGrid(e,ids));
        AdaptiveData stage;Samples data;Statistics stats;
        try{stage=EvaluateBudget(e,depths,state.counts,ids,e.seeds.data());data=Combine(stage.levels);stats=Estimate(data);
            e.lastAnalyticError=0;for(size_t r=0;r<ids.size();++r)e.lastAnalyticError=std::max(e.lastAnalyticError,stage.analytic[r]/stats.mean[r][0]);
            Require(!e.controls||e.lastAnalyticError<=e.eps/64,"analytic mean numerical budget exceeded; use complete sampled PO or refine the analytic provider");}
        catch(const std::exception& error){if(!e.DisableUnsupportedControls(error))throw;e.Calibrate(WorkingGrid(e,ids),true);std::fill(state.coarse.begin(),state.coarse.end(),false);std::fill(state.valid.begin(),state.valid.end(),false);state.streak=0;SaveAdaptive(e,state);--state.round;continue;}
        Budgets half=state.counts;std::vector<size_t> check;for(size_t id:ids)if(!state.coarse[id]&&state.counts[0][id]>=2*e.initial){half[0][id]/=2;check.push_back(id);}
        if(!check.empty()){auto previous=EvaluateBudget(e,{depths[0]},{half[0]},check,e.seeds.data());const auto pm=Mean(previous.levels[0]),cm=Mean(stage.levels[0]);
            for(size_t k=0;k<check.size();++k){const size_t r=std::lower_bound(ids.begin(),ids.end(),check[k])-ids.begin();double change=0;for(int j=0;j<16;++j)change=std::max(change,std::fabs(cm[r][j]-pm[k][j])/stats.mean[r][0]);state.coarse[check[k]]=change<=NumericalBudget(e);}}
        auto interp=InterpolateBudget(e,state.knots,ids,data);
        const bool common=e.local&&VaryingBudgets(state.counts,ids);
        e.Record(e.root+"/interpolation_production_mode.json","{\"common_prefix\":"+std::string(common?"true":"false")+"}");
        if(common){auto checked=CommonPrefixInterpolation(e,state,depths,ids,e.seeds.data());interp.errors=checked.errors;interp.matrixErrors=checked.matrixErrors;interp.ratioErrors=checked.ratioErrors;interp.maximum=checked.maximum;interp.bad=checked.bad;}
        bool positive=true;for(const auto& m:Mean(interp.output))positive=positive&&m[0]>0;

        Statistics output;if(positive)output=Estimate(interp.output);double outputCI=positive?StatisticalError(e,output):1e100;
        double change=0;bool previous=true,coarse=true;for(size_t r=0;r<ids.size();++r){const size_t id=ids[r];previous=previous&&state.valid[id];coarse=coarse&&state.coarse[id];if(state.valid[id])for(int j=0;j<16;++j)change=std::max(change,std::fabs(stats.mean[r][j]-state.previous[id][j])/stats.mean[r][0]);}
        if(!previous)change=1e100;
        PointwiseDecision budget;if(positive)budget=CheckPointwiseBudget(e,stats,output,interp,ids,"production");
        const bool pass=positive&&budget.directPass&&budget.outputPass&&interp.maximum<=e.interpolation&&change<=NumericalBudget(e)&&coarse;
        const std::string roundSnapshot=SaveRoundDiagnostics(e,state,stage,stats,ids,interp.errors,x);
        if(positive){Atomic(roundSnapshot+"/pointwise_budget.csv",Read(e.root+"/pointwise_budget_production.csv"));WriteFullAutoIntegrals(e,roundSnapshot,interp.output,stage.sources,false);}
        e.Log("round_decision","{\"round\":"+Number(state.round)+",\"statistics_pass\":"+(budget.directPass?"true":"false")
            +",\"output_statistics_pass\":"+(budget.outputPass?"true":"false")+",\"interpolation_pass\":"+(interp.maximum<=e.interpolation?"true":"false")
            +",\"stability_pass\":"+(change<=NumericalBudget(e)?"true":"false")+",\"half_prefix_pass\":"+(coarse?"true":"false")+",\"previous_available\":"+(previous?"true":"false")
            +",\"output_positive_M11\":"+(positive?"true":"false")+",\"statistical95\":"+Number(StatisticalError(e,stats))+",\"output_statistical95\":"+Number(outputCI)
            +",\"error_budget_policy\":"+Quoted(BudgetPolicy())+",\"combined_direct95\":"+Number(positive?budget.direct:1e100)+",\"combined_output95\":"+Number(positive?budget.output:1e100)+",\"stability_budget\":"+Number(NumericalBudget(e))+",\"interpolation_residual95\":"+Number(interp.maximum)+",\"change_scaled_M11\":"+Number(change)+",\"statistical_budget\":"+Number(e.statistical)+",\"interpolation_budget\":"+Number(e.interpolation)+"}");
        UpdatePassStreak(e,state,pass,stage,roundSnapshot);lastCI=std::max(stats.maxInterval,positive?output.maxInterval:1e100);lastInterpolation=interp.maximum;e.lastRatio=std::max(stats.maxRatioInterval,positive?output.maxRatioInterval:1e100);
        e.Record(e.root+"/intermediate_current.json","{\"converged\":false,\"round\":"+Number(state.round)+",\"snapshot_folder\":"+Quoted(roundSnapshot)
            +",\"matrix\":"+Quoted(roundSnapshot+"/mueller_evaluated.dat")+",\"accuracy_and_parameters\":"+Quoted(roundSnapshot+"/parameters_accuracy.csv")
            +",\"integrals\":"+(positive?Quoted(roundSnapshot+"/mueller_fullauto_integrals.tsv"):"null")
            +",\"statistical95\":"+Number(StatisticalError(e,stats))+",\"output_statistical95\":"+Number(outputCI)+",\"interpolation_residual95\":"+Number(interp.maximum)
            +",\"error_budget_policy\":"+Quoted(BudgetPolicy())+",\"combined_direct95\":"+Number(positive?budget.direct:1e100)+",\"combined_output95\":"+Number(positive?budget.output:1e100)+",\"numerical_round_passed\":"+(pass?"true":"false")+",\"pass_streak\":"+Number(state.streak)+",\"required_stable_passes\":"+Number(e.stablePasses)+"}");
        e.State("active",state.round,PeakBudgets(state.counts),ids.size(),lastCI,lastInterpolation,state.streak);
        BudgetFile(e,state,ids);Atomic(e.root+"/final_stage_sources.tsv",stage.sources);
        std::ostringstream record;record<<state.round<<' '<<PeakBudgets(state.counts)[0]<<' '<<ids.size()<<' '<<state.knots.size()<<' '<<Number(stats.maxInterval)<<' '<<Number(outputCI)<<' '<<Number(interp.maximum)<<' '<<Number(change)<<' '<<state.streak;HistoryRecord(e,state.round,record.str());
        if(state.streak>=e.stablePasses){e.State("validating",state.round,PeakBudgets(state.counts),ids.size(),lastCI,lastInterpolation,state.streak);
            if(IndependentAudit(e,state,depths,stage,ids,&interp)){
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
                e.Record(e.root+"/precision_contract.json","{\"error_budget_policy\":"+Quoted(BudgetPolicy())+",\"fixed_error_reserve\":"+Number(FixedBudgetReserve(e))+",\"statistics\":\"eight-scramble pointwise Student-t estimates\",\"normalized_ratios\":\"paired-covariance Fieller intervals for ratios of means\",\"interpolation_budget\":"+Number(e.interpolation)+",\"execution_budget\":"+Number((e.resourcePolicy=="auto"||e.config.useGpu)?e.eps/1024:0)+",\"rotation_budget\":"+Number((e.rotationBeta>1||e.rotationGamma>1)?e.eps/256:0)+",\"cutoff_budget\":"+Number(e.selectedCutoff!="requested"&&e.selectedCutoff!="off"?e.eps/8:0)+",\"selected_cutoff\":"+Quoted(e.selectedCutoff)+",\"automatic_cutoff_audit\":"+std::string(e.tuneCutoffs&&e.selectedCutoff!="off"?"true":"false")+",\"physical_depth\":"+Number(e.target)+",\"physical_accuracy_audited\":"+std::string(e.physicsAudit||e.autoDepth?"true":"false")+",\"automatic_depth\":"+std::string(e.autoDepth?"true":"false")+",\"depth_budget\":"+Number(e.depthBudget)+",\"analytic_mean_error\":\"reported internal refinement, not an absolute bound\",\"global_simultaneous_confidence\":false}");
                if(!physics){SaveAdaptive(e,state,"physics_validation_failed");e.State("physics_validation_failed",state.round,PeakBudgets(state.counts),ids.size(),lastCI,lastInterpolation,state.streak);return 4;}
                WriteFullAutoIntegrals(e,e.root,interp.output,stage.sources,true);
                std::ostringstream finalHashes;for(const char* path:{"controller_command.json","error_budget_policy.json","pointwise_budget_production.csv","mueller_fullauto.dat","mueller_evaluated.dat","all_mueller_confidence.csv","normalized_mueller_confidence.csv","theta_support.csv","theta_evaluated.csv","orientation_budget.csv","final_stage_sources.tsv","independent_final_validation.json","precision_contract.json","physics_validation.json","execution_tuning.tsv","resource_calibration.json","rotation_verification.json","cutoff_calibration.json","automatic_cutoff_validation.json"})finalHashes<<path<<' '<<Hash(Read(e.root+"/"+path))<<'\n';if(e.finalAudit)finalHashes<<"independent_orientation_budget.csv "<<Hash(Read(e.root+"/independent_orientation_budget.csv"))<<'\n';if(e.finalAudit){const auto name="pointwise_budget_independent_serial"+Number(state.auditSerial-1)+".csv";finalHashes<<name<<' '<<Hash(Read(e.root+"/"+name))<<'\n';}if(e.finalAudit)for(const char* name:{"angular_audit_plan.json","angular_audit_training.csv","angular_audit_budget.csv","angular_audit_groups.json","angular_audit_support.csv","angular_audit_sources.tsv","production_validation_receipt.tsv","independent_audit_plan.json","independent_audit_cost_gate.json","independent_audit_reservation.tsv","independent_final_sources.tsv","independent_theta.csv","interpolation_validation_mode.json","interpolation_validation_sources.tsv","interpolation_validation_budgets.csv"})if(Exists(e.root+"/"+name))finalHashes<<name<<' '<<Hash(Read(e.root+"/"+name))<<'\n';if(e.autoDepth)for(const char* name:{"automatic_depth.tsv","depth_selection.json"})finalHashes<<name<<' '<<Hash(Read(e.root+"/"+name))<<'\n';if(Exists(e.root+"/depth_startup.tsv"))for(const char* name:{"depth_startup.tsv","depth_startup.json"})finalHashes<<name<<' '<<Hash(Read(e.root+"/"+name))<<'\n';HashIntegralProducts(e.root,finalHashes);HashSchedulingReceipts(e,finalHashes);Atomic(e.root+"/final_outputs.fnv",finalHashes.str());
                SaveAdaptive(e,state,"converged");e.State("converged",state.round,PeakBudgets(state.counts),ids.size(),lastCI,lastInterpolation,state.streak);return 0;
            }
            state.streak=0;std::fill(state.valid.begin(),state.valid.end(),false);SaveAdaptive(e,state);continue;
        }
        for(size_t r=0;r<ids.size();++r){state.previous[ids[r]]=stats.mean[r];state.valid[ids[r]]=true;}
        if(!interp.bad.empty()||!positive){std::set<size_t> next(state.knots.begin(),state.knots.end());next.insert(interp.bad.begin(),interp.bad.end());if(!positive)next.insert(ids.begin(),ids.end());
            state.knots.assign(next.begin(),next.end());for(size_t id:AdaptiveGuards(state.knots,e.theta.size()))if(!std::binary_search(ids.begin(),ids.end(),id)){
                auto upper=std::upper_bound(ids.begin(),ids.end(),id);const size_t right=upper==ids.end()?ids.back():*upper,left=upper==ids.begin()?ids.front():*(upper-1);for(auto& level:state.counts)level[id]=std::max(level[left],level[right]);}state.streak=0;SaveAdaptive(e,state);continue;}
        Budgets proposed=state.counts;bool grew=false;
        for(size_t r=0;r<ids.size();++r){const size_t id=ids[r];double error=0;for(int j=0;j<16;++j){error=std::max(error,stats.interval[r][j]);if(e.ratios)error=std::max(error,stats.ratioInterval[r][j]);}
            const auto residual=ResidualAt(e,interp,ids,id);
            if(!MetricWithinBudget(e,BudgetMetric(e,stats,r,residual),residual)){auto ns=AllocateRow(e,stage,state.counts,ids,r,stats.mean,RemainingStatistics(e,residual));for(size_t l=0;l<ns.size();++l)proposed[l][id]=std::max(state.counts[l][id],std::min(e.maxCount,std::min(state.counts[l][id]>e.maxCount/4?e.maxCount:4*state.counts[l][id],ns[l])));}
            if(!state.coarse[id])proposed[0][id]=std::min(e.maxCount,std::max(proposed[0][id],state.counts[0][id]>e.maxCount/2?e.maxCount:2*state.counts[0][id]));
        }
        if(!e.local){for(size_t l=0;l<proposed.size();++l){int maximum=0;for(size_t id:ids)maximum=std::max(maximum,proposed[l][id]);for(size_t id:ids)proposed[l][id]=maximum;}}
        for(size_t l=0;l<proposed.size();++l)for(size_t id:ids)grew=grew||proposed[l][id]>state.counts[l][id];
        // A correction may dominate uncertainty while the base already has
        // millions of samples. Do not pay for another huge base before the
        // cheap correction has supplied a more informative variance estimate.
        AdditionalWork work(e,e.seeds.data());
        size_t deferred=0;double jointWork=0,stagedWork=0;std::ostringstream priorityRows;
        priorityRows<<"theta_deg,level,current_count,joint_proposal,staged_proposal,estimated_additional_worker_seconds_per_orientation,cost_profiled\n";
        for(size_t r=0;r<ids.size();++r){const size_t id=ids[r];std::vector<Samples> values;
            std::vector<double> costs;std::vector<int> current,next;std::vector<bool> profiled;
            for(size_t l=0;l<depths.size();++l){Samples row(8,Matrix(1));for(int s=0;s<8;++s)row[s][0]=stage.levels[l][s][r];
                values.push_back(row);costs.push_back(8*stage.costs[l][r]);current.push_back(state.counts[l][id]);next.push_back(proposed[l][id]);profiled.push_back(false);
                if(next.back()>current.back()){
                    const auto x=std::vector<double>{std::min(e.theta[id],e.theta[id==ids.front()?ids.back():ids.front()]),std::max(e.theta[id],e.theta[id==ids.front()?ids.back():ids.front()])};
                    auto requests=[&](int n){std::vector<Engine::WorkRequest> result={{n,0,depths[l],l&&e.depthSharing?depths[l-1]:0,x}};
                        if(l&&!e.depthSharing)result.push_back(Engine::WorkRequest{n,0,depths[l-1],0,x});return result;};
                    const auto before=work.Forecast(requests(current.back())),after=work.Forecast(requests(next.back()));
                    profiled[l]=before.profiled&&after.profiled;if(profiled[l])costs[l]=std::max(1e-18,(after.seconds-before.seconds)/(next.back()-current.back()));
                }}
            const auto selected=PrioritizeAllocation(values,costs,current,next,RemainingStatistics(e,ResidualAt(e,interp,ids,id)),e.ratios);
            for(size_t l=0;l<depths.size();++l){jointWork+=costs[l]*(next[l]-current[l]);stagedWork+=costs[l]*(selected[l]-current[l]);
                if(selected!=next)priorityRows<<Number(e.theta[id])<<','<<l<<','<<current[l]<<','<<next[l]<<','<<selected[l]<<','<<Number(costs[l])<<','<<profiled[l]<<'\n';
                deferred+=selected[l]<next[l];proposed[l][id]=selected[l];}
        }
        if(!e.local)for(auto& level:proposed){int maximum=0;for(size_t id:ids)maximum=std::max(maximum,level[id]);for(size_t id:ids)level[id]=maximum;}
        const std::string priorityFile="level_growth_priority_round"+Number(state.round)+"_"+Hash(priorityRows.str())+".csv";
        Atomic(e.root+"/"+priorityFile,priorityRows.str());
        e.Record(e.root+"/level_growth_priority.json","{\"policy\":\"cheap dominant correction first\",\"details\":"+Quoted(priorityFile)+",\"deferred_level_rows\":"+Number(deferred)
            +",\"joint_row_cost_proxy\":"+Number(jointWork)+",\"staged_row_cost_proxy\":"+Number(stagedWork)
            +",\"proxy_scope\":\"before shared tracing and local cost gate; not wall time\",\"minimum_efficiency_ratio\":8,\"accuracy_gates_unchanged\":true}");
        std::ostringstream refinement;refinement<<RefinementHeader();
        if(positive&&!budget.outputPass)for(size_t q=0;q<e.theta.size();++q){const auto residual=ResidualAt(e,interp,ids,q);const auto metric=BudgetMetric(e,output,q,residual);if(MetricWithinBudget(e,metric,residual))continue;
            RefineError(e,work,depths,state.counts,proposed,ids,stage,QueryWeights(e,state.knots,e.theta[q]),metric,output.mean[q],e.theta[q],"interpolated_output",refinement);}
        grew=false;for(size_t l=0;l<proposed.size();++l)for(size_t id:ids)grew=grew||proposed[l][id]>state.counts[l][id];
        if(!grew){ // A fresh prefix is still required for a stability pass.
            size_t worst=0;for(size_t r=1;r<ids.size();++r)if(WorstMetric(stats,r,e.ratios).error>WorstMetric(stats,worst,e.ratios).error)worst=r;
            RefineError(e,work,depths,state.counts,proposed,ids,stage,{{ids[worst],1}},WorstMetric(stats,worst,e.ratios),stats.mean[worst],e.theta[ids[worst]],"stability",refinement);}
        RefinementLog(e,state.round,"error_growth_priority",refinement);
        const auto phiPlan=ReadWorkingPhi(e);
        if(!(phiPlan.oldPhi&&phiPlan.phase==0)&&TryCorrectionProbe(e,state,depths,stage,ids,stats,proposed,lastInterpolation)){--state.round;continue;}
        if(TryWorkingPhi(e,state,depths,stage,ids,stats,lastInterpolation)){--state.round;continue;}
        GateLocalCost(e,work,proposed,depths,ids);
        if(proposed==state.counts)break;
        state.counts=proposed;SaveAdaptive(e,state);
    }
    e.State("requires_further_sampling",state.round,PeakBudgets(state.counts),lastNodes,lastCI,lastInterpolation,state.streak);return 3;
}
} // anonymous namespace
