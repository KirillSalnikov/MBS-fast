// Late phi rescue uses fresh screening and confirmation scrambles, frozen
// control coefficients, and an explicit cost for abandoning the old phi cache.
// It never accepts accuracy or combines matrices from different phi rules.
struct WorkingPhiPlan {
    int phase=0,oldPhi=0,newPhi=0,count=0;std::string sourceHash;
    Budgets counts;std::vector<size_t> rows;
};
std::string WorkingPhiPath(const Engine& e){return e.root+"/working_phi_n"+Number(e.target)+".tsv";}
void SaveWorkingPhi(Engine& e,const WorkingPhiPlan& plan){
    std::ostringstream file;file<<e.signature<<' '<<plan.phase<<' '<<plan.oldPhi<<' '<<plan.newPhi<<' '<<plan.count<<' '<<plan.sourceHash<<' '<<plan.rows.size()<<' '<<plan.counts.size()<<'\n';
    for(size_t row:plan.rows)file<<row<<' ';file<<'\n';
    for(const auto& level:plan.counts){file<<level.size();for(int n:level)file<<' '<<n;file<<'\n';}
    Atomic(WorkingPhiPath(e),Hash(file.str())+"\n"+file.str());
}
WorkingPhiPlan ReadWorkingPhi(Engine& e){
    WorkingPhiPlan p;if(!Exists(WorkingPhiPath(e)))return p;
    const auto text=Read(WorkingPhiPath(e));const auto end=text.find('\n');Require(end!=std::string::npos&&Hash(text.substr(end+1))==text.substr(0,end),"working phi receipt changed");
    std::istringstream file(text.substr(end+1));std::string signature;size_t rows,levels;
    Require(bool(file>>signature>>p.phase>>p.oldPhi>>p.newPhi>>p.count>>p.sourceHash>>rows>>levels)&&signature==e.signature&&p.phase>=0&&p.phase<=3
        &&p.oldPhi>=e.minPhi&&p.oldPhi<=e.maxPhi&&p.count>0&&p.count<=e.maxCount&&rows<=9&&levels<=8,"invalid working phi receipt");
    p.rows.resize(rows);for(auto& row:p.rows)Require(bool(file>>row)&&row<e.theta.size(),"invalid working phi angle");
    p.counts.resize(levels);for(auto& level:p.counts){size_t n;Require(bool(file>>n)&&n==e.theta.size(),"working phi budget dimensions changed");level.resize(n);
        for(auto& count:level)Require(bool(file>>count)&&count>0&&count<=e.maxCount&&(count&(count-1))==0,"invalid working phi budget");}
    std::string extra;Require(!(file>>extra),"extra working phi fields");return p;
}
void ApplyWorkingPhi(Engine& e,AdaptiveState& state,WorkingPhiPlan& plan){
    Require(plan.phase==2&&plan.newPhi>plan.oldPhi&&plan.newPhi>=e.minPhi&&plan.newPhi<=e.maxPhi&&plan.counts.size()==state.counts.size(),"invalid working phi promotion");
    Require(e.phi==plan.oldPhi||e.phi==plan.newPhi,"working phi checkpoint disagrees");
    const auto archive=e.root+"/phi_attempts/n"+Number(e.target)+"_phi"+Number(plan.oldPhi);Directory(archive);
    // Keep every old native block, matrix and log. Only the live allocation
    // changes; a future audit retains its monotonically increasing serial.
    for(const char* name:{"adaptive_checkpoint.tsv","final_stage_sources.tsv","intermediate_current.json","calibration.json","polarization_weights.tsv"})
        if(Exists(e.root+"/"+name)&&!Exists(archive+"/"+name))Atomic(archive+"/"+name,Read(e.root+"/"+name));
    e.phi=plan.newPhi;e.depthVerified=false;e.verifiedMirrorGrids.clear();
    state.counts=plan.counts;state.streak=0;state.previous.assign(e.theta.size(),Mueller{});
    std::fill(state.valid.begin(),state.valid.end(),false);std::fill(state.coarse.begin(),state.coarse.end(),false);
    e.independentPassed=false;SaveAdaptive(e,state);plan.phase=3;SaveWorkingPhi(e,plan);
    e.Log("working_phi_promoted","{\"previous_phi\":"+Number(plan.oldPhi)+",\"phi\":"+Number(plan.newPhi)+",\"old_blocks_retained\":true,\"production_restarted_as_separate_estimator\":true,\"frozen_control_coefficients_retained\":true,\"fresh_final_audit_required\":true}");
}
void RestoreWorkingPhi(Engine& e,AdaptiveState& state){
    auto plan=ReadWorkingPhi(e);if(plan.oldPhi&&plan.phase==2)ApplyWorkingPhi(e,state,plan);
}
struct WorkingPhiPilot {double error=1e100;std::vector<double> cost;std::string sources;};
WorkingPhiPilot EvaluateWorkingPhi(Engine& e,const std::vector<int>& depths,const Budgets& current,
        const std::vector<size_t>& ids,const Matrix& reference,int n,const unsigned* seeds){
    const auto x=Angles(e.theta,ids);WorkingPhiPilot result;std::vector<Engine::WorkRequest> requests;
    for(size_t l=0;l<depths.size();++l){requests.push_back({n,0,depths[l],l&&e.depthSharing?depths[l-1]:0,x});
        if(l&&!e.depthSharing)requests.push_back({n,0,depths[l-1],0,x});}
    e.Prefetch(requests,seeds);Samples projected(8,reference);std::ostringstream sources;
    for(size_t l=0;l<depths.size();++l){std::vector<Job> upper,lower;
        if(l&&e.depthSharing){auto pair=e.SharedStage(n,seeds,8,depths[l],depths[l-1],x);upper=pair.first;lower=pair.second;}
        else upper=e.Stage(n,seeds,8,depths[l],x);
        auto data=e.Data(upper);double seconds=8*n*e.Cost(upper,n);
        for(const auto& job:upper)sources<<l<<" 1 "<<n<<' '<<job.seed<<' '<<job.folder<<'\n';
        if(l){if(lower.empty())lower=e.Stage(n,seeds,8,depths[l-1],x);data=Subtract(data,e.Data(lower));seconds+=8*n*e.Cost(lower,n);
            for(const auto& job:lower)sources<<l<<" -1 "<<n<<' '<<job.seed<<' '<<job.folder<<'\n';}
        result.cost.push_back(seconds/n/ids.size());const auto mean=Mean(data);
        for(size_t row=0;row<ids.size();++row){const double scale=std::sqrt(double(n)/current[l][ids[row]]);
            for(int s=0;s<8;++s)for(int j=0;j<16;++j)projected[s][row][j]+=scale*(data[s][row][j]-mean[row][j]);}
    }
    result.error=StatisticalError(e,Estimate(projected));result.sources=sources.str();return result;
}
double WorkingPhiCost(const Engine& e,const WorkingPhiPilot& pilot,const Budgets& counts,const std::vector<size_t>& ids,bool saved){
    if(!(pilot.error<1e50))return 1e200;const double factor=std::pow(pilot.error/e.statistical,2);double cost=0;
    for(size_t l=0;l<counts.size();++l)for(size_t id:ids){const int floor=l?e.minimum:e.initial;
        const double wanted=std::min(double(e.maxCount),std::max(double(floor),counts[l][id]*factor));
        cost+=pilot.cost[l]*(saved?std::max(0.,wanted-counts[l][id]):wanted);}
    return cost;
}
bool WorkingPhiImproves(const WorkingPhiPilot& baseline,const WorkingPhiPilot& candidate,double oldRemaining,double newTotal,double probeCost){
    return baseline.error<1e50&&candidate.error<1e50&&candidate.error<=.7*baseline.error
        &&oldRemaining>0&&newTotal+probeCost<.5*oldRemaining;
}
bool TryWorkingPhi(Engine& e,AdaptiveState& state,const std::vector<int>& depths,const AdaptiveData& production,
        const std::vector<size_t>& ids,const Statistics& stats,double interpolationError){
    if(!e.coherent||!e.phiRecheck||e.fixedPhi||e.args.IsCatched("fullauto_pilot")||e.phi>=16||e.maxPhi<=e.phi)return false;
    auto plan=ReadWorkingPhi(e);if(plan.oldPhi&&plan.phase!=0)return false;
    const bool pending=plan.oldPhi!=0;
    if(!pending){if(state.round<2||PeakBudgets(state.counts)[0]<std::max(4096,8*e.initial)||StatisticalError(e,stats)<=4*e.statistical)return false;
        std::vector<std::pair<double,size_t>> worst;for(size_t row=0;row<ids.size();++row)worst.push_back({WorstMetric(stats,row,e.ratios).error,ids[row]});
        std::sort(worst.rbegin(),worst.rend());std::set<size_t> chosen;
        // The eight worst sampled locations plus a distant anchor. No final
        // holdout is inspected or used to choose the test geometry.
        const double separation=(e.theta[ids.back()]-e.theta[ids.front()])/64;
        for(const auto& row:worst){bool near=false;for(size_t id:chosen)near=near||std::fabs(e.theta[id]-e.theta[row.second])<separation;
            if(!near)chosen.insert(row.second);if(chosen.size()==8)break;}
        chosen.insert(chosen.count(ids.front())?ids.back():ids.front());plan.rows.assign(chosen.begin(),chosen.end());
        int cap=e.maxCount;if(e.config.adaptive.loadedFromFile&&e.config.adaptive.maxPilotOrientations>0)cap=std::min(cap,e.config.adaptive.maxPilotOrientations);
        plan.count=std::min(4096,Power2(std::max(1024.,4.*e.pilot)));while(plan.count>cap)plan.count/=2;
        if(plan.count<2*e.pilot)return false;
        plan.oldPhi=e.phi;plan.sourceHash=Hash(production.sources);
        SaveWorkingPhi(e,plan);
    }else Require(e.phi==plan.oldPhi&&plan.sourceHash==Hash(production.sources),"pending working phi production changed");
    std::array<unsigned,8> screen,confirm;unsigned offset=unsigned(e.target)*1000000u;
    for(;;){bool overlap=false;for(size_t s=0;s<8;++s){screen[s]=1690000000u+offset+17*s;confirm[s]=1790000000u+offset+17*s;
            overlap=overlap||std::find(e.seeds.begin(),e.seeds.end(),screen[s])!=e.seeds.end()||std::find(e.seeds.begin(),e.seeds.end(),confirm[s])!=e.seeds.end();}
        if(!overlap)break;offset+=512;}
    Matrix reference;for(size_t id:plan.rows){const auto at=std::lower_bound(ids.begin(),ids.end(),id);Require(at!=ids.end()&&*at==id,"working phi angle disappeared");reference.push_back(stats.mean[at-ids.begin()]);}
    std::vector<int> candidates={plan.oldPhi};for(int multiplier:{2,4,8,16})if(plan.oldPhi<=64/multiplier){const int p=plan.oldPhi*multiplier;if(p>=e.minPhi&&p<=e.maxPhi)candidates.push_back(p);}
    if(candidates.size()>4){if(plan.oldPhi==1)candidates.erase(candidates.begin()+1);else candidates.resize(4);} // phi=1 tests 1,4,8,16.
    const double started=Seconds();e.State("phi_recheck",state.round,PeakBudgets(state.counts),ids.size(),StatisticalError(e,stats),interpolationError,state.streak);
    WorkingPhiPilot baseline,winner;double oldCost=0,best=1e200,probeCost=0;int selected=plan.oldPhi;bool confirmed=false,budgetExhausted=false;
    std::ostringstream report;report<<"{\"count_per_seed\":"<<plan.count<<",\"theta_rows\":"<<plan.rows.size()<<",\"screening\":[";
    try{
        for(size_t c=0;c<candidates.size();++c){e.phi=candidates[c];auto value=EvaluateWorkingPhi(e,depths,state.counts,plan.rows,reference,plan.count,screen.data());
            for(double cost:value.cost)probeCost+=cost*plan.count*plan.rows.size();
            Atomic(e.root+"/working_phi_n"+Number(e.target)+"_phi"+Number(e.phi)+"_screen_sources.tsv",value.sources);
            const double cost=WorkingPhiCost(e,value,state.counts,ids,c==0);
            if(c==0){baseline=value;oldCost=cost;}else if(cost<best){best=cost;winner=value;selected=e.phi;}
            report<<(c?",":"")<<"{\"phi\":"<<e.phi<<",\"projected95\":"<<Number(value.error)<<",\"projected_worker_cost\":"<<Number(cost)<<"}";
            if(probeCost>.05*oldCost){budgetExhausted=true;break;} // Bound optional search overhead.
        }
        report<<"],\"screening_old_remaining_cost\":"<<Number(oldCost)<<",\"candidate\":"<<selected;
        bool accepted=selected!=plan.oldPhi&&WorkingPhiImproves(baseline,winner,oldCost,best,probeCost);double ratio=1;
        if(accepted){confirmed=true;e.phi=plan.oldPhi;const auto before=EvaluateWorkingPhi(e,depths,state.counts,plan.rows,reference,plan.count,confirm.data());
            e.phi=selected;const auto after=EvaluateWorkingPhi(e,depths,state.counts,plan.rows,reference,plan.count,confirm.data());
            for(double cost:before.cost)probeCost+=cost*plan.count*plan.rows.size();for(double cost:after.cost)probeCost+=cost*plan.count*plan.rows.size();
            Atomic(e.root+"/working_phi_n"+Number(e.target)+"_confirmation_old_sources.tsv",before.sources);Atomic(e.root+"/working_phi_n"+Number(e.target)+"_confirmation_new_sources.tsv",after.sources);
            const double old=WorkingPhiCost(e,before,state.counts,ids,true),next=WorkingPhiCost(e,after,state.counts,ids,false);
            accepted=WorkingPhiImproves(before,after,old,next,probeCost);
            ratio=std::max(std::pow(winner.error/baseline.error,2),std::pow(after.error/before.error,2));
            report<<",\"confirmation_old95\":"<<Number(before.error)<<",\"confirmation_new95\":"<<Number(after.error)<<",\"confirmation_old_remaining_cost\":"<<Number(old)<<",\"confirmation_new_total_cost\":"<<Number(next);
        }
        e.phi=plan.oldPhi;plan.phase=accepted?2:1;plan.newPhi=accepted?selected:plan.oldPhi;
        if(accepted){plan.counts=state.counts;ratio=std::max(1./64.,std::min(1.,ratio));
            for(size_t l=0;l<plan.counts.size();++l)for(auto& n:plan.counts[l])n=std::min(e.maxCount,Power2(std::max(double(l?e.minimum:e.initial),n*ratio)));}
        report<<",\"accepted\":"<<(accepted?"true":"false")<<",\"probe_worker_cost\":"<<Number(probeCost)<<",\"elapsed_seconds\":"<<Number(Seconds()-started)
            <<",\"cost_scope\":\"1/N variance and native worker cost forecast, not measured campaign speed\",\"independent_confirmation\":"<<(confirmed?"true":"false")
            <<",\"screening_budget_exhausted\":"<<(budgetExhausted?"true":"false")<<",\"sample_count_scale\":"<<Number(accepted?std::max(1./64.,std::min(1.,ratio)):1.)
            <<",\"screening_seeds\":[";for(size_t s=0;s<8;++s)report<<(s?",":"")<<screen[s];report<<"],\"confirmation_seeds\":[";for(size_t s=0;s<8;++s)report<<(s?",":"")<<confirm[s];
        report<<"],\"frozen_controls\":true,\"accuracy_gates_unchanged\":true}";
        e.Record(e.root+"/working_phi_n"+Number(e.target)+".json",report.str());SaveWorkingPhi(e,plan);
    }catch(const std::exception& error){e.phi=plan.oldPhi;plan.phase=1;plan.newPhi=plan.oldPhi;SaveWorkingPhi(e,plan);
        e.Log("working_phi_probe_rejected","{\"reason\":"+Quoted(error.what())+",\"original_estimator_retained\":true}");}
    // A committed promotion is a checkpoint transaction. Never convert its
    // I/O failure into an optional-pilot rejection with a half-reset state.
    if(plan.phase==2){ApplyWorkingPhi(e,state,plan);return true;}
    return false;
}
void HashSchedulingReceipts(const Engine& e,std::ostringstream& hashes){
    std::vector<std::string> names;DIR* directory=opendir(e.root.c_str());Require(directory,"cannot list scheduling receipts");
    while(auto* entry=readdir(directory)){const std::string name=entry->d_name;
        if(name=="allocation_policy.json"||name.compare(0,20,"correction_subround_")==0||name.compare(0,13,"working_phi_n")==0){
            if(name.size()<4||name.compare(name.size()-4,4,".tmp")!=0)names.push_back(name);}}
    closedir(directory);std::sort(names.begin(),names.end());for(const auto& name:names)hashes<<name<<' '<<FileHash(e.root+"/"+name)<<'\n';
}
