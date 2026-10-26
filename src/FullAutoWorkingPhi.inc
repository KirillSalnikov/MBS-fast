// Late phi rescue uses fresh screening and confirmation scrambles, frozen
// control coefficients, and an explicit cost for abandoning the old phi cache.
// It never accepts accuracy or combines matrices from different phi rules.
struct WorkingPhiPlan {
    int phase=0,oldPhi=0,newPhi=0,count=0,attempt=0,peak=0;std::string sourceHash;
    Budgets counts;std::vector<size_t> rows,retained;
};
std::string WorkingPhiPath(const Engine& e){return e.root+"/working_phi_n"+Number(e.target)+".tsv";}
void SaveWorkingPhi(Engine& e,const WorkingPhiPlan& plan){
    std::ostringstream file;file<<e.signature<<' '<<plan.phase<<' '<<plan.oldPhi<<' '<<plan.newPhi<<' '<<plan.count<<' '<<plan.sourceHash<<' '<<plan.rows.size()<<' '<<plan.counts.size()<<'\n';
    for(size_t row:plan.rows)file<<row<<' ';
    file<<'\n';
    for(const auto& level:plan.counts){file<<level.size();for(int n:level)file<<' '<<n;file<<'\n';}
    file<<"retry "<<plan.attempt<<' '<<plan.peak<<'\n';
    if(!plan.retained.empty()){file<<"retained "<<plan.retained.size();for(size_t row:plan.retained)file<<' '<<row;file<<'\n';}
    Atomic(WorkingPhiPath(e),Hash(file.str())+"\n"+file.str());
}
WorkingPhiPlan ReadWorkingPhi(Engine& e){
    WorkingPhiPlan p;if(!Exists(WorkingPhiPath(e)))return p;
    const auto text=Read(WorkingPhiPath(e));const auto end=text.find('\n');Require(end!=std::string::npos&&Hash(text.substr(end+1))==text.substr(0,end),"working phi receipt changed");
    std::istringstream file(text.substr(end+1));std::string signature;size_t rows,levels;
    Require(bool(file>>signature>>p.phase>>p.oldPhi>>p.newPhi>>p.count>>p.sourceHash>>rows>>levels)&&signature==e.signature&&p.phase>=0&&p.phase<=4
        &&p.oldPhi>=e.minPhi&&p.oldPhi<=e.maxPhi&&p.count>0&&p.count<=e.maxCount&&rows<=9&&levels<=8,"invalid working phi receipt");
    p.rows.resize(rows);for(auto& row:p.rows)Require(bool(file>>row)&&row<e.theta.size(),"invalid working phi angle");
    p.counts.resize(levels);for(auto& level:p.counts){size_t n;Require(bool(file>>n)&&n==e.theta.size(),"working phi budget dimensions changed");level.resize(n);
        for(auto& count:level)Require(bool(file>>count)&&count>0&&count<=e.maxCount&&(count&(count-1))==0,"invalid working phi budget");}
    std::string extra;
    if(file>>extra){Require(extra=="retry"&&bool(file>>p.attempt>>p.peak)&&p.attempt>=0&&p.attempt<3&&p.peak>=0,"invalid working phi retry metadata");
        if(file>>extra){size_t size;Require(extra=="retained"&&bool(file>>size)&&size<=e.theta.size(),"invalid retained phi rows");
            for(size_t i=0;i<size;++i){size_t row;Require(bool(file>>row)&&row<e.theta.size()&&(p.retained.empty()||row>p.retained.back()),"invalid retained phi index");p.retained.push_back(row);}
            Require(!(file>>extra),"extra working phi fields");}}return p;
}
void ApplyWorkingPhi(Engine& e,AdaptiveState& state,WorkingPhiPlan& plan){
    Require(plan.phase==2&&plan.newPhi>plan.oldPhi&&plan.newPhi>=e.minPhi&&plan.newPhi<=e.maxPhi&&plan.counts.size()==state.counts.size(),"invalid working phi promotion");
    Require(e.phi==plan.oldPhi||e.phi==plan.newPhi,"working phi checkpoint disagrees");
    const auto archive=e.root+"/phi_attempts/n"+Number(e.target)+"_phi"+Number(plan.oldPhi);Directory(archive);
    // Keep every old native block, matrix and log. Only the live allocation
    // changes; a future audit retains its monotonically increasing serial.
    for(const char* name:{"adaptive_checkpoint.tsv","final_stage_sources.tsv","intermediate_current.json","calibration.json","polarization_weights.tsv"})
        if(Exists(e.root+"/"+name)&&!Exists(archive+"/"+name))Atomic(archive+"/"+name,Read(e.root+"/"+name));
    e.phi=plan.newPhi;e.depthVerified=false;e.verifiedMirrorGrids.clear();e.thetaPhi.clear();
    for(size_t row:plan.retained)e.thetaPhi[row]=plan.oldPhi;
    state.counts=plan.counts;state.streak=0;
    for(size_t row=0;row<e.theta.size();++row)if(!e.thetaPhi.count(row)){state.previous[row]=Mueller{};state.valid[row]=false;state.coarse[row]=false;}
    e.independentPassed=false;SaveAdaptive(e,state);plan.phase=3;SaveWorkingPhi(e,plan);
    e.Log("working_phi_promoted","{\"previous_phi\":"+Number(plan.oldPhi)+",\"phi\":"+Number(plan.newPhi)+",\"retained_theta_rows\":"+Number(plan.retained.size())+",\"old_blocks_retained\":true,\"changed_rows_restarted_as_separate_estimators\":true,\"frozen_control_coefficients_retained\":true,\"fresh_final_audit_required\":true}");
}
void RestoreWorkingPhi(Engine& e,AdaptiveState& state){
    auto plan=ReadWorkingPhi(e);if(plan.oldPhi&&plan.phase==2)ApplyWorkingPhi(e,state,plan);
}
struct WorkingPhiPilot {double error=1e100;std::vector<double> cost;std::string sources;};
struct PhiCutoffPilotScope {
    Engine& e;bool previous;
    explicit PhiCutoffPilotScope(Engine& engine,bool corrected=true):e(engine),previous(engine.phiCutoffPilot){e.phiCutoffPilot=e.balancedCorrection&&corrected;}
    ~PhiCutoffPilotScope(){e.phiCutoffPilot=previous;}
};
WorkingPhiPilot EvaluateWorkingPhi(Engine& e,const std::vector<int>& depths,const Budgets& current,
        const std::vector<size_t>& ids,const Matrix& reference,int n,const unsigned* seeds){
    WorkingPhiPilot result;Samples projected(8,reference);std::ostringstream sources;
    for(size_t l=0;l<depths.size();++l){
        // Small correction prefixes (<512) remain direct-off in production.
        // Match that row/level too, rather than extrapolating a corrected pilot.
        std::map<bool,std::vector<size_t>> groups;
        for(size_t r=0;r<ids.size();++r){const int leaf=e.incremental?std::min(e.orientationBlock,current[l][ids[r]]):current[l][ids[r]];
            groups[e.balancedCorrection&&depths[l]>=6&&leaf>=512].push_back(r);}
        double seconds=0;
        for(const auto& group:groups){PhiCutoffPilotScope scope(e,group.first);
            Require(!group.first||std::min(n,e.incremental?e.orientationBlock:n)>=512,"phi pilot too small for the production cutoff estimator");
            std::vector<double> x;for(size_t r:group.second)x.push_back(e.theta[ids[r]]);
            if(x.size()==1){x.push_back(e.theta.front()==x[0]?e.theta.back():e.theta.front());std::sort(x.begin(),x.end());}
            std::vector<Engine::WorkRequest> requests={{n,0,depths[l],l&&e.depthSharing?depths[l-1]:0,x}};
            if(l&&!e.depthSharing)requests.push_back({n,0,depths[l-1],0,x});
            e.Prefetch(requests,seeds);
            std::vector<Job> upper,lower;
            if(l&&e.depthSharing){auto pair=e.SharedStage(n,seeds,8,depths[l],depths[l-1],x);upper=pair.first;lower=pair.second;}
            else upper=e.Stage(n,seeds,8,depths[l],x);
            auto data=e.Data(upper);seconds+=8*n*e.Cost(upper,n);
            for(const auto& job:upper)sources<<l<<" 1 "<<n<<' '<<job.seed<<' '<<job.folder<<'\n';
            if(l){if(lower.empty())lower=e.Stage(n,seeds,8,depths[l-1],x);data=Subtract(data,e.Data(lower));seconds+=8*n*e.Cost(lower,n);
                for(const auto& job:lower)sources<<l<<" -1 "<<n<<' '<<job.seed<<' '<<job.folder<<'\n';}
            const auto mean=Mean(data);
            for(size_t row:group.second){const size_t at=std::lower_bound(x.begin(),x.end(),e.theta[ids[row]])-x.begin();const double scale=std::sqrt(double(n)/current[l][ids[row]]);
                for(int s=0;s<8;++s)for(int j=0;j<16;++j)projected[s][row][j]+=scale*(data[s][at][j]-mean[at][j]);}
        }
        result.cost.push_back(seconds/n/ids.size());
    }
    result.error=StatisticalError(e,Estimate(projected));result.sources=sources.str();return result;
}

double WorkingPhiBaselineError(double sampled,double observed){
    // A short Fieller interval may be unbounded although the long completed
    // production interval is finite. Use that measured production interval as
    // the conservative scheduling baseline, never as a new accuracy verdict.
    if(!FullAutoScheduling::Resolved(observed))return sampled;
    return FullAutoScheduling::Resolved(sampled)?std::max(sampled,observed):observed;
}
double WorkingPhiCost(const Engine& e,const WorkingPhiPilot& pilot,const Budgets& counts,const std::vector<size_t>& ids,bool saved){
    if(!(pilot.error<1e50))return 1e200;
    const double factor=std::pow(pilot.error/e.statistical,2);
    double cost=0;
    for(size_t l=0;l<counts.size();++l)for(size_t id:ids){const int floor=l?e.minimum:e.initial;
        // The work forecast must describe reaching the requested error, not
        // merely reaching the launch limit. Capping it makes two rules look
        // equally cheap even when one still misses the target at maxCount.
        // Actual allocations remain bounded in ApplyWorkingPhi/the allocator.
        const double wanted=std::max(double(floor),counts[l][id]*factor);
        cost+=pilot.cost[l]*(saved?std::max(0.,wanted-counts[l][id]):wanted);}
    return cost;
}
bool WorkingPhiImproves(const WorkingPhiPilot& baseline,const WorkingPhiPilot& candidate,double oldRemaining,double newTotal,double probeCost){
    return baseline.error<1e50&&candidate.error<1e50&&candidate.error<=.7*baseline.error
        &&oldRemaining>0&&newTotal+probeCost<.5*oldRemaining;
}
// This is a scheduling diagnostic, never an accuracy acceptance test. A
// small pilot that misses a rare contribution can give a perfectly finite
// interval. Contradictory independent folds must not lock in a cheap-looking
// phi rule until another multi-million-orientation round has completed.
bool WorkingPhiDisagrees(double screen,double confirmation){
    if(!FullAutoScheduling::Resolved(screen)||!FullAutoScheduling::Resolved(confirmation))return true;
    return std::max(screen,confirmation)>4*std::max(1e-15,std::min(screen,confirmation));
}
std::vector<int> WorkingPhiCandidates(int oldPhi,int minimum,int maximum){
    std::set<int> candidates{oldPhi};
    for(int multiplier:{2,4,8,16})if(oldPhi<=64/multiplier)candidates.insert(oldPhi*multiplier);
    if(oldPhi<32)candidates.insert(32);
    if(oldPhi==1)candidates.erase(2);
    std::vector<int> result;for(int p:candidates)if(p>=minimum&&p<=maximum)result.push_back(p);
    return result;
}
bool WorkingPhiPilotDue(const Engine& e,const AdaptiveState& state,const Statistics& stats){
    const double error=StatisticalError(e,stats);
    // An exceptionally noisy rule warrants a bounded independent pilot after
    // two initial blocks, rather than another six blocks at the same rule.
    const int multiplier=error>16*e.statistical?2:8;
    return state.round>=2&&PeakBudgets(state.counts)[0]>=std::max(4096.,double(multiplier)*e.initial)&&error>4*e.statistical;
}
bool TryWorkingPhi(Engine& e,AdaptiveState& state,const std::vector<int>& depths,const AdaptiveData& production,
        const std::vector<size_t>& ids,const Statistics& stats,double interpolationError,const InterpolationResult* interpolation=nullptr){
    if(!e.coherent||!e.phiRecheck||e.fixedPhi||e.args.IsCatched("fullauto_pilot")||e.phi>=64||e.maxPhi<=e.phi)return false;
    auto plan=ReadWorkingPhi(e);int retryAttempt=0,retryCount=0;
    const auto reportPath=e.root+"/working_phi_n"+Number(e.target)+".json";
    if(e.balancedCorrection&&plan.oldPhi&&(plan.phase==1||plan.phase==4)&&Exists(reportPath)
            &&Read(reportPath).find("\"cutoff_estimator_matches_production\":true")==std::string::npos){
        const auto archive=e.root+"/phi_rechecks/unmatched_cutoff_n"+Number(e.target);Directory(archive);
        for(const char* ext:{".json",".tsv"}){const auto name="working_phi_n"+Number(e.target)+ext;
            if(Exists(e.root+"/"+name)&&!Exists(archive+"/"+name))Atomic(archive+"/"+name,Read(e.root+"/"+name));}
        plan=WorkingPhiPlan();
        e.Log("working_phi_estimator_migration","{\"legacy_diagnostic_archived\":true,\"production_blocks_retained\":true,\"reason\":\"old pilot omitted production cutoff correction\"}");
    }
    if(plan.oldPhi&&plan.phase!=0){
        const int peak=PeakBudgets(state.counts)[0];
        const bool inconclusive=plan.phase==4;
        if((plan.phase!=1&&!inconclusive)||plan.attempt>=2||plan.oldPhi!=e.phi
            ||(!inconclusive&&!(plan.peak==0||FullAutoScheduling::RetryPhi(plan.peak,peak,plan.attempt+1))))return false;
        retryAttempt=plan.attempt+1;retryCount=plan.count>e.maxCount/4?e.maxCount:plan.count*4;
        const auto archive=e.root+"/phi_rechecks/n"+Number(e.target)+"_attempt"+Number(plan.attempt);Directory(archive);
        for(const char* ext:{".json",".tsv"}){const auto name="working_phi_n"+Number(e.target)+ext;
            if(Exists(e.root+"/"+name)&&!Exists(archive+"/"+name))Atomic(archive+"/"+name,Read(e.root+"/"+name));}
        plan=WorkingPhiPlan();plan.attempt=retryAttempt;
        e.Log("working_phi_retry","{\"attempt\":"+Number(retryAttempt)+",\"production_peak\":"+Number(peak)+",\"before_production_growth\":"+(inconclusive?"true":"false")+",\"previous_receipt_retained\":true}");
    }
    const bool pending=plan.oldPhi!=0;
    if(!pending){if(!WorkingPhiPilotDue(e,state,stats))return false;
        std::vector<std::pair<double,size_t>> worst;for(size_t row=0;row<ids.size();++row)worst.push_back({WorstMetric(stats,row,e.ratios).error,ids[row]});
        std::sort(worst.rbegin(),worst.rend());std::set<size_t> chosen;
        // The eight worst sampled locations plus a distant anchor. No final
        // holdout is inspected or used to choose the test geometry.
        const double separation=(e.theta[ids.back()]-e.theta[ids.front()])/64;
        for(const auto& row:worst){bool near=false;for(size_t id:chosen)near=near||std::fabs(e.theta[id]-e.theta[row.second])<separation;
            if(!near)chosen.insert(row.second);
            if(chosen.size()==8)break;
            }
        chosen.insert(chosen.count(ids.front())?ids.back():ids.front());plan.rows.assign(chosen.begin(),chosen.end());
        int cap=e.maxCount;if(e.config.adaptive.loadedFromFile&&e.config.adaptive.maxPilotOrientations>0)cap=std::min(cap,e.config.adaptive.maxPilotOrientations);
        plan.count=retryCount?retryCount:std::min(4096,Power2(std::max(1024.,4.*e.pilot)));plan.count=std::min(plan.count,16384);while(plan.count>cap)plan.count/=2;
        if(plan.count<2*e.pilot)return false;
        plan.oldPhi=e.phi;plan.peak=PeakBudgets(state.counts)[0];plan.sourceHash=Hash(production.sources);
        SaveWorkingPhi(e,plan);
    }else Require(e.phi==plan.oldPhi&&plan.sourceHash==Hash(production.sources),"pending working phi production changed");
    std::array<unsigned,8> screen,confirm;unsigned offset=unsigned(e.target)*1000000u+unsigned(plan.attempt)*100000u+(e.balancedCorrection?50000000u:0u);
    for(;;){bool overlap=false;for(size_t s=0;s<8;++s){screen[s]=1690000000u+offset+17*s;confirm[s]=1790000000u+offset+17*s;
            overlap=overlap||std::find(e.seeds.begin(),e.seeds.end(),screen[s])!=e.seeds.end()||std::find(e.seeds.begin(),e.seeds.end(),confirm[s])!=e.seeds.end();}
        if(!overlap)break;
        offset+=512;
        }
    Matrix reference;for(size_t id:plan.rows){const auto at=std::lower_bound(ids.begin(),ids.end(),id);Require(at!=ids.end()&&*at==id,"working phi angle disappeared");reference.push_back(stats.mean[at-ids.begin()]);}
    double productionError=0;
    for(size_t id:plan.rows){const size_t row=std::lower_bound(ids.begin(),ids.end(),id)-ids.begin();productionError=std::max(productionError,WorstMetric(stats,row,e.ratios).error);}
    const auto candidates=WorkingPhiCandidates(plan.oldPhi,e.minPhi,e.maxPhi);
    const double started=Seconds();e.State("phi_recheck",state.round,PeakBudgets(state.counts),ids.size(),StatisticalError(e,stats),interpolationError,state.streak);
    WorkingPhiPilot baseline,winner;double oldCost=0,best=1e200,probeCost=0;int selected=plan.oldPhi;bool confirmed=false,budgetExhausted=false,disagrees=false;
    std::ostringstream report;report<<"{\"attempt\":"<<plan.attempt<<",\"cutoff_estimator_matches_production\":true,\"cutoff_estimator\":"<<Quoted(e.balancedCorrection?"balanced_with_paired_off_correction":"direct")<<",\"production_peak\":"<<plan.peak<<",\"count_per_seed\":"<<plan.count<<",\"theta_rows\":"<<plan.rows.size()
        <<",\"baseline_policy\":\"observed completed production interval is a floor and fallback for an unbounded short pilot; scheduling only\",\"observed_production95\":"<<Number(productionError)<<",\"screening\":[";
    try{
        for(size_t c=0;c<candidates.size();++c){e.phi=candidates[c];auto value=EvaluateWorkingPhi(e,depths,state.counts,plan.rows,reference,plan.count,screen.data());
            const double sampled=value.error;
            // Millions of actual production samples must not be replaced by
            // an accidentally quiet short pilot of the SAME estimator.
            // This floor affects scheduling only; no production CI changes.
            if(c==0)value.error=WorkingPhiBaselineError(value.error,productionError);
            for(double cost:value.cost)probeCost+=cost*plan.count*plan.rows.size();
            Atomic(e.root+"/working_phi_n"+Number(e.target)+"_phi"+Number(e.phi)+"_screen_sources.tsv",value.sources);
            const double cost=WorkingPhiCost(e,value,state.counts,ids,c==0);
            e.Log("working_phi_candidate","{\"attempt\":"+Number(plan.attempt)+",\"count_per_seed\":"+Number(plan.count)+",\"phi\":"+Number(e.phi)
                +",\"projected95\":"+Number(value.error)+",\"projected_worker_cost\":"+Number(cost)+",\"forecast_only\":true}");
            if(c==0){baseline=value;oldCost=cost;}else if(cost<best){best=cost;winner=value;selected=e.phi;}
            report<<(c?",":"")<<"{\"phi\":"<<e.phi<<",\"pilot_projected95\":"<<Number(sampled)<<",\"projected95\":"<<Number(value.error)<<",\"projected_worker_cost\":"<<Number(cost)<<",\"level_cost\":[";
            for(size_t l=0;l<value.cost.size();++l)report<<(l?",":"")<<Number(value.cost[l]);
            report<<"]}";
            if(probeCost>.05*oldCost){budgetExhausted=true;break;} // Bound optional search overhead.
        }
        report<<"],\"screening_old_remaining_cost\":"<<Number(oldCost)<<",\"candidate\":"<<selected;
        bool accepted=selected!=plan.oldPhi&&WorkingPhiImproves(baseline,winner,oldCost,best,probeCost);double ratio=1;
        if(accepted){confirmed=true;e.phi=plan.oldPhi;auto before=EvaluateWorkingPhi(e,depths,state.counts,plan.rows,reference,plan.count,confirm.data());
            const double sampledBefore=before.error;before.error=WorkingPhiBaselineError(before.error,productionError);
            e.phi=selected;const auto after=EvaluateWorkingPhi(e,depths,state.counts,plan.rows,reference,plan.count,confirm.data());
            for(double cost:before.cost)probeCost+=cost*plan.count*plan.rows.size();
            for(double cost:after.cost)probeCost+=cost*plan.count*plan.rows.size();
            Atomic(e.root+"/working_phi_n"+Number(e.target)+"_confirmation_old_sources.tsv",before.sources);Atomic(e.root+"/working_phi_n"+Number(e.target)+"_confirmation_new_sources.tsv",after.sources);
            const double old=WorkingPhiCost(e,before,state.counts,ids,true),next=WorkingPhiCost(e,after,state.counts,ids,false);
            accepted=WorkingPhiImproves(before,after,old,next,probeCost);
            disagrees=WorkingPhiDisagrees(baseline.error,before.error)||WorkingPhiDisagrees(winner.error,after.error);
            // Even a favourable confirmation does not justify a switch when
            // the two independent variance forecasts contradict each other.
            accepted=accepted&&!disagrees;
            ratio=std::max(std::pow(winner.error/baseline.error,2),std::pow(after.error/before.error,2));
            report<<",\"confirmation_old_pilot95\":"<<Number(sampledBefore)<<",\"confirmation_old95\":"<<Number(before.error)<<",\"confirmation_new95\":"<<Number(after.error)<<",\"confirmation_old_remaining_cost\":"<<Number(old)<<",\"confirmation_new_total_cost\":"<<Number(next);
        }
        int retryLimit=std::min(16384,e.maxCount);
        if(e.config.adaptive.loadedFromFile&&e.config.adaptive.maxPilotOrientations>0)retryLimit=std::min(retryLimit,e.config.adaptive.maxPilotOrientations);
        const bool retryNow=!accepted&&disagrees&&!budgetExhausted&&plan.attempt<2&&plan.count<retryLimit&&probeCost<=.05*oldCost;
        e.phi=plan.oldPhi;plan.phase=accepted?2:(retryNow?4:1);plan.newPhi=accepted?selected:plan.oldPhi;
        if(accepted){plan.counts=state.counts;ratio=std::max(1./64.,std::min(1.,ratio));
            // Retain a complete old estimator when its measured componentwise
            // budget and interpolation test pass. This is not a frozen accuracy
            // certificate: the mixed curve is rechecked by the normal gates.
            if(interpolation)for(size_t r=0;r<ids.size();++r){const size_t id=ids[r];const auto residual=ResidualAt(e,*interpolation,ids,id);
                if(state.coarse[id]&&state.valid[id]&&interpolation->errors[r]<=e.interpolation&&MetricWithinBudget(e,BudgetMetric(e,stats,r,residual),residual))plan.retained.push_back(id);}
            for(size_t l=0;l<plan.counts.size();++l)for(size_t id=0;id<e.theta.size();++id)if(!std::binary_search(plan.retained.begin(),plan.retained.end(),id)){
                auto& n=plan.counts[l][id];n=std::min(e.maxCount,Power2(std::max(double(l?e.minimum:e.initial),n*ratio)));}}
        report<<",\"accepted\":"<<(accepted?"true":"false")<<",\"probe_worker_cost\":"<<Number(probeCost)<<",\"elapsed_seconds\":"<<Number(Seconds()-started)
            <<",\"cost_scope\":\"uncapped 1/N work-to-target forecast; actual allocations remain bounded; not measured campaign speed\",\"independent_confirmation\":"<<(confirmed?"true":"false")
            <<",\"screening_budget_exhausted\":"<<(budgetExhausted?"true":"false")<<",\"sample_count_scale\":"<<Number(accepted?std::max(1./64.,std::min(1.,ratio)):1.)
            <<",\"independent_forecasts_disagree\":"<<(disagrees?"true":"false")<<",\"retry_before_production_growth\":"<<(retryNow?"true":"false")
            <<",\"screening_seeds\":[";for(size_t s=0;s<8;++s)report<<(s?",":"")<<screen[s];report<<"],\"confirmation_seeds\":[";for(size_t s=0;s<8;++s)report<<(s?",":"")<<confirm[s];
        report<<"],\"frozen_controls\":true,\"accuracy_gates_unchanged\":true}";
        e.Record(e.root+"/working_phi_n"+Number(e.target)+".json",report.str());SaveWorkingPhi(e,plan);
    }catch(const std::exception& error){e.phi=plan.oldPhi;plan.phase=1;plan.newPhi=plan.oldPhi;SaveWorkingPhi(e,plan);
        e.Log("working_phi_probe_rejected","{\"reason\":"+Quoted(error.what())+",\"original_estimator_retained\":true}");}
    // A committed promotion is a checkpoint transaction. Never convert its
    // I/O failure into an optional-pilot rejection with a half-reset state.
    if(plan.phase==2){ApplyWorkingPhi(e,state,plan);return true;}
    if(plan.phase==4)return TryWorkingPhi(e,state,depths,production,ids,stats,interpolationError,interpolation);
    return false;
}
void HashSchedulingReceipts(const Engine& e,std::ostringstream& hashes){
    std::vector<std::string> names;DIR* directory=opendir(e.root.c_str());Require(directory,"cannot list scheduling receipts");
    while(auto* entry=readdir(directory)){const std::string name=entry->d_name;
        if(name=="allocation_policy.json"||name.compare(0,20,"correction_subround_")==0||name.compare(0,13,"working_phi_n")==0){
            if(name.size()<4||name.compare(name.size()-4,4,".tmp")!=0)names.push_back(name);}}
    closedir(directory);std::sort(names.begin(),names.end());for(const auto& name:names)hashes<<name<<' '<<FileHash(e.root+"/"+name)<<'\n';
}
