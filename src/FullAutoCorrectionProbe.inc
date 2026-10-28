// A bounded production subround learns a cheap noisy correction before a
// pending expensive base extension. No audit seeds or acceptance tests change.
struct CorrectionProbeRow {size_t id,level;int before,after;};
struct CorrectionProbeReceipt {int attempt=0;bool complete=false;std::string sources;std::vector<CorrectionProbeRow> rows;};
std::string CorrectionProbePath(const Engine& e,const AdaptiveState& state){
    return e.root+"/correction_subround_n"+Number(e.target)+"_phi"+Number(e.phi)+"_round"+Number(state.round)+".tsv";
}
CorrectionProbeReceipt ReadCorrectionProbe(Engine& e,const AdaptiveState& state){
    CorrectionProbeReceipt result;const auto path=CorrectionProbePath(e,state);if(!Exists(path))return result;
    const auto text=Read(path);const auto end=text.find('\n');Require(end!=std::string::npos&&Hash(text.substr(end+1))==text.substr(0,end),"correction subround receipt changed");
    std::istringstream file(text.substr(end+1));std::string signature;size_t size;
    Require(bool(file>>signature>>result.attempt>>result.complete>>result.sources>>size)&&signature==e.signature&&result.attempt>=1&&result.attempt<=2&&size<=128,"invalid correction subround receipt");
    for(size_t i=0;i<size;++i){CorrectionProbeRow r;Require(bool(file>>r.id>>r.level>>r.before>>r.after)&&r.id<e.theta.size()&&r.level>0&&r.level<state.counts.size()
        &&r.before>0&&r.after>r.before&&r.after<=e.maxCount&&(r.after&(r.after-1))==0,"invalid correction subround allocation");result.rows.push_back(r);}
    std::string extra;Require(!(file>>extra),"extra correction subround fields");return result;
}
void SaveCorrectionProbe(Engine& e,const AdaptiveState& state,const CorrectionProbeReceipt& receipt){
    std::ostringstream file;file<<e.signature<<' '<<receipt.attempt<<' '<<receipt.complete<<' '<<receipt.sources<<' '<<receipt.rows.size()<<'\n';
    for(const auto& r:receipt.rows)file<<r.id<<' '<<r.level<<' '<<r.before<<' '<<r.after<<'\n';
    Atomic(CorrectionProbePath(e,state),Hash(file.str())+"\n"+file.str());
}
bool RestoreCorrectionProbe(Engine& e,AdaptiveState& state){
    const auto receipt=ReadCorrectionProbe(e,state);if(!receipt.complete)return false;bool changed=false;
    for(const auto& r:receipt.rows){auto& n=state.counts[r.level][r.id];Require(n==r.before||n>=r.after,"correction subround checkpoint disagrees");if(n==r.before){n=r.after;changed=true;}}
    if(changed){state.streak=0;SaveAdaptive(e,state);e.Log("correction_subround_recovered","{\"round\":"+Number(state.round)+",\"attempt\":"+Number(receipt.attempt)+"}");}return changed;
}
std::vector<Engine::WorkRequest> CorrectionProbeRequests(const Engine& e,const std::vector<int>& depths,const CorrectionProbeReceipt& receipt){
    std::map<std::pair<size_t,int>,std::vector<double>> groups;for(const auto& r:receipt.rows)groups[{r.level,r.after}].push_back(e.theta[r.id]);
    std::vector<Engine::WorkRequest> requests;
    for(auto& group:groups){const auto level=group.first.first;auto& angles=group.second;
        if(angles.size()==1&&e.thetaPhi.empty()){angles.push_back(angles[0]==e.theta.front()?e.theta.back():e.theta.front());std::sort(angles.begin(),angles.end());}
        requests.push_back({group.first.second,0,depths[level],e.depthSharing?depths[level-1]:0,angles});
        if(!e.depthSharing)requests.push_back({group.first.second,0,depths[level-1],0,angles});}
    return requests;
}
bool TryCorrectionProbe(Engine& e,AdaptiveState& state,const std::vector<int>& depths,const AdaptiveData& data,
        const std::vector<size_t>& ids,const Statistics& stats,const Budgets& proposed,double interpolationError){
    if(depths.size()<2||!e.local)return false;
    auto receipt=ReadCorrectionProbe(e,state);const bool pending=receipt.attempt&&!receipt.complete;
    if(!pending){if(receipt.attempt>=2)return false;
        struct Candidate {size_t row,level;double error;int next;};std::vector<Candidate> candidates;
        for(size_t row=0;row<ids.size();++row){const size_t id=ids[row];const auto metric=WorstMetric(stats,row,e.ratios);
            if(metric.error<=2*e.statistical||state.counts[0][id]<4096)continue;
            const double ratio=metric.ratio?stats.mean[row][metric.element]/stats.mean[row][0]:0;
            std::vector<double> variance(depths.size());double sum=0;size_t best=1;
            for(size_t l=0;l<depths.size();++l){double mean=0;for(int s=0;s<8;++s)mean+=(data.levels[l][s][row][metric.element]-ratio*data.levels[l][s][row][0])/8;
                for(int s=0;s<8;++s)variance[l]+=std::pow(data.levels[l][s][row][metric.element]-ratio*data.levels[l][s][row][0]-mean,2)/7;
                sum+=variance[l];
                if(l&&variance[l]>variance[best])best=l;}
            if(variance[best]>.75*sum&&state.counts[0][id]>=8LL*state.counts[best][id]&&Doubled(e,state.counts[best][id])>state.counts[best][id]){
                const int next=FullAutoScheduling::CorrectionProbeCount(state.counts[best][id],proposed[best][id],e.maxCount);
                candidates.push_back({row,best,metric.error,next});}
        }
        if(candidates.empty())return false;
        std::sort(candidates.begin(),candidates.end(),[](const Candidate& a,const Candidate& b){return a.error>b.error;});
        ++receipt.attempt;receipt.complete=false;receipt.sources=Hash(data.sources);receipt.rows.clear();
        AdditionalWork work(e,e.seeds.data());auto baseRequests=BudgetRequests(e,depths,proposed,ids);
        baseRequests.erase(std::remove_if(baseRequests.begin(),baseRequests.end(),[&](const Engine::WorkRequest& r){return r.depth!=depths[0];}),baseRequests.end());
        const auto base=work.Forecast(baseRequests);WorkForecast probe;bool accepted=false;
        // The former fixed 33-row, 2x probe repeatedly left the next dominant
        // angle on a tiny correction prefix. Batch more affected rows and use
        // the allocator's bounded proposal. The SAME measured 10% work gate
        // limits total cost; it is not relaxed to obtain wider packets.
        size_t width=std::min(candidates.size(),std::min(size_t(128),size_t(e.maxTheta)));
        for(;width; width=width==1?0:(width+1)/2){
            for(bool doubledOnly:{false,true}){
                receipt.rows.clear();
                for(size_t i=0;i<width;++i){const auto& c=candidates[i];const size_t id=ids[c.row];const int n=state.counts[c.level][id];
                    receipt.rows.push_back({id,c.level,n,doubledOnly?Doubled(e,n):c.next});}
                std::sort(receipt.rows.begin(),receipt.rows.end(),[](const CorrectionProbeRow& a,const CorrectionProbeRow& b){return a.id<b.id;});
                probe=work.Forecast(CorrectionProbeRequests(e,depths,receipt));
                accepted=base.profiled&&probe.profiled&&base.jobs&&base.seconds>0&&probe.seconds<=.1*base.seconds;
                if(accepted)break;
            }
            if(accepted)break;
        }
        e.Log("correction_subround_cost_gate","{\"round\":"+Number(state.round)+",\"probe_worker_seconds\":"+Number(probe.seconds)+",\"base_extension_worker_seconds\":"+Number(base.seconds)
            +",\"candidate_rows\":"+Number(candidates.size())+",\"selected_rows\":"+Number(accepted?receipt.rows.size():0)+",\"max_growth\":4"
            +",\"profiled\":"+(base.profiled&&probe.profiled?"true":"false")+",\"selected\":"+(accepted?"true":"false")+",\"scope\":\"scheduling forecast; no accuracy acceptance\"}");
        if(!accepted)return false;
        SaveCorrectionProbe(e,state,receipt);
    }else Require(receipt.sources==Hash(data.sources),"pending correction subround production changed");
    const auto requests=CorrectionProbeRequests(e,depths,receipt);const auto groups=PhiRequests(e,requests);std::ostringstream sources;
    e.State("correction_subround",state.round,PeakBudgets(state.counts),ids.size(),StatisticalError(e,stats),interpolationError,0);
    for(const auto& group:groups){PhiScope scope(e,group.first);e.Prefetch(group.second,e.seeds.data());for(const auto& request:group.second){
        if(e.depthSharing){auto pair=e.SharedStage(request.count,e.seeds.data(),8,request.depth,request.pairedDepth,request.angles);
            for(size_t s=0;s<8;++s)sources<<request.depth<<" 1 "<<request.count<<' '<<e.seeds[s]<<' '<<pair.first[s].folder<<'\n'<<request.pairedDepth<<" -1 "<<request.count<<' '<<e.seeds[s]<<' '<<pair.second[s].folder<<'\n';
        }else{auto jobs=e.Stage(request.count,e.seeds.data(),8,request.depth,request.angles);for(size_t s=0;s<8;++s)sources<<request.depth<<" 1 "<<request.count<<' '<<e.seeds[s]<<' '<<jobs[s].folder<<'\n';}
    }}
    Atomic(CorrectionProbePath(e,state)+".sources",sources.str());receipt.complete=true;SaveCorrectionProbe(e,state,receipt);
    for(const auto& r:receipt.rows)state.counts[r.level][r.id]=r.after;
    state.streak=0;
    SaveAdaptive(e,state);
    e.Log("correction_subround_complete","{\"round\":"+Number(state.round)+",\"attempt\":"+Number(receipt.attempt)+",\"theta_rows\":"+Number(receipt.rows.size())+",\"base_prefix_unchanged\":true,\"all_accuracy_gates_retained\":true}");return true;
}
