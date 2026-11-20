// Bounded, paired common-reference pilots. This code runs before production;
// published adaptive estimators retain their original physical depth on resume.
struct DepthCandidateReview {
    int depth=0,count=0;bool passed=false,rejected=false,referenceFailed=false;
    DepthResidual residual,referenceResidual;
    double forecast=1e100;std::vector<int> counts;
};
DepthCandidateReview ReviewDepthCandidate(Engine& e,int depth,const std::vector<double>& x,int limit,bool tryOtherDepth=false){
    namespace DS=FullAutoDepthSelection;
    const int referenceDepth=DS::Reference(e.depthMaximum),anchorDepth=DS::Anchor(e.depthMaximum);
    unsigned seeds[8];for(int s=0;s<8;++s){seeds[s]=1700000000u+unsigned(e.depthStart)*1024u+unsigned(s)*17u;
        Require(std::find(e.seeds.begin(),e.seeds.end(),seeds[s])==e.seeds.end(),"depth pilot/production seed overlap");}
    DepthCandidateReview result;result.depth=depth;result.counts.assign(x.size(),0);
    Samples candidate(8,Matrix(x.size())),reference=candidate,anchor=candidate;
    std::vector<size_t> active;for(size_t i=0;i<x.size();++i)active.push_back(i);
    std::vector<DepthResidual> bounds(x.size()),guardBounds(x.size());
    double workerSeconds=0,workerUnits=0;
    for(int count=e.depthPilot,attempt=0;;count=DS::DoubleCount(count,limit),++attempt){
        e.target=depth;e.State("selecting_depth",attempt,{count},active.size(),0,0,0);
        e.Record(e.root+"/depth_search_progress.json","{\"policy\":"+Quoted(DS::Policy())+",\"candidate\":"+Number(depth)
            +",\"comparison_depth\":"+Number(anchorDepth)+",\"reference_depth\":"+Number(referenceDepth)
            +",\"count_per_seed\":"+Number(count)+",\"active_theta_rows\":"+Number(active.size())+",\"attempt\":"+Number(attempt)
            +",\"status\":\"sampling paired common-reference matrices\"}");
        std::vector<size_t> sampled=active;
        // The diffraction worker needs two distinct angles. The padding row is
        // not used to replace a previously accepted row or its sample count.
        if(!e.go&&sampled.size()==1&&x.size()>1){sampled.push_back(sampled.front()?0:1);std::sort(sampled.begin(),sampled.end());}
        if(e.go){sampled.clear();for(size_t i=0;i<x.size();++i)sampled.push_back(i);}
        std::vector<double> query;for(size_t i:sampled)query.push_back(x[i]);
        Samples a,b,c;double seconds=0;
        if(e.go){a=GoSample(e,count,seeds,depth);b=GoSample(e,count,seeds,referenceDepth);c=GoSample(e,count,seeds,anchorDepth);}
        else{
            // Reference prefixes are shared across every candidate. Pairing
            // their trace is optional; pairing their orientations is mandatory.
            if(e.depthSharing&&e.depthVerified){const auto pair=e.SharedStage(count,seeds,8,anchorDepth,referenceDepth,query);b=e.Data(pair.second);c=e.Data(pair.first);}
            else{b=e.Data(e.PrefixStage(count,seeds,8,referenceDepth,query,true));c=e.Data(e.PrefixStage(count,seeds,8,anchorDepth,query,true));}
            if(depth==referenceDepth)a=b;
            else {const auto jobs=e.PrefixStage(count,seeds,8,depth,query,true);a=e.Data(jobs);for(const auto& j:jobs)seconds+=std::max(0.,j.seconds-j.controlSetup);}
        }
        workerSeconds+=seconds;workerUnits+=double(count)*8*query.size();
        const auto scale=Mean(c);std::vector<DepthResidual> nextBounds,nextGuard;
        DepthDifference(a,c,e.go,e.theta.size(),e.ratios,&scale,&nextBounds);
        DepthDifference(b,c,e.go,e.theta.size(),e.ratios,&scale,&nextGuard);
        for(size_t r=0;r<sampled.size();++r){const auto row=sampled[r];
            if(!e.go&&std::find(active.begin(),active.end(),row)==active.end())continue;
            bounds[row]=nextBounds[r];guardBounds[row]=nextGuard[r];result.counts[row]=count;
            for(int s=0;s<8;++s){candidate[s][row]=a[s][r];reference[s][row]=b[s][r];anchor[s][row]=c[s][r];}
        }
        result.count=*std::max_element(result.counts.begin(),result.counts.end());
        const auto completeScale=Mean(anchor);
        result.residual=DepthDifference(candidate,anchor,e.go,e.theta.size(),e.ratios,&completeScale);
        result.referenceResidual=DepthDifference(reference,anchor,e.go,e.theta.size(),e.ratios,&completeScale);
        active.clear();result.rejected=false;result.referenceFailed=false;
        std::ostringstream rows;rows<<"theta_deg,count_per_seed,residual_upper95,residual_lower95,reference_upper95,reference_lower95,decision\n"<<std::setprecision(17);
        for(size_t r=0;r<x.size();++r){const auto decision=DS::Classify(bounds[r].upper,bounds[r].lower,guardBounds[r].upper,guardBounds[r].lower,e.depthBudget);
            if(decision==DS::Unresolved)active.push_back(r);
            result.rejected=result.rejected||decision==DS::Fail;result.referenceFailed=result.referenceFailed||decision==DS::ReferenceFail;
            rows<<x[r]<<','<<result.counts[r]<<','<<bounds[r].upper<<','<<bounds[r].lower<<','<<guardBounds[r].upper<<','<<guardBounds[r].lower<<','
                <<(decision==DS::Pass?"pass":decision==DS::Fail?"depth_insufficient":decision==DS::ReferenceFail?"reference_unresolved":"sampling_unresolved")<<'\n';
        }
        result.passed=active.empty()&&!result.rejected&&!result.referenceFailed;
        const std::string tag="n"+Number(depth)+"_reference"+Number(anchorDepth)+"_N"+Number(count);
        const std::string snapshot=e.root+"/depth_pilot_snapshots/"+tag;Directory(snapshot);
        if(e.go){GoOutputs(e,snapshot,EstimateGO(candidate,e.theta.size(),e.ratios),count);}
        else{
            WriteMatrix(snapshot+"/mueller_evaluated.dat",x,Mean(candidate));
            WriteMatrix(snapshot+"/mueller_reference.dat",x,Mean(reference));WriteMatrix(snapshot+"/mueller_anchor.dat",x,Mean(anchor));
            SaveReplicaMatrices(e,"depth_pilot_"+tag,candidate,x,seeds);
            SaveReplicaMatrices(e,"depth_reference_"+tag,reference,x,seeds);SaveReplicaMatrices(e,"depth_anchor_"+tag,anchor,x,seeds);
        }
        // GO also exports a generic statistics CSV. Preserve the paired depth
        // decisions at the path used by the pilot receipt and live pointer.
        Atomic(snapshot+"/parameters_accuracy.csv",rows.str());
        e.Record(e.root+"/intermediate_current.json","{\"converged\":false,\"phase\":\"depth_selection\",\"policy\":"+Quoted(DS::Policy())
            +",\"candidate\":"+Number(depth)+",\"snapshot_folder\":"+Quoted(snapshot)+",\"matrix\":"+Quoted(snapshot+"/mueller_evaluated.dat")
            +",\"accuracy_and_parameters\":"+Quoted(snapshot+"/parameters_accuracy.csv")+",\"paired_depth_residual95\":"+Number(result.residual.upper)+"}");
        const bool defer=tryOtherDepth && !result.rejected && !result.referenceFailed &&
            DS::DeferSparsePilot(count,limit,active.size(),result.referenceResidual.upper,e.depthBudget);
        const bool stop=defer||result.passed||result.rejected||result.referenceFailed||count>=limit;
        const std::string action=result.passed?"accept_candidate":result.referenceFailed?"reference_unresolved":result.rejected?"depth_insufficient":defer?"try_next_depth":stop?"defer_unresolved_candidate":"double_unresolved_angles";
        std::ostringstream report;report<<"{\"policy\":"<<Quoted(DS::Policy())<<",\"automatic\":true,\"candidate\":"<<depth<<",\"comparison_depth\":"<<anchorDepth
            <<",\"reference_depth\":"<<referenceDepth<<",\"count_per_seed\":"<<count<<",\"count_scope\":\"maximum; paired per-angle counts in CSV\",\"theta_rows\":"<<x.size()
            <<",\"sampled_theta_rows\":"<<query.size()<<",\"unresolved_theta_rows\":"<<active.size()<<",\"phi\":"<<e.phi
            <<",\"mean_residual_scaled_M11\":"<<Number(result.residual.mean)<<",\"interval95_scaled_M11\":"<<Number(result.residual.interval)
            <<",\"residual_upper95\":"<<Number(result.residual.upper)<<",\"residual_lower95\":"<<Number(result.residual.lower)
            <<",\"reference_upper95\":"<<Number(result.referenceResidual.upper)<<",\"reference_budget\":"<<Number(e.depthBudget)
            <<",\"budget\":"<<Number(e.depthBudget)<<",\"passed\":"<<(result.passed?"true":"false")<<",\"action\":"<<Quoted(action)
            <<",\"row_parameters\":"<<Quoted(snapshot+"/parameters_accuracy.csv")<<",\"scope\":\"paired finite-depth pilot; not an infinite-tail bound or global confidence guarantee\"}";
        e.Record(e.root+"/depth_selection_candidate"+Number(depth)+"_reference"+Number(anchorDepth)+"_N"+Number(count)+".json",report.str());
        e.Log("depth_candidate",report.str());e.Log("depth_pilot_decision",report.str());
        if(stop)break;
    }
    if(result.passed){
        // A scheduling proxy, never an acceptance criterion. Use each row's
        // actual prefix length. Phi/control tuning and MLMC costs follow later.
        if(!e.go&&workerSeconds>0&&workerUnits>0){const auto stats=Estimate(candidate);double predicted=0;
            const double budget=std::max(e.eps/4,e.eps-e.depthBudget-e.interpolation);
            for(size_t r=0;r<x.size();++r){double error=0;for(int j=0;j<16;++j){error=std::max(error,stats.interval[r][j]);if(e.ratios)error=std::max(error,stats.ratioInterval[r][j]);}
                predicted+=result.counts[r]*std::max(1.,std::pow(error/budget,2));}
            if(DS::FiniteBound(predicted))result.forecast=workerSeconds/workerUnits*8*predicted;
        }
        e.Log("depth_candidate_cost","{\"candidate\":"+Number(depth)+",\"passed\":true,\"forecast_worker_seconds\":"+(DS::FiniteBound(result.forecast)?Number(result.forecast):"null")
            +",\"scope\":\"pilot variance and measured candidate worker cost; excludes later phi, controls and MLMC tuning\"}");
    }
    return result;
}
bool SelectAutomaticDepth(Engine& e){
    namespace DS=FullAutoDepthSelection;
    if(!e.autoDepth){e.Record(e.root+"/depth_selection.json","{\"automatic\":false,\"target_depth\":"+Number(e.target)+",\"explicit_target_preserved\":true}");return true;}
    const std::string checkpoint=e.root+"/automatic_depth.tsv",receipt=e.root+"/automatic_depth_policy.tsv";
    ResumeDepthPromotion(e);
    int searchStart=DepthPromotionFloor(e);
    if(Exists(checkpoint)){
        const auto text=Read(checkpoint);const auto split=text.find('\n');Require(split!=std::string::npos,"invalid depth checkpoint");
        const auto payload=text.substr(split+1);Require(Hash(payload)==text.substr(0,split),"depth checkpoint changed");
        std::istringstream input(payload);std::string signature,extra;int selected;
        Require(bool(input>>signature>>selected)&&!(input>>extra)&&signature==e.signature&&selected>=e.depthStart&&selected<=e.depthMaximum,"depth checkpoint inputs changed");
        const auto report=Exists(e.root+"/depth_selection.json")?CompactJson(Read(e.root+"/depth_selection.json")):std::string();
        const bool sameDepth=report.find("\"target_depth\":"+Number(selected)+",")!=std::string::npos||report.find("\"target_depth\":"+Number(selected)+"}")!=std::string::npos;
        const bool passed=selected>=searchStart&&sameDepth&&report.find("\"pilot_passed\":true")!=std::string::npos;
        const bool newPolicy=report.find(std::string("\"policy\":\"")+DS::Policy()+"\"")!=std::string::npos;
        const bool production=Exists(e.root+"/adaptive_checkpoint.tsv")||Exists(e.root+"/go_checkpoint.tsv")||Exists(e.root+"/final_outputs.fnv");
        const bool verified=newPolicy&&Exists(receipt)&&Read(receipt)==e.signature+" "+DS::Policy()+" "+Number(selected)+" "+Hash(Read(e.root+"/depth_selection.json"))+"\n";
        const bool fallback=selected==e.depthMaximum&&sameDepth
            &&report.find("\"fallback_maximum\":true")!=std::string::npos
            &&report.find("\"pilot_passed\":false")!=std::string::npos;
        if(((passed||fallback)&&verified)||(passed&&!newPolicy&&production)){
            e.depthFallback=fallback&&verified;
            e.target=selected;e.SetMirror(false);e.WarmCache();
            e.Record(e.root+"/depth_selection_resume_policy.json","{\"selected_depth\":"+Number(selected)+",\"criterion\":"+Quoted(verified?DS::Policy():"legacy-neighbor-pilot-retained")
                +",\"production_depth_preserved\":true,\"common_reference_verified\":"+((verified&&!e.depthFallback)?"true":"false")+e.DepthAccuracyFields()+"}");return true;
        }
        Require(!production,"depth receipt is missing or unverified; cannot change the depth of an existing production checkpoint");
        // A failed final audit has already ruled out the previous target.
        // Do not reselect it from its smaller cached startup pilot.
        if(report.find("promoted after fresh final depth validation")!=std::string::npos)searchStart=std::max(searchStart,selected);
        e.Log("depth_selection_resume_unresolved","{\"policy\":"+Quoted(DS::Policy())+",\"cached_worker_prefixes_retained\":true,\"production_blocked_until_pass\":true}");
    }
    std::vector<double> x;if(e.go)x=GoRows(e);
    else{const size_t size=std::min(size_t(65),std::min(e.theta.size(),size_t(e.maxTheta)));for(size_t row=0;row<size;++row)x.push_back(e.theta[row*(e.theta.size()-1)/(size-1)]);}
    WriteGrid(e.root+"/depth_pilot_theta.csv",x);
    const int limit=DS::PilotLimit(e.depthPilot,e.maxCount);
    // A marginal shallow depth can spend hours on one more doubling. Try the
    // other candidate depths first; the same bounds still decide acceptance.
    const bool joinedRosette=!e.go && e.args.IsCatched("rosette_joined");
    const int firstPass=DS::PilotFloor(e.depthPilot,e.maxCount);
    const int savedBlock=e.orientationBlock;
    // Doubling must append native prefixes, not retrace a larger monolithic
    // pilot each time. This is a startup-only default; explicit blocks win.
    if(e.incremental&&!e.args.IsCatched("fullauto_orientation_block")&&!e.config.adaptive.loadedFromFile)
        e.orientationBlock=std::min(e.orientationBlock,std::max(512,e.depthPilot));
    e.Record(e.root+"/depth_selection_policy.json","{\"policy\":"+Quoted(DS::Policy())+",\"reference_depth\":"+Number(DS::Reference(e.depthMaximum))
        +",\"comparison_depth\":"+Number(DS::Anchor(e.depthMaximum))+",\"uncertain_pilot_floor\":"+Number(firstPass)+",\"pilot_limit\":"+Number(limit)
        +",\"pilot_orientation_block\":"+Number(e.orientationBlock)+",\"depth_budget\":"+Number(e.depthBudget)+",\"reference_guard_budget\":"+Number(e.depthBudget)
        +",\"search_start\":"+Number(searchStart)+",\"search\":\"doubling coarse depths then all skipped smaller depths; no monotonicity assumption\",\"local_theta_growth\":"+(e.go?"false":"true")
        +",\"scope\":\"finite reference agreement; final checks remain a separate explicit policy\"}");
    const bool savedControls=e.controls;e.controls=false;
    std::map<int,DepthCandidateReview> reviews;DepthCandidateReview chosen;chosen.depth=e.depthStart;int smallestPass=e.depthMaximum+1;bool referenceFailed=false;double failedReferenceUpper=0;
    auto remember=[&](const DepthCandidateReview& r){reviews[r.depth]=r;referenceFailed=referenceFailed||r.referenceFailed;
        if(r.referenceFailed)failedReferenceUpper=std::max(failedReferenceUpper,r.referenceResidual.upper);
        if(r.passed){if(smallestPass>e.depthMaximum)smallestPass=r.depth;
            if(!chosen.passed||r.forecast<chosen.forecast||(r.forecast==chosen.forecast&&r.depth<chosen.depth))chosen=r;}
        else if(!chosen.passed)chosen=r;};
    try{
        if(!e.go&&e.depthSharing){e.target=DS::Anchor(e.depthMaximum);e.VerifyDepthSharing(x);}
        for(int depth:DS::Candidates(searchStart,e.depthMaximum)){
            if(referenceFailed)break;
            if(depth>smallestPass)continue;
            remember(ReviewDepthCandidate(e,depth,x,firstPass,joinedRosette && depth<e.depthMaximum));
        }
        if(!chosen.passed&&!referenceFailed){
            // Only unresolved candidates need larger pilots; confirmed failures
            // stay rejected. Review replays cached prefixes before extending the
            // remaining angles, including after an interrupted controller.
            for(int depth:DS::Candidates(searchStart,e.depthMaximum)){
                const auto old=reviews.find(depth);if(old==reviews.end()||old->second.rejected||old->second.referenceFailed||old->second.count>=limit)continue;
                if(depth>smallestPass)continue;
                remember(ReviewDepthCandidate(e,depth,x,limit));if(referenceFailed)break;
            }
        }
    }catch(...){e.orientationBlock=savedBlock;e.controls=savedControls;if(savedControls)e.depthVerified=false;throw;}
    e.orientationBlock=savedBlock;e.controls=savedControls;if(savedControls)e.depthVerified=false;
    // A subsequently resolved failure of the common reference invalidates an
    // earlier small-prefix pass. Never publish a successful receipt beside it.
    if(referenceFailed)chosen.passed=false;
    e.depthFallback=!chosen.passed;
    const int selected=DS::ProductionDepth(chosen.passed,chosen.depth,e.depthMaximum);
    SaveSelectedDepth(e,selected);
    std::ostringstream decisions;decisions<<"depth,passed,depth_insufficient,reference_failed,max_count_per_seed,residual_upper95,reference_upper95,forecast_worker_seconds,selected\n";
    for(const auto& item:reviews){const auto& r=item.second;decisions<<r.depth<<','<<r.passed<<','<<r.rejected<<','<<r.referenceFailed<<','<<r.count<<','<<Number(r.residual.upper)<<','<<Number(r.referenceResidual.upper)<<','
        <<(DS::FiniteBound(r.forecast)?Number(r.forecast):"")<<','<<(r.depth==selected)<<'\n';}
    Atomic(e.root+"/depth_candidates.csv",decisions.str());
    const std::string chosenSnapshot=e.root+"/depth_pilot_snapshots/n"+Number(chosen.depth)+"_reference"+Number(DS::Anchor(e.depthMaximum))+"_N"+Number(chosen.count);
    const std::string report="{\"policy\":"+Quoted(DS::Policy())+",\"automatic\":true,\"target_depth\":"+Number(selected)
        +",\"fallback_maximum\":"+(e.depthFallback?"true":"false")+",\"pilot_candidate_depth\":"+Number(chosen.depth)+e.DepthAccuracyFields()
        +",\"start\":"+Number(e.depthStart)+",\"maximum\":"+Number(e.depthMaximum)+",\"pilot_passed\":"+(chosen.passed?"true":"false")
        +",\"count_per_seed\":"+Number(chosen.count)+",\"count_scope\":\"maximum of local paired pilots\",\"residual_upper95\":"+Number(chosen.residual.upper)
        +",\"budget\":"+Number(e.depthBudget)+",\"reference_depth\":"+Number(DS::Reference(e.depthMaximum))+",\"comparison_depth\":"+Number(DS::Anchor(e.depthMaximum))
        +",\"reference_upper95\":"+Number(chosen.referenceResidual.upper)+",\"reference_budget\":"+Number(e.depthBudget)
        +",\"reference_failed\":"+(referenceFailed?"true":"false")+",\"failed_reference_upper95\":"+Number(failedReferenceUpper)+",\"selection\":"+Quoted(e.depthFallback?"requested maximum; pilot did not pass":"lowest forecast cost among tested passing candidates")
        +",\"row_parameters\":"+Quoted(chosenSnapshot+"/parameters_accuracy.csv")+",\"candidate_comparison\":"+Quoted(e.root+"/depth_candidates.csv")
        +",\"forecast_worker_seconds\":"+(DS::FiniteBound(chosen.forecast)?Number(chosen.forecast):"null")+",\"fresh_final_audit_required\":true}";
    e.Record(e.root+"/depth_selection.json",report);
    e.Record(e.root+"/intermediate_current.json","{\"converged\":false,\"phase\":\"depth_selection\",\"candidate\":"+Number(chosen.depth)
        +",\"snapshot_folder\":"+Quoted(chosenSnapshot)+",\"matrix\":"+Quoted(chosenSnapshot+"/mueller_evaluated.dat")
        +",\"accuracy_and_parameters\":"+Quoted(chosenSnapshot+"/parameters_accuracy.csv")+",\"paired_depth_residual95\":"+Number(chosen.residual.upper)+"}");
    Atomic(receipt,e.signature+" "+DS::Policy()+" "+Number(selected)+" "+Hash(Read(e.root+"/depth_selection.json"))+"\n");
    if(e.depthFallback){
        e.Log("depth_selection_fallback",report);
        e.State("depth_selection_fallback_maximum",0,{chosen.count},x.size(),0,0,0);
    }
    return true;
}
