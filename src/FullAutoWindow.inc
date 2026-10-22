// A scoped continuation has its own controller state and output. Native
// sampling identities stay unchanged, so complete old prefixes remain usable.
// The parent controller lock is held throughout import and execution.
namespace {
std::string RefinementToken(double value){
    std::ostringstream out;out<<std::setprecision(12)<<value;std::string token=out.str();
    for(char& c:token)if(c=='.')c='p';else if(c=='-')c='m';else if(c=='+')c='a';
    return token;
}
std::string RefinementChild(const Engine& e){
    return e.root+"/refine_theta"+RefinementToken(e.refinementLo)+"_"+RefinementToken(e.refinementHi)
        +"_eps"+RefinementToken(e.refinementEps);
}
std::vector<double> RefinementGrid(const std::vector<double>& original,double low,double high){
    Require(original.size()>=2&&low>=original.front()&&high<=original.back()&&low<high,
        "refinement window must lie within the saved requested grid");
    std::vector<double> result{low};
    for(double angle:original)if(angle>low&&angle<high)result.push_back(angle);
    result.push_back(high);return result;
}
AdaptiveState ProjectRefinementWindow(Engine& e,const AdaptiveState& original,const std::vector<double>& parentTheta){
    Require(original.previous.size()==parentTheta.size()&&original.valid.size()==parentTheta.size()
        &&original.coarse.size()==parentTheta.size()&&!original.counts.empty(),"invalid parent state for refinement");
    for(const auto& level:original.counts)Require(level.size()==parentTheta.size(),"invalid parent allocation for refinement");
    const auto oldPhi=e.thetaPhi;const auto grid=RefinementGrid(parentTheta,e.refinementLo,e.refinementHi);
    AdaptiveState result;result.round=0;result.streak=0;result.auditSerial=original.auditSerial;
    result.counts=Budgets(original.counts.size(),std::vector<int>(grid.size()));
    result.previous=Matrix(grid.size());result.valid.assign(grid.size(),false);result.coarse.assign(grid.size(),false);
    std::set<size_t> knots{0,grid.size()-1};e.thetaPhi.clear();
    auto phiAt=[&](size_t id){const auto found=oldPhi.find(id);return found==oldPhi.end()?e.phi:found->second;};
    for(size_t row=0;row<grid.size();++row){
        const auto upper=std::lower_bound(parentTheta.begin(),parentTheta.end(),grid[row]);
        const size_t hi=size_t(upper-parentTheta.begin());
        if(upper!=parentTheta.end()&&*upper==grid[row]){
            for(size_t level=0;level<result.counts.size();++level)result.counts[level][row]=original.counts[level][hi];
            result.previous[row]=original.previous[hi];result.valid[row]=original.valid[hi];
            if(std::binary_search(original.knots.begin(),original.knots.end(),hi))knots.insert(row);
            if(phiAt(hi)!=e.phi)e.thetaPhi[row]=phiAt(hi);
        }else{
            Require(hi>0&&hi<parentTheta.size(),"new refinement endpoint lacks neighboring saved angles");
            for(size_t level=0;level<result.counts.size();++level)
                result.counts[level][row]=std::max(original.counts[level][hi-1],original.counts[level][hi]);
            // This is a new estimator, never a mixture of adjacent phi rules.
            const int phi=std::max(phiAt(hi-1),phiAt(hi));if(phi!=e.phi)e.thetaPhi[row]=phi;
        }
    }
    result.knots.assign(knots.begin(),knots.end());e.theta=grid;return result;
}
void CopyRefinementFile(const std::string& source,const std::string& destination){
    // Metadata/cache copies must not share an inode: both sides may later be
    // atomically replaced and a native mean-cache writer must own its file.
    Require(source!=destination,"refinement copy must have a distinct destination");
    std::ifstream input(source.c_str(),std::ios::binary);Require(bool(input),"cannot open refinement source "+source);
    const auto temporary=destination+".import.tmp";
    std::ofstream output(temporary.c_str(),std::ios::binary|std::ios::trunc);Require(bool(output),"cannot copy refinement metadata");
    char buffer[65536];
    while(input){input.read(buffer,sizeof(buffer));const auto size=input.gcount();if(size)output.write(buffer,size);}
    Require(input.eof()&&!input.bad(),"cannot read refinement source "+source);
    output.close();Require(bool(output),"refinement copy write failed");
    int fd=open(temporary.c_str(),O_RDONLY);Require(fd>=0,"cannot flush refinement copy");
    const bool flushed=fsync(fd)==0;close(fd);Require(flushed,"refinement copy flush failed");
    Require(rename(temporary.c_str(),destination.c_str())==0,"cannot publish refinement copy");
    const auto slash=destination.find_last_of('/');const auto directory=slash==std::string::npos?".":destination.substr(0,slash);
    fd=open(directory.c_str(),O_RDONLY|O_DIRECTORY);if(fd>=0){const bool okay=fsync(fd)==0;close(fd);Require(okay,"refinement directory flush failed");}
}
int RefinementSelectedDepth(Engine& e){
    if(!e.autoDepth)return e.target;
    const auto path=e.root+"/automatic_depth.tsv";Require(Exists(path),"refinement requires completed automatic depth selection");
    const auto text=Read(path);const auto split=text.find('\n');Require(split!=std::string::npos,"invalid parent depth checkpoint");
    const auto payload=text.substr(split+1);Require(Hash(payload)==text.substr(0,split),"parent depth checkpoint changed");
    std::istringstream input(payload);std::string signature,extra;int depth=0;
    Require(bool(input>>signature>>depth)&&!(input>>extra)&&signature==e.signature&&depth>=e.depthStart&&depth<=e.depthMaximum,
        "parent depth selection is incompatible");
    Require(Exists(e.root+"/depth_selection.json")&&CompactJson(Read(e.root+"/depth_selection.json")).find("\"pilot_passed\":true")!=std::string::npos,
        "refinement requires a passed parent depth pilot");return depth;
}
void PrepareRefinementWindow(Engine& e){
    Require(e.refinementWindow&&!e.go,"scoped continuation requires physical optics");
    Require(e.refinementEps>=e.eps,"scoped continuation cannot tighten the parent physical-depth tolerance; use EPS >= the original fullauto tolerance");
    const std::string parent=e.root,child=RefinementChild(e);
    Require(Exists(parent+"/configuration.txt")&&Read(parent+"/configuration.txt")==e.signature+"\n",
        "refinement parent inputs/binary/environment changed; retain the original argv and compatible origin binary");
    Require(Exists(parent+"/adaptive_checkpoint.tsv")&&Exists(parent+"/requested_theta.csv")
        &&Exists(parent+"/execution_tuning.tsv"),"refinement requires an existing adaptive PO checkpoint and saved tuning");
    Require(!Exists(parent+"/refinement_window_policy.json"),"nested refinement windows are not supported");
    const auto parentTheta=ReadGrid(parent+"/requested_theta.csv");
    Require(parentTheta==e.theta,"parent requested angles no longer match the original command");
    const auto checkpointHash=FileHash(parent+"/adaptive_checkpoint.tsv");
    const auto tuningHash=FileHash(parent+"/execution_tuning.tsv");
    e.target=RefinementSelectedDepth(e);Require(e.RestoreTuning(),"parent execution tuning is missing");
    std::vector<int> depths;for(int depth:{std::max(1,e.target-4),std::max(1,e.target-2),e.target})
        if(depths.empty()||depths.back()!=depth)depths.push_back(depth);
    AdaptiveState original;original.counts=Budgets(depths.size(),std::vector<int>(parentTheta.size(),e.minimum));
    original.counts[0].assign(parentTheta.size(),e.initial);original.previous=Matrix(parentTheta.size());
    original.valid.assign(parentTheta.size(),false);original.coarse.assign(parentTheta.size(),false);
    Require(RestoreAdaptive(e,original,false),"parent adaptive state is unavailable");
    e.completedCheckpoint=false;
    Require(!e.controls||e.sharedMeans,"scoped continuation requires authoritative shared analytic means when controls are enabled");
    const auto meanSource=e.requestedMeanCache.empty()?parent+"/means_shared.cache":e.requestedMeanCache;
    if(e.controls){
        Require(Exists(meanSource),"parent authoritative analytic mean cache is missing");
        std::ifstream cache(meanSource.c_str());std::string version;cache>>version;
        Require(version=="MBS_PHYSICAL_FACET_MEAN_V3_COMPACT_ROWS","scoped continuation requires the V3 per-angle analytic mean cache");
    }
    const auto projected=ProjectRefinementWindow(e,original,parentTheta);
    std::vector<double> reference=parentTheta;reference.push_back(e.refinementLo);reference.push_back(e.refinementHi);
    std::sort(reference.begin(),reference.end());reference.erase(std::unique(reference.begin(),reference.end()),reference.end());
    std::ostringstream grid;for(double angle:e.theta)grid<<Number(angle)<<'\n';
    const double originalEps=e.eps;
    std::ostringstream policy;policy<<"{\"schema\":1,\"policy\":\"isolated-angular-continuation-v1\",\"parent_output\":"<<Quoted(parent)
        <<",\"parent_configuration_fnv1a64\":"<<Quoted(e.signature)<<",\"parent_checkpoint_fnv1a64\":"<<Quoted(checkpointHash)
        <<",\"parent_execution_tuning_fnv1a64\":"<<Quoted(tuningHash)
        <<",\"parent_analytic_mean_cache_fnv1a64\":"<<(e.controls?Quoted(FileHash(meanSource)):"null")
        <<",\"parent_theta_fnv1a64\":"<<Quoted(FileHash(parent+"/requested_theta.csv"))
        <<",\"parent_round\":"<<original.round<<",\"parent_target_epsilon\":"<<Number(originalEps)
        <<",\"theta_min\":"<<Number(e.refinementLo)<<",\"theta_max\":"<<Number(e.refinementHi)
        <<",\"target_epsilon\":"<<Number(e.refinementEps)<<",\"interpolation_budget\":"<<Number(e.refinementEps*.25)
        <<",\"requested_theta_fnv1a64\":"<<Quoted(Hash(grid.str()))<<",\"requested_theta_rows\":"<<e.theta.size()
        <<",\"selected_depth\":"<<e.target<<",\"child_round_starts_at\":0,\"orientation_counts_inherited\":true"
        <<",\"analytic_reference\":\"unchanged parent angles plus exact new endpoints; frozen parent coefficients\""
        <<",\"outside_window_revalidated\":false,\"parent_outputs_unchanged\":true,\"parent_jobs\":\"read-only parent cache; new jobs local to child\"}";
    const std::string policyText=policy.str();e.refinementPolicy=Hash(policyText);e.refinementParent=parent;
    Directory(child);const auto policyPath=child+"/refinement_window_policy.json";
    if(Exists(policyPath))Require(Read(policyPath)==policyText,"refinement scope or parent checkpoint changed; preserve the frozen parent");
    else Atomic(policyPath,policyText);
    e.root=child;e.meanReferenceOverride=child+"/analytic_reference_theta.csv";
    e.requestedMeanCache=child+"/means_shared.cache";
    e.eps=e.refinementEps;e.interpolation=e.eps*.25;e.depthBudget=e.autoDepth?e.eps/4:0;e.SetMirror(e.mirror);
    const auto complete=child+"/refinement_window_import_complete.json";
    if(Exists(complete)){
        Require(Read(complete)=="{\"policy_fnv1a64\":"+Quoted(e.refinementPolicy)+"}","refinement import policy changed");
        Require(Read(child+"/configuration.txt")==e.signature+"\n"&&ReadGrid(child+"/requested_theta.csv")==e.theta,
            "refinement output grid or sampling identity changed");
        Require(ReadGrid(e.meanReferenceOverride)==reference,"refinement analytic reference grid changed");
        return;
    }
    // This commit marker is written only after all projections and local metadata
    // exist. An interrupted import can safely retry before any worker runs.
    Directory(child+"/parent_snapshot");
    for(const char* file:{"configuration.txt","adaptive_checkpoint.tsv","requested_theta.csv","fullauto_parameters.json",
            "fullauto_status.json","controller_command.json","intermediate_current.json","orientation_budget.csv",
            "working_phi_n6.tsv","working_phi_n8.tsv","working_phi_n10.tsv","working_phi_n12.tsv",
            "working_phi_n14.tsv","working_phi_n16.tsv","working_phi_n18.tsv"})
        if(Exists(parent+"/"+file))CopyRefinementFile(parent+"/"+file,child+"/parent_snapshot/"+file);
    for(const char* file:{"configuration.txt","job_index_format.txt","completed_jobs.journal","execution_tuning.tsv",
            "automatic_depth.tsv","depth_selection.json","depth_startup.tsv","depth_startup.json","fullauto_parameters.json",
            "calibration.json","initial_calibration.json","polarization_weights.tsv","resource_calibration.json",
            "rotation_verification.json","cutoff_calibration.json","automatic_cutoff_validation.json","mirror_verification.json",
            "trace_backend_calibration.json","analytic_control_fallback.json"})
        if(Exists(parent+"/"+file))CopyRefinementFile(parent+"/"+file,child+"/"+file);
    if(e.controls)CopyRefinementFile(meanSource,e.requestedMeanCache);
    Directory(child+"/jobs");
    WriteGrid(e.meanReferenceOverride,reference);WriteGrid(child+"/requested_theta.csv",e.theta);
    SaveAdaptive(e,projected);
    e.Record(child+"/refinement_window_import.json","{\"parent_round\":"+Number(original.round)
        +",\"child_round\":0,\"orientation_counts_inherited\":true,\"parent_phi_rules_projected\":true"
        +",\"parent_analytic_coefficients_retained\":true,\"old_acceptance_reset\":true,\"old_audit_namespace_retained\":true"
        +",\"old_working_phi_plan\":\"archived; projected checkpoint phi rules retained\",\"parent_checkpoint_fnv1a64\":"+Quoted(checkpointHash)+"}");
    Require(FileHash(parent+"/adaptive_checkpoint.tsv")==checkpointHash,"parent checkpoint changed during import");
    Atomic(complete,"{\"policy_fnv1a64\":"+Quoted(e.refinementPolicy)+"}");
    e.Log("refinement_window_imported",policyText);
}
}
