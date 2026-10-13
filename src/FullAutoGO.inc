// GO integrates energy over finite angular bins; do not spline delta rays.
namespace {
struct GoStatistics {Statistics stats;std::vector<bool> checked;double maximum=0;};
GoStatistics EstimateGO(const Samples& samples,size_t bins,bool ratios){
    GoStatistics out;out.stats.mean=Mean(samples);const size_t rows=out.stats.mean.size();out.stats.interval.resize(rows);out.stats.ratioInterval.resize(rows);out.checked.assign(rows,false);
    for(size_t r=0;r<rows;++r){
        if(!(out.stats.mean[r][0]>0)){if(r<bins)out.maximum=1e100;continue;}
        Samples single(8,Matrix(1));for(int i=0;i<8;++i)single[i][0]=samples[i][r];auto result=Estimate(single);
        out.checked[r]=true;out.stats.interval[r]=result.interval[0];out.stats.ratioInterval[r]=result.ratioInterval[0];
        out.stats.maxInterval=std::max(out.stats.maxInterval,result.maxInterval);out.stats.maxRatioInterval=std::max(out.stats.maxRatioInterval,result.maxRatioInterval);
        out.maximum=std::max(out.maximum,ratios?std::max(result.maxInterval,result.maxRatioInterval):result.maxInterval);
    }return out;
}
std::vector<double> GoRows(const Engine& e){auto x=e.theta;x.push_back(0);x.push_back(180);return x;}
Samples GoSample(Engine& e,int count,const unsigned* seeds,int depth=-1){
    if(depth<0)depth=e.target;
    const auto x=GoRows(e);
    Samples data(8,Matrix(x.size()));
    // Dyadic contiguous windows preserve every completed prefix without
    // forcing all future work into the small initial sample size.
    int begin=0;
    while(begin<count){
        std::vector<std::pair<int,int>> windows;int end=begin;
        for(int k=0;k<8&&end<count;++k){
            const int block=e.incremental?std::min(e.orientationBlock,end?end:e.initial):count;
            const int n=std::min(block,count-end);windows.push_back({end,n});end+=n;
        }
        e.collecting=e.QueueEnabled();
        for(const auto& w:windows)e.DirectStage(w.second,seeds,8,depth,x,true,false,w.first);
        e.collecting=false;
        if(!e.pendingJobs.empty()){
            std::vector<size_t> todo(e.pendingJobs.size());for(size_t i=0;i<todo.size();++i)todo[i]=i;
            e.ExecutePrepared(e.pendingJobs,e.pendingCommands,todo);e.pendingJobs.clear();e.pendingCommands.clear();
        }
        for(const auto& w:windows){auto jobs=e.DirectStage(w.second,seeds,8,depth,x,true,false,w.first);
            for(int seed=0;seed<8;++seed)for(size_t row=0;row<x.size();++row)for(int j=0;j<16;++j)
                data[seed][row][j]+=double(w.second)/count*jobs[seed].matrix[row][j];
        }
        begin=end;e.Record(e.root+"/orientation_progress.json","{\"completed_prefix\":"+Number(begin)+",\"requested_prefix\":"+Number(count)+",\"depth\":"+Number(depth)+",\"seeds\":8,\"scope\":\"GO exact Sobol windows\"}");
    }return data;
}
double GoDifference(const Matrix& a,const Matrix& b,size_t bins,bool ratios){
    double error=0;for(size_t r=0;r<a.size();++r){
        if(!(a[r][0]>0)){if(r<bins||b[r][0]>0)error=1e100;continue;}
        for(int j=0;j<16;++j){error=std::max(error,std::fabs(a[r][j]-b[r][j])/a[r][0]);
            if(ratios&&j&&b[r][0]>0)error=std::max(error,std::fabs(a[r][j]/a[r][0]-b[r][j]/b[r][0]));}
    }return error;
}
void GoOutputs(Engine& e,const std::string& folder,const GoStatistics& result,int count){
    Directory(folder);const size_t bins=e.theta.size();Matrix density(result.stats.mean.begin(),result.stats.mean.begin()+bins);
    WriteMatrix(folder+"/mueller_evaluated.dat",e.theta,density);
    auto pole=[&](size_t index,const std::string& name){std::ostringstream text;text<<"ScAngle 2pi*dcos M11 M12 M13 M14 M21 M22 M23 M24 M31 M32 M33 M34 M41 M42 M43 M44\n"<<std::setprecision(17)<<e.theta[index==bins?0:bins-1]<<" 0";
        for(double v:result.stats.mean[index])text<<' '<<v;
        text<<'\n';Atomic(folder+"/mueller_"+name+".dat",text.str());};
    pole(bins,"forward");pole(bins+1,"backward");
    std::ostringstream precision;precision<<"kind,theta_deg,N_per_seed,positive_M11,checked95";
    for(int j=0;j<16;++j)precision<<",estimated95_M"<<j/4+1<<j%4+1<<"_scaled_M11";
    for(int j=1;j<16;++j)precision<<",estimated95_M"<<j/4+1<<j%4+1<<"_over_M11";
    precision<<'\n';
    const auto x=GoRows(e);for(size_t r=0;r<x.size();++r){precision<<(r<bins?"finite_bin":(r==bins?"forward_delta":"backward_delta"))<<','<<Number(x[r])<<','<<count<<','<<(result.stats.mean[r][0]>0)<<','<<result.checked[r];
        for(double v:result.stats.interval[r])precision<<','<<Number(result.checked[r]?v:1e100);
        for(int j=1;j<16;++j)precision<<','<<Number(result.checked[r]?result.stats.ratioInterval[r][j]:1e100);
        precision<<'\n';}
    Atomic(folder+"/parameters_accuracy.csv",precision.str());
    std::ostringstream budget;budget<<"kind,theta_deg,N_per_seed,seeds,total_orientations,depth\n";
    for(size_t row=0;row<x.size();++row)budget<<(row<bins?"finite_bin":(row==bins?"forward_delta":"backward_delta"))<<','<<Number(x[row])<<','<<count<<",8,"<<8LL*count<<','<<e.target<<'\n';
    Atomic(folder+"/orientation_budget.csv",budget.str());
    e.Record(folder+"/pole_status.json","{\"units\":\"integrated ray cross sections, not per-steradian densities\",\"forward\":"+Quoted(result.checked[bins]?"positive_mass_checked":"zero_observed_not_bounded")+",\"backward\":"+Quoted(result.checked[bins+1]?"positive_mass_checked":"zero_observed_not_bounded")+"}");
}
void GoCheckpoint(Engine& e,int round,int count,int streak){
    const std::string payload=e.signature+" "+Number(round)+" "+Number(count)+" "+Number(streak)+"\n";
    Atomic(e.root+"/go_checkpoint.tsv",Hash(payload)+"\n"+payload);e.Log("allocation_checkpoint","{\"method\":\"go\",\"round\":"+Number(round)+",\"count_per_seed\":"+Number(count)+",\"pass_streak\":"+Number(streak)+"}");
}
int RunGO(Engine& e){
    if(Exists(e.root+"/final_outputs.fnv")){
        std::istringstream files(Read(e.root+"/final_outputs.fnv"));std::string name,digest;while(files>>name>>digest)Require(FileHash(e.root+"/"+name)==digest,"completed output changed");
        e.WarmCache();
        const std::string parameters=Read(e.root+"/final_parameters.json"),key="\"selected_state\":";
        const auto at=parameters.find(key);Require(at!=std::string::npos,"missing completed GO state");
        e.Record(e.root+"/fullauto_status.json",parameters.substr(at+key.size(),parameters.size()-at-key.size()-1));
        std::cout<<"fullauto GO: completed outputs verified; no calculation repeated"<<std::endl;return 0;
    }
    e.Record(e.root+"/precision_contract.json","{\"method\":\"go\",\"angular_model\":\"fixed finite solid-angle bins; no interpolation\",\"alpha\":\"analytic scattering-plane and pole averaging; two-dimensional Haar Sobol\",\"statistics\":\"eight independent scrambles; pointwise Student and Fieller estimates\",\"statistical_budget\":"+Number(e.statistical)+",\"depth_budget\":"+Number(e.depthBudget)+",\"interpolation_budget\":0,\"zero_finite_bins\":\"require further sampling\",\"zero_poles\":\"reported as zero observed, not certified absence\",\"physical_depth\":"+Number(e.target)+",\"global_simultaneous_confidence\":false}");
    int round=0,count=e.initial,streak=0;
    if(Exists(e.root+"/go_checkpoint.tsv")){
        const auto text=Read(e.root+"/go_checkpoint.tsv");const auto line=text.find('\n');const auto payload=text.substr(line+1);
        Require(Hash(payload)==text.substr(0,line),"GO checkpoint changed");std::istringstream input(payload);std::string signature;
        Require(bool(input>>signature>>round>>count>>streak)&&signature==e.signature&&count>0&&count<=e.maxCount&&(count&(count-1))==0&&round>=0,"GO checkpoint inputs changed");
        e.WarmCache();
    }
    double lastMaximum=1e100;const size_t bins=e.theta.size();bool havePrevious=count>e.initial;Matrix previous;
    if(havePrevious)previous=Mean(GoSample(e,count/2,e.seeds.data()));
    for(;round<e.roundLimit;++round){
        GoCheckpoint(e,round,count,streak);e.State("active",round,{count},bins,0,0,streak);
        auto samples=GoSample(e,count,e.seeds.data());auto result=EstimateGO(samples,bins,e.ratios);lastMaximum=result.maximum;
        const double change=havePrevious?GoDifference(result.stats.mean,previous,bins,e.ratios):1e100;
        const bool pass=result.maximum<=e.statistical&&change<=e.statistical;streak=pass?streak+1:0;e.lastRatio=result.stats.maxRatioInterval;
        const std::string folder=e.root+"/round_snapshots/n"+Number(e.target)+"_round"+Number(round)+"_N"+Number(count);GoOutputs(e,folder,result,count);
        e.Record(e.root+"/intermediate_current.json","{\"converged\":false,\"method\":\"go\",\"round\":"+Number(round)+",\"snapshot_folder\":"+Quoted(folder)+",\"matrix\":"+Quoted(folder+"/mueller_evaluated.dat")+",\"accuracy_and_parameters\":"+Quoted(folder+"/parameters_accuracy.csv")+",\"statistical95\":"+Number(result.maximum)+",\"interpolation_residual95\":0,\"pass_streak\":"+Number(streak)+"}");
        e.Log("round_decision","{\"round\":"+Number(round)+",\"count_per_seed\":"+Number(count)+",\"statistics_pass\":"+(result.maximum<=e.statistical?"true":"false")+",\"stability_pass\":"+(change<=e.statistical?"true":"false")+",\"statistical95\":"+Number(result.maximum)+",\"change_scaled_M11\":"+Number(change)+"}");
        e.State("active",round,{count},bins,result.maximum,0,streak);
        if(streak>=e.stablePasses){
            bool validation=true;
            if(e.finalAudit){
                unsigned fresh[8];for(int i=0;i<8;++i){fresh[i]=1200000000u+unsigned(round)*1024u+unsigned(i)*17u;Require(std::find(e.seeds.begin(),e.seeds.end(),fresh[i])==e.seeds.end(),"GO validation/production seed overlap");}
                e.State("validating",round,{count},bins,result.maximum,0,streak);
                auto independent=EstimateGO(GoSample(e,count,fresh),bins,e.ratios);
                const double difference=GoDifference(result.stats.mean,independent.stats.mean,bins,e.ratios);
                const double bound=difference+result.maximum+independent.maximum;
                validation=bound<=e.statistical;
                e.Record(e.root+"/independent_final_validation.json","{\"enabled\":true,\"passed\":"+(validation?std::string("true"):"false")+",\"difference_scaled_M11_or_ratio\":"+Number(difference)+",\"conservative_bound95\":"+Number(bound)+",\"count_per_seed\":"+Number(count)+"}");
                GoOutputs(e,e.root+"/independent_snapshot_n"+Number(e.target)+"_round"+Number(round),independent,count);
            }else e.Record(e.root+"/independent_final_validation.json","{\"enabled\":false,\"passed\":false}");
            if(validation){
                bool physics=true;
                if(e.physicsAudit||e.autoDepth){
                    unsigned fresh[8];for(int i=0;i<8;++i){fresh[i]=1400000000u+unsigned(round)*1024u+unsigned(i)*17u;
                        Require(std::find(e.seeds.begin(),e.seeds.end(),fresh[i])==e.seeds.end(),"GO physics/production seed overlap");}
                    auto original=GoSample(e,count,fresh),deep=GoSample(e,count,fresh,e.target+2);
                    auto bound=[&](const Samples& other){return PhysicalDepthBound(original,other,result.stats.mean,true,bins,e.ratios);};
                    const double depthError=bound(deep);const auto savedBase=e.base;const auto savedCutoff=e.selectedCutoff;double cutoffError=0;
                    if(e.physicsAudit){e.SetCutoff("off");try{cutoffError=bound(GoSample(e,count,fresh));}
                    catch(...){e.base=savedBase;e.selectedCutoff=savedCutoff;e.SetMirror(false);throw;}
                    e.base=savedBase;e.selectedCutoff=savedCutoff;e.SetMirror(false);}
                    physics=depthError<=(e.autoDepth?e.depthBudget:e.eps)&&cutoffError<=e.eps;
                    e.Record(e.root+"/physics_validation.json","{\"enabled\":true,\"passed\":"+(physics?std::string("true"):"false")+",\"depth_plus2_residual95_scaled_M11\":"+Number(depthError)+",\"cutoffs_off_residual95_scaled_M11\":"+Number(cutoffError)+",\"depth_pass\":"+std::string(depthError<=(e.autoDepth?e.depthBudget:e.eps)?"true":"false")+",\"cutoff_pass\":"+std::string(cutoffError<=e.eps?"true":"false")+",\"cutoffs_checked\":"+std::string(e.physicsAudit?"true":"false")+",\"depth_budget\":"+Number(e.autoDepth?e.depthBudget:e.eps)+",\"count_per_seed\":"+Number(count)+",\"requested_physics_retained\":true}");
                }else e.Record(e.root+"/physics_validation.json","{\"enabled\":false,\"requested_physics_retained\":true}");
                if(!physics){e.State("physics_validation_failed",round,{count},bins,result.maximum,0,streak);return 4;}
                GoOutputs(e,e.root,result,count);Atomic(e.root+"/mueller_fullauto.dat",Read(e.root+"/mueller_evaluated.dat"));
                e.independentPassed=e.finalAudit;e.State("converged",round,{count},bins,result.maximum,0,streak);
                std::ostringstream hashes;for(const char* name:{"controller_command.json","mueller_fullauto.dat","mueller_evaluated.dat","mueller_forward.dat","mueller_backward.dat","parameters_accuracy.csv","orientation_budget.csv","pole_status.json","precision_contract.json","independent_final_validation.json","physics_validation.json","final_parameters.json"})hashes<<name<<' '<<FileHash(e.root+"/"+name)<<'\n';if(e.autoDepth)for(const char* name:{"automatic_depth.tsv","depth_selection.json"})hashes<<name<<' '<<FileHash(e.root+"/"+name)<<'\n';Atomic(e.root+"/final_outputs.fnv",hashes.str());return 0;
            }streak=0;
        }
        if(count>=e.maxCount||round+1>=e.roundLimit)break;
        previous=result.stats.mean;havePrevious=true;count=std::min(e.maxCount,count>e.maxCount/2?e.maxCount:2*count);
    }
    e.State("requires_further_sampling",round,{count},bins,lastMaximum,0,streak);return 3;
}
}
