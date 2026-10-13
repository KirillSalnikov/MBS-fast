// Shared automatic physical-depth selection. Sampling/geometry and exact
// worker commands retain the same cache identity as production.
namespace {
void SaveSelectedDepth(Engine& e,int selected){
    const std::string payload=e.signature+" "+Number(selected)+"\n";
    Atomic(e.root+"/automatic_depth.tsv",Hash(payload)+"\n"+payload);
    e.target=selected;e.SetMirror(e.mirror);
    e.Log("physical_depth_selected","{\"target_depth\":"+Number(selected)+",\"start\":"+Number(e.depthStart)+",\"maximum\":"+Number(e.depthMaximum)+",\"budget\":"+Number(e.depthBudget)+"}");
}
struct DepthResidual {double upper=0,lower=0,mean=0,interval=0;};
DepthResidual DepthDifference(const Samples& first,const Samples& second,bool go,size_t bins,bool ratios=true,const Matrix* requestedScale=nullptr){
    const auto a=Mean(first),b=Mean(second);const auto& scale=requestedScale?*requestedScale:b;
    const auto delta=Subtract(second,first);const auto mean=Mean(delta),variance=Variance(delta);
    DepthResidual result;
    auto accumulate=[&](double bias,double ci){result.upper=std::max(result.upper,bias+ci);result.lower=std::max(result.lower,std::max(0.,bias-ci));
        result.mean=std::max(result.mean,bias);result.interval=std::max(result.interval,ci);};
    for(size_t row=0;row<scale.size();++row){
        if(!(scale[row][0]>0)){
            if(!go||row<bins||mean[row][0]!=0||variance[row][0]!=0)result.upper=1e100;
            continue;
        }
        if(ratios){
            double va=0,vb=0;for(int seed=0;seed<8;++seed){va+=std::pow(first[seed][row][0]-a[row][0],2)/7;vb+=std::pow(second[seed][row][0]-b[row][0],2)/7;}
            if(t95*std::sqrt(va/8)>=a[row][0]||t95*std::sqrt(vb/8)>=b[row][0])result.upper=1e100;
        }
        for(int j=0;j<16;++j){
            accumulate(std::fabs(mean[row][j])/scale[row][0],t95*std::sqrt(variance[row][j]/8)/scale[row][0]);
            if(!ratios||j==0)continue;
            if(!(a[row][0]>0&&b[row][0]>0)){result.upper=1e100;continue;}
            const double ra=a[row][j]/a[row][0],rb=b[row][j]/b[row][0];double v=0;
            // Paired influence function for a difference of ratios of means.
            // Retains covariance both between depths and with M11.
            for(int seed=0;seed<8;++seed){
                const double influence=(second[seed][row][j]-rb*second[seed][row][0])/b[row][0]
                    -(first[seed][row][j]-ra*first[seed][row][0])/a[row][0];v+=influence*influence/7;
            }
            accumulate(std::fabs(rb-ra),t95*std::sqrt(v/8));
        }
    }return result;
}
double PhysicalDepthBound(const Samples& first,const Samples& second,const Matrix& reference,bool go,size_t bins,bool ratios){
    return DepthDifference(first,second,go,bins,ratios,&reference).upper;
}
void PrepareDepthStartup(Engine& e){
    if(!e.autoDepth||e.go||e.depthOrder=="depth-first"||(e.depthOrder=="auto"&&!e.config.useGpu))return;
    const std::string checkpoint=e.root+"/depth_startup.tsv";
    if(Exists(checkpoint)){
        const auto text=Read(checkpoint);const auto split=text.find('\n');const auto payload=text.substr(split+1);Require(Hash(payload)==text.substr(0,split),"depth startup changed");
        std::istringstream input(payload);std::string signature;Require(bool(input>>signature>>e.phi>>e.kernel.block>>e.kernel.pipeline>>e.kernel.warp>>e.controls)&&signature==e.signature,"depth startup inputs changed");
        e.startupKernel=true;e.SetMirror(false);e.WarmCache();return;
    }
    std::vector<double> x;const size_t size=std::min(size_t(9),std::min(e.theta.size(),size_t(e.maxTheta)));
    for(size_t row=0;row<size;++row)x.push_back(e.theta[row*(e.theta.size()-1)/(size-1)]);
    const int oldPilot=e.pilot,oldKernel=e.kernelCount;e.pilot=std::min(e.pilot,64);e.kernelCount=std::min(e.kernelCount,64);
    const double begin=Seconds();e.State("startup_sampling_tuning",0,{e.pilot},x.size(),0,0,0);
    try{e.Calibrate(x);}catch(const std::exception& error){
        if(!e.DisableUnsupportedControls(error)){e.pilot=oldPilot;e.kernelCount=oldKernel;throw;}e.Calibrate(x);
    }
    e.pilot=oldPilot;e.kernelCount=oldKernel;e.startupKernel=true;
    std::ostringstream payload;payload<<e.signature<<' '<<e.phi<<' '<<e.kernel.block<<' '<<e.kernel.pipeline<<' '<<e.kernel.warp<<' '<<e.controls<<'\n';
    Atomic(checkpoint,Hash(payload.str())+"\n"+payload.str());
    e.Record(e.root+"/depth_startup.json","{\"order\":\"sampling-first\",\"phi\":"+Number(e.phi)+",\"kernel_pipeline\":"+(e.kernel.pipeline?std::string("true"):"false")+",\"pilot_count\":"+Number(std::min(oldPilot,64))+",\"theta_rows\":"+Number(x.size())+",\"seconds\":"+Number(Seconds()-begin)+",\"scope\":\"cheap kernel/phi training before paired depth; final estimator is refitted at selected depth\"}");
}
void SelectAutomaticDepth(Engine& e){
    if(!e.autoDepth){e.Record(e.root+"/depth_selection.json","{\"automatic\":false,\"target_depth\":"+Number(e.target)+",\"explicit_target_preserved\":true}");return;}
    const std::string checkpoint=e.root+"/automatic_depth.tsv";
    if(Exists(checkpoint)){
        const auto text=Read(checkpoint);const auto split=text.find('\n');Require(split!=std::string::npos,"invalid depth checkpoint");
        const auto payload=text.substr(split+1);Require(Hash(payload)==text.substr(0,split),"depth checkpoint changed");
        std::istringstream input(payload);std::string signature,extra;int selected;
        Require(bool(input>>signature>>selected)&&!(input>>extra)&&signature==e.signature&&selected>=e.depthStart&&selected<=e.depthMaximum,"depth checkpoint inputs changed");
        e.target=selected;e.SetMirror(false);e.WarmCache();return;
    }
    std::vector<double> x;
    if(e.go)x=GoRows(e);
    else{const size_t size=std::min(size_t(65),std::min(e.theta.size(),size_t(e.maxTheta)));for(size_t row=0;row<size;++row)x.push_back(e.theta[row*(e.theta.size()-1)/(size-1)]);}
    WriteGrid(e.root+"/depth_pilot_theta.csv",x);
    const bool savedControls=e.controls;e.controls=false;
    int count=e.depthPilot,selected=e.depthStart;DepthResidual residual;bool passed=false;
    try{
        for(int depth=e.depthStart;;depth=std::min(e.depthMaximum,depth+2)){
            unsigned seeds[8];for(int i=0;i<8;++i){seeds[i]=1700000000u+unsigned(e.depthStart)*1024u+unsigned(i)*17u;
                Require(std::find(e.seeds.begin(),e.seeds.end(),seeds[i])==e.seeds.end(),"depth pilot/production seed overlap");}
            for(int attempt=0;;++attempt){
                e.target=depth;e.State("selecting_depth",attempt,{count},x.size(),0,0,0);
                e.Record(e.root+"/depth_search_progress.json","{\"candidate\":"+Number(depth)+",\"comparison_depth\":"+Number(depth+2)+",\"count_per_seed\":"+Number(count)+",\"attempt\":"+Number(attempt)+",\"status\":\"sampling paired depth matrices\"}");
                Samples first,second;
                if(e.go){first=GoSample(e,count,seeds,depth);second=GoSample(e,count,seeds,depth+2);}
                else{first=e.Data(e.PrefixStage(count,seeds,8,depth,x,true));second=e.Data(e.PrefixStage(count,seeds,8,depth+2,x,true));}
                residual=DepthDifference(first,second,e.go,e.theta.size(),e.ratios);passed=residual.upper<=e.depthBudget;
                const std::string snapshot=e.root+"/depth_pilot_snapshots/n"+Number(depth)+"_N"+Number(count);Directory(snapshot);
                if(e.go){GoOutputs(e,snapshot,EstimateGO(first,e.theta.size(),e.ratios),count);e.target=depth+2;
                    GoOutputs(e,snapshot+"/depth_plus2",EstimateGO(second,e.theta.size(),e.ratios),count);e.target=depth;
                }else{WriteMatrix(snapshot+"/mueller_evaluated.dat",x,Mean(first));WriteMatrix(snapshot+"/mueller_depth_plus2.dat",x,Mean(second));
                    SaveReplicaMatrices(e,"depth_pilot_n"+Number(depth)+"_N"+Number(count),first,x,seeds);
                    SaveReplicaMatrices(e,"depth_pilot_n"+Number(depth+2)+"_N"+Number(count),second,x,seeds);}
                e.Record(e.root+"/intermediate_current.json","{\"converged\":false,\"phase\":\"depth_selection\",\"candidate\":"+Number(depth)+",\"snapshot_folder\":"+Quoted(snapshot)+",\"matrix\":"+Quoted(snapshot+"/mueller_evaluated.dat")+",\"paired_depth_residual95\":"+Number(residual.upper)+"}");

                std::ostringstream report;report<<"{\"automatic\":true,\"candidate\":"<<depth<<",\"comparison_depth\":"<<depth+2
                    <<",\"count_per_seed\":"<<count<<",\"theta_rows\":"<<x.size()<<",\"phi\":"<<e.phi<<",\"mean_residual_scaled_M11\":"<<Number(residual.mean)
                    <<",\"interval95_scaled_M11\":"<<Number(residual.interval)<<",\"residual_upper95\":"<<Number(residual.upper)<<",\"residual_lower95\":"<<Number(residual.lower)
                    <<",\"budget\":"<<Number(e.depthBudget)<<",\"passed\":"<<(passed?"true":"false")<<",\"theta_file\":"<<Quoted(e.root+"/depth_pilot_theta.csv")<<",\"scope\":\"paired initial orientation pilot; representative requested angles; fresh final depth audit required\"}";
                e.Record(e.root+"/depth_selection_candidate"+Number(depth)+"_N"+Number(count)+".json",report.str());e.Log("depth_candidate",report.str());
                if(passed||residual.lower>e.depthBudget||count>=e.maxCount||attempt==4)break;
                count=std::min(e.maxCount,count>e.maxCount/2?e.maxCount:2*count);
            }
            selected=depth;if(passed||depth==e.depthMaximum)break;
        }
    }catch(...){e.controls=savedControls;throw;}
    e.controls=savedControls;SaveSelectedDepth(e,selected);
    e.Record(e.root+"/depth_selection.json","{\"automatic\":true,\"target_depth\":"+Number(selected)+",\"start\":"+Number(e.depthStart)+",\"maximum\":"+Number(e.depthMaximum)+",\"pilot_passed\":"+(passed?std::string("true"):"false")+",\"count_per_seed\":"+Number(count)+",\"residual_upper95\":"+Number(residual.upper)+",\"budget\":"+Number(e.depthBudget)+",\"fresh_final_audit_required\":true}");
}
bool PromoteAutomaticDepth(Engine& e,int code){
    if(code!=4||!e.autoDepth||e.target>=e.depthMaximum)return false;
    const auto validation=CompactJson(Read(e.root+"/physics_validation.json"));
    if(validation.find("\"depth_pass\":false")==std::string::npos||validation.find("\"cutoff_pass\":false")!=std::string::npos)return false;
    const int old=e.target,next=std::min(e.depthMaximum,old+2);const std::string archive=e.root+"/depth_attempts/n"+Number(old);Directory(archive);
    // Keep the previous physical estimator and all of its matrices. Adaptive
    // allocations/weights must be rebuilt for the new target, not restored.
    for(const char* file:{"adaptive_checkpoint.tsv","go_checkpoint.tsv","execution_tuning.tsv","physics_validation.json","fullauto_status.json","calibration.json","initial_calibration.json","depth_selection.json"}){
        if(Exists(e.root+"/"+file))Require(rename((e.root+"/"+file).c_str(),(archive+"/"+file).c_str())==0,"cannot archive previous depth attempt");
    }
    SaveSelectedDepth(e,next);e.Record(e.root+"/depth_selection.json","{\"automatic\":true,\"target_depth\":"+Number(next)+",\"reason\":\"promoted after fresh final depth validation\",\"previous\":"+Number(old)+",\"maximum\":"+Number(e.depthMaximum)+",\"budget\":"+Number(e.depthBudget)+"}");e.Log("physical_depth_promoted","{\"previous\":"+Number(old)+",\"next\":"+Number(next)+",\"reason\":\"fresh final depth audit failed\",\"archive\":"+Quoted(archive)+"}");return true;
}
}
