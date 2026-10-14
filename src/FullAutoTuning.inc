// Included inside Engine. All choices are frozen before production and restored
// before the adaptive checkpoint. Pilot caches include the complete execution
// recipe, so different cutoffs, symmetries and resources never share rows.
std::string ExecutionHash()const{
    std::ostringstream s;s<<workerThreads<<' '<<workerChunk<<'\n';for(const auto& gpu:gpus)s<<gpu<<'\n';for(const auto& v:base)s<<v<<'\n';return Hash(s.str());
}
void RemoveOption(const std::string& option,int count=1){
    for(size_t i=0;i<base.size();)if(base[i]==option){Require(i+size_t(count)<base.size(),"incomplete worker option");base.erase(base.begin()+i,base.begin()+i+count+1);}else ++i;
}
void SetRotation(int beta,int gamma){
    RemoveOption("--verified-rotation-symmetry",0);RemoveOption("--symmetry",2);
    rotationBeta=beta;rotationGamma=gamma;SetMirror(mirror);
    base.push_back("--symmetry");base.push_back(Number(beta));base.push_back(Number(gamma));
}
void SetCutoff(const std::string& profile){
    // Hard tree limits/retries are accuracy safeguards, not approximations.
    for(const char* name:{"--cutoff-profile","--beam-cutoff","--beam-cutoff-jones","--beam-cutoff-area","--beam-cutoff-importance",
            "--trace-cutoff","--trace-cutoff-jones","--trace-cutoff-area","--trace-cutoff-importance","--beam-area-ratio"})RemoveOption(name);
    base.push_back("--cutoff-profile");base.push_back(profile);selectedCutoff=profile;SetMirror(mirror);
}
double PairedBias(const Samples& candidate,const Samples& reference)const{
    const Matrix mean=Mean(reference),other=Mean(candidate),delta=Mean(Subtract(candidate,reference)),variance=Variance(Subtract(candidate,reference));
    for(size_t r=0;r<mean.size();++r)if(!(mean[r][0]>0&&other[r][0]>0))return 1e100;
    const auto a=Estimate(candidate),b=Estimate(reference);
    if(ratios&&(a.maxRatioInterval>1e50||b.maxRatioInterval>1e50))return 1e100;
    double maximum=0;
    for(size_t r=0;r<mean.size();++r)for(int j=0;j<16;++j){maximum=std::max(maximum,(std::fabs(delta[r][j])+t95*std::sqrt(variance[r][j]/8))/mean[r][0]);
        if(ratios&&j){const double ra=other[r][j]/other[r][0],rb=mean[r][j]/mean[r][0];double v=0;
            // Paired influence functions for ratios of means; no mean of ratios.
            for(int s=0;s<8;++s){const double z=(candidate[s][r][j]-ra*candidate[s][r][0])/other[r][0]
                -(reference[s][r][j]-rb*reference[s][r][0])/mean[r][0];v+=z*z/7;}
            maximum=std::max(maximum,std::fabs(ra-rb)+t95*std::sqrt(v/8));}
    }return maximum;
}
void SaveTuning(){
    std::ostringstream s;s<<std::setprecision(17)<<signature<<'\n'<<workerThreads<<' '<<workerChunk<<' '<<orientationBlock<<' '
        <<rotationBeta<<' '<<rotationGamma<<' '<<rotationError<<' '<<selectedCutoff<<' '<<cutoffError<<' '<<tuneCutoffs<<' '<<depthSharing<<' '<<hostBudgetMb<<'\n';
    s<<gpus.size()<<'\n';for(const auto& v:gpus)s<<v<<'\n';s<<base.size()<<'\n';for(const auto& v:base){Require(v.find('\n')==std::string::npos,"newline in worker argument");s<<v<<'\n';}
    Record(root+"/execution_tuning.tsv","1 "+Hash(s.str())+"\n"+s.str());
}
bool RestoreTuning(){
    if(!Exists(root+"/execution_tuning.tsv"))return false;
    const auto text=Read(root+"/execution_tuning.tsv");const auto newline=text.find('\n');std::istringstream header(text.substr(0,newline));int version;std::string digest;
    Require(bool(header>>version>>digest)&&version==1&&Hash(text.substr(newline+1))==digest,"execution tuning changed");
    std::istringstream s(text.substr(newline+1));std::string id;std::getline(s,id);Require(id==signature,"execution tuning inputs changed");
    Require(bool(s>>workerThreads>>workerChunk>>orientationBlock>>rotationBeta>>rotationGamma>>rotationError>>selectedCutoff>>cutoffError>>tuneCutoffs>>depthSharing>>hostBudgetMb),"invalid execution tuning");
    size_t n;std::string line;s>>n;std::getline(s,line);gpus.clear();for(size_t i=0;i<n;++i){std::getline(s,line);
        Require(std::find(gpuCandidates.begin(),gpuCandidates.end(),line)!=gpuCandidates.end(),"tuned GPU outside allowed set");gpus.push_back(line);}
    s>>n;std::getline(s,line);base.clear();for(size_t i=0;i<n;++i){Require(bool(std::getline(s,line)),"truncated execution recipe");base.push_back(line);}
    Require(workerThreads>0&&workerChunk>=0&&orientationBlock>0&&!gpus.empty(),"invalid tuned resources");SetMirror(false);return true;
}
void TuneResources(const std::vector<double>& x){
    const unsigned seedSet[2]={8001,8023};const int n=std::min(kernelCount,maxCount);
    std::ostringstream report;report<<"{\"policy\":"<<Quoted(resourcePolicy)<<",\"gpu_candidates\":[";
    std::vector<std::pair<double,std::string>> speed;Matrix reference;bool first=true;
    const auto allowed=gpus;
    if(autoGpus){for(const auto& gpu:allowed){
            const int device=int(std::find(gpuCandidates.begin(),gpuCandidates.end(),gpu)-gpuCandidates.begin());long long freeBytes=0,totalBytes=0;std::string memoryError;
            if(!QueryVisibleGpuMemory(device,freeBytes,totalBytes,memoryError)||freeBytes<512LL*1024*1024){
                report<<(first?"":",")<<"{\"device\":"<<Quoted(gpu)<<",\"skipped\":true,\"free_bytes\":"<<freeBytes<<",\"reason\":"<<Quoted(memoryError.empty()?"less than 512 MiB free":memoryError)<<"}";first=false;continue;}
            gpus={gpu};auto jobs=DirectStage(n,seedSet,2,target,x,false);double time=jobs[0].seconds+jobs[1].seconds,error=0;
            if(reference.empty())reference=jobs[0].matrix;else error=Difference(jobs[0].matrix,reference,reference);
            report<<(first?"":",")<<"{\"device\":"<<Quoted(gpu)<<",\"seconds\":"<<Number(time)<<",\"all16_difference\":"<<Number(error)<<"}";first=false;
            if(error<std::min(1e-7,eps/1024))speed.emplace_back(time,gpu);}
        Require(!speed.empty(),"no equivalent GPU pilot");std::sort(speed.begin(),speed.end());gpus.clear();
        // Very slow devices lengthen eight-seed stage completion instead of
        // helping it. The dynamic queue assigns the next seed to a free device.
        for(const auto& entry:speed)if(entry.first<=2*speed.front().first)gpus.push_back(entry.second);
    }else gpus=allowed;
    const auto selected=gpus;const int baselineThreads=workerThreads,baselineChunk=workerChunk;
    const int cap=std::max(1,QueryRuntimeResourceSnapshot().physicalCores/int(selected.size()));
    gpus={selected.front()};auto baseline=DirectStage(n,seedSet,2,target,x,false);
    reference=baseline[0].matrix;double best=baseline[0].seconds+baseline[1].seconds;int chosenThreads=workerThreads,chosenChunk=workerChunk;
    report<<"],\"resource_candidates\":[{\"threads\":"<<workerThreads<<",\"chunk\":"<<workerChunk<<",\"seconds\":"<<Number(best)<<",\"all16_difference\":0}";
    double nativePilot=0;for(const auto& job:baseline)nativePilot+=job.NativeWall();
    const double forecastNative=nativePilot*double(initial+4LL*minimum)/n*4/std::max(size_t(1),selected.size());
    const double predictedSearch=best*10;
    const bool amortized=forecastNative*.25>=predictedSearch;
    if(resourcePolicy=="auto"&&amortized){
        std::set<int> threads;if(!args.IsCatched("threads"))for(int t:{1,2,4,8,16,32})if(t<=cap)threads.insert(t);
        threads.insert(workerThreads);
        for(int t:threads){workerThreads=t;auto jobs=DirectStage(n,seedSet,2,target,x,false);const double seconds=jobs[0].seconds+jobs[1].seconds;
            const double error=std::max(Difference(jobs[0].matrix,baseline[0].matrix,reference),Difference(jobs[1].matrix,baseline[1].matrix,reference));
            Log("resource_candidate","{\"threads\":"+Number(workerThreads)+",\"chunk\":"+Number(workerChunk)+",\"seconds\":"+Number(seconds)+",\"all16_difference_scaled_M11\":"+Number(error)+",\"budget\":"+Number(std::min(1e-7,eps/1024))+",\"numerically_accepted\":"+(error<std::min(1e-7,eps/1024)?"true":"false")+"}");
            report<<",{\"threads\":"<<t<<",\"chunk\":"<<workerChunk<<",\"seconds\":"<<Number(seconds)<<",\"all16_difference\":"<<Number(error)<<"}";
            if(error<std::min(1e-7,eps/1024)&&seconds<best*.95){best=seconds;chosenThreads=t;}}
        workerThreads=chosenThreads;
        if(!args.IsCatched("chunk"))for(int chunk:{64,256,1024,4096})if(chunk<=n||chunk==64){workerChunk=chunk;
            auto jobs=DirectStage(n,seedSet,2,target,x,false);const double seconds=jobs[0].seconds+jobs[1].seconds;
            const double error=std::max(Difference(jobs[0].matrix,baseline[0].matrix,reference),Difference(jobs[1].matrix,baseline[1].matrix,reference));
            Log("resource_candidate","{\"threads\":"+Number(workerThreads)+",\"chunk\":"+Number(workerChunk)+",\"seconds\":"+Number(seconds)+",\"all16_difference_scaled_M11\":"+Number(error)+",\"budget\":"+Number(std::min(1e-7,eps/1024))+",\"numerically_accepted\":"+(error<std::min(1e-7,eps/1024)?"true":"false")+"}");
            report<<",{\"threads\":"<<workerThreads<<",\"chunk\":"<<chunk<<",\"seconds\":"<<Number(seconds)<<",\"all16_difference\":"<<Number(error)<<"}";
            if(error<std::min(1e-7,eps/1024)&&seconds<best*.95){best=seconds;chosenChunk=chunk;}}
    }
    workerThreads=chosenThreads;workerChunk=chosenChunk;auto chosen=DirectStage(n,seedSet,2,target,x,false);
    double overhead=0,perOrientation=0;bool processingMeasured=true;for(const auto& job:chosen){
        overhead+=job.OutsideProcessing();perOrientation+=job.NativeWall()/n;processingMeasured=processingMeasured&&job.processingSeconds>=0;}
    overhead/=2;perOrientation/=2;
    if(resourcePolicy=="auto"&&!args.IsCatched("fullauto_orientation_block")&&perOrientation>0){
        const int lower=std::min(4096,maxCount),upper=std::min(1048576,maxCount);
        int proposed=Power2(std::max(double(lower),overhead/(.05*perOrientation)));
        while(proposed>upper)proposed/=2;
        while(proposed>lower&&proposed*perOrientation>120)proposed/=2;
        orientationBlock=proposed;
    }
    gpus=selected;report<<"],\"selected_gpus\":[";for(size_t i=0;i<gpus.size();++i)report<<(i?",":"")<<Quoted(gpus[i]);
    report<<"],\"host_memory_budget_mib\":"<<hostBudgetMb<<",\"selected_threads\":"<<workerThreads<<",\"selected_chunk\":"<<workerChunk<<",\"selected_orientation_block\":"<<orientationBlock
        <<",\"search_amortized\":"<<(amortized?"true":"false")<<",\"forecast_native_production_seconds\":"<<Number(forecastNative)<<",\"forecast_search_seconds\":"<<Number(predictedSearch)
        <<",\"explicit_threads_preserved\":"<<(args.IsCatched("threads")?"true":"false")<<",\"explicit_chunk_preserved\":"<<(args.IsCatched("chunk")?"true":"false")
        <<",\"outside_processing_seconds_per_worker\":"<<Number(overhead)<<",\"processing_seconds_per_orientation\":"<<Number(perOrientation)
        <<",\"processing_wall_measured\":"<<(processingMeasured?"true":"false")
        <<",\"block_model\":\"worker wall minus processing wall (overlapping phases are not summed); aim <=5% outside-processing overhead and <=120s processing work, bounded by 4096..1048576 and max orientations\",\"baseline_threads\":"<<baselineThreads<<",\"baseline_chunk\":"<<baselineChunk<<"}";
    Record(root+"/resource_calibration.json",report.str());
}
void TuneRotations(const std::vector<double>& x){
    if(args.IsCatched("sym")){rotationBeta=args.GetIntValue("sym",0);rotationGamma=args.GetIntValue("sym",1);Record(root+"/rotation_verification.json","{\"explicit_override\":true}");return;}
    if(rotationPolicy=="off"){SetRotation(1,1);Record(root+"/rotation_verification.json","{\"enabled\":false,\"reason\":\"disabled\"}");return;}
    RemoveOption("--symmetry",2);base.push_back("--verified-rotation-symmetry");const unsigned probe=8101;
    auto job=DirectStage(16,&probe,1,target,x,false);const std::string log=Read(job[0].folder+"/stdout.log"),marker="Verified body rotation factors: ";const auto pos=log.find(marker);
    Require(pos!=std::string::npos,"missing verified rotation factors");std::istringstream factors(log.substr(pos+marker.size()));int beta,gamma;Require(bool(factors>>beta>>gamma),"invalid rotation factors");SetRotation(beta,gamma);
    const unsigned seeds[8]={8101,8111,8117,8123,8147,8161,8167,8171};std::ostringstream report;report<<"{\"geometry_beta_factor\":"<<beta<<",\"geometry_gamma_factor\":"<<gamma<<",\"checks\":[";bool first=true;
    for(const char* axis:{"gamma","beta"})if((std::string(axis)=="gamma"?rotationGamma:rotationBeta)>1){double worst=0;
        std::set<int> depths={std::max(1,target-4),target};for(int depth:depths){auto reduced=DirectStage(std::min(pilot,128),seeds,8,depth,x,false);
            base.push_back("--haar-rotation-audit");base.push_back(axis);auto paired=DirectStage(std::min(pilot,128),seeds,8,depth,x,false);RemoveOption("--haar-rotation-audit");
            worst=std::max(worst,PairedBias(Data(reduced),Data(paired)));}
        const double bound=5*worst;report<<(first?"":",")<<"{\"axis\":"<<Quoted(axis)<<",\"paired_bias95\":"<<Number(bound)<<"}";first=false;
        if(bound>eps/512)SetRotation(std::string(axis)=="beta"?1:rotationBeta,std::string(axis)=="gamma"?1:rotationGamma);
        else rotationError+=bound;
    }
    report<<"],\"selected_beta_factor\":"<<rotationBeta<<",\"selected_gamma_factor\":"<<rotationGamma<<",\"estimated_error\":"<<Number(rotationError)<<",\"budget\":"<<Number(eps/256)
        <<",\"scope\":\"oriented-facet bijection around body origin plus equivalent-rotation paired pilot; fixed axes only\"}";Record(root+"/rotation_verification.json",report.str());
}
Samples ProfilePilot(const std::vector<double>& x,int n,const unsigned* seeds,double& cost){
    std::vector<int> depths;for(int d:{std::max(1,target-4),std::max(1,target-2),target})if(depths.empty()||depths.back()!=d)depths.push_back(d);
    std::vector<Job> final;cost=0;
    for(size_t l=0;l<depths.size();++l){std::vector<Job> upper;
        if(l&&depthSharing&&selectedCutoff=="off"){
            try{upper=DirectStage(n,seeds,8,depths[l],x,false,false,0,depths[l-1]);}
            catch(const std::exception& error){const std::string message=error.what();const auto pos=message.find("see ");
                if(pos==std::string::npos||Read(message.substr(pos+4)).find("paired-depth cache unsupported")==std::string::npos)throw;
                depthSharing=false;}
        }
        if(upper.empty()){upper=DirectStage(n,seeds,8,depths[l],x,false);if(l){const auto lower=DirectStage(n,seeds,8,depths[l-1],x,false);for(const auto& job:lower)cost+=job.seconds;}}
        for(const auto& job:upper)cost+=job.seconds;
        final=upper;
    }return Data(final);
}
void IncludeBackscatter(std::vector<double>& x)const{
    for(double angle:{170.,175.,178.,179.,179.5,180.})if(angle>=theta.front()&&angle<=theta.back())x.push_back(angle);
    std::sort(x.begin(),x.end());x.erase(std::unique(x.begin(),x.end()),x.end());
}
void TuneCutoffs(const std::vector<double>& input){
    std::vector<double> x=input;IncludeBackscatter(x);
    bool explicitCutoff=args.IsCatched("r");for(const char* key:{"cutoff_profile","beam_cutoff","beam_cutoff_j","beam_cutoff_area","beam_cutoff_importance","trace_cutoff","trace_cutoff_j","trace_cutoff_area","trace_cutoff_importance"})explicitCutoff=explicitCutoff||args.IsCatched(key);
    if(cutoffPolicy=="off"||explicitCutoff){Record(root+"/cutoff_calibration.json","{\"enabled\":false,\"explicit_settings_preserved\":"+std::string(explicitCutoff?"true":"false")+"}");return;}
    if(x.size()>size_t(maxTheta)){
        tuneCutoffs=true;SetCutoff("off");
        Record(root+"/cutoff_calibration.json","{\"enabled\":true,\"selected\":\"off\",\"reason\":\"explicit theta limit cannot fit cutoff guards\"}");return;
    }
    tuneCutoffs=true;SetCutoff("off");const bool sharing=depthSharing;const unsigned seeds[8]={9001,9007,9011,9013,9029,9041,9043,9049};
    double best=0;auto reference=ProfilePilot(x,std::min(std::min(pilot,2048),maxCount),seeds,best);const bool supportedSharing=depthSharing;std::string winner="off";double winnerBias=0;
    double nativePilot=0;for(const auto& entry:cache){const auto& job=entry.second;
        if(job.executionHash!=ExecutionHash()||job.assembled||job.count!=std::min(std::min(pilot,2048),maxCount)||std::find(seeds,seeds+8,job.seed)==seeds+8)continue;
        nativePilot+=job.NativeWall();}
    const double forecastGain=.25*nativePilot*double(initial)/std::min(std::min(pilot,2048),maxCount);
    if(best<30&&forecastGain<3*best){
        Record(root+"/cutoff_calibration.json","{\"enabled\":true,\"selected\":\"off\",\"reference\":\"off\",\"search_amortized\":false,\"baseline_seconds\":"+Number(best)+",\"forecast_gain_seconds\":"+Number(forecastGain)+",\"reason\":\"extra preset pilots cost more than forecast savings on this short workload\"}");
        depthSharing=supportedSharing&&sharing;depthVerified=false;return;
    }
    std::ostringstream report;report<<"{\"enabled\":true,\"reference\":\"off\",\"budget\":"<<Number(eps/8)<<",\"candidates\":[{\"profile\":\"off\",\"seconds\":"<<Number(best)<<",\"bias95\":0,\"shared_depths\":"<<(supportedSharing?"true":"false")<<"}";
    for(const char* profile:{"safe","balanced","fast"}){SetCutoff(profile);depthSharing=false;double seconds=0;auto candidate=ProfilePilot(x,std::min(std::min(pilot,2048),maxCount),seeds,seconds);
        double bias=PairedBias(candidate,reference);
        Log("cutoff_candidate","{\"profile\":"+Quoted(profile)+",\"seconds\":"+Number(seconds)+",\"bias95\":"+Number(bias)+",\"budget\":"+Number(eps/8)+",\"within_bias_budget\":"+(bias<=eps/8?"true":"false")+",\"ten_percent_faster\":"+(seconds<best*.90?"true":"false")+"}");report<<",{\"profile\":"<<Quoted(profile)<<",\"seconds\":"<<Number(seconds)<<",\"bias95\":"<<Number(bias)<<",\"shared_depths\":false}";
        if(bias<=eps/8&&seconds<best*.90){best=seconds;winner=profile;winnerBias=bias;}}
    SetCutoff(winner);cutoffError=winnerBias;depthSharing=winner=="off"?supportedSharing&&sharing:false;depthVerified=false;
    report<<"],\"selected\":"<<Quoted(winner)<<",\"final_fresh_seed_audit_required\":"<<(winner=="off"?"false":"true")<<",\"scope\":\"all16 raw and ratios; measured multilevel workload including shared-depth eligibility; representative theta\"}";
    Record(root+"/cutoff_calibration.json",report.str());
}
void PrepareTuning(const std::vector<double>& angles){
    if(RestoreTuning())return;
    Record(root+"/automatic_cutoff_validation.json","{\"enabled\":false}");
    std::vector<double> x;const size_t rows=std::min(size_t(17),angles.size());for(size_t i=0;i<rows;++i)x.push_back(angles[i*(angles.size()-1)/(rows-1)]);
    const bool oldControls=controls;controls=false;SetMirror(false);phi=fixedPhi?fixedPhi:4;
    if(resourcePolicy=="auto"||autoGpus)TuneResources(x);else Record(root+"/resource_calibration.json","{\"policy\":\"off\"}");
    TuneRotations(x);
    TuneCutoffs(x);controls=oldControls;SetMirror(false);SaveTuning();
}
bool AuditAutomaticCutoffs(const std::vector<double>& angles,int n,int round){
    if(!tuneCutoffs||selectedCutoff=="off"||selectedCutoff=="requested")return true;
    std::vector<double> x;const size_t rows=std::min(size_t(65),angles.size());for(size_t i=0;i<rows;++i)x.push_back(angles[i*(angles.size()-1)/(rows-1)]);
    IncludeBackscatter(x);
    unsigned seedSet[8];for(int s=0;s<8;++s){seedSet[s]=1600000000u+unsigned(round)*1024u+17u*s;
        Require(std::find(seeds.begin(),seeds.end(),seedSet[s])==seeds.end(),"cutoff audit/production seed overlap");}
    n=std::min(n,std::max(minimum,std::min(pilot,4096)));
    const bool oldControls=controls;controls=false;const std::string profile=selectedCutoff;double error=1e100;
    for(int attempt=0;attempt<3;++attempt){auto requested=DirectStage(n,seedSet,8,target,x,false);SetCutoff("off");auto reference=DirectStage(n,seedSet,8,target,x,false);SetCutoff(profile);
        error=PairedBias(Data(requested),Data(reference));if(error<=eps/8||n>=maxCount||attempt==2)break;n=std::min(maxCount,n>maxCount/2?maxCount:2*n);}
    controls=oldControls;SetMirror(mirror);cutoffError=error;const bool pass=error<=eps/8;
    Record(root+"/automatic_cutoff_validation.json","{\"passed\":"+std::string(pass?"true":"false")+",\"profile\":"+Quoted(profile)+",\"bias95\":"+Number(error)+",\"budget\":"+Number(eps/8)+",\"count_per_seed\":"+Number(n)+",\"theta_count\":"+Number(x.size())+",\"fresh_seeds\":true,\"fallback_to_off\":"+(pass?"false":"true")+"}");
    if(!pass){SetCutoff("off");cutoffError=0;depthSharing=coherent&&On("fullauto_depth_share");depthVerified=false;weights.clear();polarWeights.clear();weightTheta.clear();SetMirror(mirror);SaveTuning();}
    SaveTuning();return pass;
}
