// Included inside Engine. All choices are frozen before production and restored
// before the adaptive checkpoint. Pilot caches include the complete execution
// recipe, so different cutoffs, symmetries and resources never share rows.
// Controller-owned, bounded to one recipe. Check the actual inputs instead of
// relying on setters: pilots and checkpoint restore also replace base/gpus
// directly. Preserve the serialized recipe so old jobs remain compatible.
mutable bool executionRecipeValid=false;
mutable int executionRecipeThreads=0,executionRecipeChunk=0;
mutable bool executionRecipePhiPilot=false;
mutable std::vector<std::string> executionRecipeGpus,executionRecipeBase;
mutable std::string executionRecipeAnglePilot,executionRecipeHash;
const std::string& ExecutionHash()const{
    if(executionRecipeValid&&executionRecipeThreads==workerThreads&&executionRecipeChunk==workerChunk
            &&executionRecipePhiPilot==phiCutoffPilot&&executionRecipeAnglePilot==angleBatchPilotRecipe
            &&executionRecipeGpus==gpus&&executionRecipeBase==base)return executionRecipeHash;
    std::ostringstream s;s<<workerThreads<<' '<<workerChunk<<'\n';for(const auto& gpu:gpus)s<<gpu<<'\n';for(const auto& v:base)s<<v<<'\n';
    // A variance pilot of B_N+(O_m-B_m) is a different estimator from a
    // full-off pilot with the same seed. Never reuse the latter's cache.
    if(phiCutoffPilot)s<<"matched-cutoff-phi-pilot-v1\n";
    if(!angleBatchPilotRecipe.empty())s<<"angle-batch-pilot "<<angleBatchPilotRecipe<<'\n';
    executionRecipeThreads=workerThreads;executionRecipeChunk=workerChunk;executionRecipePhiPilot=phiCutoffPilot;
    executionRecipeGpus=gpus;executionRecipeBase=base;executionRecipeAnglePilot=angleBatchPilotRecipe;
    executionRecipeHash=Hash(s.str());executionRecipeValid=true;return executionRecipeHash;
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
        <<rotationBeta<<' '<<rotationGamma<<' '<<rotationError<<' '<<selectedCutoff<<' '<<cutoffError<<' '<<tuneCutoffs<<' '<<depthSharing<<' '<<hostBudgetMb<<' '<<analyticMeanThreads<<'\n';
    s<<gpus.size()<<'\n';for(const auto& v:gpus)s<<v<<'\n';s<<base.size()<<'\n';for(const auto& v:base){Require(v.find('\n')==std::string::npos,"newline in worker argument");s<<v<<'\n';}
    Record(root+"/execution_tuning.tsv","2 "+Hash(s.str())+"\n"+s.str());
}
bool RestoreTuning(){
    const long long dispatchHostBudget=hostBudgetMb;
    if(!Exists(root+"/execution_tuning.tsv"))return false;
    const auto text=Read(root+"/execution_tuning.tsv");const auto newline=text.find('\n');std::istringstream header(text.substr(0,newline));int version;std::string digest;
    Require(bool(header>>version>>digest)&&version==2&&Hash(text.substr(newline+1))==digest,"execution tuning changed");
    std::istringstream s(text.substr(newline+1));std::string id;std::getline(s,id);Require(id==signature,"execution tuning inputs changed");
    Require(bool(s>>workerThreads>>workerChunk>>orientationBlock>>rotationBeta>>rotationGamma>>rotationError>>selectedCutoff>>cutoffError>>tuneCutoffs>>depthSharing>>hostBudgetMb>>analyticMeanThreads),"invalid execution tuning");
    size_t n;std::string line;s>>n;std::getline(s,line);gpus.clear();for(size_t i=0;i<n;++i){std::getline(s,line);
        Require(std::find(gpuCandidates.begin(),gpuCandidates.end(),line)!=gpuCandidates.end(),"tuned GPU outside allowed set");gpus.push_back(line);}
    s>>n;std::getline(s,line);base.clear();for(size_t i=0;i<n;++i){Require(bool(std::getline(s,line)),"truncated execution recipe");base.push_back(line);}
    Require(workerThreads>0&&analyticMeanThreads>0&&workerChunk>=0&&orientationBlock>0&&!gpus.empty(),"invalid tuned resources");
    // Saved tuning describes the previous worker layout. A compatible dispatch
    // override must retain the freshly divided host budget, not multiply the
    // old per-worker allocation by the new number of slots.
    if(!dispatchGpus.empty()){hostBudgetMb=std::min(hostBudgetMb,dispatchHostBudget);ValidateDispatchResources();}
    const auto backend=std::find(base.begin(),base.end(),"--po-trace-backend");if(backend!=base.end()&&std::next(backend)!=base.end()&&*std::next(backend)=="hybrid")kernel.pipeline=true;
    SetMirror(false);return true;
}
void TuneResources(const std::vector<double>& x){
    const unsigned seedSet[2]={8001,8023};const int n=std::min(kernelCount,maxCount);
    std::ostringstream report;report<<"{\"policy\":"<<Quoted(resourcePolicy)<<",\"gpu_candidates\":[";
    std::vector<std::pair<double,std::string>> speed;Matrix reference,scale;bool first=true;
    auto pilotScale=[&](const Job& job){Matrix result=job.matrix;
        if(controls)for(size_t row=0;row<result.size();++row)result[row][0]=job.controls.at(row)[0];
        return result;};
    const auto allowed=gpus;
    if(autoGpus){for(const auto& gpu:allowed){
            const int device=int(std::find(gpuCandidates.begin(),gpuCandidates.end(),gpu)-gpuCandidates.begin());long long freeBytes=0,totalBytes=0;std::string memoryError;
            if(!QueryVisibleGpuMemory(device,freeBytes,totalBytes,memoryError)||freeBytes<512LL*1024*1024){
                report<<(first?"":",")<<"{\"device\":"<<Quoted(gpu)<<",\"skipped\":true,\"free_bytes\":"<<freeBytes<<",\"reason\":"<<Quoted(memoryError.empty()?"less than 512 MiB free":memoryError)<<"}";first=false;continue;}
            gpus={gpu};auto jobs=DirectStage(n,seedSet,2,target,x,controls);double time=jobs[0].seconds+jobs[1].seconds,error=0;
            if(reference.empty()){reference=jobs[0].matrix;scale=pilotScale(jobs[0]);}else error=Difference(jobs[0].matrix,reference,scale);
            report<<(first?"":",")<<"{\"device\":"<<Quoted(gpu)<<",\"seconds\":"<<Number(time)<<",\"all16_difference\":"<<Number(error)<<"}";first=false;
            if(error<std::min(1e-7,eps/1024))speed.emplace_back(time,gpu);}
        Require(!speed.empty(),"no equivalent GPU pilot");std::sort(speed.begin(),speed.end());gpus.clear();
        // Very slow devices lengthen eight-seed stage completion instead of
        // helping it. The dynamic queue assigns the next seed to a free device.
        for(const auto& entry:speed)if(entry.first<=2*speed.front().first)gpus.push_back(entry.second);
    }else gpus=allowed;
    const auto selected=gpus;const int baselineThreads=workerThreads,baselineChunk=workerChunk;
    const int cap=std::max(1,QueryRuntimeResourceSnapshot().physicalCores/int(selected.size()));
    gpus={selected.front()};auto baseline=DirectStage(n,seedSet,2,target,x,controls);
    reference=baseline[0].matrix;scale=pilotScale(baseline[0]);double best=baseline[0].seconds+baseline[1].seconds;int chosenThreads=workerThreads,chosenChunk=workerChunk;
    report<<"],\"resource_candidates\":[{\"threads\":"<<workerThreads<<",\"chunk\":"<<workerChunk<<",\"seconds\":"<<Number(best)<<",\"all16_difference\":0}";
    double nativePilot=0;for(const auto& job:baseline)nativePilot+=job.NativeWall();
    const double forecastNative=nativePilot*double(initial+4LL*minimum)/n*4/std::max(size_t(1),selected.size());
    const double predictedSearch=best*10;
    const bool amortized=forecastNative*.25>=predictedSearch;
    if(resourcePolicy=="auto"&&amortized){
        std::set<int> threads;if(!args.IsCatched("threads"))for(int t:{1,2,4,8,16,32})if(t<=cap)threads.insert(t);
        threads.insert(workerThreads);
        for(int t:threads){workerThreads=t;auto jobs=DirectStage(n,seedSet,2,target,x,controls);const double seconds=jobs[0].seconds+jobs[1].seconds;
            const double error=std::max(Difference(jobs[0].matrix,baseline[0].matrix,scale),Difference(jobs[1].matrix,baseline[1].matrix,scale));
            Log("resource_candidate","{\"threads\":"+Number(workerThreads)+",\"chunk\":"+Number(workerChunk)+",\"seconds\":"+Number(seconds)+",\"all16_difference_scaled_M11\":"+Number(error)+",\"budget\":"+Number(std::min(1e-7,eps/1024))+",\"numerically_accepted\":"+(error<std::min(1e-7,eps/1024)?"true":"false")+"}");
            report<<",{\"threads\":"<<t<<",\"chunk\":"<<workerChunk<<",\"seconds\":"<<Number(seconds)<<",\"all16_difference\":"<<Number(error)<<"}";
            if(error<std::min(1e-7,eps/1024)&&seconds<best*.95){best=seconds;chosenThreads=t;}}
        workerThreads=chosenThreads;
        if(!args.IsCatched("chunk"))for(int chunk:{64,256,1024,4096})if(chunk<=n||chunk==64){workerChunk=chunk;
            auto jobs=DirectStage(n,seedSet,2,target,x,controls);const double seconds=jobs[0].seconds+jobs[1].seconds;
            const double error=std::max(Difference(jobs[0].matrix,baseline[0].matrix,scale),Difference(jobs[1].matrix,baseline[1].matrix,scale));
            Log("resource_candidate","{\"threads\":"+Number(workerThreads)+",\"chunk\":"+Number(workerChunk)+",\"seconds\":"+Number(seconds)+",\"all16_difference_scaled_M11\":"+Number(error)+",\"budget\":"+Number(std::min(1e-7,eps/1024))+",\"numerically_accepted\":"+(error<std::min(1e-7,eps/1024)?"true":"false")+"}");
            report<<",{\"threads\":"<<workerThreads<<",\"chunk\":"<<chunk<<",\"seconds\":"<<Number(seconds)<<",\"all16_difference\":"<<Number(error)<<"}";
            if(error<std::min(1e-7,eps/1024)&&seconds<best*.95){best=seconds;chosenChunk=chunk;}}
    }
    workerThreads=chosenThreads;workerChunk=chosenChunk;auto chosen=DirectStage(n,seedSet,2,target,x,controls);
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
    report<<"],\"phi\":"<<phi<<",\"depth\":"<<target<<",\"theta_rows\":"<<x.size()<<",\"controls\":"<<(controls?"true":"false")
        <<",\"balanced_correction\":"<<(balancedCorrection?"true":"false")
        <<",\"host_memory_budget_mib\":"<<hostBudgetMb<<",\"analytic_mean_threads\":"<<analyticMeanThreads<<",\"selected_threads\":"<<workerThreads<<",\"selected_chunk\":"<<workerChunk<<",\"selected_orientation_block\":"<<orientationBlock
        <<",\"search_amortized\":"<<(amortized?"true":"false")<<",\"forecast_native_production_seconds\":"<<Number(forecastNative)<<",\"forecast_search_seconds\":"<<Number(predictedSearch)
        <<",\"explicit_threads_preserved\":"<<(args.IsCatched("threads")?"true":"false")<<",\"explicit_chunk_preserved\":"<<(args.IsCatched("chunk")?"true":"false")
        <<",\"outside_processing_seconds_per_worker\":"<<Number(overhead)<<",\"processing_seconds_per_orientation\":"<<Number(perOrientation)
        <<",\"processing_wall_measured\":"<<(processingMeasured?"true":"false")
        <<",\"block_model\":\"worker wall minus processing wall (overlapping phases are not summed); aim <=5% outside-processing overhead and <=120s processing work, bounded by 4096..1048576 and max orientations\",\"baseline_threads\":"<<baselineThreads<<",\"baseline_chunk\":"<<baselineChunk<<"}";
    Record(root+"/resource_calibration.json",report.str());
}
void TuneFinalResources(const std::vector<double>& angles){
    if(resourcePolicy!="auto"&&!autoGpus)return;
    std::vector<double> x;const size_t rows=std::min(size_t(17),angles.size());
    for(size_t i=0;i<rows;++i)x.push_back(angles[i*(angles.size()-1)/std::max(size_t(1),rows-1)]);
    // Analytic shared means are defined on the requested grid. Use its nearest
    // actual nodes around backscatter instead of introducing new coordinates.
    for(double angle:{170.,175.,178.,179.,179.5,180.})if(angle>=theta.front()&&angle<=theta.back()){
        auto at=std::lower_bound(theta.begin(),theta.end(),angle);
        if(at==theta.end())--at;
        if(at!=theta.begin()&&std::fabs(*std::prev(at)-angle)<std::fabs(*at-angle))--at;
        x.push_back(*at);}
    std::sort(x.begin(),x.end());x.erase(std::unique(x.begin(),x.end()),x.end());
    Log("final_resource_calibration","{\"phi\":"+Number(phi)+",\"depth\":"+Number(target)+",\"theta_rows\":"+Number(x.size())
        +",\"controls\":"+(controls?"true":"false")+",\"balanced_correction\":"+(balancedCorrection?"true":"false")
        +",\"complete_phi_selection_finished\":true,\"explicit_options_preserved\":true}");
    TuneResources(x);SaveTuning();
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
            worst=std::max(worst,RotationEquivalenceBound(Data(reduced),Data(paired),ratios));}
        const double bound=5*worst;report<<(first?"":",")<<"{\"axis\":"<<Quoted(axis)<<",\"paired_equivalence_bound\":"<<Number(bound)<<",\"safety_multiplier\":5}";first=false;
        if(bound>eps/512)SetRotation(std::string(axis)=="beta"?1:rotationBeta,std::string(axis)=="gamma"?1:rotationGamma);
        else rotationError+=bound;
    }
    report<<"],\"selected_beta_factor\":"<<rotationBeta<<",\"selected_gamma_factor\":"<<rotationGamma<<",\"estimated_error\":"<<Number(rotationError)<<",\"budget\":"<<Number(eps/256)
        <<",\"scope\":\"oriented-facet bijection plus all16 per-block relative equivalence and ratio perturbation bound; no production-convergence gate; fixed axes only\"}";Record(root+"/rotation_verification.json",report.str());
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
    // Fixed production estimator. Its off-balanced correction removes cutoff
    // bias; uncertainty is measured from the complete eight combined replicas.
    tuneCutoffs=false;SetCutoff("off");cutoffError=0;
    balancedCorrection=!go;
    Record(root+"/cutoff_calibration.json","{\"enabled\":true,\"selected\":\"balanced_with_paired_off_correction\",\"physical_target\":\"off\",\"profile_search\":false,\"correction_fraction\":0.125,\"minimum_block\":512,\"supported\":"+std::string(balancedCorrection?"true":"false")+"}");
}
void TuneTracing(const std::vector<double>& x){
#ifdef MBS_GPU_FP32
    const double traceTolerance=std::min(1e-5,eps/1024);
#else
    const double traceTolerance=std::min(1e-7,eps/1024);
#endif
    const bool explicitEnvironment=getenv("MBS_GPU_TRACE_CONVEX_QUEUE")||getenv("MBS_GPU_TRACE_CONVEX_DIRECT")||getenv("MBS_GPU_TRACE_HYBRID");
    if(tracePolicy!="auto"||explicitEnvironment||!config.useGpu||go||controls||config.refractiveImag!=0||config.geometry==GeometryClassification::Nonconvex)return;
    if(std::find(base.begin(),base.end(),"--po-trace-backend")!=base.end())return;
    // A short startup orientation pilot says nothing about the cost of the
    // later adaptive prefixes, especially for incoherent PO. Compare actual
    // native packets after phi/theta calibration even when initial is 512.
    const int n=std::min(maxCount,std::max(512,std::min(2048,initial)));const bool previousPipeline=kernel.pipeline;kernel.pipeline=true;
    const auto originalBase=base;
    if(!args.IsCatched("trace_cache_mb")){RemoveOption("--trace-cache-mib");base.push_back("--trace-cache-mib");base.push_back("0");}
    const unsigned seeds[8]={28001,28019,28027,28031,28051,28057,28069,28081};
    auto select=[&](const std::string& mode){RemoveOption("--po-trace-backend");base.push_back("--po-trace-backend");base.push_back(mode);};
    // Warm the persistent CUDA contexts before comparing execution paths.
    select("cpu");DirectStage(std::min(16,n),seeds,std::min(size_t(8),gpus.size()),target,x,false);
    const double started=Seconds();auto reference=DirectStage(n,seeds,8,target,x,false);const double cpuSeconds=Seconds()-started;
    double best=cpuSeconds;std::string chosen="cpu";std::ostringstream report;
    report<<"{\"count_per_seed\":"<<n<<",\"theta_rows\":"<<x.size()<<",\"phi\":"<<phi<<",\"depth\":"<<target
          <<",\"tolerance\":"<<Number(traceTolerance)<<",\"candidates\":[{\"backend\":\"cpu\",\"wall_seconds\":"<<Number(cpuSeconds)<<",\"error\":0}";
    for(const char* mode:{"resident","hybrid"}){
        select(mode);const double begin=Seconds();auto candidate=DirectStage(n,seeds,8,target,x,false);const double wall=Seconds()-begin;
        double error=0;bool acceptedGpu=false;
        for(size_t i=0;i<candidate.size();++i){error=std::max(error,Difference(candidate[i].matrix,reference[i].matrix,reference[i].matrix));
            for(size_t t=0;t<reference[i].matrix.size();++t){
                const auto& a=candidate[i].matrix[t];const auto& b=reference[i].matrix[t];
                if(!(a[0]>0)){error=std::numeric_limits<double>::infinity();continue;}
                for(int j=1;j<16;++j){const double ar=a[j]/a[0],br=b[j]/b[0];
                    error=std::max(error,std::fabs(ar-br)/std::max(1.,std::max(std::fabs(ar),std::fabs(br))));}
            }
            acceptedGpu=acceptedGpu||FullAutoTraceEvidence::Executed(candidate[i].folder);}
        const bool accepted=acceptedGpu&&error<traceTolerance;
        report<<",{\"backend\":"<<Quoted(mode)<<",\"wall_seconds\":"<<Number(wall)<<",\"error\":"<<Number(error)<<",\"gpu_executed\":"<<(acceptedGpu?"true":"false")<<",\"accepted\":"<<(accepted?"true":"false")<<"}";
        if(accepted&&wall<best*.95){best=wall;chosen=mode;}
    }
    base=originalBase;select(chosen);if(chosen=="cpu")kernel.pipeline=previousPipeline;
    report<<"],\"selected\":"<<Quoted(chosen)<<",\"speedup\":"<<Number(cpuSeconds/best)<<",\"sampling_unchanged\":true}";
    Record(root+"/trace_backend_calibration.json",report.str());SaveTuning();
}

void PrepareTuning(const std::vector<double>& angles){
    if(RestoreTuning()){
        bool explicitCutoff=args.IsCatched("r");for(const char* key:{"cutoff_profile","beam_cutoff","beam_cutoff_j","beam_cutoff_area","beam_cutoff_importance","trace_cutoff","trace_cutoff_j","trace_cutoff_area","trace_cutoff_importance"})explicitCutoff=explicitCutoff||args.IsCatched(key);
        // Old incoherent checkpoints used direct-off. Retain that estimator;
        // only a saved successful opt-in may restore the new corrected path.
        const bool savedIncoherentCorrection=!coherent&&Exists(root+"/cutoff_calibration.json")
            &&CompactJson(Read(root+"/cutoff_calibration.json")).find("\"supported\":true")!=std::string::npos;
        balancedCorrection=(coherent||savedIncoherentCorrection)&&!go&&cutoffPolicy!="off"&&!explicitCutoff&&selectedCutoff=="off";
        Log("cutoff_estimator_restore","{\"balanced_with_paired_off_correction\":"+std::string(balancedCorrection?"true":"false")+",\"saved_physical_target\":"+Quoted(selectedCutoff)+",\"old_completed_blocks_preserved\":true,\"profile_search\":false}");
        return;
    }
    Record(root+"/automatic_cutoff_validation.json","{\"enabled\":false}");
    std::vector<double> x;const size_t rows=std::min(size_t(17),angles.size());for(size_t i=0;i<rows;++i)x.push_back(angles[i*(angles.size()-1)/(rows-1)]);
    const bool oldControls=controls;controls=false;SetMirror(false);phi=fixedPhi?fixedPhi:std::max(minPhi,std::min(maxPhi,4));
    if(autoGpus&&!go){
        // Keep the original resource probe's availability guard before any
        // final phi workers run, without benchmarking the wrong phi recipe.
        std::vector<std::string> available;
        for(const auto& gpu:gpus){const int device=int(std::find(gpuCandidates.begin(),gpuCandidates.end(),gpu)-gpuCandidates.begin());
            long long freeBytes=0,totalBytes=0;std::string reason;
            if(QueryVisibleGpuMemory(device,freeBytes,totalBytes,reason)&&freeBytes>=512LL*1024*1024)available.push_back(gpu);
            else Log("gpu_availability_skip","{\"device\":"+Quoted(gpu)+",\"free_bytes\":"+Number(freeBytes)+",\"reason\":"+Quoted(reason.empty()?"less than 512 MiB free":reason)+"}");}
        Require(!available.empty(),"no GPU with enough free memory for calibration");gpus=available;
    }
    if((resourcePolicy=="auto"||autoGpus)&&go)TuneResources(x);
    else if(resourcePolicy=="auto"||autoGpus)Record(root+"/resource_calibration.json","{\"policy\":"+Quoted(resourcePolicy)+",\"deferred_until_final_phi\":true}");
    else Record(root+"/resource_calibration.json","{\"policy\":\"off\"}");
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
    if(!pass){SetCutoff("off");balancedCorrection=!go;cutoffError=0;depthSharing=!go&&On("fullauto_depth_share");depthVerified=false;weights.clear();polarWeights.clear();weightTheta.clear();SetMirror(mirror);SaveTuning();}
    SaveTuning();return pass;
}
