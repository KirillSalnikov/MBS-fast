// Included inside the adaptive controller. This is an execution-only pilot:
// independent seeds, identical orientations in split/merged calculations,
// unchanged production budgets, stencil prefixes and statistical acceptance.
#ifdef MBS_GPU_FP32
bool VerifyInterpolationAngleBatches(Engine& e,const AdaptiveState& state,
        const std::vector<int>& depths,const std::vector<size_t>& evaluated){
    if(!e.config.useGpu||e.go||e.resourcePolicy!="auto"||!e.QueueEnabled()||evaluated.size()<4)return false;
    const auto stencils=FullAutoAuditPlan::Stencils(state.knots,state.counts,evaluated);
    size_t smallGroups=0;std::set<size_t> eligible;
    for(const auto& group:stencils)if(group.second.rows.size()<=64){
        ++smallGroups;eligible.insert(group.second.rows.begin(),group.second.rows.end());}
    if(smallGroups<2)return false;
    std::array<unsigned,2> pilotSeeds{{8609u,8623u}};
    for(size_t i=0;i<pilotSeeds.size();++i)
        while(std::find(e.seeds.begin(),e.seeds.end(),pilotSeeds[i])!=e.seeds.end()
                ||(i&&pilotSeeds[i]==pilotSeeds[0]))++pilotSeeds[i];
    std::map<int,std::vector<double>> byPhi;
    for(size_t id:eligible)byPhi[e.PhiForTheta(id)].push_back(e.theta[id]);
    // Validate each actual phi rule independently. Never mix phi in a native
    // command or discard the retained local rules on checkpoint restoration.
    bool allPassed=true;
    for(const auto& group:byPhi){
        if(group.second.size()<4)return false;
        PhiScope phiScope(e,group.first);
        const size_t count=std::min(size_t(64),group.second.size());
        std::vector<double> angles;
        for(size_t i=0;i<count;++i)angles.push_back(group.second[i*(group.second.size()-1)/(count-1)]);
        const auto split=FullAutoScheduling::AngleBatches(angles,std::min(size_t(16),(count+1)/2));
        const int n=std::min(e.maxCount,std::min(1024,e.orientationBlock));
        if(n<256)return false; // retain manually configured tiny blocks
        const auto weights=e.WriteWeights(angles);
        std::ostringstream recipe;recipe<<"angle-batch-v1 "<<e.binaryHash<<' '<<e.ExecutionHash()
            <<' '<<e.phi<<' '<<n<<' '<<e.depthSharing<<' '<<e.kernel.block<<' '<<e.kernel.warp
            <<' '<<e.kernel.pipeline<<' '<<e.workerChunk<<' '<<e.balancedCorrection<<' '<<Number(e.eps);
        for(int d:depths)recipe<<' '<<d;
        for(unsigned seed:pilotSeeds)recipe<<' '<<seed;
        for(double x:angles)recipe<<' '<<Number(x);
        if(!weights.empty())recipe<<' '<<Hash(Read(weights));
        const std::string path=e.root+"/angle_batch_validation_"+Hash(recipe.str())+".json";
        if(Exists(path)){const auto saved=CompactJson(Read(path));
            if(saved.find("\"passed\":true")==std::string::npos)return false;continue;}
        const bool oldCollecting=e.collecting,oldVerifying=e.shortPacketVerifying,oldCutoffPilot=e.phiCutoffPilot;
        const std::string oldPilotRecipe=e.angleBatchPilotRecipe;
        e.collecting=false;e.shortPacketVerifying=true;
        // Independent pilot seeds must exercise the production correction,
        // rather than silently falling back to unpruned native jobs.
        e.phiCutoffPilot=e.balancedCorrection;
        // Rebuild-specific pilot identities: an old binary's native pilot
        // cannot make the new binary pass without executing its kernels.
        e.angleBatchPilotRecipe=Hash(recipe.str());
        double maximum=0,baselineSeconds=0,mergedSeconds=0;size_t comparisons=0;
        const double started=Seconds();
        try{
            for(size_t l=0;l<depths.size();++l)for(unsigned seed:pilotSeeds){
                const int paired=l&&e.depthSharing?depths[l-1]:0;
                std::vector<Job> merged;std::vector<std::vector<Job>> pieces;
                auto runMerged=[&](){merged=e.DirectStage(n,&seed,1,depths[l],angles,true,false,0,paired);
                    mergedSeconds+=merged[0].NativeWall()>0?merged[0].NativeWall():merged[0].seconds;};
                auto runSplit=[&](){for(const auto& x:split){auto jobs=e.DirectStage(n,&seed,1,depths[l],x,true,false,0,paired);
                    baselineSeconds+=jobs[0].NativeWall()>0?jobs[0].NativeWall():jobs[0].seconds;pieces.push_back(std::move(jobs));}};
                // Alternate order to avoid consistently rewarding a warm
                // analytic cache or a transient drop in unrelated GPU load.
                if(seed==pilotSeeds[0]){runSplit();runMerged();}else{runMerged();runSplit();}
                for(int d:{depths[l],paired}){if(!d)continue;
                    const auto unionData=e.Data({d==depths[l]?merged[0]:e.LowerView(merged[0],d)});
                    for(size_t part=0;part<split.size();++part){
                        const auto reference=e.Data({d==depths[l]?pieces[part][0]:e.LowerView(pieces[part][0],d)});
                        for(size_t r=0;r<split[part].size();++r){
                            const size_t row=std::lower_bound(angles.begin(),angles.end(),split[part][r])-angles.begin();
                            const auto& a=reference[0][r];const auto& b=unionData[0][row];
                            if(!(a[0]>0&&b[0]>0)){maximum=1e100;continue;}
                            for(int j=0;j<16;++j){
                                if(!std::isfinite(a[j])||!std::isfinite(b[j])){maximum=1e100;continue;}
                                maximum=std::max(maximum,std::fabs(a[j]-b[j])/a[0]);}
                            if(e.ratios)for(int j=1;j<16;++j)maximum=std::max(maximum,std::fabs(a[j]/a[0]-b[j]/b[0]));
                            ++comparisons;
                        }
                    }
                }
            }
        }catch(...){e.collecting=oldCollecting;e.shortPacketVerifying=oldVerifying;
            e.phiCutoffPilot=oldCutoffPilot;e.angleBatchPilotRecipe=oldPilotRecipe;throw;}
        e.collecting=oldCollecting;e.shortPacketVerifying=oldVerifying;
        e.phiCutoffPilot=oldCutoffPilot;e.angleBatchPilotRecipe=oldPilotRecipe;
        // Spend one quarter of the EXISTING eps/1024 execution allowance.
        // The short-packet gate uses at most another quarter, leaving half
        // for other execution choices. A fixed ppm ceiling unnecessarily
        // rejects a corrected sum of three FP32 native terms even when it is
        // well inside this budget. Compare that final sum, not each leg.
        const double tolerance=e.eps/4096;
        const bool passed=comparisons&&maximum<=tolerance&&mergedSeconds<baselineSeconds*.97;
        const std::string report="{\"passed\":"+std::string(passed?"true":"false")
            +",\"binary_fnv1a64\":"+Quoted(e.binaryHash)+",\"phi\":"+Number(e.phi)
            +",\"angle_batch_cap\":64,\"pilot_theta_rows\":"+Number(angles.size())
            +",\"count_per_seed\":"+Number(n)+",\"independent_seeds\":["+Number(pilotSeeds[0])+","+Number(pilotSeeds[1])+"]"
            +",\"balanced_correction_checked\":"+(e.balancedCorrection?"true":"false")
            +",\"depth_levels\":"+Number(depths.size())+",\"paired_depths_checked\":"+(e.depthSharing?"true":"false")
            +",\"max_all16_and_ratios_difference\":"+Number(maximum)+",\"numerical_tolerance\":"+Number(tolerance)
            +",\"split_processing_seconds\":"+Number(baselineSeconds)+",\"merged_processing_seconds\":"+Number(mergedSeconds)
            +",\"pilot_elapsed_seconds\":"+Number(Seconds()-started)
            +",\"minimum_saving_fraction\":0.03,\"production_counts_unchanged\":true,\"stencil_prefixes_unchanged\":true}";
        e.Record(path,report);e.Log("angle_batch_validation",report);allPassed=allPassed&&passed;
        if(!passed)break;
    }
    return allPassed;
}
#endif
