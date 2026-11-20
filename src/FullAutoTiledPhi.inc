bool measuringPhiMemory=false;
double phiPeakKiB=0;

void ObservePhiMemory(const Job& job){
    if(job.pairedDepth)return;
    const double kib=FullAutoPhiTileMemory::Peak(Read(job.folder+"/stdout.log"));
    phiPeakKiB=std::max(phiPeakKiB,kib);
}

// Scope is the azimuth pilot only. Final production keeps its own seeds and evaluator.
struct PhiPilotTable {
    std::map<int,std::vector<Job>> fit, validation;
    std::map<int,double> standaloneSeconds;
};

bool CanTilePhi(const std::vector<int>& rules) const {
    return coherent && !go && config.useGpu && !phiCutoffPilot
        && args.IsCatched("rosette_joined") && rules.size()>=3 && !fixedPhi && pilot>=16
        && WorkerGpus().size()==1 && workerPolicy=="persistent" && memoryMode && checkpointInterval==0
        && !args.IsCatched("trace_cache_mb") && !args.IsCatched("chunk") && workerChunk==0 && !mirror && (!controls||sharedMeans);
}

void FillTiledPhi(std::map<int,std::vector<Job>>& table,const unsigned* seedSet,int first,int end,
                  const std::vector<int>& rules,int depth,const std::vector<double>& x) {
    // Restrict the temporary cache to the existing host memory allowance. The worker
    // also clamps the actual cache to the live per-process budget and available RAM.
    const char* cap=std::getenv("MBS_HOST_MEM_BUDGET_MB");
    const long long allowance=cap?std::min(hostBudgetMb,std::atoll(cap)):hostBudgetMb;
    const auto limits=TraceCachePolicy::Resolve(true,true,512,QueryRuntimeMemorySnapshot().ramAvailableKb,true,allowance,false);
    const int initialTile=phiPeakKiB>0?FullAutoPhiTileMemory::Fit(limits.hostMiB,phiPeakKiB,16):4;
    int tileSize=initialTile;
    unsigned long long hitCount=0,missCount=0;
    const int savedPhi=phi;
    struct RestorePhi {int& phi;int old;~RestorePhi(){phi=old;}} restorePhi{phi,savedPhi};
    KernelCalibrationScope saved(*this);
    RemoveOption("--trace-cache-mib");
    base.push_back("--trace-cache-mib");base.push_back("512");
    const std::string session=Number(getpid())+" tile "+Number(++kernelTrialSequence);
    size_t serial=0;
    for(int seedIndex=first;seedIndex<end;++seedIndex){
        for(int rule:rules){
            auto& result=table[rule];result.resize(8);
            Job& out=result[seedIndex];
            out.seed=seedSet[seedIndex];out.count=pilot;out.depth=depth;out.phi=rule;
            out.angles=x;out.inputAngles=x;out.matrix.assign(x.size(),Mueller{});out.hasControls=controls;
            if(controls)out.controls.assign(x.size(),std::array<double,5>{});
        }
        for(int offset=0;offset<pilot;){
            const int length=std::min(tileSize,pilot-offset);
            const double weight=double(length)/pilot;
            unsigned tileHits=0,tileMisses=0;
            for(int rule:rules){
                phi=rule;kernelTrialNonce=Hash(session+" "+Number(++serial));
                Job part=DirectStage(length,seedSet+seedIndex,1,depth,x,false,false,offset).front();
                if(rule!=rules.front()){
                    const std::string log=Read(part.folder+"/stdout.log");
                    const std::string mark="Prepared trace reuse: hits=";
                    const auto at=log.rfind(mark);unsigned hits=0,misses=0;
                    if(at!=std::string::npos && std::sscanf(log.c_str()+at+mark.size(),"%u misses=%u",&hits,&misses)==2){
                        tileHits+=hits;tileMisses+=misses;
                    }
                }
                RemoveKernelTrial(part.folder);
                Job& out=table[rule][seedIndex];
                Require(part.seed==out.seed && part.offset==offset && part.phi==rule &&
                        part.count==length && part.depth==depth && part.angles.size()==x.size() &&
                        part.matrix.size()==x.size() && part.hasControls==controls &&
                        (!controls||part.controls.size()==x.size()),"incompatible native azimuth tile");
                for(size_t row=0;row<x.size();++row){
                    Require(std::fabs(part.angles[row]-x[row])<1e-10,"azimuth tile theta differs");
                    for(size_t k=0;k<16;++k)out.matrix[row][k]+=weight*part.matrix[row][k];
                    if(controls)for(size_t k=0;k<5;++k)out.controls[row][k]+=weight*part.controls[row][k];
                }
                out.seconds+=part.seconds;
                out.controlSetup=std::max(out.controlSetup,part.controlSetup);
                out.controlRefinement=std::max(out.controlRefinement,part.controlRefinement);
            }
            hitCount+=tileHits;missCount+=tileMisses;
            if(tileHits<3*tileMisses && tileSize>1)tileSize/=2;
            offset+=length;
        }
    }
    Log("phi_tile_reuse","{\"seed_start\":"+Number(first)+",\"hits\":"+Number(hitCount)+",\"misses\":"+Number(missCount)
        +",\"safe_kib_per_orientation\":"+Number(phiPeakKiB)+",\"cache_mib\":"+Number(limits.hostMiB)
        +",\"first_tile\":"+Number(initialTile)+",\"final_tile\":"+Number(tileSize)+"}");
}

bool BuildTiledPhi(const std::vector<int>& rules,int depth,const std::vector<double>& x,PhiPilotTable& t) {
    const int savedPhi=phi;
    struct RestorePhi {int& phi;int old;~RestorePhi(){phi=old;}} restore{phi,savedPhi};
    const double start=Seconds();
    struct MeasureScope {
        bool& active;
        explicit MeasureScope(bool& flag):active(flag){active=true;}
        ~MeasureScope(){active=false;}
    };
    phiPeakKiB=0;
    // The first two validation scrambles supply their matrices and a cold
    // standalone timing for each rule. Warm timings cannot price a later cold solve.
    const auto submittedBefore=nativeSubmitted;
    double baselineSeconds=0;
    {
        MeasureScope measure(measuringPhiMemory);
        for(int rule:rules){
            phi=rule;
            auto refs=Stage(pilot,validation,2,depth,x,false);
            double seconds=0;
            for(const auto& j:refs)seconds+=std::max(1e-6,j.seconds-j.controlSetup);
            baselineSeconds+=refs.back().seconds;
            t.standaloneSeconds[rule]=4*seconds;
            t.validation[rule]=std::move(refs);
        }
    }
    if(nativeSubmitted==submittedBefore){
        Log("phi_reuse_decision","{\"use_tiled_phi\":false,\"reason\":\"cold calibration samples already cached\"}");
        return false;
    }
    // A private replay on the SAME scramble measures the actual cost of small
    // jobs, including packing and controller I/O. Check the optical result
    // before choosing this schedule for the six remaining validation seeds.
    std::map<int,std::vector<Job>> replay;
    const double replayStart=Seconds();
    FillTiledPhi(replay,validation,1,2,rules,depth,x);
    const double replaySeconds=Seconds()-replayStart;
    bool consistent=true;
    double difference=0,ratioDifference=0;
    for(int rule:rules){
        const Job& a=t.validation.at(rule).at(1);
        const Job& b=replay.at(rule).at(1);
        for(size_t row=0;row<x.size();++row){
            const double scale=a.matrix[row][0],replayScale=b.matrix[row][0];
            if(!(scale>0) || !(replayScale>0) || !std::isfinite(scale) || !std::isfinite(replayScale)){
                consistent=false;continue;
            }
            auto compare=[&](double left,double right){
                const double delta=std::fabs(left-right)/scale;
                if(!std::isfinite(delta))consistent=false;
                else difference=std::max(difference,delta);
            };
            for(int j=0;j<16;++j){
                compare(a.matrix[row][j],b.matrix[row][j]);
                const double delta=std::fabs(a.matrix[row][j]/scale-b.matrix[row][j]/replayScale);
                if(!std::isfinite(delta))consistent=false;
                else ratioDifference=std::max(ratioDifference,delta);
            }
            if(controls)for(int j=0;j<5;++j)compare(a.controls[row][j],b.controls[row][j]);
        }
    }
    const bool use=consistent && std::max(difference,ratioDifference)<=eps/100 && 7*replaySeconds<6*baselineSeconds;
    Log("phi_reuse_decision","{\"depth\":"+Number(depth)+",\"orientations\":"+Number(pilot)
        +",\"baseline_seconds\":"+Number(baselineSeconds)+",\"replay_seconds\":"+Number(replaySeconds)
        +",\"finite\":"+(consistent?"true":"false")+",\"max_difference_over_M11\":"+Number(difference)
        +",\"max_normalized_difference\":"+Number(ratioDifference)+",\"use_tiled_phi\":"+(use?"true":"false")+"}");
    if(!use)return false;
    FillTiledPhi(t.validation,validation,2,8,rules,depth,x);
    // With score-floor screening, fit only rules that survive the screen.
    if(controls&&!phiScoreFloor)FillTiledPhi(t.fit,training,0,8,rules,depth,x);
    Log("tiled_phi_calibration","{\"rules\":"+Number(rules.size())+",\"depth\":"+Number(depth)
        +",\"count_per_seed\":"+Number(pilot)+",\"standalone_timing_seeds\":2,\"validation_seeds\":8"
        +",\"training_seeds\":"+Number(t.fit.empty()?0:8)+",\"wall_seconds\":"+Number(Seconds()-start)
        +",\"accuracy_accepted\":false}");
    return true;
}
