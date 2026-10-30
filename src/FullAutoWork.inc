// Included inside FullAutoAdaptive.inc. Scheduling forecasts only: never an
// accuracy criterion. Times are summed native worker seconds, not elapsed time.
struct WorkForecast {double seconds=0;size_t jobs=0;bool profiled=true;};
class AdditionalWork {
    Engine& e;std::array<unsigned,8> seeds;bool ignoreCoverage,productionOnly,uniformPhi;
    struct Profile {double trace=0,diffraction=0,setup=0,overlap=0;size_t observations=0;};
    // Keep native timing recipes separate from physical cache identities.
    // Corrected parents are assembled; their actual workers have child hashes.
    // role: 0 physical/off command, 1 balanced child, 2 off correction child.
    using Depth=std::tuple<int,int,int>;
    using Slot=std::array<uint64_t,5>; // N, offset, depth, paired depth, seed
    std::map<Depth,Profile> profiles;
    struct Coverage {std::vector<const Job*> candidates;std::unique_ptr<double[]> angles;size_t size=0;bool indexed=false;};
    struct CoverageBudget {size_t bytes=0,limit;explicit CoverageBudget(size_t maximum):limit(maximum){}};
    std::shared_ptr<CoverageBudget> coverageBudget;
    std::map<Slot,Coverage> coverage;
    struct RequestLess {
        bool operator()(const std::vector<Engine::WorkRequest>& a,const std::vector<Engine::WorkRequest>& b)const{
            for(size_t i=0;i<std::min(a.size(),b.size());++i){
                const auto x=std::tie(a[i].count,a[i].offset,a[i].depth,a[i].pairedDepth);
                const auto y=std::tie(b[i].count,b[i].offset,b[i].depth,b[i].pairedDepth);
                if(x!=y)return x<y;
                if(a[i].angles!=b[i].angles)return a[i].angles<b[i].angles;
            }return a.size()<b.size();
        }
    };
    // These keys are process-local. Avoid formatting every theta as text for
    // each candidate allocation; keep the original request/angle order.
    std::map<std::vector<Engine::WorkRequest>,WorkForecast,RequestLess> forecasts;
    struct LeadingCoverage {int count=0;bool stopped=false;};
    std::map<std::tuple<int,int,int,double>,LeadingCoverage> leadingCoverage;
    std::map<int,std::unique_ptr<AdditionalWork>> phiWork;
    void NativeCost(WorkForecast& result,int count,int depth,int paired,size_t rows,int role)const{
        auto profile=profiles.find(Depth{depth,paired,role});
        // The off child runs the same physical command. Reuse its timings,
        // never its samples, for a small uncorrected production block.
        if(profile==profiles.end()&&role==0)profile=profiles.find(Depth{depth,paired,2});
        ++result.jobs;
        if(profile==profiles.end()){result.profiled=false;return;}
        const auto& p=profile->second;const double trace=count*p.trace,diffraction=count*p.diffraction*rows;
        result.seconds+=p.setup+trace+diffraction-p.overlap*std::min(trace,diffraction);
    }
    void CorrectedCost(WorkForecast& result,int count,int depth,int paired,size_t rows)const{
        const int m=std::max(64,count/8);
        // Conservatively charge all three legs, including B_m when a verified
        // embedded prefix could save it. No native pilot or filesystem scan
        // is allowed inside a metadata-only forecast. Legacy paired parents
        // used five legs; retain that estimate if paired balanced timings are
        // absent instead of silently treating their base work as free.
        for(int n:{count,m}){
            if(paired&&!profiles.count(Depth{depth,paired,1})){
                NativeCost(result,n,depth,0,rows,1);NativeCost(result,n,paired,0,rows,1);
            }else NativeCost(result,n,depth,paired,rows,1);
        }
        NativeCost(result,m,depth,paired,rows,2);
    }
    Slot Key(int n,int offset,int depth,int paired,unsigned seed)const{return {{uint64_t(n),uint64_t(offset),uint64_t(depth),uint64_t(paired),seed}};}
    bool Covered(int n,int offset,int depth,int paired,unsigned seed,double angle){
        if(ignoreCoverage)return false;
        const auto key=Key(n,offset,depth,paired,seed);auto found=coverage.find(key);
        if(found==coverage.end()){
            Coverage value;size_t rawAngles=0;bool fits=true;
            const size_t available=(coverageBudget->limit-coverageBudget->bytes)/sizeof(double);
            for(const auto& item:e.SamplingCandidates(n,seed,depth,offset)){
                const Job& job=*item.second;
                if(e.SameSampling(job,n,seed,depth,offset)&&(!paired||job.pairedDepth==paired)){
                    value.candidates.push_back(&job);
                    if(fits){if(job.inputAngles.size()>available-rawAngles)fits=false;else rawAngles+=job.inputAngles.size();}
                }
            }
            // Prebound the unmerged input, including duplicates, before any
            // angle allocation. A full budget retains the original exact
            // candidate lookup; neither case changes saved-work coverage.
            if(fits){
                if(rawAngles){
                    value.angles.reset(new double[rawAngles]);double* next=value.angles.get();
                    for(const Job* job:value.candidates)next=std::copy(job->inputAngles.begin(),job->inputAngles.end(),next);
                    std::sort(value.angles.get(),next);
                    value.size=std::unique(value.angles.get(),next)-value.angles.get();
                }
                coverageBudget->bytes+=rawAngles*sizeof(double);value.indexed=true;
                std::vector<const Job*>().swap(value.candidates);
            }
            found=coverage.emplace(key,std::move(value)).first;
        }
        if(found->second.indexed)return found->second.size&&std::binary_search(found->second.angles.get(),found->second.angles.get()+found->second.size,angle);
        for(const Job* job:found->second.candidates)if(std::binary_search(job->inputAngles.begin(),job->inputAngles.end(),angle))return true;
        return false;
    }
    int LeadingCovered(const Engine::WorkRequest& r,int block){
        if(ignoreCoverage)return 0;
        int common=r.count-r.count%block;
        if(!common)return 0;
        for(double angle:r.angles){
            const auto key=std::make_tuple(block,r.depth,r.pairedDepth,angle);
            auto found=leadingCoverage.find(key);
            // Bounded accelerator only; falling back scans the same leaves.
            if(found==leadingCoverage.end()){
                // Charge a conservative map-node bound to the same pool as
                // angular coverage (also shared by phi-specific planners).
                if(leadingCoverage.size()>=32768||coverageBudget->limit-coverageBudget->bytes<128)return 0;
                bool complete=true;for(unsigned seed:seeds)complete=complete&&Covered(block,0,r.depth,r.pairedDepth,seed,angle);
                if(!complete||coverageBudget->limit-coverageBudget->bytes<128)return 0;
                LeadingCoverage initial;initial.count=block;
                found=leadingCoverage.emplace(key,initial).first;
                coverageBudget->bytes+=128;
            }
            auto& saved=found->second;
            while(!saved.stopped&&saved.count<common){
                bool complete=true;for(unsigned seed:seeds)complete=complete&&Covered(block,saved.count,r.depth,r.pairedDepth,seed,angle);
                if(complete)saved.count+=block;else saved.stopped=true;
            }
            common=std::min(common,saved.count);if(!common)break;
        }return common;
    }

public:
    AdditionalWork(Engine& engine,const unsigned* seedSet,bool ignoreSaved=false,bool productionProfilesOnly=false,bool uniform=false,size_t coverageLimit=8*1024*1024):e(engine),ignoreCoverage(ignoreSaved),productionOnly(productionProfilesOnly),uniformPhi(uniform),coverageBudget(new CoverageBudget(coverageLimit)){
        std::copy(seedSet,seedSet+8,seeds.begin());const auto recipe=e.ExecutionHash();
        const auto balancedRecipe=Hash(recipe+" balanced-child-v1"),offRecipe=Hash(recipe+" off-child-v1");
        for(const auto& item:e.cache){const Job& j=item.second;
            if(productionProfilesOnly&&std::find(e.seeds.begin(),e.seeds.end(),j.seed)==e.seeds.end())continue;
            const int role=j.executionHash==recipe?0:j.executionHash==balancedRecipe?1:j.executionHash==offRecipe?2:-1;
            if(j.assembled||j.paired||j.count<=0||j.inputAngles.empty()||role<0||j.phi!=e.phi
                ||j.block!=e.kernel.block||j.pipeline!=e.kernel.pipeline||j.warp!=e.kernel.warp||j.mirrored!=e.mirror||j.hasControls!=e.controls)continue;
            if(j.seconds<=0||!std::isfinite(j.seconds))continue;
            auto& p=profiles[{j.depth,j.pairedDepth,role}];
            // CPU/older logs can lack a split. Their measured total remains a
            // useful conservative trace profile, with no invented theta gain.
            if(j.traceSeconds+j.diffractionSeconds>0){p.trace+=j.traceSeconds/j.count;p.diffraction+=j.diffractionSeconds/j.count/j.inputAngles.size();
                const double overlap=j.OverlapFraction();p.overlap+=overlap;
                p.setup+=std::max(0.,j.seconds-(j.traceSeconds+j.diffractionSeconds-overlap*std::min(j.traceSeconds,j.diffractionSeconds)));}
            else p.trace+=j.seconds/j.count;
            ++p.observations;
        }
        for(auto& item:profiles){auto& p=item.second;p.trace/=p.observations;p.diffraction/=p.observations;p.setup/=p.observations;p.overlap/=p.observations;}
    }
    size_t CoverageIndexBytes()const{return coverageBudget->bytes;}
    size_t CoverageIndexLimit()const{return coverageBudget->limit;}
    int CompletePrefix(int maximum,int depth,int paired,double angle){
        if(!uniformPhi&&!e.thetaPhi.empty()){
            const auto at=std::lower_bound(e.theta.begin(),e.theta.end(),angle);
            Require(at!=e.theta.end()&&*at==angle,"prefix coverage requires a requested theta");
            const int phi=e.PhiForTheta(at-e.theta.begin());PhiScope scope(e,phi);auto& work=phiWork[phi];
            if(!work){work.reset(new AdditionalWork(e,seeds.data(),ignoreCoverage,productionOnly,true));work->coverageBudget=coverageBudget;}
            return work->CompletePrefix(maximum,depth,paired,angle);
        }
        if(ignoreCoverage||!e.incremental||!e.angularCache||(e.controls&&!e.sharedMeans))return 0;
        const int block=e.orientationBlock;if(maximum<block)return 0;
        // Same exact coverage predicate as Forecast, including all eight seeds.
        // Unlike the optional accelerator, restoration must also work with a
        // full metadata budget. No payload reads, no completion from children.
        int complete=LeadingCovered({maximum,0,depth,paired,{angle}},block);
        for(int offset=complete;offset<=maximum-block;offset+=block){
            bool covered=true;for(unsigned seed:seeds)covered=covered&&Covered(block,offset,depth,paired,seed,angle);
            if(!covered)break;
            complete=offset+block;
        }
        int power=0;for(int n=1;n<=complete;){power=n;if(n>complete/2)break;n*=2;}return power;
    }
    WorkForecast Forecast(const std::vector<Engine::WorkRequest>& requests){
        if(!uniformPhi&&!e.thetaPhi.empty()){
            WorkForecast total;const auto groups=PhiRequests(e,requests);
            for(const auto& group:groups){PhiScope scope(e,group.first);auto& work=phiWork[group.first];if(!work){work.reset(new AdditionalWork(e,seeds.data(),ignoreCoverage,productionOnly,true));work->coverageBudget=coverageBudget;}const auto part=work->Forecast(group.second);
                total.seconds+=part.seconds;total.jobs+=part.jobs;total.profiled=total.profiled&&part.profiled;}return total;
        }
        auto saved=forecasts.find(requests);if(saved!=forecasts.end())return saved->second;
        WorkForecast result;std::map<Slot,std::vector<double>> planned;
        auto covered=[&](const Engine::WorkRequest& r,unsigned seed,double angle){
            if(Covered(r.count,r.offset,r.depth,r.pairedDepth,seed,angle))return true;
            const auto it=planned.find(Key(r.count,r.offset,r.depth,r.pairedDepth,seed));
            return it!=planned.end()&&std::binary_search(it->second.begin(),it->second.end(),angle);
        };
        auto missingFor=[&](const Engine::WorkRequest& r,bool requireAny){
            bool any=true;std::vector<double> missing;
            if(requireAny)for(unsigned seed:seeds){bool found=false;for(double q:r.angles)found=found||covered(r,seed,q);any=any&&found;}
            for(double q:r.angles){bool complete=true;for(unsigned seed:seeds)complete=complete&&covered(r,seed,q);if(!complete)missing.push_back(q);}
            if(requireAny&&!any)missing=r.angles;
            if(missing.size()==1&&r.angles.size()>1){missing.push_back(r.angles.front()==missing.front()?r.angles.back():r.angles.front());std::sort(missing.begin(),missing.end());}return missing;
        };
        auto remember=[&](int n,int offset,int depth,int paired,unsigned seed,const std::vector<double>& angles){
            auto& rows=planned[Key(n,offset,depth,paired,seed)];rows.insert(rows.end(),angles.begin(),angles.end());std::sort(rows.begin(),rows.end());rows.erase(std::unique(rows.begin(),rows.end()),rows.end());
        };
        auto leaf=[&](const Engine::WorkRequest& r){
            auto missing=missingFor(r,!e.QueueEnabled()&&!r.pairedDepth);if(missing.empty())return;
            for(unsigned seed:seeds){
                // DirectStage can reuse a completed native command for a seed
                // even when another seed of the batch was interrupted.
                bool exact=false;if(!ignoreCoverage&&!e.controls)for(const auto& item:e.SamplingCandidates(r.count,seed,r.depth,r.offset)){
                    const Job& j=*item.second;exact=exact||(!j.assembled&&j.pairedDepth==r.pairedDepth&&j.inputAngles==missing&&e.SameSampling(j,r.count,seed,r.depth,r.offset));}
                if(!exact){
                    const bool corrected=e.balancedCorrection&&r.count>=512&&r.depth>=6
                        &&(e.phiCutoffPilot||std::find(e.seeds.begin(),e.seeds.end(),seed)!=e.seeds.end());
                    if(corrected)CorrectedCost(result,r.count,r.depth,r.pairedDepth,missing.size());
                    else NativeCost(result,r.count,r.depth,r.pairedDepth,missing.size(),0);
                }
                // Prefetch prepares paired and standalone jobs separately;
                // LowerView becomes available only during the later Stage.
                if(!e.QueueEnabled()){
                    for(int depth:{r.depth,r.pairedDepth})if(depth)remember(r.count,r.offset,depth,0,seed,missing);
                    remember(r.count,r.offset,r.depth,r.pairedDepth,seed,missing);
                }
            }
        };
        int maximum=0;for(const auto& r:requests)maximum=std::max(maximum,r.count);
        if(e.QueueEnabled()){
            const int block=e.incremental?e.orientationBlock:maximum;
            // A complete leading prefix contributes no new jobs or seconds.
            // Prove coverage for ALL eight seeds at the native block size;
            // then omit only those zero-cost leaves. Keep the remaining leaf
            // groups, padding rules and summation order exactly as before.
            int first=maximum;
            if(e.incremental){for(const auto& r:requests){const int covered=LeadingCovered(r,block);if(covered<r.count)first=std::min(first,covered);}}
            else first=0;
            for(int offset=first;offset<maximum;offset+=block){std::map<std::string,Engine::WorkRequest> leaves;
                // Retain lexical leaf ordering: it fixes floating-point cost
                // accumulation order. All three fields are exact integers.
                for(const auto& r:requests)if(offset<r.count){const int n=std::min(block,r.count-offset);const auto key=std::to_string(n)+" "+std::to_string(r.depth)+" "+std::to_string(r.pairedDepth);
                    auto it=leaves.find(key);if(it==leaves.end())it=leaves.emplace(key,Engine::WorkRequest{n,offset,r.depth,r.pairedDepth,{}}).first;
                    it->second.angles.insert(it->second.angles.end(),r.angles.begin(),r.angles.end());}
                for(auto& item:leaves){auto& r=item.second;std::sort(r.angles.begin(),r.angles.end());r.angles.erase(std::unique(r.angles.begin(),r.angles.end()),r.angles.end());leaf(r);}}
        }else{
            for(const auto& request:requests){auto r=request;
                if(!r.pairedDepth&&e.angularCache){r.angles=missingFor(r,true);if(r.angles.empty())continue;}
                // Without the angular cache only exact native/prefix commands
                // are reusable; decline decisions rather than assume row reuse.
                if(!e.angularCache){result.profiled=false;}
                const int block=e.incremental?e.orientationBlock:r.count;
                for(int offset=0;offset<r.count;offset+=block)leaf(Engine::WorkRequest{std::min(block,r.count-offset),offset,r.depth,r.pairedDepth,r.angles});
            }
        }
        forecasts.emplace(requests,result);return result;
    }
};
std::vector<Engine::WorkRequest> BudgetRequests(const Engine& e,const std::vector<int>& depths,const Budgets& counts,const std::vector<size_t>& ids,bool uniform=false){
    std::vector<Engine::WorkRequest> requests;
    for(size_t l=0;l<depths.size();++l){std::map<int,std::vector<double>> groups;
        for(size_t id:ids)groups[counts[l][id]].push_back(e.theta[id]);
        for(auto& g:groups){auto& x=g.second;if(x.size()==1&&(uniform||e.thetaPhi.empty())){x.push_back(ids.size()>1?(e.theta[ids.front()]==x[0]?e.theta[ids.back()]:e.theta[ids.front()]):(e.theta.front()==x[0]?e.theta.back():e.theta.front()));std::sort(x.begin(),x.end());}
            requests.push_back(Engine::WorkRequest{g.first,0,depths[l],l&&e.depthSharing?depths[l-1]:0,x});
            if(l&&!e.depthSharing)requests.push_back(Engine::WorkRequest{g.first,0,depths[l-1],0,x});}}
    return requests;
}
