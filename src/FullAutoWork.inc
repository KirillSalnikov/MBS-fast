// Included inside FullAutoAdaptive.inc. Scheduling forecasts only: never an
// accuracy criterion. Times are summed native worker seconds, not elapsed time.
struct WorkForecast {double seconds=0;size_t jobs=0;bool profiled=true;};
class AdditionalWork {
    const Engine& e;std::array<unsigned,8> seeds;
    struct Profile {double trace=0,diffraction=0,setup=0,overlap=0;size_t observations=0;};
    using Depth=std::pair<int,int>;
    using Slot=std::array<uint64_t,5>; // N, offset, depth, paired depth, seed
    std::map<Depth,Profile> profiles;
    std::map<Slot,std::vector<const Job*>> coverage;
    std::map<std::string,WorkForecast> forecasts;
    Slot Key(int n,int offset,int depth,int paired,unsigned seed)const{return {{uint64_t(n),uint64_t(offset),uint64_t(depth),uint64_t(paired),seed}};}
    bool Covered(int n,int offset,int depth,int paired,unsigned seed,double angle){
        const auto key=Key(n,offset,depth,paired,seed);auto found=coverage.find(key);
        if(found==coverage.end()){std::vector<const Job*> candidates;for(const auto& item:e.SamplingCandidates(n,seed,depth,offset)){
            const Job& job=*item.second;if(e.SameSampling(job,n,seed,depth,offset)&&(!paired||job.pairedDepth==paired))candidates.push_back(&job);}
            found=coverage.emplace(key,candidates).first;}
        for(const Job* job:found->second)if(std::binary_search(job->inputAngles.begin(),job->inputAngles.end(),angle))return true;return false;
    }

public:
    AdditionalWork(const Engine& engine,const unsigned* seedSet):e(engine){
        std::copy(seedSet,seedSet+8,seeds.begin());const auto recipe=e.ExecutionHash();
        for(const auto& item:e.cache){const Job& j=item.second;
            if(j.assembled||j.paired||j.count<=0||j.inputAngles.empty()||j.executionHash!=recipe||j.phi!=e.phi
                ||j.block!=e.kernel.block||j.pipeline!=e.kernel.pipeline||j.warp!=e.kernel.warp||j.mirrored!=e.mirror||j.hasControls!=e.controls)continue;
            if(j.seconds<=0||!std::isfinite(j.seconds))continue;
            auto& p=profiles[{j.depth,j.pairedDepth}];
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
    WorkForecast Forecast(const std::vector<Engine::WorkRequest>& requests){
        std::ostringstream name;for(const auto& r:requests){name<<r.count<<':'<<r.depth<<':'<<r.pairedDepth<<':';for(double q:r.angles)name<<Number(q)<<',';name<<';';}
        auto saved=forecasts.find(name.str());if(saved!=forecasts.end())return saved->second;
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
            auto profile=profiles.find({r.depth,r.pairedDepth});if(profile==profiles.end())result.profiled=false;
            for(unsigned seed:seeds){
                // DirectStage can reuse a completed native command for a seed
                // even when another seed of the batch was interrupted.
                bool exact=false;if(!e.controls)for(const auto& item:e.SamplingCandidates(r.count,seed,r.depth,r.offset)){
                    const Job& j=*item.second;exact=exact||(!j.assembled&&j.pairedDepth==r.pairedDepth&&j.inputAngles==missing&&e.SameSampling(j,r.count,seed,r.depth,r.offset));}
                if(!exact){++result.jobs;if(profile!=profiles.end()){const auto& p=profile->second;const double trace=r.count*p.trace,diffraction=r.count*p.diffraction*missing.size();
                    result.seconds+=p.setup+trace+diffraction-p.overlap*std::min(trace,diffraction);}}
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
            for(int offset=0;offset<maximum;offset+=block){std::map<std::string,Engine::WorkRequest> leaves;
                for(const auto& r:requests)if(offset<r.count){const int n=std::min(block,r.count-offset);const auto key=Number(n)+" "+Number(r.depth)+" "+Number(r.pairedDepth);
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
        forecasts[name.str()]=result;return result;
    }
};
std::vector<Engine::WorkRequest> BudgetRequests(const Engine& e,const std::vector<int>& depths,const Budgets& counts,const std::vector<size_t>& ids){
    std::vector<Engine::WorkRequest> requests;
    for(size_t l=0;l<depths.size();++l){std::map<int,std::vector<double>> groups;
        for(size_t id:ids)groups[counts[l][id]].push_back(e.theta[id]);
        for(auto& g:groups){auto& x=g.second;if(x.size()==1){x.push_back(ids.size()>1?(e.theta[ids.front()]==x[0]?e.theta[ids.back()]:e.theta[ids.front()]):(e.theta.front()==x[0]?e.theta.back():e.theta.front()));std::sort(x.begin(),x.end());}
            requests.push_back(Engine::WorkRequest{g.first,0,depths[l],l&&e.depthSharing?depths[l-1]:0,x});
            if(l&&!e.depthSharing)requests.push_back(Engine::WorkRequest{g.first,0,depths[l-1],0,x});}}
    return requests;
}
