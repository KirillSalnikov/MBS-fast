// Metadata-only reuse of complete production prefixes. Sampling choices are
// determined by available work, never by a favourable matrix or error estimate.
// Stage remains responsible for validating selected payloads and dependencies.
bool CompleteCachedPrefix(Engine& e,int count,int depth,const std::vector<double>& angles,const unsigned* seeds){
    if(!e.angularCache||(e.controls&&!e.sharedMeans)||angles.empty()||count<=0)return false;
    for(size_t s=0;s<8;++s){
        std::vector<bool> covered(angles.size(),false);size_t remaining=angles.size();
        for(const auto& entry:e.SamplingCandidates(count,seeds[s],depth,0)){
            const Job& job=*entry.second;if(!e.SameSampling(job,count,seeds[s],depth,0))continue;
            for(size_t r=0;r<angles.size();++r)if(!covered[r]
                    &&std::binary_search(job.inputAngles.begin(),job.inputAngles.end(),angles[r])){
                covered[r]=true;--remaining;
            }
            if(!remaining)break;
        }
        if(remaining)return false;
    }
    return true;
}

bool CachedPrefixPadding(Engine& e,int count,int depth,int lowerDepth,double angle,const unsigned* seeds,double& padding){
    if(!e.angularCache||(e.controls&&!e.sharedMeans)||count<=0)return false;
    // Keep only the smallest candidates if an exceptionally wide angular cache
    // exceeds this metadata budget. A returned row is still the smallest fully
    // covered candidate; exhausting the bounded set merely declines reuse.
    const size_t limit=65536;std::vector<double> candidates;
    for(const auto& entry:e.SamplingCandidates(count,seeds[0],depth,0)){
        const Job& job=*entry.second;if(!e.SameSampling(job,count,seeds[0],depth,0))continue;
        std::vector<double> merged;merged.reserve(std::min(limit,candidates.size()+job.inputAngles.size()));
        size_t a=0,b=0;
        while(merged.size()<limit&&(a<candidates.size()||b<job.inputAngles.size())){
            double value;
            if(b==job.inputAngles.size()||(a<candidates.size()&&candidates[a]<job.inputAngles[b]))value=candidates[a++];
            else if(a==candidates.size()||job.inputAngles[b]<candidates[a])value=job.inputAngles[b++];
            else {value=candidates[a++];++b;}
            if(value!=angle&&(merged.empty()||merged.back()!=value))merged.push_back(value);
        }
        candidates.swap(merged);
    }
    if(candidates.empty())return false;
    // Intersect complete row coverage across all eight seeds and both legs.
    // Sorted merge scans avoid asking CompleteCachedPrefix once per candidate,
    // which would repeatedly scan a large bucket for every possible angle.
    const std::vector<int> legs=lowerDepth?std::vector<int>{depth,lowerDepth}:std::vector<int>{depth};
    for(int leg:legs)for(size_t s=0;s<8;++s){
        std::vector<unsigned char> covered(candidates.size(),0);bool target=false;size_t remaining=candidates.size();
        for(const auto& entry:e.SamplingCandidates(count,seeds[s],leg,0)){
            const Job& job=*entry.second;if(!e.SameSampling(job,count,seeds[s],leg,0))continue;
            target=target||std::binary_search(job.inputAngles.begin(),job.inputAngles.end(),angle);
            size_t a=0,b=std::lower_bound(job.inputAngles.begin(),job.inputAngles.end(),candidates.front())-job.inputAngles.begin();
            while(a<candidates.size()&&b<job.inputAngles.size()){
                if(candidates[a]<job.inputAngles[b])++a;
                else if(job.inputAngles[b]<candidates[a])++b;
                else {if(!covered[a]){covered[a]=1;--remaining;}++a;++b;}
            }
            if(target&&!remaining)break;
        }
        if(!target)return false;
        size_t keep=0;for(size_t r=0;r<candidates.size();++r)if(covered[r])candidates[keep++]=candidates[r];
        candidates.resize(keep);if(candidates.empty())return false;
    }
    padding=candidates.front();return true;
}

Budgets CachedCompletedBudgets(Engine& e,const std::vector<int>& depths,const Budgets& current,const std::vector<size_t>& ids){
    if(!e.angularCache||(e.controls&&!e.sharedMeans)||ids.empty())return current;
    Require(depths.size()==current.size()&&e.seeds.size()==8,"invalid cached production budget dimensions");
    // A singleton's auxiliary row is checked separately, so the candidate
    // index needs only the actual production rows.
    std::set<size_t> wanted(ids.begin(),ids.end());
    const std::vector<size_t> indexed(wanted.begin(),wanted.end());
    std::map<double,size_t> positions;
    for(size_t r=0;r<indexed.size();++r)positions.emplace(e.theta[indexed[r]],r);
    using LevelPhi=std::pair<size_t,int>;
    std::map<LevelPhi,int> minimum;
    std::map<int,std::vector<size_t>> phiRows;
    for(size_t id:ids){const int phi=e.PhiForTheta(id);phiRows[phi].push_back(id);
        for(size_t l=0;l<current.size();++l){const LevelPhi key{l,phi};const auto at=minimum.find(key);
            if(at==minimum.end())minimum.emplace(key,current[l][id]);else at->second=std::min(at->second,current[l][id]);}}
    // There are only power-of-two prefix counts, but a heterogeneous phi map
    // can still have many combinations. Bound masks independently of payloads.
    using Key=std::tuple<size_t,int,int>; // level, phi, count
    std::map<Key,std::vector<uint16_t>> coverage;
    const size_t limit=32*1024*1024,rowBytes=indexed.size()*sizeof(uint16_t);
    size_t bytes=0;bool bounded=false;
    for(const auto& entry:e.cache){const Job& job=entry.second;
        if(job.offset||job.count<=0||job.count>e.maxCount||(job.count&(job.count-1))||!phiRows.count(job.phi))continue;
        const auto seed=std::find(e.seeds.begin(),e.seeds.end(),job.seed);if(seed==e.seeds.end())continue;
        PhiScope scope(e,job.phi);
        if(!e.SameSampling(job,job.count,job.seed,job.depth,0))continue;
        const size_t seedIndex=seed-e.seeds.begin();
        for(size_t l=0;l<depths.size();++l){
            const bool upper=job.depth==depths[l],lower=l&&job.depth==depths[l-1];
            if(!upper&&!lower)continue;
            const auto threshold=minimum.find({l,job.phi});
            if(threshold==minimum.end()||job.count<=threshold->second)continue;
            const Key key{l,job.phi,job.count};auto found=coverage.find(key);
            if(found==coverage.end()){
                if(rowBytes>limit-bytes){bounded=true;continue;}
                found=coverage.emplace(key,std::vector<uint16_t>(indexed.size(),0)).first;bytes+=rowBytes;
            }
            const uint16_t bit=uint16_t(1u<<(seedIndex+(lower?8:0)));
            for(double angle:job.inputAngles){const auto row=positions.find(angle);
                if(row!=positions.end())found->second[row->second]|=bit;}
        }
    }
    auto complete=[&](size_t level,int phi,int count,size_t id){
        const auto key=coverage.find(Key{level,phi,count});if(key==coverage.end())return false;
        const auto row=positions.find(e.theta[id]);if(row==positions.end())return false;
        const uint16_t expected=level?uint16_t(0xffff):uint16_t(0xff);
        return key->second[row->second]==expected;
    };
    Budgets result=current;
    for(const auto& entry:coverage){const size_t l=std::get<0>(entry.first);const int phi=std::get<1>(entry.first),count=std::get<2>(entry.first);
        for(size_t id:phiRows.at(phi))if(count>result[l][id]&&complete(l,phi,count,id))result[l][id]=count;
    }
    // Match EvaluateUniformBudget's cached singleton padding. A promotion must
    // not turn a free cache read into an extra optical calculation. Reverting
    // one group can change other group sizes, so repeat until no row changes.
    bool changed=true;
    while(changed){changed=false;
        for(size_t l=0;l<result.size();++l)for(const auto& phiGroup:phiRows){
            const int phi=phiGroup.first;const auto& rows=phiGroup.second;
            std::map<int,std::vector<size_t>> groups;
            for(size_t id:rows)groups[result[l][id]].push_back(id);
            for(const auto& group:groups){if(group.second.size()!=1)continue;
                const size_t id=group.second.front();if(result[l][id]==current[l][id])continue;
                PhiScope scope(e,phi);double padding=0;
                if(!CachedPrefixPadding(e,result[l][id],depths[l],l?depths[l-1]:0,e.theta[id],e.seeds.data(),padding)){
                    result[l][id]=current[l][id];changed=true;
                }
            }
        }
    }
    if(bounded)e.Log("cached_prefix_promotion_index_limit","{\"limit_bytes\":"+Number(limit)+",\"indexed_bytes\":"+Number(bytes)+",\"unindexed_candidates_retained_at_current_budget\":true}");
    return result;
}
