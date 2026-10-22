// Included inside Engine. Stream a saved orientation block through several
// pending prefixes while its wide matrices are resident. Each prefix keeps
// its own source selection, weights and increasing-offset summation order.
struct PreparedPrefix {
    int count=0,depth=0;std::vector<double> angles;std::string weightHash;
    std::vector<Job> values;std::vector<std::ostringstream> sources;
    std::vector<std::array<double,2>> coefficients;size_t blocks=0,bytes=0;
};
std::vector<WorkRequest> prefixAssemblyPlans;
std::map<std::string,std::unique_ptr<PreparedPrefix>> preparedPrefixes;
std::set<std::string> attemptedPrefixBatches;
size_t preparedPrefixBytes=0;
size_t prefixLeafGeneration=0,preparedPrefixGeneration=0;
std::array<unsigned,8> prefixPlanSeeds{{0,0,0,0,0,0,0,0}};
std::string PrefixBatchKey(int n,int depth,const std::vector<double>& x)const{
    std::ostringstream key;key<<n<<' '<<depth<<' '<<phi;for(double a:x)key<<' '<<Number(a);return key.str();
}
void ClearPrefixAssemblyPlan(){prefixAssemblyPlans.clear();preparedPrefixes.clear();attemptedPrefixBatches.clear();preparedPrefixBytes=0;}
void SetPrefixAssemblyPlan(const std::vector<WorkRequest>& requests,const unsigned* seeds){
    ClearPrefixAssemblyPlan();if(!incremental||!angularCache||(controls&&!sharedMeans)||!thetaPhi.empty())return;
    std::copy(seeds,seeds+8,prefixPlanSeeds.begin());
    std::set<std::string> unique;
    for(const auto& request:requests)for(int depth:{request.depth,request.pairedDepth}){
        if(!depth||request.offset||request.count<=orientationBlock)continue;
        // Reproduce Stage's missing-row/padding decision from metadata only.
        bool eachSeedHasRows=true;std::vector<std::set<double>> rows(8);
        for(size_t s=0;s<8;++s){for(const auto& item:SamplingCandidates(request.count,seeds[s],depth,0)){
            const Job& j=*item.second;if(!SameSampling(j,request.count,seeds[s],depth,0))continue;
            for(double a:request.angles)if(std::binary_search(j.inputAngles.begin(),j.inputAngles.end(),a))rows[s].insert(a);
        }eachSeedHasRows=eachSeedHasRows&&!rows[s].empty();}
        std::vector<double> missing;
        for(double a:request.angles){bool present=eachSeedHasRows;for(const auto& r:rows)present=present&&r.count(a);if(!present)missing.push_back(a);}
        if(missing.empty())continue;
        if(missing.size()==1){for(double a:request.angles)if(a!=missing[0]){missing.push_back(a);break;}std::sort(missing.begin(),missing.end());}
        if(missing.size()<2)continue;
        const auto key=PrefixBatchKey(request.count,depth,missing);
        if(unique.insert(key).second)prefixAssemblyPlans.push_back({request.count,0,depth,0,missing});
    }
}
std::unique_ptr<PreparedPrefix> TakePreparedPrefix(int n,const unsigned* seeds,size_t ns,int depth,const std::vector<double>& x){
    if(ns!=8||!std::equal(seeds,seeds+8,prefixPlanSeeds.begin()))return {};
    if(preparedPrefixGeneration!=prefixLeafGeneration){preparedPrefixes.clear();preparedPrefixBytes=0;preparedPrefixGeneration=prefixLeafGeneration;}
    const auto wanted=PrefixBatchKey(n,depth,x);auto saved=preparedPrefixes.find(wanted);
    if(saved!=preparedPrefixes.end()){auto result=std::move(saved->second);preparedPrefixBytes-=result->bytes;preparedPrefixes.erase(saved);return result;}
    if(prefixAssemblyPlans.empty()||!attemptedPrefixBatches.insert(wanted).second)return {};
    const size_t limit=128*1024*1024;const double began=Seconds();
    std::vector<std::unique_ptr<PreparedPrefix>> plans;std::set<std::string> selected;
    size_t estimate=preparedPrefixBytes;
    auto add=[&](const std::vector<double>& angles){
        const auto key=PrefixBatchKey(n,depth,angles);
        if(selected.count(key)||preparedPrefixes.count(key))return;
        const size_t cost=size_t((n-1)/orientationBlock+1)*ns*std::min(size_t(4),angles.size())*640;
        if(cost>limit-estimate||plans.size()>=8)return;
        estimate+=cost;selected.insert(key);std::unique_ptr<PreparedPrefix> p(new PreparedPrefix);
        p->count=n;p->depth=depth;p->angles=angles;p->values.resize(ns);p->sources.resize(ns);
        p->coefficients.assign(angles.size(),{{1.,1.}});const auto wf=WriteWeights(angles);p->weightHash=Hash(wf.empty()?"":Read(wf));
        if(!wf.empty()){std::istringstream input(Read(wf));std::string header;std::getline(input,header);
            for(auto& b:p->coefficients){double angle;Require(bool(input>>angle>>b[0]>>b[1]),"invalid batched prefix weights");}}
        for(auto& j:p->values){j.matrix.assign(angles.size(),Mueller{});if(controls)j.controls.assign(angles.size(),std::array<double,5>{{0,0,0,0,0}});}
        plans.push_back(std::move(p));
    };
    add(x);for(const auto& r:prefixAssemblyPlans)if(r.count==n&&r.depth==depth)add(r.angles);
    if(plans.size()<2)return {};
    std::set<std::string> verified;bool bounded=true;
    for(int offset=0;offset<n;offset+=orientationBlock){const int count=std::min(orientationBlock,n-offset);
        CheckHostMemory();
        for(auto& p:plans)if(p){
            if(!AccumulateCachedPrefixBlock(count,offset,n,seeds,ns,depth,p->angles,p->weightHash,p->coefficients,p->values,p->sources,verified)){p.reset();continue;}
            ++p->blocks;p->bytes=sizeof(PreparedPrefix)+ns*p->angles.size()*(sizeof(Mueller)+sizeof(std::array<double,5>));
            for(auto& source:p->sources)p->bytes+=size_t(source.tellp())*2; // include stream allocation slack
        }
        size_t used=preparedPrefixBytes;for(const auto& p:plans)if(p)used+=p->bytes;
        // The source set is shared across all plans and bounded as well.
        used+=verified.size()*(root.size()+512);
        if(used>limit){bounded=false;break;}
        if(std::none_of(plans.begin(),plans.end(),[](const std::unique_ptr<PreparedPrefix>& p){return bool(p);}))break;
    }
    size_t completed=0;
    if(bounded)for(auto& p:plans)if(p&&p->blocks==size_t((n-1)/orientationBlock+1)){
        const auto key=PrefixBatchKey(n,depth,p->angles);preparedPrefixBytes+=p->bytes;preparedPrefixes.emplace(key,std::move(p));++completed;
    }
    Log("prefix_assembly_batch","{\"count_per_seed\":"+Number(n)+",\"depth\":"+Number(depth)+",\"completed_prefixes\":"+Number(completed)
        +",\"retained_bytes\":"+Number(preparedPrefixBytes)+",\"memory_limit_bytes\":"+Number(limit)+",\"bounded\":"+(bounded?"true":"false")
        +",\"seconds\":"+Number(Seconds()-began)+",\"native_work\":0,\"source_selection_and_sum_order_preserved\":true}");
    saved=preparedPrefixes.find(wanted);if(saved==preparedPrefixes.end())return {};
    auto result=std::move(saved->second);preparedPrefixBytes-=result->bytes;preparedPrefixes.erase(saved);return result;
}
