// Included inside Engine. Fold already verified angular rows straight into a
// prefix, avoiding an on-disk angular view for every block and every seed.
// Selection order and the sequence of floating point additions match Stage.
using PrefixSource=std::pair<const Job*,uint64_t>;
bool AccumulateCachedPrefixBlock(int count,int offset,int total,const unsigned* seedSet,size_t ns,int depth,
        const std::vector<double>& x,const std::string& weightHash,const std::vector<std::array<double,2>>& coefficients,
        std::vector<Job>& result,std::vector<std::ostringstream>& sources,std::set<std::string>& verified,
        std::vector<std::vector<PrefixSource>>* references=nullptr){
    if(!angularCache)return false;
    using SourceRow=std::pair<const Job*,size_t>;
    std::vector<std::map<double,SourceRow>> rows(ns);std::vector<const Job*> exact(ns,nullptr);
    bool allExact=true;
    for(size_t s=0;s<ns;++s){
        for(const auto& item:SamplingCandidates(count,seedSet[s],depth,offset)){
            const Job& job=*item.second;if(!SameSampling(job,count,seedSet[s],depth,offset))continue;
            // Stage retains the last exact match, but the first identity for
            // each overlapping row when the eight-seed stage is not exact.
            if(job.inputAngles==x&&job.weightHash==weightHash)exact[s]=&job;
            MatchingAngularRows(job.inputAngles,x,[&](size_t r){rows[s].insert({job.inputAngles[r],{&job,r}});});
        }
        if(rows[s].size()!=x.size())return false;
        allExact=allExact&&exact[s];
    }
    // A non-shared analytic mean requires a separate full-grid native pilot.
    // Keep that established path, rather than change the control estimator.
    if(!allExact&&controls&&!sharedMeans)return false;
    const double fraction=double(count)/total;
    for(size_t s=0;s<ns;++s){
        Job& destination=result[s];double represented=0,setup=0,refinement=0;
        std::map<std::string,const Job*> used;
        for(size_t r=0;r<x.size();++r){
            const SourceRow source=allExact?SourceRow{exact[s],r}:rows[s].at(x[r]);
            const Job& job=*source.first;const size_t index=source.second;const auto data=Payload(job);Mueller matrix=data->matrix.at(index);
            used[job.folder]=&job;
            if(controls){const auto& c=data->controls.at(index);const auto& b=coefficients[r];
                if(!allExact)matrix[0]=c[0]+b[0]*(c[1]-c[3])+b[1]*(c[2]-c[4]);
                for(size_t k=0;k<5;++k)destination.controls[r][k]+=fraction*c[k];
                refinement=std::max(refinement,data->controlRefinement);setup=data->controlSetup;
            }
            for(size_t k=0;k<16;++k)destination.matrix[r][k]+=fraction*matrix[k];
        }
        // Retain direct dependencies in the final durable prefix. Integrity
        // validation and integral traversal still visit every selected source.
        for(const auto& entry:used){ValidateStoredJob(entry.first,verified);
            const auto marker=FileHash(entry.first+"/complete.txt");
            if(references)(*references)[s].push_back({entry.second,std::stoull(marker,nullptr,16)});
            else sources[s]<<entry.first<<'\t'<<marker<<'\n';
            represented+=entry.second->seconds;}
        Require(represented>0,"missing cached prefix source cost");
        destination.seconds+=allExact?exact[s]->seconds:ForecastFullGridCost(count,seedSet[s],depth,x.size(),setup,represented);
        destination.controlRefinement=std::max(destination.controlRefinement,refinement);destination.controlSetup+=setup;
    }
    return true;
}
