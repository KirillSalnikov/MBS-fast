// Included inside Engine. Group nearby dyadic azimuth rules while retaining
// independent wall timings on a separate validation group at every candidate.
bool phiTimingWave=false;
size_t PhiTimingSeeds()const{
    // A partial wave on one GPU can underestimate the cost of wider rules.
    return phiTimingWave?std::min<size_t>(8,std::max<size_t>(2,WorkerGpus().size())):2;
}
struct NestedPhiSweep {
    std::map<int,std::pair<int,int>> groups;
    std::map<std::tuple<int,int,bool,size_t>,std::vector<Job>> samples;
    explicit NestedPhiSweep(const std::vector<int>& rules,bool completeSweep=false) {
        for(size_t first=0;first<rules.size();){
            size_t last=first;
            while(last+1<rules.size()){
                const int next=rules[last+1],lower=rules[first];
                // A full pass already computes every azimuth in a dyadic
                // subset, including a subset separated by an unrequested rule.
                if(completeSweep){
                    if(lower<=0||next<=rules[last]||next>64||next%lower)break;
                    const int ratio=next/lower;
                    if(ratio&(ratio-1))break;
                }else if(next!=2*rules[last]||next>8*lower)break;
                ++last;
            }
            for(size_t i=first;i<=last;++i)groups[rules[i]]={rules[first],rules[last]};
            first=last+1;
        }
    }
    bool Grouped(int rule)const{return groups.at(rule).first!=groups.at(rule).second;}
};
void AppendNestedPhiIntegrity(const Job& job,std::ostringstream& manifest)const {
    if(!job.nestedPhiMin)return;
    for(int count=job.phi;count>job.nestedPhiMin;){
        Require(count%2==0,"invalid nested phi division");
        count/=2;
        const std::string stem=job.folder+"/result/result_phi"+Number(count);
        for(const char* suffix:{".dat","_integrals.tsv","_analytic_facets.tsv"}){
            if(!job.hasControls&&std::string(suffix)=="_analytic_facets.tsv")continue;
            const std::string file=stem+suffix;manifest<<file<<' '<<FileHash(file)<<'\n';
        }
    }
}
std::vector<Job> NestedPhiJobs(NestedPhiSweep& sweep,int candidate,int n,int depth,const std::vector<double>& x,bool train,bool deferStandalone=false) {
    const auto limits=sweep.groups.at(candidate);
    const int upper=limits.second;
    const size_t timed=train||deferStandalone?0:PhiTimingSeeds();
    const unsigned* scrambles=train?training:validation;
    std::vector<Job> rows;
    if(timed)rows=DirectStage(n,scrambles,timed,depth,x,false);
    if(timed==8)return rows;
    auto& sources=sweep.samples[std::make_tuple(upper,n,train,timed)];
    if(sources.empty()){
        phi=upper;
        try {sources=DirectStage(n,scrambles+timed,8-timed,depth,x,false,false,0,0,limits.first);}
        catch(...){phi=candidate;throw;}
        phi=candidate;
        std::set<std::string> verified;
        for(const auto& s:sources)ValidateStoredJob(s.folder,verified);
        Log("nested_phi_source","{\"upper\":"+Number(upper)+",\"lower\":"+Number(limits.first)
            +",\"orientations_per_seed\":"+Number(n)+",\"depth\":"+Number(depth)+",\"theta_rows\":"+Number(x.size())
            +",\"train\":"+(train?"true":"false")+",\"scrambles\":"+Number(sources.size())+"}");
    }
    for(const auto& source:sources){
        if(candidate==upper){rows.push_back(source);continue;}
        Job view=source;
        const std::string stem=source.folder+"/result/result_phi"+Number(candidate);
        view.phi=candidate;view.data=stem+".dat";view.summary=stem+"_analytic_facets.tsv";
        view.matrix.clear();view.controls.clear();view.angles.clear();
        // Only Calibrate consumes this view. The parent's measured wall is
        // not a timing for this azimuth rule; never index it for production.
        view=Load(view,false);
        Require(view.angles==source.angles,"nested phi changed theta coordinates");
        rows.push_back(std::move(view));
    }
    return rows;
}
