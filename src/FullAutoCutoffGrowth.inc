// Kept separate from production sample identities. A parent with a larger N
// can reuse an earlier fast, balanced, or off child only for the SAME profile.
struct CutoffGrowth {
    Job target, prefix;
    std::vector<std::string> command;
};
std::vector<std::string> CutoffWindowPhysics(const std::vector<std::string>& args)const{
    std::vector<std::string> physics;
    for(size_t i=0;i<args.size();++i){
        // An embedded prefix is a side output. Only full native matrices are
        // candidates here; targets requesting that side output are excluded below.
        if(args[i]=="--orientation-prefix-output"||args[i]=="--output"||args[i]=="--theta-grid-file"||args[i]=="--sobol-offset"||args[i]=="--analytic-control-weights"){++i;continue;}
        if(args[i]=="--no-text-reports")continue;
        if(args[i]=="--sobol-seed"){i+=2;continue;}
        physics.push_back(args[i]);
    }
    return physics;
}
using CutoffWeights=std::vector<std::array<double,2>>;
CutoffWeights CutoffWindowWeights(const std::vector<std::string>& cmd,const std::vector<double>& angles,const std::string& expected)const{
    CutoffWeights out(angles.size(),{{1.,1.}});
    const auto option=std::find(cmd.begin(),cmd.end(),"--analytic-control-weights");
    if(option==cmd.end()){Require(expected==Hash(""),"missing signed control recipe");return out;}
    Require(std::next(option)!=cmd.end(),"missing cutoff weights path");
    const auto text=Read(*std::next(option));Require(Hash(text)==expected,"signed control weights changed");
    std::istringstream file(text);std::string header;std::getline(file,header);
    double x,a,b;size_t row=0;
    while(row<angles.size()&&file>>x>>a>>b){
        if(x<angles[row])continue;
        Require(x==angles[row]&&std::isfinite(a)&&std::isfinite(b),"inconsistent signed control weights");
        out[row++]={{a,b}};
    }
    Require(row==angles.size(),"missing signed control weights");return out;
}
const Job* CutoffWindowSource(const Job& target,int count,const std::vector<std::string>& cmd){
    const auto original=CutoffWindowPhysics(cmd);
    auto actual=cmd;Job dispatch=target;
    // The worker records its final thread and pipeline arguments after its
    // optical request is assembled. Reproduce those choices without launching
    // a timing pilot or weakening comparison of the physical CLI parameters.
    ApplyDispatchThreads(dispatch,actual);
    ApplyAutomaticPreparation(dispatch,actual);
    const auto physics=CutoffWindowPhysics(actual);
    auto match=[&](const Job& s,bool project)->const Job*{
        // LowerView registers its derived result in the same-depth index. Always
        // enter through the native paired parent, including after restoration,
        // so an old view cannot bypass verification or the role of its source.
        if(s.folder==root+"/jobs/d"+s.identity)return nullptr;
        if(s.executionHash!=target.executionHash||s.pairedDepth!=(project?target.depth:target.pairedDepth)||s.hasControls!=target.hasControls||s.paired||
           s.block!=target.block||s.phi!=target.phi||s.warp!=target.warp||s.pipeline!=target.pipeline||s.mirrored!=target.mirrored||
           (!target.hasControls&&s.weightHash!=target.weightHash)||!std::includes(s.inputAngles.begin(),s.inputAngles.end(),target.inputAngles.begin(),target.inputAngles.end()))return nullptr;
        // An embedded prefix view's replay command is reconstructed for standalone
        // execution. Compare the original native invocation which wrote the view.
        // Confirm its unique parent and the exact count encoded in the view ID.
        std::string commandFolder=s.folder;
        if(s.assembled&&s.folder==root+"/jobs/v"+s.identity){
            std::istringstream parents(Read(s.folder+"/angular_sources.tsv"));
            std::string parent,extra;
            if(!std::getline(parents,parent)||std::getline(parents,extra))return nullptr;
            const auto tab=parent.find('\t');
            if(tab==std::string::npos)return nullptr;
            const auto path=parent.substr(0,tab);
            if(path.compare(0,root.size()+7,root+"/jobs/c")!=0)return nullptr;
            const auto id=path.substr(root.size()+7);
            if(Hash(id+"\nembedded cutoff prefix v1 "+Number(count))!=s.identity)return nullptr;
            commandFolder=path;
        }
        std::istringstream input(Read(commandFolder+"/command.args"));std::string name;
        if(!std::getline(input,name)||name!=exe)return nullptr;
        std::vector<std::string> saved;while(std::getline(input,name))saved.push_back(name);
        if(project){
            // Only project a native paired child or its authenticated native prefix.
            // Other assembled jobs can carry replay commands from production.
            if(s.folder!=root+"/jobs/c"+s.identity&&s.folder!=root+"/jobs/v"+s.identity)return nullptr;
            if(s.assembled&&s.folder!=root+"/jobs/v"+s.identity)return nullptr;
            const auto stem=s.folder+"/result/result_depth"+Number(target.depth);
            if(!Exists(stem+"_integrals_revision.txt")||Read(stem+"_integrals_revision.txt")!="paired_extinction_v1\n")return nullptr;
            auto pair=std::find(saved.begin(),saved.end(),"--paired-depth");
            if(pair==saved.end()||std::next(pair)==saved.end()||*std::next(pair)!=Number(target.depth))return nullptr;
            saved.erase(pair,pair+2);
            auto maxDepth=std::find(saved.begin(),saved.end(),"--max-reflections");
            auto sobol=std::find(saved.begin(),saved.end(),"--sobol-seed");
            if(maxDepth==saved.end()||std::next(maxDepth)==saved.end()||*std::next(maxDepth)!=Number(s.depth)
                ||sobol==saved.end()||saved.end()-sobol<3||sobol[2]!=Number(s.seed))return nullptr;
            maxDepth[1]=Number(target.depth);sobol[1]=Number(count);
            auto prefix=std::find(saved.begin(),saved.end(),"--orientation-prefix-output");
            if(prefix!=saved.end()){
                if(std::next(prefix)==saved.end())return nullptr;
                saved.erase(prefix,prefix+2);
            }
        }
        const auto previous=CutoffWindowPhysics(saved);
        if(previous!=physics&&!(s.assembled&&previous==original))return nullptr;
        if(target.hasControls&&CutoffWindowWeights(saved,target.inputAngles,s.weightHash)!=CutoffWindowWeights(cmd,target.inputAngles,target.weightHash))return nullptr;
        if(!project)return &s;
        std::set<std::string> verified;ValidateStoredJob(s.folder,verified);
        const auto view=LowerView(s,target.depth,false,&saved);
        return &cache.at(view.identity);
    };
    for(const auto& item:SamplingCandidates(count,target.seed,target.depth,target.offset))
        if(const auto* found=match(*item.second,false))return found;
    if(!reuseCutoffDepth||!depthSharing||!depthVerified||target.pairedDepth||target.paired)return nullptr;
    // Skip whole depth buckets, regardless of the number of saved offsets.
    for(auto it=samplingIndex.upper_bound(SamplingSlot(count,target.seed,target.depth,INT_MAX));
        it!=samplingIndex.end()&&it->first[0]==uint64_t(count)&&it->first[1]==target.seed;){
        const int depth=int(it->first[2]);
        for(const auto& item:SamplingCandidates(count,target.seed,depth,target.offset))
            if(const auto* found=match(*item.second,true))return found;
        it=samplingIndex.upper_bound(SamplingSlot(count,target.seed,depth,INT_MAX));
    }
    return nullptr;
}
Job CutoffWindowRows(const Job& source,const std::vector<double>& x){
    std::set<std::string> visited;ValidateStoredJob(source.folder,visited);
    auto s=Materialize(source);if(s.angles==x)return s;
    Matrix selected;selected.reserve(x.size());std::vector<std::array<double,5>> cv;
    for(double t:x){auto at=std::lower_bound(s.inputAngles.begin(),s.inputAngles.end(),t);
        Require(at!=s.inputAngles.end()&&*at==t,"cutoff child missing requested angle");const size_t index=at-s.inputAngles.begin();selected.push_back(s.matrix[index]);if(s.hasControls)cv.push_back(s.controls.at(index));}
    s.matrix.swap(selected);s.controls.swap(cv);s.angles=x;return s;
}
void PlanCutoffGrowth(std::vector<Job>& jobs,std::vector<std::vector<std::string>>& commands,
                      std::vector<size_t>& launch,std::map<size_t,CutoffGrowth>& growth){
    if(!extendCutoffChild||!coherent||(controls&&!sharedMeans)||go)return;
    std::vector<size_t> pending;size_t aliases=0,savedSamples=0;
    for(size_t i:launch){auto& j=jobs[i];auto& cmd=commands[i];
        // Only raw correction legs. A standalone direct, paired output, or
        // embedded prefix has a different output contract.
        if(j.folder.compare(0,root.size()+7,root+"/jobs/c")!=0||j.paired||
            std::find(cmd.begin(),cmd.end(),"--orientation-prefix-output")!=cmd.end()){
            pending.push_back(i);continue;
        }
        const Job* equal=CutoffWindowSource(j,j.count,cmd);
        if(equal&&(!j.hasControls||equal->inputAngles==j.inputAngles)){
            if(equal->inputAngles==j.inputAngles)j=CutoffWindowRows(*equal,j.inputAngles);
            else j=CommitCutoffGrowth(CutoffGrowth{j,*equal,cmd},nullptr);
            ++aliases;savedSamples+=j.count;continue;
        }
        const int half=j.count/2;
        const Job* source=half>=64&&half*2==j.count?CutoffWindowSource(j,half,cmd):nullptr;
        if(!source){pending.push_back(i);continue;}
        growth.emplace(i,CutoffGrowth{j,*source,cmd});
        j.count=half;j.offset+=half;j.identity=Hash(growth.at(i).target.identity+" exact native second half v1");
        j.folder=root+"/jobs/h"+j.identity;j.data=j.folder+"/result/result.dat";j.summary=j.folder+"/result/result_analytic_facets.tsv";
        for(size_t k=0;k<cmd.size();++k){
            if(cmd[k]=="--sobol-seed")cmd[k+1]=Number(half);
            if(cmd[k]=="--output")cmd[k+1]=j.folder+"/result";
            if(cmd[k]=="--sobol-offset"){cmd.erase(cmd.begin()+k,cmd.begin()+k+2);--k;}
        }
        cmd.push_back("--sobol-offset");cmd.push_back(Number(j.offset));savedSamples+=half;
        if(cache.count(j.identity)){j=Materialize(cache.at(j.identity));continue;}
        if(Exists(j.folder+"/complete.txt")){std::set<std::string> seen;ValidateStoredJob(j.folder,seen);
            std::istringstream marker(Read(j.folder+"/complete.txt"));std::string key;marker>>key>>j.seconds;
            Require(key==j.identity,"cutoff suffix identity changed");j=Load(j);CacheJob(j);continue;}
        Directory(j.folder);
        if(Exists(j.folder+"/worker.pid"))Require(!WorkerAlive(j.folder+"/worker.pid"),"cutoff suffix still running");
        if(Exists(j.folder+"/result"))Require(RenamePath(j.folder+"/result",j.folder+"/result_incomplete_"+Number(getpid())+"_"+Number(i))==0,"cannot preserve incomplete cutoff suffix");
        pending.push_back(i);
    }
    launch.swap(pending);
    if(savedSamples)Log("cutoff_child_growth","{\"exact_windows_reused\":"+Number(aliases)+",\"half_windows_reused\":"+Number(growth.size())+",\"saved_native_orientations\":"+Number(savedSamples)+"}");
}
Job CommitCutoffGrowth(const CutoffGrowth& g,const Job* suffix){
    Job j=g.target;
    // A null suffix is a projection of a whole source onto fewer directions.
    // Both depths and the integral provenance then refer to that same source.
    const bool extend=suffix!=nullptr;
    const auto& prefix=g.prefix;
    Require(prefix.hasControls==j.hasControls&&prefix.pairedDepth==j.pairedDepth&&prefix.offset==j.offset&&
            prefix.seed==j.seed&&prefix.executionHash==j.executionHash&&
            prefix.count==(extend?j.count/2:j.count),"incompatible native cutoff prefix");
    if(extend)Require(prefix.count*2==j.count&&suffix->count==prefix.count&&
            suffix->pairedDepth==j.pairedDepth&&suffix->offset==j.offset+prefix.count&&
            suffix->seed==j.seed&&suffix->executionHash==j.executionHash&&
            suffix->inputAngles==j.inputAngles,"incompatible native cutoff suffix");
    auto meanAt=[&](int depth){
        const auto view=[&](const Job& src){return CutoffWindowRows(depth==j.depth?src:LowerView(src,depth),j.inputAngles);};
        Job left=view(prefix);
        if(!extend)return left;
        const auto right=view(*suffix);
        for(size_t r=0;r<left.matrix.size();++r){
            for(int k=0;k<16;++k)left.matrix[r][k]=.5*left.matrix[r][k]+.5*right.matrix[r][k];
            if(j.hasControls)for(int k=0;k<5;++k)left.controls[r][k]=.5*left.controls[r][k]+.5*right.controls[r][k];
        }
        // The native suffix evaluated exactly the target theta rows. The shared
        // analytic reference defines the bound as their maximum, independent
        // of the orientation window. A larger source grid can have a larger maximum.
        left.controlRefinement=right.controlRefinement;
        left.controlSetup+=right.controlSetup;
        return left;
    };
    const auto upper=meanAt(j.depth);Job lower;if(j.pairedDepth)lower=meanAt(j.pairedDepth);
    j.angles=j.inputAngles;j.assembled=true;j.seconds=prefix.seconds+(extend?suffix->seconds:0.);
    const auto weights=CutoffWindowWeights(g.command,j.inputAngles,j.weightHash);
    auto publish=[&](const Job& view,const std::string& stem){
        WriteMatrix(stem+".dat",j.angles,view.matrix);
        if(!j.hasControls)return;
        std::ostringstream os;os<<std::setprecision(17)<<"theta_deg\traw_M11\treflection_mean\tshadow_mean\tsampled_reflection\tsampled_shadow\tresidual_mean\thybrid_M11\trefinement\tsetup_seconds\treflection_weight\tshadow_weight\n";
        for(size_t r=0;r<j.angles.size();++r){
            const auto& c=view.controls[r];const auto& w=weights[r];os<<j.angles[r];
            for(double v:c)os<<'\t'<<v;
            os<<'\t'<<c[0]-w[0]*c[3]-w[1]*c[4]<<'\t'<<view.matrix[r][0]<<'\t'<<view.controlRefinement<<'\t'<<view.controlSetup<<'\t'<<w[0]<<'\t'<<w[1]<<'\n';
        }
        Atomic(stem+"_analytic_facets.tsv",os.str());
    };
    Directory(j.folder+"/result");
    {DerivedFileBatch files(j.folder,PrepareDerivedPublication(),false);
        publish(upper,j.folder+"/result/result");
        std::ostringstream command;command<<exe<<'\n';for(const auto& arg:g.command)command<<arg<<'\n';Atomic(j.folder+"/command.args",command.str());
        Atomic(j.folder+"/stdout.log",std::string(extend?"Exact mean of two disjoint native cutoff windows.":"Exact direction subset of a native cutoff window.")+"\nHard tree-limit hits: 0\nOrientations still incomplete: 0\n");
        const auto origin=[](const Job& part){return part.folder+"\t"+FileHash(part.folder+"/complete.txt")+"\n";};
        const auto sources=origin(prefix)+(extend?origin(*suffix):"");
        Atomic(j.folder+"/angular_sources.tsv",sources);
        std::string lowerIntegrity;
        if(j.pairedDepth){
            const auto path=j.folder+"/result/result_depth"+Number(j.pairedDepth)+".dat";
            const auto stem=path.substr(0,path.size()-4);
            publish(lower,stem);
            std::string manifest=path+" "+FileHash(path)+"\n";
            if(j.hasControls){const auto cv=stem+"_analytic_facets.tsv";manifest+=cv+" "+FileHash(cv)+"\n";}
            Atomic(j.folder+"/integral_outputs.fnv",manifest);
            lowerIntegrity=" "+Hash(manifest);
        }
        PublishDerivedComplete(j,j.identity+" "+Number(j.seconds)+" "+FileHash(j.data)+" "+(j.hasControls?FileHash(j.summary):"none")+" "+Hash(sources)+lowerIntegrity+"\n",files);
    }
    std::set<std::string> seen;ValidateStoredJob(j.folder,seen);j.matrix.clear();j.angles.clear();j.controls.clear();j=Load(j);CacheJob(j);return j;
}
