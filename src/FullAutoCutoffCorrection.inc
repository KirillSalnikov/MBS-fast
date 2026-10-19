// Included inside Engine. A fixed control-variate estimator for beam pruning:
// B_N + (O_m - B_m), m=max(64,N/8), on the SAME scramble and prefix offset.
// Its expectation is the unpruned target. Never select m from these samples.
// Each combined scramble remains one observation in all16/Fieller statistics.
bool CorrectCutoffJob(const Job& j)const{
    return balancedCorrection&&j.count>=512&&j.depth>=6&&!j.paired
        &&(phiCutoffPilot||std::find(seeds.begin(),seeds.end(),j.seed)!=seeds.end())
        &&j.executionHash==ExecutionHash();
}
std::vector<std::string> CutoffChildCommand(const std::vector<std::string>& original,int n,int depth,int pairedDepth,bool balanced,const std::string& output){
    auto c=original;
    auto remove=[&](const std::string& name){for(size_t i=0;i<c.size();)if(c[i]==name){Require(i+1<c.size(),"truncated cutoff child command");c.erase(c.begin()+i,c.begin()+i+2);}else ++i;};
    auto set=[&](const std::string& name,const std::string& value){remove(name);c.push_back(name);c.push_back(value);};
    for(size_t i=0;i+2<c.size();++i)if(c[i]=="--sobol-seed")c[i+1]=Number(n);
    set("--max-reflections",Number(depth));remove("--paired-depth");if(pairedDepth)set("--paired-depth",Number(pairedDepth));
    if(balanced){for(const char* name:{"--cutoff-profile","--beam-cutoff","--beam-cutoff-jones","--beam-cutoff-area","--beam-cutoff-importance","--trace-cutoff","--trace-cutoff-jones","--trace-cutoff-area","--trace-cutoff-importance","--beam-area-ratio"})remove(name);
        set("--cutoff-profile","balanced");}
    set("--output",output);return c;
}
void WriteCutoffMatrix(Job& job,const std::vector<Job>& parts,const std::vector<std::string>& command,int depth,const std::string& stem){
    Require(parts.size()==3,"cutoff correction needs three terms");
    const auto& x=job.inputAngles;Matrix matrix(x.size(),Mueller{});std::vector<std::array<double,5>> cv(job.hasControls?x.size():0,std::array<double,5>{});
    double refinement=0,setup=0;
    for(size_t i=0;i<3;++i){const auto& p=parts[i];const int sign=i==2?-1:1;
        Require(p.seed==job.seed&&p.offset==job.offset&&p.depth==depth&&p.angles.size()==x.size()&&p.hasControls==job.hasControls,"incompatible cutoff correction source");
        for(size_t r=0;r<x.size();++r){Require(std::fabs(p.angles[r]-x[r])<1e-10,"cutoff source theta differs");for(int j=0;j<16;++j)matrix[r][j]+=sign*p.matrix[r][j];
            if(job.hasControls)for(int j=0;j<5;++j)cv[r][j]+=sign*p.controls[r][j];}
        refinement+=p.controlRefinement;setup+=p.controlSetup;
    }
    WriteMatrix(stem+".dat",x,matrix);
    if(job.hasControls){std::vector<std::array<double,2>> coefficients(x.size(),{{1,1}});
        for(size_t i=0;i+1<command.size();++i)if(command[i]=="--analytic-control-weights"){
            std::istringstream input(Read(command[i+1]));std::string line;std::getline(input,line);
            for(size_t r=0;r<x.size();++r){double angle;Require(bool(input>>angle>>coefficients[r][0]>>coefficients[r][1])&&std::fabs(angle-x[r])<1e-10,"cutoff control weights differ");}}
        std::ostringstream out;out<<std::setprecision(17)<<"theta_deg\traw_M11\treflection_mean\tshadow_mean\tsampled_reflection\tsampled_shadow\tresidual_mean\thybrid_M11\trefinement\tsetup_seconds\treflection_weight\tshadow_weight\n";
        for(size_t r=0;r<x.size();++r){const auto& c=cv[r];const auto& b=coefficients[r];out<<x[r];for(double v:c)out<<'\t'<<v;
            out<<'\t'<<c[0]-b[0]*c[3]-b[1]*c[4]<<'\t'<<matrix[r][0]<<'\t'<<refinement<<'\t'<<setup<<'\t'<<b[0]<<'\t'<<b[1]<<'\n';}
        Atomic(stem+"_analytic_facets.tsv",out.str());}
}
bool ExecuteCorrected(std::vector<Job>& result,std::vector<std::vector<std::string>>& commands,const std::vector<size_t>& todo,const std::function<void(size_t)>& onComplete){
    bool any=false;for(size_t i:todo)any=any||CorrectCutoffJob(result[i]);if(!any)return false;
    struct Recipe {size_t parent;std::vector<size_t> children;};std::vector<Recipe> recipes;
    std::vector<Job> native;std::vector<std::vector<std::string>> nativeCommands;std::vector<size_t> launch;
    std::map<size_t,size_t> direct;
    auto add=[&](const Job& parent,const std::vector<std::string>& command,int n,int depth,int pairedDepth,bool balanced)->size_t{
        Job j=parent;j.count=n;j.depth=depth;j.pairedDepth=pairedDepth;j.assembled=false;
        // Separate child identities from production targets: a corrected target
        // must never recursively consume its own correction source.
        j.identity=Hash(parent.identity+"\ncutoff-child-v1 "+Number(n)+" "+Number(depth)+" "+Number(pairedDepth)+" "+Number(balanced));
        j.executionHash=Hash(parent.executionHash+(balanced?" balanced-child-v1":" off-child-v1"));
        j.folder=root+"/jobs/c"+j.identity;j.data=j.folder+"/result/result.dat";j.summary=j.folder+"/result/result_analytic_facets.tsv";
        const size_t index=native.size();nativeCommands.push_back(CutoffChildCommand(command,n,depth,pairedDepth,balanced,j.folder+"/result"));
        if(cache.count(j.identity))j=Materialize(cache.at(j.identity));
        else if(Exists(j.folder+"/complete.txt")){std::set<std::string> seen;ValidateStoredJob(j.folder,seen);std::istringstream metadata(Read(j.folder+"/complete.txt"));std::string identity;metadata>>identity>>j.seconds;Require(identity==j.identity,"cutoff child identity changed");j=Load(j);CacheJob(j);}
        else{Directory(j.folder);if(Exists(j.folder+"/worker.pid"))Require(!WorkerAlive(j.folder+"/worker.pid"),"cutoff child still running");
            if(Exists(j.folder+"/result"))Require(rename((j.folder+"/result").c_str(),(j.folder+"/result_incomplete_"+Number(getpid())+"_"+Number(index)).c_str())==0,"cannot preserve cutoff child partial output");
            launch.push_back(index);}
        native.push_back(j);return index;
    };
    for(size_t i:todo){const auto& j=result[i];
        if(!CorrectCutoffJob(j)){direct[i]=native.size();launch.push_back(native.size());native.push_back(j);nativeCommands.push_back(commands[i]);continue;}
        Recipe recipe;recipe.parent=i;const int m=std::max(64,j.count/8);
        bool pairedBalanced=j.pairedDepth!=0;
        // Resume any existing five-leg transaction with its original children.
        // New transactions share both balanced depth traces and need three legs.
        if(pairedBalanced)for(int n:{j.count,m})for(int depth:{j.depth,j.pairedDepth}){
            const auto identity=Hash(j.identity+"\ncutoff-child-v1 "+Number(n)+" "+Number(depth)+" 0 1");
            if(cache.count(identity)||Exists(root+"/jobs/c"+identity))pairedBalanced=false;
        }
        recipe.children.push_back(add(j,commands[i],j.count,j.depth,pairedBalanced?j.pairedDepth:0,true));
        recipe.children.push_back(add(j,commands[i],m,j.depth,j.pairedDepth,false));
        recipe.children.push_back(add(j,commands[i],m,j.depth,pairedBalanced?j.pairedDepth:0,true));
        if(j.pairedDepth&&!pairedBalanced){recipe.children.push_back(add(j,commands[i],j.count,j.pairedDepth,0,true));recipe.children.push_back(add(j,commands[i],m,j.pairedDepth,0,true));}
        recipes.push_back(recipe);
        Log("cutoff_correction_start","{\"folder\":"+Quoted(j.folder)+",\"profile\":\"balanced\",\"physical_target\":\"off\",\"count\":"+Number(j.count)+",\"correction_count\":"+Number(m)+",\"seed\":"+Number(j.seed)+",\"offset\":"+Number(j.offset)+",\"depth\":"+Number(j.depth)+",\"paired_depth\":"+Number(j.pairedDepth)+",\"native_trace_jobs\":"+Number(recipe.children.size())+",\"paired_balanced_depths\":"+(pairedBalanced?"true":"false")+",\"fraction_selected_before_samples\":true}");
    }
    // Each parent is committed immediately after its last native source. The
    // executor has already fed the next GPU request before calling this code.
    auto assemble=[&](const Recipe& recipe){Job& j=result[recipe.parent];const auto& k=recipe.children;
        const double assemblyBegin=Seconds();
        std::vector<Job> upper{native[k[0]],native[k[1]],native[k[2]]};
        Directory(j.folder+"/result");WriteCutoffMatrix(j,upper,commands[recipe.parent],j.depth,j.folder+"/result/result");
        std::ostringstream sources,signedSources,integrity;
        for(size_t i=0;i<3;++i)signedSources<<j.depth<<' '<<(i==2?-1:1)<<' '<<upper[i].folder<<'\n';
        if(j.pairedDepth){Job low=LowerView(native[k[1]],j.pairedDepth);
            std::vector<Job> lower=k.size()==3?std::vector<Job>{LowerView(native[k[0]],j.pairedDepth),low,LowerView(native[k[2]],j.pairedDepth)}
                :std::vector<Job>{native[k[3]],low,native[k[4]]};const auto stem=j.folder+"/result/result_depth"+Number(j.pairedDepth);
            WriteCutoffMatrix(j,lower,commands[recipe.parent],j.pairedDepth,stem);
            for(size_t i=0;i<3;++i){signedSources<<j.pairedDepth<<' '<<(i==2?-1:1)<<' '<<lower[i].folder<<'\n';
                // Views are durable inputs too; validate their provenance on resume.
                sources<<lower[i].folder<<'\t'<<FileHash(lower[i].folder+"/complete.txt")<<'\n';}
            std::ostringstream pair;pair<<stem+".dat"<<' '<<FileHash(stem+".dat")<<'\n';if(j.hasControls)pair<<stem+"_analytic_facets.tsv"<<' '<<FileHash(stem+"_analytic_facets.tsv")<<'\n';
            Atomic(j.folder+"/paired_depth.tsv",pair.str());integrity<<j.folder+"/paired_depth.tsv"<<' '<<Hash(pair.str())<<'\n';}
        j.seconds=0;for(size_t index:k){const auto& part=native[index];j.seconds+=part.seconds;sources<<part.folder<<'\t'<<FileHash(part.folder+"/complete.txt")<<'\n';}
        Atomic(j.folder+"/cutoff_sources.tsv",signedSources.str());integrity<<j.folder+"/cutoff_sources.tsv"<<' '<<Hash(signedSources.str())<<'\n';
        Atomic(j.folder+"/integral_outputs.fnv",integrity.str());Atomic(j.folder+"/angular_sources.tsv",sources.str());
        std::ostringstream replay;replay<<exe<<'\n';for(const auto& c:commands[recipe.parent])replay<<c<<'\n';Atomic(j.folder+"/command.args",replay.str());
        Atomic(j.folder+"/stdout.log","Balanced control-variate assembly: B_N + O_m - B_m. command.args replays the physical unpruned target; cutoff_sources.tsv records the actual signed estimator.\nHard tree-limit hits: 0\nOrientations still incomplete: 0\n");
        j.assembled=true;Atomic(j.folder+"/complete.txt",j.identity+" "+Number(j.seconds)+" "+FileHash(j.data)+" "+(j.hasControls?FileHash(j.summary):"none")+" "+Hash(sources.str())+" "+Hash(integrity.str())+"\n");
        std::set<std::string> verified;ValidateStoredJob(j.folder,verified);j=Load(j);CacheJob(j);
        Log("cutoff_correction_complete","{\"folder\":"+Quoted(j.folder)+",\"seconds\":"+Number(j.seconds)+",\"assembly_seconds\":"+Number(Seconds()-assemblyBegin)+",\"streaming_assembly\":true,\"combined_replica_enters_accuracy_estimate\":true,\"old_off_blocks_compatible\":true}");
    };
    std::vector<bool> ready(native.size(),true),assembled(recipes.size(),false);
    for(size_t i:launch)ready[i]=false;
    std::vector<std::vector<size_t>> dependents(native.size());
    for(size_t i=0;i<recipes.size();++i)for(size_t child:recipes[i].children)dependents[child].push_back(i);
    auto finish=[&](size_t parent){
        if(assembled[parent])return;
        for(size_t child:recipes[parent].children)if(!ready[child])return;
        assemble(recipes[parent]);assembled[parent]=true;
    };
    for(size_t i=0;i<recipes.size();++i)finish(i); // interrupted assembly: children already durable
    const bool wasCollecting=collecting;collecting=false;correctionExecuting=true;
    try{ExecutePrepared(native,nativeCommands,launch,[&](size_t index){
        ready[index]=true;for(size_t parent:dependents[index])finish(parent);
        // This is a progress pump, not a parent-completion notification. It may
        // prepare the next bounded batch while current native work runs.
        if(onComplete)onComplete(index);
    });}catch(...){correctionExecuting=false;collecting=wasCollecting;throw;}
    correctionExecuting=false;collecting=wasCollecting;
    for(const auto& item:direct)result[item.first]=native[item.second];
    for(size_t i=0;i<recipes.size();++i)Require(assembled[i],"cutoff parent missing committed sources");
    return true;
}
