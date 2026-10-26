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
    remove("--orientation-prefix-output");
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
// Prefix products are part of the native transaction, including scalar
// diagnostics. A completion marker never blesses an unchecked prefix.
void AppendNativePrefixIntegrity(const Job& j,std::ostringstream& manifest)const{
    std::istringstream command(Read(j.folder+"/command.args"));std::string argument;int prefix=0;
    while(std::getline(command,argument))if(argument=="--orientation-prefix-output"){
        Require(bool(std::getline(command,argument)),"missing native prefix count");prefix=std::stoi(argument);}
    if(!prefix)return;
    Require(prefix>0&&prefix<j.count,"invalid native prefix count");
    for(int depth:{0,j.pairedDepth}){
        const auto stem=j.folder+"/result/result_prefix"+Number(prefix)+(depth?"_depth"+Number(depth):"");
        for(const char* suffix:{".dat","_integrals.tsv"}){const auto path=stem+suffix;manifest<<path<<' '<<FileHash(path)<<'\n';}
        if(j.hasControls){const auto path=stem+"_analytic_facets.tsv";manifest<<path<<' '<<FileHash(path)<<'\n';}
        if(depth){const auto path=stem+"_integrals_revision.txt";manifest<<path<<' '<<FileHash(path)<<'\n';}
        if(!j.pairedDepth)break;
    }
}
Job CutoffPrefixView(const Job& source,int count){
    Job job=source;job.count=count;job.assembled=true;
    job.identity=Hash(source.identity+"\nembedded cutoff prefix v1 "+Number(count));job.folder=root+"/jobs/v"+job.identity;
    job.data=job.folder+"/result/result.dat";job.summary=job.folder+"/result/result_analytic_facets.tsv";
    job.angles.clear();job.matrix.clear();job.controls.clear();job.seconds=std::max(1e-9,source.seconds*count/source.count);
    const auto found=cache.find(job.identity);if(found!=cache.end())return Materialize(found->second);
    if(!Exists(job.folder+"/complete.txt")){
        std::set<std::string> verified;ValidateStoredJob(source.folder,verified);Directory(job.folder+"/result");
        AtomicFileBatch outputFiles(job.folder);
        std::ostringstream integrity;
        for(int depth:{0,job.pairedDepth}){
            const std::string suffix=depth?"_depth"+Number(depth):"";
            const auto from=source.folder+"/result/result_prefix"+Number(count)+suffix,to=job.folder+"/result/result"+suffix;
            for(const char* extension:{".dat","_integrals.tsv"}){Atomic(to+extension,Read(from+extension));integrity<<to+extension<<' '<<FileHash(to+extension)<<'\n';}
            if(job.hasControls){const auto extension="_analytic_facets.tsv";Atomic(to+extension,Read(from+extension));integrity<<to+extension<<' '<<FileHash(to+extension)<<'\n';}
            if(depth){const auto extension="_integrals_revision.txt";Atomic(to+extension,Read(from+extension));integrity<<to+extension<<' '<<FileHash(to+extension)<<'\n';}
            if(!job.pairedDepth)break;
        }
        std::istringstream sourceCommand(Read(source.folder+"/command.args"));std::string argument;std::getline(sourceCommand,argument);
        std::vector<std::string> original;while(std::getline(sourceCommand,argument))original.push_back(argument);
        const auto replay=CutoffChildCommand(original,count,job.depth,job.pairedDepth,true,job.folder+"/result");
        std::ostringstream recipe;recipe<<exe<<'\n';for(const auto& part:replay)recipe<<part<<'\n';Atomic(job.folder+"/command.args",recipe.str());
        Atomic(job.folder+"/stdout.log","First orientations saved during balanced native pass; no repeated trace or diffraction.\nHard tree-limit hits: 0\nOrientations still incomplete: 0\n");
        const auto sources=source.folder+"\t"+FileHash(source.folder+"/complete.txt")+"\n";
        Atomic(job.folder+"/angular_sources.tsv",sources);Atomic(job.folder+"/integral_outputs.fnv",integrity.str());
        outputFiles.Commit();
        PublishComplete(job,job.identity+" "+Number(job.seconds)+" "+FileHash(job.data)+" "+(job.hasControls?FileHash(job.summary):"none")+" "+Hash(sources)+" "+Hash(integrity.str())+"\n");
    }
    std::set<std::string> verified;ValidateStoredJob(job.folder,verified);job=Load(job);CacheJob(job);return job;
}
// FP32 accumulation is order-sensitive when the prefix boundary splits a
// prepared chunk. Spend one quarter of the existing execution reserve on an
// independent, cached numerical check of the FINAL corrected estimator.
std::map<std::string,bool> cutoffPrefixChecks;
bool VerifyCutoffPrefix(const Job& job,const std::vector<std::string>& actualCommand){
    if(!NativeWorkerCapabilities::OrientationPrefixOutput())return false;
#ifndef MBS_GPU_FP32
    (void)job;(void)actualCommand;return true;
#else
    if(!config.useGpu)return true;
    const int n=std::min(job.count,1024),m=std::max(64,n/8);
    std::ostringstream profile;profile<<"cutoff-prefix-v1 "<<binaryHash<<' '<<job.executionHash<<' '<<job.count
        <<' '<<job.depth<<' '<<job.pairedDepth<<' '<<job.phi<<' '<<job.block<<' '<<job.warp<<' '<<job.pipeline
        <<' '<<job.mirrored<<' '<<job.hasControls<<' '<<workerChunk<<' '<<Number(eps)
        <<' '<<FullAutoPreparation::Bucket(job.inputAngles.size())<<' '<<int(job.inputAngles.front()/15)<<' '<<int(job.inputAngles.back()/15);
    for(const auto& row:weights)for(double value:row)profile<<' '<<Number(value);
    for(const auto& row:polarWeights)for(const auto& component:row)for(double value:component)profile<<' '<<Number(value);
    const auto key=Hash(profile.str()),path=root+"/cutoff_prefix_validation_"+key+".json";
    const auto known=cutoffPrefixChecks.find(key);if(known!=cutoffPrefixChecks.end())return known->second;
    if(Exists(path)){const auto saved=CompactJson(Read(path));const bool passed=saved.find("\"passed\":true")!=std::string::npos;
        cutoffPrefixChecks[key]=passed;return passed;}
    const bool oldCollecting=collecting,oldCorrecting=correctionExecuting,oldVerifying=shortPacketVerifying;
    collecting=false;correctionExecuting=true;shortPacketVerifying=true;
    double error=0,baselineSeconds=0,candidateSeconds=0;size_t comparisons=0;const double started=Seconds();
    auto pilot=[&](unsigned seed,const std::string& role,int count,bool balanced,bool embedded){
        Job j=job;j.count=count;j.offset=0;j.seed=seed;j.assembled=false;j.pid=-1;j.seconds=0;j.matrix.clear();j.controls.clear();j.angles.clear();
        j.identity=Hash(key+" "+role+" "+Number(seed));j.executionHash=Hash(job.executionHash+" cutoff-prefix-pilot "+key+" "+role);
        j.folder=root+"/jobs/v"+j.identity;j.data=j.folder+"/result/result.dat";j.summary=j.folder+"/result/result_analytic_facets.tsv";
        const auto found=cache.find(j.identity);if(found!=cache.end())return Materialize(found->second);
        if(Exists(j.folder+"/complete.txt")){std::set<std::string> visited;ValidateStoredJob(j.folder,visited);
            std::istringstream saved(Read(j.folder+"/complete.txt"));std::string identity;saved>>identity>>j.seconds;Require(identity==j.identity,"prefix pilot identity differs");j=Load(j);CacheJob(j);return j;}
        Directory(j.folder);
        if(Exists(j.folder+"/worker.pid"))Require(!WorkerAlive(j.folder+"/worker.pid"),"prefix pilot still running");
        if(Exists(j.folder+"/result"))Require(rename((j.folder+"/result").c_str(),(j.folder+"/result_incomplete_"+Number(getpid())).c_str())==0,"cannot preserve prefix pilot output");
        auto command=CutoffChildCommand(actualCommand,count,j.depth,j.pairedDepth,balanced,j.folder+"/result");
        for(size_t i=0;i+2<command.size();++i)if(command[i]=="--sobol-seed")command[i+2]=Number(seed);
        for(size_t i=0;i+1<command.size();)if(command[i]=="--sobol-offset")command.erase(command.begin()+i,command.begin()+i+2);else ++i;
        if(embedded){command.push_back("--orientation-prefix-output");command.push_back(Number(m));}
        std::vector<Job> jobs{j};std::vector<std::vector<std::string>> commands{command};ExecutePrepared(jobs,commands,{0});return jobs[0];
    };
    std::array<unsigned,2> pilotSeeds{{8803u,8819u}};
    for(size_t i=0;i<pilotSeeds.size();++i)while(std::find(seeds.begin(),seeds.end(),pilotSeeds[i])!=seeds.end()||(i&&pilotSeeds[i]==pilotSeeds[0]))++pilotSeeds[i];
    try{
        for(unsigned seed:pilotSeeds){
            Job full,small,off,candidate;
            if(seed==pilotSeeds[0]){full=pilot(seed,"full",n,true,false);small=pilot(seed,"small",m,true,false);off=pilot(seed,"off",m,false,false);candidate=pilot(seed,"embedded",n,true,true);}
            else {candidate=pilot(seed,"embedded",n,true,true);off=pilot(seed,"off",m,false,false);small=pilot(seed,"small",m,true,false);full=pilot(seed,"full",n,true,false);}
            const Job embedded=CutoffPrefixView(candidate,m);
            baselineSeconds+=full.NativeWall()+small.NativeWall();candidateSeconds+=candidate.NativeWall();
            for(int depth:{job.depth,job.pairedDepth}){if(!depth)continue;
                auto data=[&](const Job& source){return Data({depth==job.depth?source:LowerView(source,depth)})[0];};
                const auto a=data(full),b=data(small),o=data(off),c=data(candidate),p=data(embedded);
                Require(a.size()==c.size()&&a.size()==o.size(),"prefix pilot theta count changed");
                for(size_t row=0;row<a.size();++row){Mueller reference{},result{};
                    for(int k=0;k<16;++k){reference[k]=a[row][k]+o[row][k]-b[row][k];result[k]=c[row][k]+o[row][k]-p[row][k];}
                    if(!(reference[0]>0&&result[0]>0)){error=1e100;continue;}
                    for(int k=0;k<16;++k){if(!std::isfinite(reference[k])||!std::isfinite(result[k])){error=1e100;continue;}
                        error=std::max(error,std::fabs(reference[k]-result[k])/reference[0]);
                        error=std::max(error,std::fabs(reference[k]/reference[0]-result[k]/result[0]));}
                    ++comparisons;
                }
            }
        }
    }catch(...){collecting=oldCollecting;correctionExecuting=oldCorrecting;shortPacketVerifying=oldVerifying;throw;}
    collecting=oldCollecting;correctionExecuting=oldCorrecting;shortPacketVerifying=oldVerifying;
    const double tolerance=eps/4096;
    const bool passed=comparisons&&error<=tolerance&&candidateSeconds<baselineSeconds*.97;
    const auto report="{\"passed\":"+std::string(passed?"true":"false")+",\"binary_fnv1a64\":"+Quoted(binaryHash)
        +",\"execution_profile\":"+Quoted(key)+",\"count_per_seed\":"+Number(n)+",\"independent_seeds\":["+Number(pilotSeeds[0])+","+Number(pilotSeeds[1])+"]"
        +",\"theta_rows\":"+Number(job.inputAngles.size())+",\"phi\":"+Number(job.phi)+",\"depth\":"+Number(job.depth)+",\"paired_depth\":"+Number(job.pairedDepth)
        +",\"max_all16_and_ratios_difference\":"+Number(error)+",\"numerical_tolerance\":"+Number(tolerance)
        +",\"baseline_processing_seconds\":"+Number(baselineSeconds)+",\"prefix_processing_seconds\":"+Number(candidateSeconds)
        +",\"pilot_elapsed_seconds\":"+Number(Seconds()-started)+",\"fallback\":\"original_three_leg_estimator\",\"production_sampling_unchanged\":true}";
    Record(path,report);Log("cutoff_prefix_validation",report);cutoffPrefixChecks[key]=passed;return passed;
#endif
}
// A completed native source can still be read after a controller-only upgrade
// removes the ability to produce new embedded prefixes. Never reinterpret or
// overwrite such a source; validate all prefix products before retaining it.
bool CompleteEmbeddedCutoffSource(const Job& parent,int prefix){
    Job source=parent;
    source.identity=Hash(parent.identity+"\ncutoff-child-v1 "+Number(parent.count)+" "+Number(parent.depth)+" "+Number(parent.pairedDepth)+" 1 prefix-v1 "+Number(prefix));
    source.folder=root+"/jobs/c"+source.identity;
    const auto found=cache.find(source.identity);if(found!=cache.end())source=found->second;
    if(!Exists(source.folder+"/complete.txt")){
        Require(found==cache.end(),"indexed embedded cutoff source lost its completion marker");return false;
    }
    std::set<std::string> visited;ValidateStoredJob(source.folder,visited);
    std::istringstream marker(Read(source.folder+"/complete.txt"));std::string identity;marker>>identity;
    Require(identity==source.identity,"embedded cutoff source identity changed");
    const auto command=Read(source.folder+"/command.args");
    Require(command.find("\n--orientation-prefix-output\n"+Number(prefix)+"\n")!=std::string::npos,
        "completed embedded cutoff source lacks its prefix recipe");
    std::ostringstream products;AppendNativePrefixIntegrity(source,products);return true;
}
bool RestoreCutoffPrefixRecipe(const Job& parent,int prefix){
    const auto path=parent.folder+"/cutoff_prefix_recipe.txt",recipe=Read(path);
    if(recipe=="three-leg-v1\n")return false;
    Require(recipe=="embedded-prefix-v1\n","unsupported cutoff prefix recipe");
    if(NativeWorkerCapabilities::OrientationPrefixOutput()||CompleteEmbeddedCutoffSource(parent,prefix))return true;
    Require(!Exists(parent.folder+"/complete.txt")&&!cache.count(parent.identity),
        "cannot migrate a completed cutoff estimator");
    // Commit the explanation and the original recipe first. If interrupted,
    // retrying this migration is safe and all partial native files remain.
    const auto saved=parent.folder+"/cutoff_prefix_recipe.before_fallback.txt";
    if(Exists(saved))Require(Read(saved)==recipe,"saved cutoff prefix recipe changed");else Atomic(saved,recipe);
    const auto receipt="{\"policy\":\"unavailable-native-prefix-three-leg-v1\",\"parent_identity\":"+Quoted(parent.identity)
        +",\"native_prefix_output\":false,\"old_recipe\":\"embedded-prefix-v1\",\"new_recipe\":\"three-leg-v1\""
        +",\"estimator\":\"B_N + O_m - B_m\",\"count\":"+Number(parent.count)+",\"prefix_count\":"+Number(prefix)
        +",\"completed_embedded_source\":false,\"partial_sources_preserved\":true}";
    Record(parent.folder+"/cutoff_prefix_fallback.json",receipt);
    Atomic(path,"three-leg-v1\n");Log("cutoff_prefix_recipe_fallback",receipt);return false;
}
bool ExecuteCorrected(std::vector<Job>& result,std::vector<std::vector<std::string>>& commands,const std::vector<size_t>& todo,const std::function<void(size_t)>& onComplete){
    bool any=false;for(size_t i:todo)any=any||CorrectCutoffJob(result[i]);if(!any)return false;
    IndexAppendBatch indexWrites(*this);
    struct Recipe {size_t parent;bool embeddedPrefix=false;std::vector<size_t> children;};std::vector<Recipe> recipes;
    std::vector<Job> native;std::vector<std::vector<std::string>> nativeCommands;std::vector<size_t> launch;
    std::map<size_t,size_t> direct;
    auto add=[&](const Job& parent,const std::vector<std::string>& command,int n,int depth,int pairedDepth,bool balanced,int prefix)->size_t{
        Job j=parent;j.count=n;j.depth=depth;j.pairedDepth=pairedDepth;j.assembled=false;
        // Separate child identities from production targets: a corrected target
        // must never recursively consume its own correction source.
        j.identity=Hash(parent.identity+"\ncutoff-child-v1 "+Number(n)+" "+Number(depth)+" "+Number(pairedDepth)+" "+Number(balanced)+(prefix?" prefix-v1 "+Number(prefix):""));
        j.executionHash=Hash(parent.executionHash+(balanced?" balanced-child-v1":" off-child-v1"));
        j.folder=root+"/jobs/c"+j.identity;j.data=j.folder+"/result/result.dat";j.summary=j.folder+"/result/result_analytic_facets.tsv";
        const size_t index=native.size();nativeCommands.push_back(CutoffChildCommand(command,n,depth,pairedDepth,balanced,j.folder+"/result"));
        if(prefix){nativeCommands.back().push_back("--orientation-prefix-output");nativeCommands.back().push_back(Number(prefix));}
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
        bool prefixEligible=NativeWorkerCapabilities::OrientationPrefixOutput();
        for(const char* flag:{"--no-shadow-output","--analytic-backscatter","--analytic-facet-samples","--analytic-azimuth-gaussian"})
            if(std::find(commands[i].begin(),commands[i].end(),flag)!=commands[i].end())prefixEligible=false;
        // Existing interrupted transactions retain their exact native recipe.
        // The durable choice is written before dispatching any of its children.
        const auto prefixRecipe=j.folder+"/cutoff_prefix_recipe.txt";
        if(Exists(prefixRecipe))recipe.embeddedPrefix=RestoreCutoffPrefixRecipe(j,m);
        else if(prefixEligible&&(!j.pairedDepth||pairedBalanced)){
            bool existing=false;
            for(int n:{j.count,m})for(int depth:{j.depth,j.pairedDepth})if(depth)for(bool balanced:{false,true})for(int pair:{0,j.pairedDepth}){
                const auto identity=Hash(j.identity+"\ncutoff-child-v1 "+Number(n)+" "+Number(depth)+" "+Number(pair)+" "+Number(balanced));
                existing=existing||cache.count(identity)||Exists(root+"/jobs/c"+identity);
            }
            if(!existing&&VerifyCutoffPrefix(j,commands[i])){Atomic(prefixRecipe,"embedded-prefix-v1\n");recipe.embeddedPrefix=true;}
        }
        // An embedded recipe was created only with paired balanced depths.
        // Unrelated incomplete ordinary children must not change that identity.
        if(recipe.embeddedPrefix)pairedBalanced=j.pairedDepth!=0;
        recipe.children.push_back(add(j,commands[i],j.count,j.depth,pairedBalanced?j.pairedDepth:0,true,recipe.embeddedPrefix?m:0));
        recipe.children.push_back(add(j,commands[i],m,j.depth,j.pairedDepth,false,0));
        if(recipe.embeddedPrefix){
            // The view becomes ready only after its native source has a durable
            // completion marker. It consumes no worker slot.
            recipe.children.push_back(native.size());Job view=j;view.count=m;native.push_back(view);nativeCommands.emplace_back();
        }else recipe.children.push_back(add(j,commands[i],m,j.depth,pairedBalanced?j.pairedDepth:0,true,0));
        if(j.pairedDepth&&!pairedBalanced){recipe.children.push_back(add(j,commands[i],j.count,j.pairedDepth,0,true,0));recipe.children.push_back(add(j,commands[i],m,j.pairedDepth,0,true,0));}
        recipes.push_back(recipe);
        Log("cutoff_correction_start","{\"folder\":"+Quoted(j.folder)+",\"profile\":\"balanced\",\"physical_target\":\"off\",\"count\":"+Number(j.count)+",\"correction_count\":"+Number(m)+",\"seed\":"+Number(j.seed)+",\"offset\":"+Number(j.offset)+",\"depth\":"+Number(j.depth)+",\"paired_depth\":"+Number(j.pairedDepth)+",\"native_trace_jobs\":"+Number(recipe.children.size()-(recipe.embeddedPrefix?1:0))+",\"native_prefix_output_capability\":"+(NativeWorkerCapabilities::OrientationPrefixOutput()?"true":"false")+",\"embedded_balanced_prefix\":"+(recipe.embeddedPrefix?"true":"false")+",\"paired_balanced_depths\":"+(pairedBalanced?"true":"false")+",\"fraction_selected_before_samples\":true}");
    }
    // Each parent is committed immediately after its last native source. The
    // executor has already fed the next GPU request before calling this code.
    auto assemble=[&](const Recipe& recipe){Job& j=result[recipe.parent];const auto& k=recipe.children;
        const double assemblyBegin=Seconds();
        std::vector<Job> upper{native[k[0]],native[k[1]],native[k[2]]};
        // Lower views can publish their own transactions. Finish those before
        // opening the parent's private write batch; never nest transactions.
        std::vector<Job> lower;
        if(j.pairedDepth){Job low=LowerView(native[k[1]],j.pairedDepth);
            lower=k.size()==3?std::vector<Job>{LowerView(native[k[0]],j.pairedDepth),low,LowerView(native[k[2]],j.pairedDepth)}
                :std::vector<Job>{native[k[3]],low,native[k[4]]};}
        Directory(j.folder+"/result");AtomicFileBatch outputFiles(j.folder);
        WriteCutoffMatrix(j,upper,commands[recipe.parent],j.depth,j.folder+"/result/result");
        std::ostringstream sources,signedSources,integrity;
        for(size_t i=0;i<3;++i)signedSources<<j.depth<<' '<<(i==2?-1:1)<<' '<<upper[i].folder<<'\n';
        if(j.pairedDepth){const auto stem=j.folder+"/result/result_depth"+Number(j.pairedDepth);
            WriteCutoffMatrix(j,lower,commands[recipe.parent],j.pairedDepth,stem);
            for(size_t i=0;i<3;++i){signedSources<<j.pairedDepth<<' '<<(i==2?-1:1)<<' '<<lower[i].folder<<'\n';
                // Views are durable inputs too; validate their provenance on resume.
                sources<<lower[i].folder<<'\t'<<FileHash(lower[i].folder+"/complete.txt")<<'\n';}
            std::ostringstream pair;pair<<stem+".dat"<<' '<<FileHash(stem+".dat")<<'\n';if(j.hasControls)pair<<stem+"_analytic_facets.tsv"<<' '<<FileHash(stem+"_analytic_facets.tsv")<<'\n';
            Atomic(j.folder+"/paired_depth.tsv",pair.str());integrity<<j.folder+"/paired_depth.tsv"<<' '<<Hash(pair.str())<<'\n';}
        j.seconds=0;for(size_t index:k){const auto& part=native[index];if(!recipe.embeddedPrefix||index!=k[2])j.seconds+=part.seconds;sources<<part.folder<<'\t'<<FileHash(part.folder+"/complete.txt")<<'\n';}
        Atomic(j.folder+"/cutoff_sources.tsv",signedSources.str());integrity<<j.folder+"/cutoff_sources.tsv"<<' '<<Hash(signedSources.str())<<'\n';
        Atomic(j.folder+"/integral_outputs.fnv",integrity.str());Atomic(j.folder+"/angular_sources.tsv",sources.str());
        std::ostringstream replay;replay<<exe<<'\n';for(const auto& c:commands[recipe.parent])replay<<c<<'\n';Atomic(j.folder+"/command.args",replay.str());
        Atomic(j.folder+"/stdout.log","Balanced control-variate assembly: B_N + O_m - B_m. command.args replays the physical unpruned target; cutoff_sources.tsv records the actual signed estimator.\nHard tree-limit hits: 0\nOrientations still incomplete: 0\n");
        j.assembled=true;outputFiles.Commit();
        PublishComplete(j,j.identity+" "+Number(j.seconds)+" "+FileHash(j.data)+" "+(j.hasControls?FileHash(j.summary):"none")+" "+Hash(sources.str())+" "+Hash(integrity.str())+"\n");
        std::set<std::string> verified;ValidateStoredJob(j.folder,verified);j=Load(j);CacheJob(j);
        Log("cutoff_correction_complete","{\"folder\":"+Quoted(j.folder)+",\"seconds\":"+Number(j.seconds)+",\"assembly_seconds\":"+Number(Seconds()-assemblyBegin)+",\"streaming_assembly\":true,\"combined_replica_enters_accuracy_estimate\":true,\"old_off_blocks_compatible\":true}");
    };
    std::vector<bool> ready(native.size(),true),assembled(recipes.size(),false);
    for(size_t i:launch)ready[i]=false;
    for(const auto& recipe:recipes)if(recipe.embeddedPrefix)ready[recipe.children[2]]=false;
    std::vector<std::vector<size_t>> dependents(native.size());
    for(size_t i=0;i<recipes.size();++i)for(size_t child:recipes[i].children)dependents[child].push_back(i);
    auto finish=[&](size_t parent){
        if(assembled[parent])return;
        const auto& recipe=recipes[parent];
        if(recipe.embeddedPrefix&&ready[recipe.children[0]]&&!ready[recipe.children[2]]){
            native[recipe.children[2]]=CutoffPrefixView(native[recipe.children[0]],std::max(64,native[recipe.children[0]].count/8));
            ready[recipe.children[2]]=true;
        }
        for(size_t child:recipe.children)if(!ready[child])return;
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
    indexWrites.Commit();
    return true;
}
