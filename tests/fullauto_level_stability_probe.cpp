// Saved real 8-seed matrices reproduce the backward-scattering false pass.
// See data/fullauto_level_stability.json for provenance and limitations.
#include "../src/FullAuto.cpp"
using namespace FullAuto;

void CloseStability(double actual,double expected,double scale,const char* why){
    Require(std::fabs(actual-expected)<=1e-12*std::max(1.,scale),why);
}

void CompatibilityReceipt(const std::string& root){
    Directory(root);const auto binary=root+"/solver",evidence=root+"/comparison.json";
    const std::string old="old-exact-hash",current="new-exact-hash";
    Require(!VerifiedCompatibleExecutable(binary,old,current),"uncertified binary pair accepted");
    Atomic(evidence,"{\"status\":\"passed\"}\n");
    const auto body="FULLAUTO_NATIVE_COMPATIBILITY_V1\n"+old+" "+current+" "+evidence+" "+FileHash(evidence)+"\n";
    Atomic(binary+".native-compatibility.tsv",Hash(body)+"\n"+body);
    Require(VerifiedCompatibleExecutable(binary,old,current),"certified exact binary pair rejected");
    Require(!VerifiedCompatibleExecutable(binary,"other",current)&&!VerifiedCompatibleExecutable(binary,old,"other"),
            "compatibility certificate leaked to a different binary");
    Atomic(evidence,"changed evidence\n");bool rejected=false;
    try{VerifiedCompatibleExecutable(binary,old,current);}catch(const std::exception&){rejected=true;}
    Require(rejected,"changed compatibility evidence accepted");
    Atomic(binary+".native-compatibility.tsv",Hash(body)+"\n"+body+"tamper");rejected=false;
    try{VerifiedCompatibleExecutable(binary,old,current);}catch(const std::exception&){rejected=true;}
    Require(rejected,"tampered compatibility receipt accepted");
}

void AdversarialLevelChanges(){
    Mueller total{},current{},alternate{};total[0]=100;
    // Each Mueller entry must independently prevent a pass. The test uses a
    // zero-M11 correction, whose normalized ratios are undefined; normalization
    // belongs to the complete physical estimator, never the correction alone.
    for(int j=0;j<16;++j){
        alternate.fill(0);alternate[j]=3;
        Require(LevelPrefixChange(total,current,alternate,false)>.02,
                "a Mueller correction entry was omitted from stability");
    }
    alternate.fill(0);alternate[1]=.1;
    CloseStability(LevelPrefixChange(total,current,alternate,true),.001,1,
                   "zero-intensity correction was incorrectly normalized alone");
    // A small raw change can move a normalized ratio by more than tolerance.
    // Large ratio is intentional: isolate the contract from physical fixtures.
    total[1]=1000;alternate.fill(0);alternate[0]=1;
    Require(LevelPrefixChange(total,current,alternate,false)<.02,
            "raw-matrix control for ratio-only fixture failed");
    Require(LevelPrefixChange(total,current,alternate,true)>.02,
            "full-estimator ratio instability was omitted");
    total[1]=0;alternate.fill(0);alternate[0]=-101;
    Require(LevelPrefixChange(total,current,alternate,true)>.02,
            "nonpositive counterfactual M11 accepted");
    alternate.fill(0);alternate[7]=std::numeric_limits<double>::quiet_NaN();
    Require(LevelPrefixChange(total,current,alternate,true)>.02,
            "nonfinite correction accepted");
    alternate.fill(0);alternate[0]=3;
    const double positive=LevelPrefixChange(total,current,alternate,true);
    alternate[0]=-3;
    Require(positive>.02&&LevelPrefixChange(total,current,alternate,true)>.02,
            "opposite unstable levels could cancel instead of being checked separately");
}

void JointLevelInstability(){
    ArgPP args;RunConfig config;Engine e(args,config);e.ratios=true;
    e.statistical=.00951171875;e.interpolation=.005;
    Samples total(8,Matrix(1)),joint(8,Matrix(1));
    for(size_t s=0;s<8;++s)total[s][0][0]=joint[s][0][0]=100;
    for(size_t l=0;l<3;++l){
        Samples current(8,Matrix(1)),alternate(8,Matrix(1));
        for(size_t s=0;s<8;++s){current[s][0][0]=l?0:100;alternate[s][0][0]=current[s][0][0]+.8;joint[s][0][0]+=.8;}
        Require(PairedLevelPrefixChange(e,total,current,alternate,0)<NumericalBudget(e),
                "joint fixture should have three individually small changes");
    }
    const double sum=PairedLevelPrefixChange(e,total,total,joint,0);
    CloseStability(sum,.024,1,"joint all-level change lost aligned corrections");
    Require(sum>NumericalBudget(e),"three aligned sub-tolerance changes were accepted together");
}

void CacheStabilityJob(Engine& e,int n,int depth,unsigned seed,const std::vector<double>& x,int offset=0){
    Job job;e.SetSampling(job,n,seed,depth,x,false);job.offset=offset;
    job.identity="stability-fixture-"+Number(e.cache.size());job.folder="/nonexistent/"+job.identity;
    // No payload, completion marker, process or GPU is available to this test.
    // Prefix promotion must depend only on complete compatible metadata.
    e.CacheJob(job);
}

void CachedPrefixPromotion(){
    for(int scenario=0;scenario<9;++scenario){
        ArgPP args;RunConfig config;Engine e(args,config);
        e.controls=false;e.compactIndex=false;e.forbidNativeWork=true;e.phi=32;
        e.theta={170,171,172};e.seeds.assign(production,production+8);e.maxCount=4096;
        const std::vector<int> depths={8,10};const Budgets current(2,std::vector<int>(3,512));
        std::vector<size_t> ids={0,1};
        if(scenario==4)ids={1};
        if(scenario==5)ids={0,1,2};
        if(scenario==6)e.thetaPhi[1]=64;
        if(scenario==7)e.controls=true;
        if(scenario==3)e.phi=64;
        const std::vector<double> angles=scenario==4?std::vector<double>{171}
            :scenario==5?std::vector<double>{170,171}:e.theta;
        for(int depth:depths)for(size_t seed=0;seed<8;++seed){
            if(scenario==1&&depth==8)continue;
            if(scenario==2&&depth==8&&seed==7)continue;
            CacheStabilityJob(e,2048,depth,e.seeds[seed],angles);
        }
        if(scenario==3)e.phi=32;
        if(scenario==8)e.angularCache=false;
        const size_t jobs=e.cache.size();const auto promoted=CachedCompletedBudgets(e,depths,current,ids);
        Require(e.cache.size()==jobs,"metadata promotion created optical work");
        if(scenario==0){
            Require(CompleteCachedPrefix(e,2048,8,angles,e.seeds.data()),"complete eight-seed prefix not detected");
            for(size_t l=0;l<2;++l)for(size_t id:ids)Require(promoted[l][id]==2048,"complete cached correction not promoted");
        }else if(scenario==5){
            for(size_t l=0;l<2;++l){Require(promoted[l][0]==2048&&promoted[l][1]==2048,
                "two-row cached group unnecessarily blocked");Require(promoted[l][2]==512,"uncached angle promoted");}
        }else if(scenario==6){
            for(size_t l=0;l<2;++l)Require(promoted[l][0]==2048&&promoted[l][1]==512,
                "cached promotion mixed different phi estimators");
        }else Require(promoted==current,"incomplete, incompatible or unavailable cached prefix promoted");
        if(scenario==2)Require(!CompleteCachedPrefix(e,2048,8,angles,e.seeds.data()),
                              "seven seeds mistaken for a complete estimator");
        if(scenario==4)Require(CompleteCachedPrefix(e,2048,8,angles,e.seeds.data()),
                              "singleton fixture has no prefix to test native padding");
    }
}

void CachedPaddingSelection(){
    for(int scenario=0;scenario<6;++scenario){
        ArgPP args;RunConfig config;Engine e(args,config);e.controls=scenario==5;
        e.compactIndex=false;e.forbidNativeWork=true;e.theta={170,171,172,180};
        e.phi=32;e.seeds.assign(production,production+8);
        if(scenario==3)e.phi=64;
        for(int depth:{8,10})for(size_t s=0;s<8;++s){
            if((scenario==1||scenario==4)&&depth==8)continue;
            const std::vector<double> angles=scenario==2&&depth==8&&s==7
                ?std::vector<double>{171}:std::vector<double>{171,172};
            CacheStabilityJob(e,2048,depth,e.seeds[s],angles);
        }
        // A smaller angle available to only one seed cannot pad all eight.
        for(int depth:{8,10})CacheStabilityJob(e,2048,depth,e.seeds[0],{170,171});
        if(scenario==3)e.phi=32;
        double padding=-1;
        const bool available=CachedPrefixPadding(e,2048,10,scenario==4?0:8,171,e.seeds.data(),padding);
        if(scenario==0||scenario==4)Require(available&&padding==172,
            "singleton padding was not selected from complete paired coverage");
        else Require(!available,"singleton padding accepted incomplete or incompatible replicas");
    }
}

void PrefixPolicyMigration(const std::string& root){
    ArgPP args;RunConfig config;Engine e(args,config);e.root=root;e.signature="level-prefix-policy-test";
    Directory(root);e.theta={170,171,172};e.initial=4096;e.minimum=256;e.maxCount=1<<20;
    Require(PrefixMinimumCount(e,0)==8192&&PrefixMinimumCount(e,1)==4096,
            "correction can accept its first quiet pilot or inherited an excessive base count");
    e.initial=512;Require(PrefixMinimumCount(e,1)==1024,"small-particle correction pilot floor wrong");
    e.initial=1<<24;Require(PrefixMinimumCount(e,2)==4096,"large-particle base count forced onto correction");
    e.initial=4096;
    AdaptiveState state;state.round=4;state.streak=2;state.auditSerial=7;state.knots={0,2};
    state.counts=Budgets(3,std::vector<int>(3,512));state.previous=Matrix(3);
    state.previous[0][0]=123;state.valid.assign(3,true);state.coarse.assign(3,true);
    e.independentPassed=true;const auto counts=state.counts;const auto previous=state.previous;
    Atomic(root+"/prefix_stability_policy.json","{\"policy\":\"legacy-base-only\"}");
    PreparePrefixStabilityPolicy(e,state,true);
    Require(state.streak==0&&!e.independentPassed&&state.coarse==std::vector<bool>(3,false),
            "legacy correction acceptance survived policy migration");
    Require(state.counts==counts&&state.previous==previous&&state.valid==std::vector<bool>(3,true)
            &&state.round==4&&state.auditSerial==7,"policy migration discarded computed work");
    Require(Exists(root+"/adaptive_checkpoint.tsv"),"policy migration was not persisted");
    state.streak=1;state.coarse.assign(3,true);PreparePrefixStabilityPolicy(e,state,true);
    Require(state.streak==1&&state.coarse==std::vector<bool>(3,true),"same policy repeatedly reset progress");
    LevelPrefixDecision decision;decision.required=counts;decision.required[1][1]=1024;
    Budgets proposed=counts;RequirePrefixGrowth(decision,{0,1,2},proposed);
    Require(proposed[1][1]==1024&&proposed[0]==counts[0]&&proposed[2]==counts[2]
            &&proposed[1][0]==512&&proposed[1][2]==512,"failed correction caused blanket refinement");
    e.forbidNativeWork=true;e.seeds.assign(production,production+8);
    AdaptiveData quiet;quiet.levels.assign(3,Samples(8,Matrix(3)));
    for(auto& sample:quiet.levels[0])for(auto& row:sample)row[0]=100;
    const auto gate=CheckLevelPrefixes(e,{8,10,12},state,quiet,{0,1,2});
    Require(!gate.passed&&gate.failed>=9&&state.coarse==std::vector<bool>(3,false),
            "quiet under-sampled levels inherited a stale accepted flag");
    for(const auto& level:gate.required)for(int count:level)Require(count==1024,
            "quiet first-pilot correction was not scheduled for growth");
    e.maxCount=512;const auto capped=CheckLevelPrefixes(e,{8,10,12},state,quiet,{0,1,2});
    Require(!capped.passed&&capped.required==state.counts,
            "exhausted orientation cap accepted a quiet pilot or exceeded the explicit cap");
}

void CacheCompleteMatrixRows(Engine& e,int n,int depth,const Samples& values){
    for(size_t s=0;s<8;++s){
        Job job;e.SetSampling(job,n,e.seeds[s],depth,e.theta,false);
        job.identity="complete-fixture-"+Number(e.cache.size());job.folder=e.root+"/jobs/"+job.identity;
        job.data=job.folder+"/result/result.dat";job.angles=e.theta;job.matrix=values[s];
        job.weightHash=Hash("");job.seconds=1;Directory(job.folder+"/result");WriteMatrix(job.data,e.theta,job.matrix);
        Atomic(job.folder+"/complete.txt",job.identity+" 1 "+FileHash(job.data)+" none\n");e.CacheJob(job);
        Atomic(job.folder+"/stdout.log","Hard tree-limit hits: 0\nOrientations still incomplete: 0\n");
    }
}

void PartialAngularPrevalidation(const std::string& root,const std::vector<double>& theta,const Samples& values){
    for(int corrupt=0;corrupt<2;++corrupt){
        ArgPP args;RunConfig config;Engine e(args,config);e.root=root+"/case"+Number(corrupt);
        e.signature="partial-angular-validation";Directory(e.root);Directory(e.root+"/jobs");
        e.exe="/bin/false";e.binaryHash=FileHash(e.exe);e.controls=false;e.compactIndex=false;
        e.forbidNativeWork=true;e.gpus={""};e.workerThreads=4;e.orientationBlock=64;
        e.seeds.assign(production,production+8);
        for(size_t first:{size_t(0),size_t(2)}){
            e.theta={theta[first],theta[first+1]};Samples part(8,Matrix(2));
            for(size_t s=0;s<8;++s)for(size_t r=0;r<2;++r)part[s][r]=values[s][first+r];
            CacheCompleteMatrixRows(e,2048,8,part);
        }
        e.theta=theta;e.payloads.Clear();
        if(corrupt){const auto& job=e.cache.begin()->second;Atomic(job.data,Read(job.data)+"corrupt\n");}
        bool failed=false;std::vector<Job> jobs;
        try{jobs=e.Stage(2048,e.seeds.data(),8,8,theta);}
        catch(const std::exception& error){failed=std::string(error.what()).find("changed")!=std::string::npos;if(!failed)throw;}
        Require(failed==bool(corrupt),"partial angular validation missed corrupted source or rejected intact data");
        if(!corrupt){
            Require(e.Data(jobs)==values,"parallel partial-source validation changed all16 replica values");
            cpu_set_t allowed;CPU_ZERO(&allowed);
            if(sched_getaffinity(0,sizeof(allowed),&allowed)==0&&CPU_COUNT(&allowed)>1)
                Require(Read(e.root+"/fullauto_history.jsonl").find("partial angular source rows")!=std::string::npos,
                        "partial angular source graphs were not submitted to parallel validation");
        }
        const auto history=e.root+"/fullauto_history.jsonl";
        Require(e.forbidNativeWork&&(!Exists(history)||Read(history).find("worker_start")==std::string::npos),
                "partial source validation launched optical work");
    }
}

void CachedPaddingRuntime(const std::string& root){
    ArgPP args;RunConfig config;Engine e(args,config);e.root=root;e.signature="cached-padding-runtime";
    Directory(root);Directory(root+"/jobs");e.exe="/bin/false";e.binaryHash=FileHash(e.exe);
    e.controls=false;e.compactIndex=false;e.forbidNativeWork=true;e.depthSharing=true;
    e.gpus={""};e.workerThreads=1;e.seeds.assign(production,production+8);
    e.theta={171,172};Samples lower(8,Matrix(2)),upper(8,Matrix(2));
    for(size_t s=0;s<8;++s)for(size_t r=0;r<2;++r){
        lower[s][r][0]=100+.1*s+r;lower[s][r][5]=20+.2*s;
        upper[s][r]=lower[s][r];upper[s][r][0]+=2;upper[s][r][5]-=.5;
    }
    CacheCompleteMatrixRows(e,2048,8,lower);CacheCompleteMatrixRows(e,2048,10,upper);
    e.theta={170,171,172,180};const Budgets counts(2,std::vector<int>(4,2048));
    const size_t jobs=e.cache.size();
    const auto data=EvaluateBudget(e,{8,10},counts,{1},e.seeds.data(),true);
    for(size_t s=0;s<8;++s)for(int j=0;j<16;++j){
        CloseStability(data.levels[0][s][0][j],lower[s][0][j],100,"cached singleton changed base data");
        CloseStability(data.levels[1][s][0][j],upper[s][0][j]-lower[s][0][j],100,"cached singleton changed paired correction");
    }
    Require(e.cache.size()==jobs&&e.forbidNativeWork,
            "singleton prefetch or paired stage launched optics despite available cached padding");
    Atomic(root+"/result.json","{\"target_theta\":171,\"cached_padding_theta\":172,\"absent_default_endpoints\":[170,180],\"native_work_forbidden\":true,\"depth_sharing_enabled\":true,\"all16_exact\":true}\n");
}

void JointPrefixGate(const std::string& root){
    for(int scenario=0;scenario<2;++scenario){
        ArgPP args;RunConfig config;Engine e(args,config);e.root=root+"/case"+Number(scenario);
        e.signature="joint-prefix-case"+Number(scenario);Directory(e.root);Directory(e.root+"/jobs");
        e.exe="/bin/false";e.binaryHash=FileHash(e.exe);e.controls=false;e.compactIndex=false;
        e.forbidNativeWork=true;e.depthSharing=false;e.gpus={""};e.workerThreads=1;
        e.theta={170,171};e.seeds.assign(production,production+8);e.initial=512;e.minimum=64;e.maxCount=4096;
        e.statistical=.00951171875;e.interpolation=.005;
        const std::vector<int> depths={8,10,12};
        for(size_t l=0;l<3;++l){Samples half(8,Matrix(2));
            for(auto& sample:half)for(auto& row:sample)row[0]=scenario?(l?100:103):100+.8*(l+1);
            CacheCompleteMatrixRows(e,512,depths[l],half);
        }
        AdaptiveState state;state.counts=Budgets(3,std::vector<int>(2,1024));
        state.coarse.assign(2,true);AdaptiveData stage;stage.levels.assign(3,Samples(8,Matrix(2)));
        for(auto& sample:stage.levels[0])for(auto& row:sample)row[0]=100;
        stage.costs.assign(3,std::vector<double>(2,1));stage.analytic.assign(2,0);
        const size_t jobs=e.cache.size();const auto decision=CheckLevelPrefixes(e,depths,state,stage,{0,1});
        Require(!decision.passed&&state.coarse==std::vector<bool>(2,false),
                scenario?"opposite failing levels canceled the prefix gate":"aligned small corrections escaped the joint prefix gate");
        bool growth=false;for(const auto& level:decision.required)for(int n:level)growth=growth||n>1024;
        Require(growth,"failed joint prefix gate did not schedule any correction");
        for(size_t l=0;l<3;++l)for(size_t r=0;r<2;++r){
            if(!scenario)Require(decision.errors[l][r]<NumericalBudget(e),"aligned fixture is not individually stable");
            else if(l<2)Require(decision.errors[l][r]>NumericalBudget(e),"cancellation fixture lacks failing individual levels");
        }
        Require(e.forbidNativeWork&&e.cache.size()==jobs,"joint gate launched optics despite complete half prefixes");
    }
}

void CachedRealCorrectionGate(const std::string& root,const std::vector<double>& theta,
        const std::vector<Budgets>& counts,const std::vector<std::vector<Samples>>& levels){
    ArgPP args;RunConfig config;Engine e(args,config);e.root=root;e.signature="real-correction-gate";
    Directory(root);Directory(root+"/jobs");e.exe="/bin/false";e.binaryHash=FileHash(e.exe);
    e.controls=false;e.compactIndex=false;e.forbidNativeWork=true;e.depthSharing=false;
    e.phi=32;e.seeds.assign(production,production+8);e.theta=theta;e.maxCount=1<<22;
    e.initial=512;e.minimum=256;
    e.statistical=.00951171875;e.interpolation=.005;e.workerThreads=1;
    const std::vector<int> depths={8,10,12};const std::vector<size_t> ids={0,1,2,3};
    // Fixture values already include analytic polarization controls. Recreate
    // cached upper/lower wrappers with exactly those observed differences.
    // The common baseline is a surrogate, not a new optical measurement.
    Samples baseline=levels[1][0];CacheCompleteMatrixRows(e,counts[1][0][0],depths[0],baseline);
    for(size_t l=1;l<3;++l){CacheCompleteMatrixRows(e,counts[1][l][0],depths[l-1],baseline);
        for(size_t s=0;s<8;++s)for(size_t r=0;r<theta.size();++r)for(int j=0;j<16;++j)
            baseline[s][r][j]+=levels[1][l][s][r][j];
        CacheCompleteMatrixRows(e,counts[1][l][0],depths[l],baseline);
    }
    AdaptiveState state;state.counts=counts[0];state.knots={0,2,3};state.previous=Mean(Combine(levels[0]));
    state.valid.assign(theta.size(),true);state.coarse.assign(theta.size(),true);state.streak=1;
    AdaptiveData stage;stage.levels=levels[0];stage.costs.assign(3,std::vector<double>(theta.size(),1));
    stage.analytic.assign(theta.size(),0);
    const auto before=state.counts;const size_t jobs=e.cache.size();
    Require(CheckCachedLevelPrefixes(e,depths,state,stage,ids,false),"saved outliers not invalidated by complete cached evidence");
    for(size_t l=0;l<3;++l)for(size_t r=0;r<theta.size();++r){
        const int expected=l==1&&r<2?counts[1][l][r]:before[l][r];
        Require(state.counts[l][r]==expected,"cached consistency check refined an unrelated angle/level");
    }
    Require(e.forbidNativeWork&&e.cache.size()==jobs&&state.streak==0
            &&!state.coarse[0]&&!state.coarse[1]&&state.coarse[2]&&state.coarse[3],
            "cached check launched work, lost caller prohibition, or invalidated unrelated acceptance");
    const std::string initial=Read(root+"/cached_prefix_consistency.csv");std::string durable;
    DIR* directory=opendir(root.c_str());Require(directory,"cannot inspect saved consistency evidence");
    while(auto* entry=readdir(directory)){const std::string name=entry->d_name;
        if(name.find("cached_prefix_consistency_round")==0&&name.size()>4&&name.substr(name.size()-4)==".csv"
                &&Read(root+"/"+name)==initial)durable=name;}
    closedir(directory);Require(!durable.empty(),"cached promotion has no immutable decision receipt");
    const std::string sources=durable.substr(0,durable.size()-4)+"_sources.tsv";
    Require(Read(root+"/"+sources)==Read(root+"/cached_prefix_consistency_sources.tsv")
            &&Read(root+"/fullauto_history.jsonl").find(Quoted(durable))!=std::string::npos,
            "immutable promotion sources or history link missing");
    Atomic(root+"/cached_prefix_initial_decision.csv",initial);
    // A second pass after correcting only the failed level must be quiet and
    // restore the eight-replica receipt rather than re-run optical workers.
    for(size_t r=0;r<2;++r)for(size_t s=0;s<8;++s)stage.levels[1][s][r]=levels[1][1][s][r];
    Require(!CheckCachedLevelPrefixes(e,depths,state,stage,ids,false),"stable corrected estimator repeatedly promoted");
    Require(e.forbidNativeWork&&e.cache.size()==jobs,"receipt reuse created new optical work");
    Require(Read(root+"/"+durable)==initial,"quiet comparison overwrote the triggering outlier evidence");
    // A genuinely quiet correction still needs the minimum prefix check.
    // Reuse already completed larger samples instead of launching floor work.
    e.initial=4096;state.counts=before;state.coarse.assign(theta.size(),true);state.streak=1;
    stage.levels=levels[1];state.previous=Mean(Combine(stage.levels));
    Require(CheckCachedLevelPrefixes(e,depths,state,stage,ids,false),"complete cached correction did not satisfy its minimum count for free");
    for(size_t l=0;l<3;++l)for(size_t r=0;r<theta.size();++r)
        Require(state.counts[l][r]==(l==1&&r<2?counts[1][l][r]:before[l][r]),
                "free minimum-count promotion changed unrelated budgets");
    Atomic(root+"/cached_prefix_floor_decision.csv",Read(root+"/cached_prefix_consistency.csv"));
    Require(e.forbidNativeWork&&e.cache.size()==jobs,"free prefix-floor promotion launched optical work");
    Require(Read(root+"/"+durable)==initial,"later floor promotion overwrote the triggering outlier evidence");

    // A large tolerance used to leave small angular teeth in the output even
    // though a longer, complete eight-seed estimate was already available.
    // Output promotion must depend on coverage, not observed significance.
    e.initial=512;e.statistical=100;e.interpolation=100;
    state.counts=before;state.coarse.assign(theta.size(),true);state.streak=1;
    stage.levels=levels[0];
    const auto expected=CachedCompletedBudgets(e,depths,state.counts,ids);
    Require(expected!=before,"quiet cache-promotion fixture has no longer samples");
    Require(!CheckCachedLevelPrefixes(e,depths,state,stage,ids,false),"quiet fixture unexpectedly fails the legacy discrepancy gate");
    const int round=state.round;
    Require(CheckCachedLevelPrefixes(e,depths,state,stage,ids,true),"sub-tolerance cached estimates were not adopted");
    Require(state.counts==expected&&state.round==round&&state.streak==1,
            "in-place output promotion restarted a round or discarded a prior pass");
    Require(Read(root+"/output_cache_policy.json").find(OutputCachePolicy())!=std::string::npos,
            "promoted output has no versioned policy receipt");
    const auto reference=EvaluateBudget(e,depths,expected,ids,e.seeds.data(),false);
    Require(stage.levels==reference.levels&&stage.sources==reference.sources,
            "promoted output did not use the complete corrected eight-seed/all16 estimator");
    const auto actual=Estimate(Combine(stage.levels)),checked=Estimate(Combine(reference.levels));
    Require(actual.mean==checked.mean&&actual.interval==checked.interval&&actual.ratioInterval==checked.ratioInterval,
            "promoted matrix and confidence/ratio intervals belong to different estimators");
    Require(!CheckCachedLevelPrefixes(e,depths,state,stage,ids,true),"output promotion loops without new cached work");
    Require(e.forbidNativeWork,"output promotion relaxed native-work prohibition");
}

int main(int argc,char** argv){try{
    Require(argc==3,"saved level-stability fixture and output directory required");
    std::ifstream in(argv[1]);size_t rows,levelCount,seeds,elements;
    Require(bool(in>>rows>>levelCount>>seeds>>elements)&&rows==4&&levelCount==3&&seeds==8&&elements==16,
            "unexpected saved level-stability fixture shape");
    const double tolerance=.01451171875;
    std::vector<double> realTheta;
    std::vector<Budgets> allCounts(2,Budgets(levelCount,std::vector<int>(rows)));
    std::vector<std::vector<Samples>> allLevels(2,std::vector<Samples>(levelCount,Samples(seeds,Matrix(rows))));
    // Independent NumPy replay: max over all matrix entries and all paired
    // ratio influence functions of |change| + t7 * standard_error(change).
    const double expectedPaired[4][3]={
        {.010366995419229285,.03781281838449419,.00950935076428766},
        {.010525353294982132,.0376207145048629,.009546665336958762},
        {.006040193856992987,.004151405438993786,.009242182029123901},
        {0,0,.006554939783129589}};
    std::cout<<"theta,low_stat95,base_change,correction1_change,correction2_change,unstable_levels\n";
    for(size_t r=0;r<rows;++r){
        double theta;Require(bool(in>>theta),"missing saved angle");
        realTheta.push_back(theta);
        std::vector<std::vector<int>> counts(2,std::vector<int>(levelCount));
        std::vector<std::vector<Samples>> levels(2,std::vector<Samples>(levelCount,Samples(seeds,Matrix(1))));
        std::vector<Statistics> stats(2);
        for(size_t kind=0;kind<2;++kind){
            for(int& n:counts[kind])Require(bool(in>>n)&&n>0,"invalid saved count");
            Samples combined(seeds,Matrix(1));
            for(size_t l=0;l<levelCount;++l)for(size_t s=0;s<seeds;++s)for(size_t j=0;j<elements;++j){
                double& value=levels[kind][l][s][0][j];Require(bool(in>>value),"incomplete saved replicas");
                combined[s][0][j]+=value;
            }
            stats[kind]=Estimate(combined);
            for(size_t l=0;l<levelCount;++l){allCounts[kind][l][r]=counts[kind][l];
                for(size_t s=0;s<seeds;++s)allLevels[kind][l][s][r]=levels[kind][l][s][0];}
            for(const Matrix* values:{&stats[kind].mean,&stats[kind].interval,&stats[kind].ratioInterval}){
                for(size_t j=0;j<elements;++j){double expected;Require(bool(in>>expected),"missing independent expectation");
                    CloseStability((*values)[0][j],expected,values==&stats[kind].mean?stats[kind].mean[0][0]:1,
                        "saved all-element/Fieller replay does not reproduce archived estimator");}
            }
        }
        // These saved 2% decisions were below the remaining numerical budget,
        // but the low-N correction at the first two angles moved by >3%.
        const double interval=std::max(stats[0].maxInterval,stats[0].maxRatioInterval);
        Require(interval<tolerance,"fixture no longer reproduces optimistic statistical acceptance");
        ArgPP args;RunConfig config;Engine e(args,config);e.ratios=true;
        Samples reference(seeds,Matrix(1));for(size_t l=0;l<levelCount;++l)
            for(size_t s=0;s<seeds;++s)for(size_t j=0;j<elements;++j)reference[s][0][j]+=levels[0][l][s][0][j];
        std::vector<double> change(levelCount);size_t unstable=0,pairedUnstable=0;
        for(size_t l=0;l<levelCount;++l){
            change[l]=LevelPrefixChange(stats[0].mean[0],Mean(levels[0][l])[0],Mean(levels[1][l])[0],true);
            unstable+=change[l]>tolerance;
            Samples alternate=reference;
            for(size_t s=0;s<seeds;++s)for(size_t j=0;j<elements;++j)
                alternate[s][0][j]+=levels[1][l][s][0][j]-levels[0][l][s][0][j];
            const double paired=e.PairedBias(alternate,reference);
            CloseStability(paired,expectedPaired[r][l],1,"paired full-estimator stability differs from independent NumPy replay");
            CloseStability(PairedLevelPrefixChange(e,reference,levels[0][l],levels[1][l],0),paired,1,
                           "per-level stability helper did not compare the complete physical estimator");
            pairedUnstable+=paired>tolerance;
        }
        Require(pairedUnstable==unstable,"paired-noise criterion masked bad correction or forced stable neighbors");
        if(r<2){
            Require(unstable==1&&change[1]>.03,"bad small-N depth correction was not isolated");
            Require(change[0]<.005&&change[2]<.005,"stable levels unnecessarily rejected");
        }else Require(unstable==0,"well-sampled neighboring angle unnecessarily rejected");
        std::cout<<Number(theta)<<','<<Number(interval);
        for(double value:change)std::cout<<','<<Number(value);
        std::cout<<','<<unstable<<'\n';
    }
    std::string trailing;Require(!(in>>trailing),"unparsed fixture data");
    CompatibilityReceipt(std::string(argv[2])+"/binary_compatibility");
    AdversarialLevelChanges();
    JointLevelInstability();
    CachedPrefixPromotion();
    CachedPaddingSelection();
    PrefixPolicyMigration(argv[2]);
    CachedPaddingRuntime(std::string(argv[2])+"/padding_runtime");
    JointPrefixGate(std::string(argv[2])+"/joint");
    CachedRealCorrectionGate(std::string(argv[2])+"/cached_real",realTheta,allCounts,allLevels);
    PartialAngularPrevalidation(std::string(argv[2])+"/partial_angular",realTheta,allLevels[1][0]);
    std::cout<<"PASS: four real angles; all16/Student-t/Fieller replay; selective cached outlier correction and free floor; joint aligned/canceling levels; policy migration, stale acceptance, ratios, dark corrections, nonfinite; eight-seed/depth/phi coverage and singleton padding; no optical work\n";
    return 0;
}catch(const std::exception& error){std::cerr<<error.what()<<'\n';return 1;}}
