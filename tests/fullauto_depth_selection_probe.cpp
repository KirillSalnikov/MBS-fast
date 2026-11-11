// No native worker is launched: exercise the production paired-residual
// implementation directly on known eight-replica matrices.
#include "../src/FullAuto.cpp"
#include "../src/FullAutoDepthSelection.h"
using namespace FullAuto;
namespace D=FullAutoDepthSelection;
namespace {
void Close(double a,double b,const char* why){Require(std::isfinite(a)&&std::isfinite(b)&&std::fabs(a-b)<=1e-12*std::max(1.,std::fabs(b)),why);}
Samples Constant(size_t rows,double m11=100,double m22=20){
    Samples x(8,Matrix(rows));for(auto& sample:x)for(auto& row:sample){row[0]=m11;row[5]=m22;}return x;
}
void LocalResiduals(){
    Samples a=Constant(3),b=a;
    for(size_t s=0;s<8;++s){b[s][0][5]+=.1;b[s][1][5]+=.3+(s<4?-.2:.2);b[s][2][5]+=.6;}
    std::vector<DepthResidual> rows(4);const auto total=DepthDifference(a,b,false,3,true,nullptr,&rows);
    Require(rows.size()==3,"per-row output not reset/exactly theta-sized");
    const double ci=t95*.2/std::sqrt(7.);
    Close(rows[0].upper,.001,"constant row acquired another row's noise");
    Close(rows[1].upper,(.3+ci)/100,"local paired standard error incorrect");
    Close(rows[1].lower,(.3-ci)/100,"local paired lower bound incorrect");
    Close(rows[2].upper,.006,"second constant row acquired global noise");
    Close(total.upper,.006,"aggregate upper not worst local upper");
    Close(total.interval,ci/100,"aggregate interval not worst local interval");
    double maximum=0;for(const auto& row:rows)maximum=std::max(maximum,row.upper);
    Close(total.upper,maximum,"per-row aggregation differs from production result");
    Matrix scale=Mean(b);scale[1][0]=10;rows.clear();
    const auto normalized=DepthDifference(a,b,false,3,false,&scale,&rows);
    Close(rows[1].upper,(.3+ci)/10,"common reference normalization ignored");
    Close(normalized.upper,rows[1].upper,"common scale did not alter worst row");
}
void PairingAndRatios(){
    Samples a=Constant(1),b=a;
    for(size_t s=0;s<8;++s){const double shared=s<4?-20:20;a[s][0][0]+=shared;b[s][0][0]+=shared;a[s][0][5]+=shared;b[s][0][5]+=shared+.5;}
    const auto pair=DepthDifference(a,b,false,1,false);
    Close(pair.upper,.005,"common-noise cancellation lost between depths");
    Close(pair.interval,0,"independent rather than paired depth variance used");
    a=Constant(1,10,5);b=Constant(1,20,15);Matrix scale=Mean(b);scale[0][0]=1000;
    Close(DepthDifference(a,b,false,1,false,&scale).upper,.01,"raw scaling incorrect");
    Close(DepthDifference(a,b,false,1,true,&scale).upper,.25,"difference of normalized ratios omitted");
    // Same positive per-seed values at both depths: zero paired difference,
    // but M11's own interval crosses zero and ratios must remain unresolved.
    a=Constant(2,1,0);for(size_t s=0;s<8;++s)a[s][1][0]=s==0?1000:1;b=a;
    std::vector<DepthResidual> rows;
    const auto unresolved=DepthDifference(a,b,false,2,true,nullptr,&rows);
    Require(unresolved.upper>=1e50&&rows[1].upper>=1e50&&rows[0].upper==0,
            "unresolved M11 denominator accepted or spread to unrelated theta");
}
void NonfiniteAndEmptyRows(){
    for(double bad:{std::numeric_limits<double>::quiet_NaN(),std::numeric_limits<double>::infinity()}){
        auto a=Constant(2),b=a;b[0][1][5]=bad;std::vector<DepthResidual> rows;
        bool rejected=false;
        try{const auto total=DepthDifference(a,b,false,2,false,nullptr,&rows);rejected=total.upper>=1e50;}
        catch(const std::exception&){rejected=true;} // Mean() validates samples before residuals.
        Require(rejected,"nonfinite Mueller component silently accepted");
    }
    auto finiteA=Constant(1),finiteB=finiteA;for(auto& sample:finiteB)sample[0][5]+=10;
    for(double bad:{std::numeric_limits<double>::quiet_NaN(),std::numeric_limits<double>::infinity()}){
        Matrix scale=Mean(finiteB);scale[0][0]=bad;bool rejected=false;
        try{rejected=DepthDifference(finiteA,finiteB,false,1,false,&scale).upper>=1e50;}
        catch(const std::exception&){rejected=true;}
        Require(rejected,"invalid requested common scale erased a real residual");
    }
    auto a=Constant(2),b=a;for(auto& sample:b)sample[1][0]=0;
    Require(DepthDifference(a,b,false,2,true).upper>=1e50,"nonpositive PO M11 accepted");
    a=Constant(2);for(auto& sample:a)sample[1].fill(0);b=a;
    std::vector<DepthResidual> rows;
    Require(DepthDifference(a,b,true,1,true,nullptr,&rows).upper==0&&rows.size()==2,
            "zero unobserved GO pole changed row count or caused false failure");
    Require(DepthDifference(a,b,false,2,true).upper>=1e50,
            "GO empty-pole exception leaked into PO");
}
void DeepTailRegression(){
    auto n6=Constant(1,1,.2),n8=Constant(1,1,.201);
    auto n20=Constant(1,1,.4),n22=Constant(1,1,.401);
    const auto adjacent=DepthDifference(n6,n8,false,1,true);
    const auto anchor=DepthDifference(n6,n22,false,1,true);
    const auto guard=DepthDifference(n20,n22,false,1,true);
    Require(adjacent.upper<.01,"adjacent-pair shadow fixture invalid");
    Require(D::Classify(anchor.upper,anchor.lower,guard.upper,guard.lower,.01)==D::Fail,
            "small6/8 pair hid substantial6/22 tail");
    Close(guard.upper,.001,"20/22 guard fixture incorrect");
    auto unstable=Constant(1,1,.415);const auto badGuard=DepthDifference(n20,unstable,false,1,true);
    Require(D::Classify(0,0,badGuard.upper,badGuard.lower,.01)==D::ReferenceFail,
            "unsettled deep reference accepted a candidate");
}
void SetupDepth(Engine& e,const std::string& root,bool go=false){
    e.root=root;e.signature="depth-policy-resume-probe";Directory(root);Directory(root+"/jobs");
    e.exe="/bin/true";e.binaryHash=FileHash(e.exe);e.controls=false;e.sharedMeans=false;
    e.gpus={""};e.seeds.assign(production,production+8);e.workerThreads=1;e.analyticMeanThreads=1;
    e.hostBudgetMb=512;e.phi=1;e.target=6;e.theta={150,155,160};e.initial=8;e.maxCount=8;
    e.orientationBlock=8;e.depthPilot=8;e.maxTheta=3;e.depthStart=6;e.depthMaximum=18;
    e.autoDepth=true;e.depthSharing=false;e.forbidNativeWork=true;e.go=go;e.requestedOptions="{}";
    e.workerPolicy="persistent";e.resourcePolicy="off";e.workQueuePolicy="on";
    e.base={"--method",go?"go":"po","--backend","cpu","--particle","1",".2",".1","--ri","1.3116","0","--wavelength-um",".532","--symmetry","1","1","--cutoff-profile","off"};
    WriteGrid(root+"/requested_theta.csv",e.theta);
}
void PublishFixture(Engine& e,int selected,int reported,bool modern,bool productionState){
    SaveSelectedDepth(e,selected);
    const std::string report="{\"target_depth\":"+Number(reported)+",\"pilot_passed\":true"
        +(modern?",\"policy\":"+Quoted(D::Policy()):std::string())+"}";
    Atomic(e.root+"/depth_selection.json",report);
    if(modern)Atomic(e.root+"/automatic_depth_policy.tsv",e.signature+" "+D::Policy()+" "+Number(selected)+" "+Hash(report)+"\n");
    if(productionState)Atomic(e.root+(e.go?"/go_checkpoint.tsv":"/adaptive_checkpoint.tsv"),"preserved production checkpoint\n");
}
void ResumeRegressions(const std::string& root){
    for(bool go:{false,true})for(bool modern:{false,true}){
        ArgPP args;RunConfig config;Engine e(args,config);SetupDepth(e,root+"/resume"+Number(go)+Number(modern),go);
        PublishFixture(e,10,10,modern,true);
        const auto checkpoint=Read(e.root+(go?"/go_checkpoint.tsv":"/adaptive_checkpoint.tsv"));
        Require(SelectAutomaticDepth(e)&&e.target==10,"verified/legacy production depth not preserved");
        Require(e.forbidNativeWork&&e.cache.empty(),"resume unexpectedly computed worker data");
        Require(Read(e.root+(go?"/go_checkpoint.tsv":"/adaptive_checkpoint.tsv"))==checkpoint,"resume changed production checkpoint");
        const auto policy=Read(e.root+"/depth_selection_resume_policy.json");
        Require(policy.find(modern?D::Policy():"legacy-neighbor-pilot-retained")!=std::string::npos,"resume mislabeled depth criterion");
    }
    for(int variant=0;variant<4;++variant){
        ArgPP args;RunConfig config;Engine e(args,config);SetupDepth(e,root+"/mismatch"+Number(variant));
        PublishFixture(e,12,variant==0?10:12,true,true);
        if(variant==1)Atomic(e.root+"/automatic_depth_policy.tsv",e.signature+" "+D::Policy()+" 10 "+Hash(Read(e.root+"/depth_selection.json"))+"\n");
        if(variant==2)Require(unlink((e.root+"/automatic_depth_policy.tsv").c_str())==0,"fixture receipt removal failed");
        if(variant==3)Atomic(e.root+"/depth_selection.json",Read(e.root+"/depth_selection.json")+" ");
        bool blocked=false;try{SelectAutomaticDepth(e);}catch(const std::exception& error){blocked=std::string(error.what()).find("depth receipt is missing or unverified")!=std::string::npos;}
        Require(blocked,"mismatched/missing receipt entered existing production");
        Require(e.forbidNativeWork&&e.cache.empty()&&!Exists(e.root+"/depth_search_progress.json"),"invalid receipt began replacement pilot");
    }
    {
        ArgPP args;RunConfig config;Engine e(args,config);SetupDepth(e,root+"/manual");e.autoDepth=false;e.target=9;
        Atomic(e.root+"/automatic_depth.tsv","irrelevant malformed old automatic checkpoint");
        Require(SelectAutomaticDepth(e)&&e.target==9&&!Exists(e.root+"/depth_search_progress.json"),"manual depth did not bypass automatic search");
    }
    for(bool promoted:{false,true}){
        ArgPP args;RunConfig config;Engine e(args,config);SetupDepth(e,root+(promoted?"/promoted":"/legacy_unstarted"));
        PublishFixture(e,14,14,false,false);
        if(promoted)Atomic(e.root+"/depth_selection.json","{\"automatic\":true,\"target_depth\":14,\"reason\":\"promoted after fresh final depth validation\"}");
        bool guarded=false;try{SelectAutomaticDepth(e);}catch(const std::exception& error){guarded=std::string(error.what()).find("read-only probe would launch native optical work")!=std::string::npos;}
        Require(guarded,"unverified unstarted pilot did not require fresh scientific evidence");
        const auto progress=Read(e.root+"/depth_search_progress.json");
        Require(progress.find("\"candidate\":"+Number(promoted?14:6)+",")!=std::string::npos,
                "promotion reset below failed target or legacy pilot incorrectly trusted");
    }
}
void MaximumFallbackResume(const std::string& root){
    for(bool go:{false,true})for(int variant=0;variant<3;++variant){
        ArgPP args;RunConfig config;Engine e(args,config);SetupDepth(e,root+"/fallback"+Number(go)+Number(variant),go);
        const int selected=variant==2?e.depthMaximum-2:e.depthMaximum;
        SaveSelectedDepth(e,selected);
        const std::string report="{\"policy\":"+Quoted(D::Policy())+",\"target_depth\":"+Number(selected)+",\"pilot_passed\":false,\"fallback_maximum\":true}";
        Atomic(e.root+"/depth_selection.json",report);
        Atomic(e.root+"/automatic_depth_policy.tsv",e.signature+" "+D::Policy()+" "+Number(selected)+" "+Hash(report)+"\n");
        const std::string checkpoint=e.root+(go?"/go_checkpoint.tsv":"/adaptive_checkpoint.tsv");
        Atomic(checkpoint,"existing estimator must survive\n");
        if(variant==1)Atomic(e.root+"/depth_selection.json",report+" ");
        bool rejected=false;try{
            Require(SelectAutomaticDepth(e),"valid maximum fallback did not continue");
        }catch(const std::exception&){rejected=true;}
        Require(rejected==(variant!=0),"fallback accepted invalid receipt or wrong maximum");
        Require(Read(checkpoint)=="existing estimator must survive\n","fallback modified estimator");
        if(variant==0){
            Require(e.depthFallback&&e.target==e.depthMaximum,"fallback lost fixed depth marker");
            Require(e.DepthAccuracyFields().find("\"accuracy_conditional_on_fixed_depth\":true")!=std::string::npos,"fallback advertised full depth accuracy");
            const auto policy=Read(e.root+"/depth_selection_resume_policy.json");
            Require(policy.find("\"common_reference_verified\":false")!=std::string::npos,"fallback advertised verified reference");
        }
    }
}
void PromotionTransactions(const std::string& root){
    for(bool go:{false,true})for(int stop=0;stop<6;++stop){
        ArgPP args;RunConfig config;Engine e(args,config);SetupDepth(e,root+"/crash"+Number(go)+Number(stop),go);
        PublishFixture(e,10,10,true,true);Atomic(e.root+"/physics_validation.json","{\"depth_pass\":false,\"cutoff_pass\":true}");
        Atomic(e.root+"/jobs/preserved_native.txt","native bytes must survive all promotion phases\n");
        BeginDepthPromotion(e,10,12);auto plan=ReadDepthPromotionPlan(e);
        const auto archive=e.root+"/depth_attempts/n10";
        if(stop>0){const size_t count=stop==1?1:plan.files.size();for(size_t i=0;i<count;++i)
            Require(rename((e.root+"/"+plan.files[i].first).c_str(),(archive+"/"+plan.files[i].first).c_str())==0,"fixture partial rename failed");}
        if(stop>=3){plan.phase="publish";Atomic(e.root+"/depth_promotion_pending.tsv",DepthPromotionText(plan));}
        if(stop>=4){const auto payload=e.signature+" 12\n";Atomic(e.root+"/depth_promotion_floor.tsv",Hash(payload)+"\n"+payload);SaveSelectedDepth(e,12);}
        if(stop>=5)Atomic(e.root+"/depth_selection.json","{\"automatic\":true,\"target_depth\":12,\"reason\":\"promoted after fresh final depth validation\"}");
        Require(ResumeDepthPromotion(e)&&e.target==12&&DepthPromotionFloor(e)==12,"interrupted promotion did not finish at the new floor");
        Require(!Exists(e.root+"/depth_promotion_pending.tsv")&&!ResumeDepthPromotion(e),"completed promotion was replayed twice");
        for(const auto& file:plan.files)Require(FileHash(archive+"/"+file.first)==file.second,"promotion replaced an original archived file");
        Require(Read(e.root+"/jobs/preserved_native.txt")=="native bytes must survive all promotion phases\n","promotion changed native data");
        Require(!Exists(e.root+(go?"/go_checkpoint.tsv":"/adaptive_checkpoint.tsv")),"old physical estimator remained active");
        // The reason/report is optional for the monotone search floor now.
        Require(unlink((e.root+"/depth_selection.json").c_str())==0,"fixture report removal failed");
        if(!go){bool blocked=false;try{SelectAutomaticDepth(e);}catch(const std::exception& error){blocked=std::string(error.what()).find("read-only probe would launch native optical work")!=std::string::npos;}
            Require(blocked&&Read(e.root+"/depth_search_progress.json").find("\"candidate\":12,")!=std::string::npos,"lost report reset promotion search below its durable floor");}
    }
    for(int conflict=0;conflict<3;++conflict){
        ArgPP args;RunConfig config;Engine e(args,config);SetupDepth(e,root+"/conflict"+Number(conflict));
        PublishFixture(e,10,10,true,true);Atomic(e.root+"/physics_validation.json","{\"depth_pass\":false,\"cutoff_pass\":true}");
        BeginDepthPromotion(e,10,12);const auto plan=ReadDepthPromotionPlan(e);const auto first=e.root+"/"+plan.files.front().first;
        if(conflict==0)Atomic(first,"unrelated changed content\n");
        if(conflict==1)Atomic(e.root+"/depth_attempts/n10/"+plan.files.front().first,Read(first));
        if(conflict==2)Atomic(e.root+"/depth_promotion_pending.tsv",Read(e.root+"/depth_promotion_pending.tsv")+"x");
        const auto before=Read(first);bool rejected=false;try{ResumeDepthPromotion(e);}catch(const std::exception&){rejected=true;}
        Require(rejected&&Read(first)==before&&Exists(e.root+"/depth_promotion_pending.tsv"),"conflicting transaction mutated or discarded data");
        for(const auto& file:plan.files)Require(Exists(e.root+"/"+file.first),"failed preflight partially archived later files");
    }
    {
        ArgPP args;RunConfig config;Engine e(args,config);SetupDepth(e,root+"/normal");
        PublishFixture(e,10,10,true,true);e.target=10;
        Atomic(e.root+"/physics_validation.json","{\"depth_pass\":false,\"cutoff_pass\":true}");
        Require(PromoteAutomaticDepth(e,4)&&e.target==12&&DepthPromotionFloor(e)==12,"normal failed-audit promotion failed");
    }
}
}
int main(int argc,char** argv){try{
    LocalResiduals();PairingAndRatios();NonfiniteAndEmptyRows();DeepTailRegression();
    Require(argc==2,"private fixture output path required");ResumeRegressions(argv[1]);MaximumFallbackResume(std::string(argv[1])+"/maximum");PromotionTransactions(std::string(argv[1])+"/promotion");
    std::cout<<"PASS production DepthDifference: per-theta bounds, exact common reference scaling, paired covariance, normalized ratios, nonfinite/denominator guards, GO poles, hidden6/8 versus6/22 tail, manual/legacy/new receipt resume and interrupted publication guards, restartable promotion phases and conflict preservation\n";
    return 0;
}catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}}
