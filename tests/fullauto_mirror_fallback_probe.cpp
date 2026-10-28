// Native regression: simulate an analytic-budget exception, then exercise the
// actual controller recovery and real ordinary-PO mirror verification.
#include "../src/FullAuto.cpp"
using namespace FullAuto;

void SetupMirrorFallback(Engine& e,const std::string& root,const std::string& binary){
    e.root=root;e.signature="mirror-fallback-regression-v1";Directory(root);Directory(root+"/jobs");
    e.exe=binary;e.binaryHash=FileHash(binary);e.compactIndex=true;
    e.controls=true;e.controlPolicy="auto";e.sharedMeans=false;e.depthSharing=false;
    e.gpus={""};e.workerThreads=1;e.analyticMeanThreads=1;e.workerPolicy="process";
    e.resourcePolicy="off";e.hostBudgetMb=512;e.ratios=true;
    e.eps=.04;e.interpolation=.01;e.depthBudget=.01;
    e.phi=4;e.fixedPhi=4;e.minPhi=1;e.maxPhi=64;e.pilot=16;e.maxCount=64;
    e.target=6;e.initial=64;e.minimum=16;e.orientationBlock=64;
    e.theta={.5,45,90,170,179.5};e.seeds.assign(production,production+8);
    e.base={"--method","po","--backend","cpu","--particle","4","1",
        "--ri","1.3116","0","--wavelength-um",".532","--symmetry","2","6","--cutoff-profile","off"};
    WriteGrid(root+"/requested_theta.csv",e.theta);e.SetMirror(true);
}

const std::runtime_error analyticFailure("fullauto: analytic mean numerical budget exceeded");

void VerificationKeys(const std::string& root,const std::string& binary){
    ArgPP args;RunConfig config;Engine e(args,config);SetupMirrorFallback(e,root,binary);
    const auto analytic=e.MirrorVerificationKey(e.theta);e.controls=false;
    const auto plain=e.MirrorVerificationKey(e.theta);Require(plain!=analytic,"analytic proof authorized plain PO");
    e.phi=8;Require(plain!=e.MirrorVerificationKey(e.theta),"phi change reused mirror proof");e.phi=4;
    e.target=8;Require(plain!=e.MirrorVerificationKey(e.theta),"depth change reused mirror proof");e.target=6;
    e.binaryHash="different-native";Require(plain!=e.MirrorVerificationKey(e.theta),"native change reused mirror proof");e.binaryHash=FileHash(binary);
    e.controls=true;e.weightTheta=e.theta;e.weights.assign(e.theta.size(),{{1,1}});
    const auto weighted=e.MirrorVerificationKey(e.theta);e.weights[0][0]=.5;
    Require(weighted!=e.MirrorVerificationKey(e.theta),"changed analytic coefficient reused proof");
    e.controls=true;e.controlPolicy="on";Require(!e.RecoverUnsupportedControls(analyticFailure,e.theta,true)&&e.controls,"explicit analytic control overridden");
    e.controlPolicy="auto";Require(!e.RecoverUnsupportedControls(std::runtime_error("GPU worker failed"),e.theta,true)&&e.controls,"unrelated native failure swallowed");
    std::cout<<"PASS proof context and explicit/error priority\n";
}

void NativeRecovery(const std::string& root,const std::string& binary,const std::string& policy,bool numericalRejection){
    ArgPP args;RunConfig config;Engine e(args,config);SetupMirrorFallback(e,root,binary);e.mirrorPolicy=policy;
    if(numericalRejection)e.eps=1e-20; // Deliberately impossible tolerance tests fail-closed behavior.
    if(policy=="off")e.SetMirror(false);
    e.verifiedMirrorGrids.insert(e.MirrorVerificationKey(e.theta));
    Require(e.RecoverUnsupportedControls(analyticFailure,e.theta,true),"analytic fallback not handled");
    Require(!e.controls,"analytic controls remained enabled");
    const bool expected=policy!="off"&&!numericalRejection;
    Require(e.mirror==expected,"mirror recovery ignored the fresh verification or explicit off");
    Require(Read(root+"/analytic_fallback_mirror.json").find(expected?"\"enabled\":true":"\"enabled\":false")!=std::string::npos,"missing recovery decision");
    if(expected){
        Require(e.mirrorError<=e.eps/256,"plain PO mirror exceeded unchanged budget");
        Require(Read(root+"/mirror_verification.json").find("\"analytic_controls\":false")!=std::string::npos,"recovery reused analytic verification");
        // The restored recipe must keep its round and sample allocations and
        // prove equivalence entirely from the existing immutable native jobs.
        AdaptiveState saved;saved.round=7;saved.knots={0,2,4};saved.counts=Budgets(3,std::vector<int>(e.theta.size(),16));
        saved.previous=Matrix(e.theta.size());saved.valid.assign(e.theta.size(),false);saved.coarse.assign(e.theta.size(),false);
        SaveAdaptive(e,saved);const auto checkpoint=Read(root+"/adaptive_checkpoint.tsv");
        Engine resumed(args,config);SetupMirrorFallback(resumed,root,binary);resumed.mirrorPolicy=policy;
        AdaptiveState restored;restored.counts=Budgets(3,std::vector<int>(e.theta.size(),16));restored.previous=Matrix(e.theta.size());
        restored.valid.assign(e.theta.size(),false);restored.coarse.assign(e.theta.size(),false);
        Require(RestoreAdaptive(resumed,restored)&&restored.round==7&&restored.counts==saved.counts&&!resumed.controls&&resumed.mirror,"recovery state changed on resume");
        resumed.forbidNativeWork=true;Require(resumed.VerifyMirror(resumed.theta),"resumed plain PO failed cached mirror audit");
        Require(Read(root+"/adaptive_checkpoint.tsv")==checkpoint,"mirror audit rewrote adaptive checkpoint");
    }
    std::cout<<"PASS native droxtal recovery policy="<<policy<<" reject="<<numericalRejection<<" mirror="<<e.mirror<<"\n";
}

void AsymmetricRecovery(const std::string& root,const std::string& binary,const std::string& particle){
    ArgPP args;RunConfig config;Engine e(args,config);SetupMirrorFallback(e,root,binary);
    e.base={"--method","po","--backend","cpu","--particle-file",particle,
        "--ri","1.3116","0","--wavelength-um",".532","--symmetry","1","1","--cutoff-profile","off"};
    Require(e.RecoverUnsupportedControls(analyticFailure,e.theta,true)&&!e.mirror&&!e.controls,"asymmetric particle acquired a mirror shortcut");
    Require(Read(root+"/mirror_verification.json").find("\"geometry_verified\":false")!=std::string::npos,"asymmetric geometry not checked");
    std::cout<<"PASS asymmetric geometry keeps full gamma\n";
}

void ExplicitRejection(const std::string& root,const std::string& binary){
    ArgPP args;RunConfig config;Engine e(args,config);SetupMirrorFallback(e,root,binary);
    e.mirrorPolicy="on";e.eps=1e-20;bool rejected=false;
    try{e.RecoverUnsupportedControls(analyticFailure,e.theta,true);}
    catch(const std::exception& error){rejected=std::string(error.what()).find("paired Haar mirror verification failed")!=std::string::npos;}
    Require(rejected&&!Exists(root+"/analytic_fallback_mirror.json"),"explicit mirror failure silently fell back");
    std::cout<<"PASS explicit mirror rejection remains an error\n";
}

int main(int argc,char** argv){try{
    Require(argc==4,"OUTPUT NATIVE_BINARY ASYMMETRIC_PARTICLE");
    VerificationKeys(std::string(argv[1])+"/keys",argv[2]);
    NativeRecovery(std::string(argv[1])+"/auto",argv[2],"auto",false);
    NativeRecovery(std::string(argv[1])+"/off",argv[2],"off",false);
    NativeRecovery(std::string(argv[1])+"/rejected",argv[2],"auto",true);
    NativeRecovery(std::string(argv[1])+"/explicit_on",argv[2],"on",false);
    AsymmetricRecovery(std::string(argv[1])+"/asymmetric",argv[2],argv[3]);
    ExplicitRejection(std::string(argv[1])+"/explicit_rejected",argv[2]);
    return 0;
}catch(const std::exception& error){std::cerr<<error.what()<<'\n';return 1;}}
