// Small deterministic graph fixture: no optics, no existing outputs modified.
#include "../src/FullAuto.cpp"
using namespace FullAuto;
void SetFixture(Engine& e,const std::string& root){
    e.root=root;e.signature="imported-parent-resume-fixture-v1";e.exe="/bin/false";
    e.theta={150,152.5,155,157.5,160};e.initial=512;e.minimum=256;e.maxCount=67108864;
    e.minPhi=1;e.maxPhi=64;e.phi=16;e.eps=.04;e.interpolation=.01;e.autoDepth=false;e.target=6;
    e.controls=false;e.sharedMeans=false;e.mirror=false;e.workerThreads=1;e.workerChunk=256;
    e.orientationBlock=4096;e.analyticMeanThreads=1;e.hostBudgetMb=256;e.gpus={"0"};e.gpuCandidates={"0"};
    e.refinementWindow=true;e.refinementLo=150;e.refinementHi=160;e.refinementEps=.04;
    e.forbidNativeWork=true;e.seeds.assign(production,production+8);
}
AdaptiveState EmptyState(const Engine& e){AdaptiveState s;s.counts=Budgets(3,std::vector<int>(e.theta.size(),256));
    s.counts[0].assign(e.theta.size(),512);s.previous=Matrix(e.theta.size());s.valid.assign(e.theta.size(),true);s.coarse.assign(e.theta.size(),false);s.knots={0,2,4};return s;}
bool Throws(const std::function<void()>& f){try{f();return false;}catch(const std::exception&){return true;}}
void PublishFinal(const std::string& root,const std::string& sources,const std::string& integral=""){
    Atomic(root+"/published.dat","150 1 0 0 0 1\n160 2 0 0 0 1\n");Atomic(root+"/final_stage_sources.tsv",sources);
    Atomic(root+"/integral_sources.fnv",integral);
    Atomic(root+"/final_outputs.fnv","published.dat "+FileHash(root+"/published.dat")+"\nfinal_stage_sources.tsv "+FileHash(root+"/final_stage_sources.tsv")+"\nintegral_sources.fnv "+FileHash(root+"/integral_sources.fnv")+"\n");
}
int main(int argc,char**argv){try{
    Require(argc==2,"private test output required");const std::string root=argv[1],parent=root+"/parent";Directory(parent);Directory(parent+"/jobs");
    ArgPP args;RunConfig config;Engine original(args,config);SetFixture(original,parent);
    Atomic(parent+"/configuration.txt",original.signature+"\n");WriteGrid(parent+"/requested_theta.csv",original.theta);original.SaveTuning();
    auto state=EmptyState(original);state.round=9;state.streak=1;state.auditSerial=8;
    for(size_t i=0;i<state.previous.size();++i){state.previous[i][0]=2+i;state.previous[i][5]=.25+i;}
    SaveAdaptive(original,state,"converged");
    std::ostringstream sources,integrals;Job victim;const size_t leaves=512;
    for(size_t i=0;i<leaves;++i){Job j;original.SetSampling(j,512,11,6,original.theta,false);j.identity="leaf"+Number(i);j.folder=parent+"/jobs/"+j.identity;j.data=j.folder+"/result/result.dat";j.summary=j.folder+"/result/result_analytic_facets.tsv";j.seconds=1;j.angles=original.theta;j.matrix=state.previous;
        Directory(j.folder+"/result");WriteMatrix(j.data,j.angles,j.matrix);Atomic(j.folder+"/complete.txt",j.identity+" 1 "+FileHash(j.data)+" none\n");
        sources<<"0 1 512 11 "<<j.folder<<'\n';integrals<<j.data<<'\t'<<FileHash(j.data)<<'\n';if(!i)victim=j;
    }
    PublishFinal(parent,sources.str(),integrals.str());
    const std::string originalCheckpoint=Read(parent+"/adaptive_checkpoint.tsv"),originalTuning=Read(parent+"/execution_tuning.tsv");
    size_t checkedImport=0;double importSeconds=0;
    {ValidationBatch scope;auto& batch=*activeValidationBatch;const size_t before=batch.checked;const double t=Seconds();Engine e(args,config);SetFixture(e,parent);PrepareRefinementWindow(e);checkedImport=batch.checked-before;importSeconds=Seconds()-t;}
    Require(checkedImport==leaves,"fresh import skipped parent dependencies");const auto child=parent+"/refine_theta150_160_eps0p04";
    const std::string childCheckpoint=Read(child+"/adaptive_checkpoint.tsv"),policy=Read(child+"/refinement_window_policy.json"),commit=Read(child+"/refinement_window_import_complete.json");
    size_t checkedResume=0;double resumeSeconds=0;
    {ValidationBatch scope;auto& batch=*activeValidationBatch;const size_t before=batch.checked;const double t=Seconds();Engine e(args,config);SetFixture(e,parent);PrepareRefinementWindow(e);checkedResume=batch.checked-before;resumeSeconds=Seconds()-t;
        Require(e.theta==original.theta&&e.phi==original.phi,"metadata-only preparation changed recipe");}
    Require(checkedResume==0,"committed child still walked unused parent dependencies");
    Require(Read(child+"/adaptive_checkpoint.tsv")==childCheckpoint,"child checkpoint overwritten");
    // Matched warmed ABBA: the same completed parent, parser and published
    // file checks; only the deep graph pass differs. All saved values agree.
    double fullSeconds=0,metadataSeconds=0;size_t fullRuns=0,metadataRuns=0;
    for(bool full:{true,false,false,true,true,false,false,true}){
        ValidationBatch scope;Engine e(args,config);SetFixture(e,parent);auto restored=EmptyState(e);const double t=Seconds();
        Require(RestoreAdaptive(e,restored,false,full?CompletedRestoreValidation::FullDependencies:CompletedRestoreValidation::PublishedFilesOnly),"matched parent restore missing");
        const double seconds=Seconds()-t;if(full){fullSeconds+=seconds;++fullRuns;}else{metadataSeconds+=seconds;++metadataRuns;}
        Require(restored.previous==state.previous&&restored.counts==state.counts&&restored.knots==state.knots
            &&restored.round==state.round&&restored.auditSerial==state.auditSerial,"all16 or estimator identity changed in metadata restore");
    }
    Require(unlink((child+"/refinement_window_policy.json").c_str())==0,"remove child policy");
    size_t missingPolicyChecked=0;
    {ValidationBatch scope;auto& batch=*activeValidationBatch;const auto before=batch.checked;Engine e(args,config);SetFixture(e,parent);PrepareRefinementWindow(e);missingPolicyChecked=batch.checked-before;}
    Require(missingPolicyChecked==leaves&&Read(child+"/refinement_window_policy.json")==policy
        &&Read(child+"/adaptive_checkpoint.tsv")==childCheckpoint,"missing policy did not use full fallback without resetting child");
    // Source corruption is unrelated to metadata import, but must still be
    // rejected when its values are actually consumed by the child.
    const auto native=Read(victim.data);Atomic(victim.data,native+"corrupt\n");
    {Engine e(args,config);SetFixture(e,parent);PrepareRefinementWindow(e);Engine::DropPayload(victim);Require(Throws([&]{e.Materialize(victim);}),"child consumed corrupt inherited native values");}
    {Engine e(args,config);SetFixture(e,parent);auto s=EmptyState(e);Require(Throws([&]{RestoreAdaptive(e,s,false);}),"ordinary parent resume skipped dependency corruption");}
    // A missing import marker takes the full original path, even when a child
    // happens to have an apparently valid checkpoint already.
    Require(unlink((child+"/refinement_window_import_complete.json").c_str())==0,"remove fixture marker");
    {Engine e(args,config);SetFixture(e,parent);Require(Throws([&]{PrepareRefinementWindow(e);}),"incomplete import skipped corrupt parent dependency");}
    Atomic(child+"/refinement_window_import_complete.json",commit);Atomic(victim.data,native);
    size_t rejected=0;
    auto rejectFile=[&](const std::string& path,const std::string& bad){const auto saved=Read(path);Atomic(path,bad);Engine e(args,config);SetFixture(e,parent);Require(Throws([&]{PrepareRefinementWindow(e);}),"invalid metadata accepted: "+path);Atomic(path,saved);++rejected;};
    rejectFile(parent+"/published.dat","corrupt published matrix\n");
    rejectFile(child+"/refinement_window_policy.json",policy+" ");
    rejectFile(child+"/refinement_window_import_complete.json","{\"policy_fnv1a64\":\"wrong\"}");
    rejectFile(child+"/adaptive_checkpoint.tsv",childCheckpoint+"trailing corruption\n");
    {const auto end=childCheckpoint.find('\n');std::string body=childCheckpoint.substr(end+1);const auto at=body.find("active");Require(at!=std::string::npos,"fixture phase");body.replace(at,6,"active");body+="invalid_extra_field\n";rejectFile(child+"/adaptive_checkpoint.tsv","3 "+Hash(body)+"\n"+body);}
    Require(unlink((child+"/adaptive_checkpoint.tsv").c_str())==0,"remove child checkpoint");
    {Engine e(args,config);SetFixture(e,parent);Require(Throws([&]{PrepareRefinementWindow(e);}),"missing child checkpoint accepted");++rejected;}Atomic(child+"/adaptive_checkpoint.tsv",childCheckpoint);
    // Completed child's final dependency graph is still verified by ordinary
    // child RestoreAdaptive, even after optimized parent preparation.
    {Engine e(args,config);SetFixture(e,parent);PrepareRefinementWindow(e);auto s=EmptyState(e);Require(RestoreAdaptive(e,s,false),"child state missing");SaveAdaptive(e,s,"converged");PublishFinal(child,sources.str());}
    Atomic(victim.data,native+"corrupt again\n");
    {Engine e(args,config);SetFixture(e,parent);PrepareRefinementWindow(e);auto s=EmptyState(e);Require(Throws([&]{RestoreAdaptive(e,s,false);}),"completed child dependency validation bypassed");}
    Atomic(victim.data,native);
    Require(Read(parent+"/adaptive_checkpoint.tsv")==originalCheckpoint&&Read(parent+"/execution_tuning.tsv")==originalTuning,"parent solver state modified");
    std::cout<<"{\"status\":\"passed\",\"parent_native_dependencies\":"<<leaves<<",\"fresh_import_checked_nodes\":"<<checkedImport<<",\"committed_resume_checked_nodes\":"<<checkedResume<<",\"fresh_import_seconds\":"<<importSeconds<<",\"committed_resume_seconds\":"<<resumeSeconds<<",\"matched_full_parent_restore_seconds\":"<<fullSeconds/fullRuns<<",\"matched_metadata_parent_restore_seconds\":"<<metadataSeconds/metadataRuns<<",\"missing_policy_fallback_checked_nodes\":"<<missingPolicyChecked<<",\"all16_metadata_restore_exact\":true,\"invalid_metadata_cases_rejected\":"<<rejected<<",\"ordinary_parent_corruption_rejected\":true,\"actual_child_native_corruption_rejected\":true,\"completed_child_corruption_rejected\":true,\"interrupted_import_full_validation\":true,\"parent_checkpoint_and_tuning_unchanged\":true,\"native_optical_jobs\":0}\n";return 0;
}catch(const std::exception&e){std::cerr<<e.what()<<'\n';return 1;}}
