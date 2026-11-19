#define main cutoff_existing_probe_main
#include "fullauto_cutoff_correction_probe.cpp"
#undef main
int main(int argc,char** argv){try{
 Require(argc==3,"OUTPUT BINARY");ArgPP args;RunConfig config;
 for(bool incoherent:{false,true}){
  const auto root=std::string(argv[1])+(incoherent?"/incoherent":"/coherent");Engine old(args,config);Setup(old,root,argv[2],false,false,incoherent);old.gpuCandidates=old.gpus;old.balancedCorrection=true;old.SaveTuning();
  const auto hash=old.ExecutionHash();
  { Job dispatch; std::vector<std::string> command={"binary","--threads","1"};
    old.dispatchPipeline=true;old.ApplyDispatchThreads(dispatch,command);
    Require(std::count(command.begin(),command.end(),"--orientation-pipeline")==1,"dispatch pipeline missing");
    old.ApplyDispatchThreads(dispatch,command);Require(std::count(command.begin(),command.end(),"--orientation-pipeline")==1,"dispatch pipeline duplicated");
    Require(old.ExecutionHash()==hash,"dispatch pipeline changed physical identity");old.dispatchPipeline=false; }
  auto before=old.DirectStage(512,old.seeds.data(),8,8,old.theta);const auto marker=FileHash(before[0].folder+"/complete.txt");
  old.phiCutoffPilot=true;const auto pilotHash=old.ExecutionHash();old.phiCutoffPilot=false;
  old.collecting=true;old.DirectStage(4096,old.seeds.data(),1,8,old.theta,true,false,8192);old.collecting=false;
  const Job partial=old.pendingJobs.at(0);Atomic(partial.folder+"/cutoff_prefix_recipe.txt","embedded-prefix-v1\n");
  const auto oldChild=old.CutoffChildIdentity(partial,4096,8,0,true,512);
  Engine e(args,config);Setup(e,root,argv[2],false,false,incoherent);e.gpuCandidates=e.gpus;e.cutoffUpgrade=true;e.balancedCorrection=true;e.WarmCache();e.PrepareTuning(e.theta);
  Require(e.ExecutionHash()==hash,"upgrade changed physical target identity");e.phiCutoffPilot=true;Require(e.ExecutionHash()!=pilotHash,"upgrade reused old variance pilot");e.phiCutoffPilot=false;
  e.forbidNativeWork=true;auto replay=e.DirectStage(512,e.seeds.data(),8,8,e.theta);Require(replay[0].matrix==before[0].matrix&&FileHash(replay[0].folder+"/complete.txt")==marker,"upgrade reran or changed completed blocks");e.forbidNativeWork=false;
  Require(e.CutoffChildIdentity(partial,4096,8,0,true,512)!=oldChild,"upgrade mixed native pruning profiles");
  const auto added=e.DirectStage(512,e.seeds.data(),8,8,e.theta,true,false,512);
  for(const auto& j:added){std::istringstream in(Read(j.folder+"/cutoff_sources.tsv"));int depth,sign,terms=0;std::string path;while(in>>depth>>sign>>path)++terms;Require(terms==3&&e.StagedCutoff(),"equal staged prefixes did not cancel");}
  auto mixed=e.PrefixStage(1024,e.seeds.data(),8,8,e.theta);for(size_t i=0;i<8;++i)for(size_t r=0;r<e.theta.size();++r)for(int k=0;k<16;++k)Close(mixed[i].matrix[r][k],.5*(before[i].matrix[r][k]+added[i].matrix[r][k]),"mixed old/new prefix changed expectation weights");
  auto recovered=e.DirectStage(4096,e.seeds.data(),1,8,e.theta,true,false,8192);Require(recovered[0].folder==partial.folder,"physical parent recipe changed");
  Require(Read(partial.folder+"/cutoff_prefix_recipe.txt")=="embedded-prefix-v1\n"&&Exists(e.CutoffPrefixRecipePath(partial)),"old partial prefix recipe overwritten");
  e.forbidNativeWork=true;Require(e.PrefixStage(1024,e.seeds.data(),8,8,e.theta)[0].matrix==mixed[0].matrix,"mixed prefix replay changed");
 }
 std::cout<<"PASS staged upgrade: completed blocks retained, pilots isolated, mixed prefixes and old partial recipes safe\n";return 0;
}catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}}
