#include "../src/FullAuto.cpp"
using namespace FullAuto;
void Setup(Engine& e,const std::string& root,const std::string& binary,bool controls,bool pack,bool prefix,bool incoherent){
 e.root=root;e.signature="native-private-optics-v1";Directory(root+"/jobs");e.exe=binary;e.binaryHash=FileHash(binary);e.compactIndex=true;e.packJobs=pack;
 e.controls=controls;e.sharedMeans=controls;e.coherent=!incoherent;e.gpus={""};e.seeds={11};e.workerThreads=1;e.analyticMeanThreads=1;e.hostBudgetMb=512;
 e.phi=4;e.target=10;e.theta={150,153.5};e.initial=64;e.orientationBlock=64;e.workerPolicy="persistent";e.resourcePolicy="off";e.workQueuePolicy="on";
 e.base={"--method","po","--backend","cpu","--particle","1",".2",".1","--ri","1.3116","0","--wavelength-um",".532","--symmetry","1","1","--cutoff-profile","balanced"};
 if(prefix){e.base.push_back("--orientation-prefix-output");e.base.push_back("16");}if(incoherent)e.base.push_back("--incoherent");
 WriteGrid(root+"/requested_theta.csv",e.theta);
}
int main(int argc,char**argv){try{Require(argc==3,"fresh scratch path and worker required");int cases=0;const std::string root=argv[1];
 for(bool controls:{false,true})for(bool pair:{false,true})for(bool prefix:{false,true}){
  Matrix reference;std::vector<std::array<double,5>> referenceControls;std::map<std::string,std::string> outputs;
  for(bool packed:{false,true}){ArgPP args;RegisterCliOptions(args);RunConfig config;Engine e(args,config);const auto dir=root+"/c"+Number(cases)+(packed?"/pack":"/raw");Setup(e,dir,argv[2],controls,packed,prefix,true);
   const auto result=e.DirectStage(64,e.seeds.data(),1,pair?10:6,e.theta,true,false,0,pair?8:0);Require(result.size()==1,"native result size");const auto& j=result[0];
   std::set<std::string> verified;ValidateStoredJob(j.folder,verified);
   if(!packed){reference=j.matrix;referenceControls=j.controls;}else Require(reference==j.matrix&&referenceControls==j.controls,"packed real optics all16/controls changed");
   std::vector<std::string> names={"result.dat","result_integrals.tsv"};if(pair){names.push_back("result_depth8.dat");names.push_back("result_depth8_integrals.tsv");}if(prefix){names.push_back("result_prefix16.dat");names.push_back("result_prefix16_integrals.tsv");if(pair){names.push_back("result_prefix16_depth8.dat");names.push_back("result_prefix16_depth8_integrals.tsv");}}
   for(const auto& name:names){const auto text=Read(j.folder+"/result/"+name);if(!packed)outputs[name]=text;else Require(outputs.at(name)==text,"packed product changed: "+name);}
   for(auto& service:e.services)service.reset();Engine resumed(args,config);Setup(resumed,dir,argv[2],controls,packed,prefix,true);resumed.WarmCache();resumed.forbidNativeWork=true;const size_t before=resumed.cache.size();const auto restored=resumed.DirectStage(64,resumed.seeds.data(),1,pair?10:6,resumed.theta,true,false,0,pair?8:0);Require(restored[0].matrix==reference&&restored[0].controls==referenceControls&&resumed.cache.size()==before,"resume recalculated or changed estimator");
  }std::cout<<"PASS real incoherent controls="<<controls<<" pair="<<pair<<" prefix="<<prefix<<" all16/controls/products exact and no-work resume\n";++cases;
 }
 Matrix goReference;std::string goIntegrals;
 for(bool packed:{false,true}){ArgPP args;RegisterCliOptions(args);RunConfig config;Engine e(args,config);const auto dir=root+(packed?"/go_pack":"/go_raw");Setup(e,dir,argv[2],false,packed,false,false);e.go=true;e.phi=1;e.theta={0,90,180};e.base[1]="go";e.base.back()="off";
  const auto rows=GoRows(e);const auto jobs=e.DirectStage(64,e.seeds.data(),1,6,rows);const auto& j=jobs[0];Require(j.matrix.size()==5&&j.angles[3]==0&&j.angles[4]==180,"GO bins/poles missing");
  const auto ints=Read(j.folder+"/result/result_integrals.tsv");if(!packed){goReference=j.matrix;goIntegrals=ints;}else Require(j.matrix==goReference&&ints==goIntegrals,"packed GO all16/poles/integrals changed");
  for(auto& service:e.services)service.reset();Engine resumed(args,config);Setup(resumed,dir,argv[2],false,packed,false,false);resumed.go=true;resumed.phi=1;resumed.theta=e.theta;resumed.base=e.base;resumed.WarmCache();resumed.forbidNativeWork=true;
  Require(resumed.DirectStage(64,resumed.seeds.data(),1,6,rows)[0].matrix==goReference,"GO resume changed bins/poles");
 }
 ++cases;std::cout<<"PASS real GO all16 bins/poles/integrals exact and no-work resume\n";
 std::cout<<"PASS "<<cases<<" real native optics cases\n";return 0;
 }catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}}
