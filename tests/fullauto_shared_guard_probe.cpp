#include "../src/FullAuto.cpp"
using namespace FullAuto;
int main(int argc,char**argv){try{
 Require(argc==3,"OUTPUT BINARY");ArgPP args;RunConfig config;Engine e(args,config);
 e.root=argv[1];e.signature="shared-guards-v1";Directory(e.root);Directory(e.root+"/jobs");e.exe=argv[2];e.binaryHash=FileHash(e.exe);
 e.controls=true;e.sharedMeans=true;e.gpus={""};e.seeds.assign(production,production+8);e.workerThreads=1;e.analyticMeanThreads=1;e.hostBudgetMb=512;
 e.phi=1;e.target=6;e.initial=128;e.orientationBlock=32;e.workerPolicy="persistent";e.resourcePolicy="off";e.workQueuePolicy="on";
 e.base={"--method","po","--backend","cpu","--particle","1",".2",".1","--ri","1.3116","0","--wavelength-um",".532","--symmetry","1","1","--cutoff-profile","off"};
 std::vector<size_t> ids;AdaptiveState state;state.counts=Budgets(1,std::vector<int>(17,64));
 for(size_t i=0;i<17;++i){ids.push_back(i);e.theta.push_back(180.*i/16);if(i%2==0)state.knots.push_back(i);if(i>=10)state.counts[0][i]=128;}
 WriteGrid(e.root+"/requested_theta.csv",e.theta);
 std::vector<Job> compoundBlocks;
 const std::vector<double> left(e.theta.begin(),e.theta.begin()+9),right(e.theta.begin()+8,e.theta.end());
 for(int offset=0;offset<128;offset+=32){
  e.DirectStage(32,e.seeds.data(),8,6,left,true,false,offset);
  e.DirectStage(32,e.seeds.data(),8,6,right,true,false,offset);
  const auto views=e.Stage(32,e.seeds.data(),8,6,e.theta,true,false,offset);
  for(const auto& j:views){Require(j.assembled,"fixture needs compound blocks");compoundBlocks.push_back(j);}
 }
 e.cache.clear();e.samplingIndex.clear();e.folderIndex.clear();e.payloads.Clear();
 for(const auto& j:compoundBlocks)e.CacheJob(j);
 e.PrefixStage(128,e.seeds.data(),8,6,e.theta);e.payloads.SetLimit(1);e.forbidNativeWork=true;
 Require(FullAutoAuditPlan::Stencils(state.knots,state.counts,ids).size()>1,"fixture needs multiple groups");
 const auto separate=CommonPrefixInterpolationImpl(e,state,{6},ids,e.seeds.data());size_t reused=0;
 {
  ValidationBatch outer("shared-guard regression");const auto shared=CommonPrefixInterpolation(e,state,{6},ids,e.seeds.data());
  Require(separate.errors==shared.errors&&separate.matrixErrors==shared.matrixErrors&&separate.ratioErrors==shared.ratioErrors&&separate.bad==shared.bad&&separate.maximum==shared.maximum,"guard/ratio results changed");reused=outer.reused;Require(reused>0,"source graphs not reused across guards");
 }
 Require(activeValidationBatch==nullptr,"interpolation validation scope leaked");
 const auto standalone=CommonPrefixInterpolation(e,state,{6},ids,e.seeds.data());Require(standalone.errors==separate.errors,"standalone scope changed errors");
 Job victim;std::istringstream sourceRows(Read(compoundBlocks[0].folder+"/angular_sources.tsv"));std::string leaf,marker;Require(bool(sourceRows>>leaf>>marker),"no native source");victim.data=leaf+"/result/result.dat";
 const auto bytes=Read(victim.data);Atomic(victim.data,bytes+"corrupt\n");bool rejected=false;try{CommonPrefixInterpolation(e,state,{6},ids,e.seeds.data());}catch(const std::exception&){rejected=true;}Atomic(victim.data,bytes);Require(rejected,"changed descendant reused in new interpolation verification");Require(activeValidationBatch==nullptr,"failed verification leaked scope");
 std::cout<<"PASS shared guard graphs: all16 and normalized guard bounds bitwise equal, no new optical work, fresh-pass corruption rejected, reused="<<reused<<'\n';return 0;
}catch(const std::exception&e){std::cerr<<e.what()<<std::endl;return 1;}}
