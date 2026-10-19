#include "../src/FullAuto.cpp"
using namespace FullAuto;
struct GuardRun {InterpolationResult result;size_t native=0;double seconds=0;std::string budgets;std::map<std::string,Mueller> samples;};
void Setup(Engine& e,const std::string& root,const std::string& binary,int mode){
 e.root=root;e.signature="guard-prefetch-regression-v1";Directory(root);Directory(root+"/jobs");e.exe=binary;e.binaryHash=FileHash(binary);
 e.controls=mode==1;e.depthSharing=true;e.sharedMeans=true;e.gpus={""};e.seeds.assign(production,production+8);e.workerThreads=1;e.analyticMeanThreads=1;e.hostBudgetMb=512;
 e.phi=1;e.target=8;e.initial=mode==3?512:64;e.orientationBlock=mode==3?512:32;e.workerPolicy="persistent";e.resourcePolicy="off";e.workQueuePolicy="on";
 e.balancedCorrection=mode==3;
 e.base={"--method","po","--backend","cpu","--particle","1",".2",".1","--ri","1.3116",mode==1?".001":"0","--wavelength-um",".532","--symmetry","1","1","--cutoff-profile","off"};
 for(size_t i=0;i<17;++i){e.theta.push_back(170.+10.*i/16);if(mode==2)e.thetaPhi[i]=i<8?1:4;}
 WriteGrid(root+"/requested_theta.csv",e.theta);
}
size_t Native(const std::string& root){std::istringstream in(Read(root+"/fullauto_history.jsonl"));std::string line;size_t n=0;while(std::getline(in,line))if(line.find("\"event\":\"worker_complete\"")!=std::string::npos)++n;return n;}
GuardRun Compute(const std::string& root,const std::string& binary,int mode,bool merge,bool restore=false,bool partial=false){
 ArgPP args;RunConfig config;Engine e(args,config);Setup(e,root,binary,mode);AdaptiveState state;std::vector<size_t> ids;
 state.counts=Budgets(2,std::vector<int>(17,e.initial));
 for(size_t i=0;i<17;++i){ids.push_back(i);if(i%2==0)state.knots.push_back(i);if(i>=10)state.counts[0][i]*=2;}
 Require(FullAutoAuditPlan::Stencils(state.knots,state.counts,ids).size()>1,"requires distinct stencil budgets");
 if(restore){e.WarmCache();e.forbidNativeWork=true;}
 std::map<std::string,std::string> retained;
 if(partial){
  const auto groups=FullAutoAuditPlan::Stencils(state.knots,state.counts,ids);const auto& group=*groups.begin();
  const std::vector<size_t> subset(group.second.rows.begin(),group.second.rows.end());
  EvaluateBudget(e,{6,8},StencilBudgets(e,state.counts,group),subset,e.seeds.data());
  for(const auto& item:e.cache)retained[item.second.folder+"/complete.txt"]=Read(item.second.folder+"/complete.txt");
  e.cache.clear();e.samplingIndex.clear();e.folderIndex.clear();e.payloads.Clear();e.WarmCache();
 }
 const double start=Seconds();GuardRun r;
 {ValidationBatch scope;r.result=CommonPrefixInterpolationImpl(e,state,{6,8},ids,e.seeds.data(),merge);}
 for(const auto& item:retained)Require(Read(item.first)==item.second,"partial checkpoint overwritten");
 r.seconds=Seconds()-start;r.native=Native(root);r.budgets=Read(root+"/interpolation_production_budgets.csv");
 std::istringstream sources(Read(root+"/interpolation_production_sources.tsv"));int group,level,sign,n;unsigned seed;std::string folder;
 while(sources>>group>>level>>sign>>n>>seed>>folder){const auto job=e.Materialize(*e.folderIndex.at(folder));const auto data=e.Data({job});
  for(size_t row=0;row<job.angles.size();++row){
   double angle=job.angles[row];
#ifdef MBS_TEST_FP32
   angle=std::round(angle*1e8)/1e8; // DAT degree conversion, not an estimator change
#endif
   auto& sum=r.samples[Number(group)+":"+Number(level)+":"+Number(seed)+":"+Number(angle)];for(int j=0;j<16;++j)sum[j]+=sign*data[0][row][j];}}
 return r;
}
void Equal(const GuardRun& a,const GuardRun& b){
 Require(a.budgets==b.budgets,"stencil budgets changed");
 Require(a.result.bad==b.result.bad,"guard acceptance changed");
#ifdef MBS_TEST_FP32
 // Float reductions need not be bitwise identical when angular packets change.
 // Require four times less than the 2% run's eps/1024 execution budget; scale
 // matrix entries by the production mean M11. The estimator is a ratio of
 // means, not a ratio from an individual control-corrected scramble.
 const double tolerance=.02/4096;
#else
 const double tolerance=2e-11;
#endif
 auto close=[&](double x,double y){Require(std::fabs(x-y)<=tolerance*std::max(1.,std::max(std::fabs(x),std::fabs(y))),"guard estimate changed");};
 for(size_t r=0;r<a.result.errors.size();++r){close(a.result.errors[r],b.result.errors[r]);for(int j=0;j<16;++j){close(a.result.matrixErrors[r][j],b.result.matrixErrors[r][j]);close(a.result.ratioErrors[r][j],b.result.ratioErrors[r][j]);}}
 close(a.result.maximum,b.result.maximum);
 Require(a.samples.size()==b.samples.size(),"sample coordinates changed");
#ifdef MBS_TEST_FP32
 auto coordinate=[](const std::string& key){const size_t first=key.find(':'),second=key.find(':',first+1),third=key.find(':',second+1);return key.substr(0,first)+":"+key.substr(third+1);};
 auto levelCoordinate=[](const std::string& key){const size_t first=key.find(':'),second=key.find(':',first+1),third=key.find(':',second+1);return key.substr(0,second)+":"+key.substr(third+1);};
 std::map<std::string,Mueller> meanA,meanB;std::map<std::string,size_t> counts;
 for(const auto& item:a.samples){const auto found=b.samples.find(item.first);Require(found!=b.samples.end(),"sample seed/depth/angle changed");
  auto& x=meanA[coordinate(item.first)];auto& y=meanB[coordinate(item.first)];++counts[levelCoordinate(item.first)];
  for(int j=0;j<16;++j){x[j]+=item.second[j]/8;y[j]+=found->second[j]/8;}}
 for(const auto& item:counts)Require(item.second==8,"not eight independent scramble samples");
 for(const auto& item:meanA){const auto& x=item.second;const auto& y=meanB.at(item.first);Require(x[0]>0&&y[0]>0,"nonpositive mean M11");
  for(int j=0;j<16;++j)Require(std::fabs(x[j]-y[j])<=tolerance*x[0],"float mean matrix budget exceeded");
  for(int j=1;j<16;++j)Require(std::fabs(x[j]/x[0]-y[j]/y[0])<=tolerance,"float ratio-of-means budget exceeded");}
#endif
 for(const auto& item:a.samples){const auto found=b.samples.find(item.first);Require(found!=b.samples.end(),"sample seed/depth/angle changed");
#ifdef MBS_TEST_FP32
  const double scale=meanA.at(coordinate(item.first))[0];
  for(int j=0;j<16;++j)Require(std::fabs(item.second[j]-found->second[j])<=tolerance*scale,"all16 float budget exceeded");
#else
  for(int j=0;j<16;++j)close(item.second[j],found->second[j]);
#endif
 }
}
int main(int argc,char**argv){try{Require(argc==3,"OUTPUT BINARY");
#ifdef MBS_TEST_FP32
 // A control-corrected scramble can have zero M11 while its positive mean
 // and ratio-of-means remain well defined. Still reject an excessive raw
 // perturbation, even when its mean could cancel across independent seeds.
 {GuardRun a,b;for(int s=0;s<8;++s){Mueller m{};m[0]=s?8./7:0;m[1]=.1;
   a.samples["0:0:"+Number(s)+":170"]=m;}b=a;Equal(a,b);
  b.samples.begin()->second[1]+=.02/512;bool rejected=false;try{Equal(a,b);}catch(const std::exception&){rejected=true;}
  Require(rejected,"raw scramble perturbation escaped float gate");}
#endif
 for(int mode=0;mode<4;++mode){const std::string root=std::string(argv[1])+"/mode"+Number(mode);
  const auto old=Compute(root+"/separate",argv[2],mode,false),merged=Compute(root+"/merged",argv[2],mode,true);Equal(old,merged);
  Require(merged.native<old.native,"prefetch did not reduce native trace jobs");
  const auto resumed=Compute(root+"/merged",argv[2],mode,true,true);Equal(merged,resumed);Require(resumed.native==merged.native,"resume launched optical work");
  if(mode==3){const auto partial=Compute(root+"/partial",argv[2],mode,true,false,true);Equal(old,partial);Require(partial.native<=old.native,"partial warmup repeated optical work");std::cout<<"PASS partial checkpoint: saved blocks retained, all16 and ratios equal\n";}
  std::cout<<"PASS mode="<<mode<<" separate_jobs="<<old.native<<" merged_jobs="<<merged.native<<" separate_seconds="<<old.seconds<<" merged_seconds="<<merged.seconds<<" all16 ratios budgets resume\n";
 }
 return 0;}catch(const std::exception&e){std::cerr<<e.what()<<std::endl;return 1;}}
