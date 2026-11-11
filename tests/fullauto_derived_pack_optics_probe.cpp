#include <string>
#include <map>
#include <set>
#include <sys/stat.h>
static void PackEvent(const char*,const std::string&);
#define MBS_PACK_EVENT(event,folder) PackEvent(event,folder)
#include "../src/FullAuto.cpp"
#undef MBS_PACK_EVENT
using namespace FullAuto;
static bool preparedExpected=false;
static std::set<char> packedKinds;
static bool Raw(const std::string& p){struct stat st{};return lstat(p.c_str(),&st)==0;}
static void PackEvent(const char* event,const std::string& folder){
 if(std::string(event)=="after_collect")packedKinds.insert(folder[folder.find_last_of('/')+1]);
 if(preparedExpected)Require(!Raw(folder+"/complete.txt"),"prepared path exposed loose completion");
}
void Setup(Engine& e,const std::string& root,const std::string& binary,bool controls,bool pack,bool compact){
 e.root=root;e.signature="derived-optics-v1";Directory(root+"/jobs");e.exe=binary;e.binaryHash=FileHash(binary);e.compactIndex=compact;e.packJobs=pack;
 e.controls=controls;e.sharedMeans=controls;e.coherent=false;e.gpus={""};e.seeds={11};e.workerThreads=1;e.analyticMeanThreads=1;e.hostBudgetMb=512;
 e.phi=4;e.target=10;e.theta={150,153.5,160};e.initial=512;e.orientationBlock=512;e.workerPolicy="persistent";e.resourcePolicy="off";e.workQueuePolicy="on";
 e.balancedCorrection=true;e.selectedCutoff="off";
 e.base={"--method","po","--backend","cpu","--particle","1",".2",".1","--ri","1.3116","0","--wavelength-um",".532","--symmetry","1","1","--cutoff-profile","off","--incoherent"};
 WriteGrid(root+"/requested_theta.csv",e.theta);
}
std::map<std::string,Job> Run(Engine& e){
 const unsigned* seeds=e.seeds.data();const auto& x=e.theta;std::map<std::string,Job> out;
 const auto one=e.DirectStage(512,seeds,1,10,x,true,false,0,8);
 const auto two=e.DirectStage(512,seeds,1,10,x,true,false,512,8);
 out["balanced_upper"]=one[0];out["paired_lower"]=e.LowerView(one[0],8);
 out["paired_parts"]=e.AssemblePairParts({one,two},1024,x,"proof upper")[0];
 out["prefix"]=e.PrefixStage(1024,seeds,1,10,x)[0];
 out["angular"]=e.Stage(512,seeds,1,10,{150,160})[0];
 for(const auto& item:e.cache){const auto& j=item.second;if(j.folder[j.folder.find_last_of('/')+1]=='v'&&j.offset==0){
   out["embedded_prefix"]=e.Materialize(j);break;}}
 Require(out.count("embedded_prefix"),"embedded prefix callsite was not exercised");
 const std::map<std::string,char> kinds={{"balanced_upper",'j'},{"paired_lower",'d'},{"paired_parts",'q'},{"prefix",'p'},{"angular",'a'},{"embedded_prefix",'v'}};
 for(const auto& it:kinds)Require(out.at(it.first).folder[out.at(it.first).folder.find_last_of('/')+1]==it.second,"wrong derived callsite "+it.first);
 for(const auto& item:out){std::set<std::string> visited;ValidateStoredJob(item.second.folder,visited);Require(visited.size()>1,"derived graph not validated");}
 return out;
}
static double maximumDifference=0;
static bool EqualNumber(double a,double b){const double delta=std::fabs(a-b)/std::max(1e-30,std::max(std::fabs(a),std::fabs(b)));maximumDifference=std::max(maximumDifference,delta);return delta<2e-12;}
static void Compare(const Job& a,const Job& b,const std::string& label){
 Require(a.matrix.size()==b.matrix.size()&&a.controls.size()==b.controls.size(),"row count differs "+label);
 double maximum=0;
 for(size_t r=0;r<a.matrix.size();++r)for(size_t k=0;k<16;++k){const double delta=std::fabs(a.matrix[r][k]-b.matrix[r][k])/std::max(1e-30,std::fabs(a.matrix[r][0]));maximum=std::max(maximum,delta);}
 for(size_t r=0;r<a.controls.size();++r)for(size_t k=0;k<5;++k)Require(EqualNumber(a.controls[r][k],b.controls[r][k]),"control differs "+label);
 std::cout<<"matrix_comparison "<<label<<" max_abs_delta_over_M11="<<std::setprecision(17)<<maximum<<std::endl;
 Require(maximum<2e-12,"all16 changed beyond FP64 roundoff "+label);maximumDifference=std::max(maximumDifference,maximum);
}
int main(int argc,char**argv){try{
 Require(argc==3,"fresh scratch root and frozen CPU worker required");const std::string root=argv[1];
 for(bool controls:{false,true}){
  std::map<std::string,Job> reference;std::map<std::string,std::pair<double,double>> scalars;
  for(bool compact:{true,false})for(bool packed:{false,true}){
   ArgPP args;RegisterCliOptions(args);RunConfig config;Engine e(args,config);
   const auto dir=root+"/control"+Number(controls)+"_compact"+Number(compact)+"_pack"+Number(packed);
   Setup(e,dir,argv[2],controls,packed,compact);preparedExpected=packed&&compact;packedKinds.clear();
   const auto results=Run(e);const bool first=compact&&!packed;IntegralReader integrals;
   for(const auto& item:results){const auto& j=item.second;const auto v=integrals.Resolve(j.folder,j.depth);Require(std::isfinite(v.area)&&std::isfinite(v.ext),"nonfinite integral");
    if(first){reference[item.first]=j;scalars[item.first]={v.area,v.ext};}else{
     const auto& r=reference.at(item.first);Compare(r,j,item.first);
     Require(EqualNumber(scalars.at(item.first).first,v.area)&&EqualNumber(scalars.at(item.first).second,v.ext),"signed/paired/prefix integral changed "+item.first);}
    if(packed)Require(!Raw(j.folder+"/complete.txt"),"packed derived left loose marker");else Require(Raw(j.folder+"/complete.txt"),"raw derived marker missing");
   }
   if(packed)for(char k:std::string("jvdqpa"))Require(packedKinds.count(k),"uncovered publication kind "+std::string(1,k));
   for(auto& service:e.services)service.reset();Engine resumed(args,config);Setup(resumed,dir,argv[2],controls,packed,compact);resumed.WarmCache();resumed.forbidNativeWork=true;
   const size_t before=resumed.cache.size();const auto restored=Run(resumed);Require(before==resumed.cache.size(),"resume added job");
   for(const auto& item:restored)Require(item.second.matrix==results.at(item.first).matrix&&item.second.controls==results.at(item.first).controls,"resume changed estimator");
   std::cout<<"PASS controls="<<controls<<" compact="<<compact<<" packed="<<packed<<" all6 callsites all16/controls/integrals within FP64 rounding; graph valid; no-work resume; jobs="<<before<<std::endl;
  }
 }
 std::cout<<"PASS eight derived optics configurations; maximum scaled difference="<<std::setprecision(17)<<maximumDifference<<"\n";return 0;
}catch(const std::exception& ex){std::cerr<<ex.what()<<'\n';return 1;}}
