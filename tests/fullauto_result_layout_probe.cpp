#include "../src/FullAuto.cpp"
using namespace FullAuto;
int main(int argc,char** argv){try{
 Require(argc==2,"fixture path required");
 for(bool go:{false,true}){
  std::string parent=std::string(argv[1])+(go?"/go":"/po");Directory(parent);
  ArgPP args;RegisterCliOptions(args);const char* flags[]={"test","--fullauto-layout","compact"};args.Parse(3,flags);
  RunConfig cfg;Engine e(args,cfg);e.root=parent;e.go=go;e.ConfigureResultLayout();
  Require(e.root==parent+"/work","work routing wrong");
  Atomic(e.root+"/fullauto_status.json","{\"status\":\"evaluating\",\"round\":14}");
  Atomic(e.root+"/fullauto_parameters.json","{\"requested_options\":{},\"selected_state\":{\"target_depth\":8}}");
  e.ExportCompactResults(true);Require(!Exists(parent+"/latest"),"empty matrix published");
  const auto f=e.root+"/round14";Directory(f);const std::vector<double> x={10,20};Matrix m(2);m[0][0]=12;m[1][0]=34;
  WriteMatrix(f+(go?"/mueller_evaluated.dat":"/mueller_fullauto.dat"),x,m);
  Atomic(f+"/mueller_fullauto_integrals.tsv","quantity\tvalue\nCsca\t123\n");
  if(go){std::vector<double> xx={0};Matrix p(1);p[0][0]=56;WriteMatrix(f+"/mueller_forward.dat",xx,p);xx[0]=180;p[0][0]=78;WriteMatrix(f+"/mueller_backward.dat",xx,p);}
  Atomic(e.root+"/intermediate_current.json","{\"snapshot_folder\":"+Quoted(f)+",\"phase\":\"depth_selection\"}");e.ExportCompactResults(true);Require(!Exists(parent+"/latest"),"pilot published as production");
  Atomic(e.root+"/intermediate_current.json","{\"snapshot_folder\":"+Quoted(f)+",\"round\":14,\"statistical95\":0.03}");e.ExportCompactResults(true);
  const auto old=e.lastCompactGeneration;Require(Read(parent+"/latest/mueller_bins.dat")==Read(f+(go?"/mueller_evaluated.dat":"/mueller_fullauto.dat")),"matrix changed");
  const auto pole=Read(parent+"/latest/pole_contributions.dat");Require(go?pole.find("180 78")!=std::string::npos:pole.find("not_applicable")!=std::string::npos,"pole semantics wrong");
  Require(Read(parent+"/latest/accuracy.json").find("\"converged\":false")!=std::string::npos,"false completion");
  Atomic(e.root+"/fullauto_status.json","{\"status\":\"evaluating\",\"round\":15}");e.compactPublishedAt=-1;e.ExportCompactResults(false);Require(e.lastCompactGeneration==old,"observation republished matrix");
  Atomic(e.root+"/fullauto_status.json","{\"status\":\"converged\",\"round\":14}");Atomic(e.root+"/final_outputs.fnv","fixture\n");
  Atomic(e.root+(go?"/mueller_evaluated.dat":"/mueller_fullauto.dat"),Read(f+(go?"/mueller_evaluated.dat":"/mueller_fullauto.dat")));
  if(go)for(const char* n:{"mueller_forward.dat","mueller_backward.dat"})Atomic(e.root+"/"+n,Read(f+"/"+n));
  e.ExportCompactResults(true);Require(e.lastCompactGeneration!=old,"completion not published");Require(Read(parent+"/latest/accuracy.json").find("\"converged\":true")!=std::string::npos,"missing accepted completion");
  Engine resume(args,cfg);resume.root=parent;resume.ConfigureResultLayout();Require(resume.root==e.root,"resume path differs");
  Atomic(parent+"/latest/mueller_bins.dat","corrupt");e.lastCompactGeneration.clear();bool rejected=false;try{e.ExportCompactResults(true);}catch(const std::exception&){rejected=true;}Require(rejected,"corrupted immutable result accepted");
 }
 std::cout<<"PASS compact PO/GO: layout, no pilot publication, exact matrix, pole semantics, immutable generations, status throttling, resume path, corruption rejection\n";
}catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}}
