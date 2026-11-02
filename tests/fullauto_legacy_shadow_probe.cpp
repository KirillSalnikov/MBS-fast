// Identity migration guards; optical equivalence is certified separately.
#include "../src/FullAuto.cpp"
using namespace FullAuto;
int main(int argc,const char* argv[]){try{
    if(argc==2&&std::string(argv[1])=="--version"){
        std::cout<<"MBS-fast shadow-probe\nbackend-build: CPU double\n";return 0;
    }
    Require(argc==2,"expected private test directory");
    const std::string root=argv[1];Directory(root);
    const std::string origin=root+"/old-version",evidence=root+"/evidence.json";
    Atomic(origin,"#!/bin/sh\nprintf 'MBS-fast 8ae9e85\\nbackend-build: CPU double\\n'\n");
    Require(chmod(origin.c_str(),0755)==0,"chmod old fixture");
    Atomic(evidence,"{\"fixture\":\"identity guards only\"}\n");
    char self[4096];const ssize_t n=readlink("/proc/self/exe",self,sizeof(self)-1);Require(n>0,"self");self[n]=0;
    const std::string receipt=std::string(self)+".native-compatibility.tsv";
    const std::string body="FULLAUTO_NATIVE_COMPATIBILITY_V1\n"+FileHash(origin)+" "+FileHash(self)+" "+evidence+" "+FileHash(evidence)+"\n";
    Atomic(receipt,Hash(body)+"\n"+body);
    std::vector<std::string> command={self,"--method","po","--backend","cpu","--particle","1","2","1",
        "--ri","1.3116","0","--wavelength-um",".532","--symmetry","2","6","--threads","1",
        "--fullauto",".04","--output",root,"--close","--fullauto-resume-compatible",origin};
    std::ostringstream old;old<<FileHash(origin)<<'\n';
    for(size_t i=1;i<command.size();++i){if(command[i]=="--fullauto-resume-compatible"){++i;continue;}old<<command[i]<<'\n';}
    old<<FormatRuntimeEnvironmentReport(QueryRuntimeEnvironmentSnapshot());
    char host[256]={};Require(gethostname(host,sizeof(host)-1)==0,"hostname");
    old<<"host="<<host<<"\nthreads=1\nanalytic_mean_threads=1\nvisible="<<(getenv("CUDA_VISIBLE_DEVICES")?getenv("CUDA_VISIBLE_DEVICES"):"")<<"\ngpu=\n";
    const std::string legacy=Hash(old.str());
    Atomic(root+"/configuration.txt",legacy+"\n");Atomic(root+"/controller_command.json","{}\n");Atomic(root+"/fullauto_parameters.json","{}\n");
    auto check=[&](const std::vector<std::string>& values,bool expectMigration){
        std::vector<const char*> a;for(const auto& s:values)a.push_back(s.c_str());
        ArgPP args;RegisterCliOptions(args);args.Parse(a.size(),a.data());
        RunConfig config=RunConfig::FromCommandLine(args,false,false);Engine e(args,config);e.Init(a.size(),a.data());
        Require(e.legacyProjectedShadowIdentity==expectMigration,"migration decision differs");
        Require((e.signature==legacy)==expectMigration,"legacy identity acceptance differs");
        const bool explicitRecipe=std::find(e.base.begin(),e.base.end(),"--po-shadow-model")!=e.base.end();
        Require(explicitRecipe!=expectMigration,"worker recipe identity changed");
    };
    check(command,true);
    auto changed=command;changed.insert(changed.end(),{"--po-shadow-model","projected"});check(changed,false);
    changed=command;changed.insert(changed.end(),{"--po-shadow-model","surface"});check(changed,false);
    changed=command;changed[std::find(changed.begin(),changed.end(),"--fullauto")-changed.begin()+1]=".05";check(changed,false);
    changed=command;changed[std::find(changed.begin(),changed.end(),"--ri")-changed.begin()+1]="1.4";check(changed,false);
    Require(unlink(receipt.c_str())==0,"remove test certificate");
    bool rejected=false;try{check(command,true);}catch(const std::exception& e){rejected=std::string(e.what()).find("numerical certification")!=std::string::npos;}
    Require(rejected,"uncertified legacy resume accepted");
    std::cout<<"legacy shadow: preserved identity/recipe; explicit models, changed accuracy/material and missing certificate rejected\n";
    return 0;
}catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}}
