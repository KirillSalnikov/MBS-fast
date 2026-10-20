// Regression for production prefixes assembled from legacy paired lower views.
// --sources reads an existing production source graph without modifying it.
#include "../src/FullAuto.cpp"
using namespace FullAuto;
namespace {
std::string Native(const std::string& root,const std::string& id,int n,int offset,unsigned seed,
                   int depth,double area,double ext,bool paired=false,bool revision=false,bool explicitArea=false){
    const auto f=root+"/"+id;Directory(f+"/result");
    Atomic(f+"/sampling.index","fixture "+id+" "+Number(n)+" "+Number(offset)+" "+Number(seed)+" "+Number(depth)+"\n");
    Atomic(f+"/result/result.dat","verified matrix fixture\n");
    auto scalar=[&](const std::string& path,double a,double c){Atomic(path,
        "Cext\tQext"+std::string(explicitArea?"\tprojected_area":"")+"\n"+Number(c)+"\t"+Number(c/a)+(explicitArea?"\t"+Number(a):"")+"\n");};
    scalar(f+"/result/result_integrals.tsv",area,ext);
    if(paired){scalar(f+"/result/result_depth4_integrals.tsv",999,888);
        if(revision)Atomic(f+"/result/result_depth4_integrals_revision.txt","paired_extinction_v1\n");}
    Atomic(f+"/complete.txt",id+" 1 "+FileHash(f+"/result/result.dat")+" none\n");return f;
}
std::string View(const std::string& root,const std::string& id,int n,int offset,unsigned seed,int depth,const std::vector<std::string>& inputs){
    const auto f=root+"/"+id;Directory(f+"/result");std::ostringstream sources;
    for(const auto& p:inputs)sources<<p<<'\t'<<FileHash(p+"/complete.txt")<<'\n';
    Atomic(f+"/sampling.index","fixture "+id+" "+Number(n)+" "+Number(offset)+" "+Number(seed)+" "+Number(depth)+"\n");
    Atomic(f+"/result/result.dat","assembled fixture\n");Atomic(f+"/angular_sources.tsv",sources.str());
    Atomic(f+"/complete.txt",id+" 1 "+FileHash(f+"/result/result.dat")+" none "+Hash(sources.str())+"\n");return f;
}
void Check(bool condition,const std::string& message){Require(condition,"legacy integral regression: "+message);}
}
int main(int argc,char** argv){try{
    Require(argc>=2,"output directory or --sources SOURCE_FILE required");const std::vector<unsigned> seeds(production,production+8);
    if(std::string(argv[1])=="--sources"){
        Require(argc==3,"--sources requires a file");const auto result=IntegralProductionSources(Read(argv[2]),seeds);
        std::cout<<"fallback\t"<<result.directTargetFallback<<"\nseed\tprojected_area\tCext\n"<<std::setprecision(17);
        for(size_t i=0;i<8;++i)std::cout<<seeds[i]<<'\t'<<result.samples[i].area<<'\t'<<result.samples[i].ext<<'\n';
        std::cout<<"source_files\t"<<result.files.size()<<'\n';return 0;
    }
    const std::string root=argv[1];Directory(root);
    for(bool explicitArea:{false,true})for(bool allLegacy:{false,true}){
        const auto dir=root+"/"+Number(explicitArea)+Number(allLegacy);Directory(dir);std::ostringstream sources;
        std::string pairedNative,pairedView;
        for(size_t i=0;i<8;++i){const unsigned seed=seeds[i];const auto id=Number(seed);
            const bool legacy=allLegacy||i==2; // Fault can first occur at seed 37.
            const auto first=Native(dir,"native"+id,2,0,seed,legacy?6:4,10+i,20+i,legacy,false,explicitArea);
            const auto lower=legacy?View(dir,"lower"+id,2,0,seed,4,{first}):first;
            const auto tail=Native(dir,"tail"+id,6,2,seed,4,30+i,40+i,false,false,explicitArea);
            const auto prefix=View(dir,"prefix"+id,8,0,seed,4,{lower,lower,tail});
            const auto duplicate=View(dir,"duplicate"+id,8,0,seed,4,{prefix});
            const auto upper=Native(dir,"upper"+id,2,0,seed,8,50+i,100+i,false,false,explicitArea);
            sources<<"0 1 8 "<<seed<<' '<<prefix<<"\n0 1 8 "<<seed<<' '<<duplicate
                   <<"\n1 1 2 "<<seed<<' '<<upper<<"\n1 -1 2 "<<seed<<' '<<lower<<'\n';
            if(legacy){pairedNative=first;pairedView=lower;}
        }
        Atomic(dir+"/sources.tsv",sources.str());const auto result=IntegralProductionSources(sources.str(),seeds);
        const auto parallel=IntegralProductionSources(sources.str(),seeds,4);
        Check(parallel.files==result.files&&parallel.selectedSources==result.selectedSources&&parallel.directTargetFallback==result.directTargetFallback,"parallel provenance/fallback differs");
        for(size_t i=0;i<8;++i)Check(parallel.samples[i].area==result.samples[i].area&&parallel.samples[i].ext==result.samples[i].ext,"parallel scalar differs");
        Check(result.directTargetFallback,"fallback not recorded");
        for(size_t i=0;i<8;++i){Check(std::fabs(result.samples[i].area-(25+i))<1e-12,"wrong window-weighted base area");Check(result.samples[i].ext==100+i,"legacy extinction used or duplicated");}
        for(const auto& entry:result.files)Check(entry.first.find("_depth4_integrals.tsv")==std::string::npos,"fallback reads untrusted lower scalar");
        IntegralReader reader;bool rejected=false;try{reader.Resolve(pairedView,4);}catch(const LegacyPairedIntegral&){rejected=true;}
        Check(rejected,"legacy extinction silently accepted");
        // Quantity-sensitive memoization: reading area must not authorize extinction.
        IntegralReader quantities;quantities.Resolve(pairedView,4,true);rejected=false;
        try{quantities.Resolve(pairedView,4);}catch(const LegacyPairedIntegral&){rejected=true;}
        Check(rejected,"area cache masks invalid extinction");
        Atomic(pairedNative+"/result/result.dat","corrupt matrix\n");rejected=false;
        try{IntegralProductionSources(sources.str(),seeds,4);}catch(const std::exception&){rejected=true;}
        Check(rejected,"fallback bypasses source integrity");
    }
    std::cout<<"PASS: mixed/all legacy base views; eight replicas; explicit/derived area; nested prefixes; duplicate angles; weighted windows; quantity cache; corruption rejection\n";
    return 0;
}catch(const std::exception& error){std::cerr<<error.what()<<'\n';return 1;}}
