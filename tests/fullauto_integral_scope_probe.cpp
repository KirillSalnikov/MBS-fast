// The fixtures exercise numerical/provenance and failure semantics only.
// Production performance evidence is a separate replay of real saved graphs.
#define main fullauto_legacy_scope_fixture_main
#include "fullauto_legacy_integrals_probe.cpp"
#undef main
namespace {
void SameIntegral(const IntegralProduction& a,const IntegralProduction& b){
    Require(a.files==b.files&&a.selectedSources==b.selectedSources&&a.directTargetFallback==b.directTargetFallback,"scoped integral provenance changed");
    for(size_t i=0;i<8;++i)Require(a.samples[i].area==b.samples[i].area&&a.samples[i].ext==b.samples[i].ext,"scoped integral scalar changed");
}
std::string ProtectedNative(const std::string& root,const std::string& id,int n,int offset,unsigned seed,int depth,double area,double ext){
    const auto folder=Native(root,id,n,offset,seed,depth,area,ext,false,false,true);
    const auto manifest=folder+"/result/result_integrals.tsv "+FileHash(folder+"/result/result_integrals.tsv")+"\n";
    Atomic(folder+"/integral_outputs.fnv",manifest);
    Atomic(folder+"/complete.txt",id+" 1 "+FileHash(folder+"/result/result.dat")+" none none "+Hash(manifest)+"\n");return folder;
}
std::string MemoDigest(const ValidationBatch& batch){std::map<std::string,std::string> values;batch.compounds.LeastRecentFirst([&](const std::string& key,std::shared_ptr<const std::string> value){values[key]=*value;});std::ostringstream out;for(const auto& value:values)out<<value.first<<' '<<value.second<<'\n';return Hash(out.str());}
}
int main(int argc,char** argv){try{
    Require(argc==2,"OUTPUT_DIRECTORY");const std::string root=argv[1];Directory(root);const std::vector<unsigned>seeds(production,production+8);
    std::vector<std::string> roots;std::ostringstream sources;std::string victim;
    for(size_t i=0;i<8;++i){const auto id=Number(seeds[i]);
        const auto first=ProtectedNative(root,"first"+id,2,0,seeds[i],4,10+i,20+i);
        const auto tail=ProtectedNative(root,"tail"+id,6,2,seeds[i],4,30+i,40+i);
        const auto prefix=View(root,"prefix"+id,8,0,seeds[i],4,{first,first,tail});
        const auto duplicate=View(root,"duplicate"+id,8,0,seeds[i],4,{prefix});
        const auto upper=ProtectedNative(root,"upper"+id,2,0,seeds[i],8,50+i,100+i);
        sources<<"0 1 8 "<<seeds[i]<<' '<<prefix<<"\n0 1 8 "<<seeds[i]<<' '<<duplicate<<"\n1 1 2 "<<seeds[i]<<' '<<upper<<"\n1 -1 2 "<<seeds[i]<<' '<<first<<'\n';
        for(const auto& path:{prefix,duplicate,upper,first})roots.push_back(path);if(i==2)victim=first;
    }
    const auto source=sources.str();const auto serial=IntegralProductionSources(source,seeds,1);
    for(size_t i=0;i<8;++i){Require(serial.samples[i].area==25+i&&serial.samples[i].ext==115+i,"window weighting or signed depth correction changed");}
    Require(!serial.directTargetFallback,"new scalar fixture unexpectedly needs legacy fallback");
    Require(activeValidationBatch==nullptr,"serial integral scope leaked");
    for(int threads:{1,4}){
        ValidationBatch parent("preceding production round");ValidateStoredRoots(roots,4);
        const auto before=parent.reused;const auto result=IntegralProductionSources(source,seeds,threads);SameIntegral(serial,result);
        Require(activeValidationBatch==&parent,"integral replaced parent scope");
        Require(result.validationReusedGraphs>=16&&parent.reused-before==result.validationReusedGraphs,"not all eight replicas inherited parent graph memo");
        Require(result.validationCheckedNodes>0&&result.validationTemporaryMemoBytes<=integralValidationMemoLimit,"integral validation counters/budget invalid");
        const auto unchanged=MemoDigest(parent);const auto checked=parent.checked,reused=parent.reused,bytes=parent.compounds.Bytes();
        const auto file=victim+"/result/result.dat",original=Read(file);std::string bad=original;bad[0]=bad[0]=='x'?'y':'x';struct stat stamp;Require(stat(file.c_str(),&stamp)==0,"fault fixture stat");
        {std::fstream edit(file.c_str(),std::ios::binary|std::ios::in|std::ios::out);edit.write(bad.data(),std::streamsize(bad.size()));edit.flush();Require(bool(edit),"fault fixture write");}
        timespec times[2]={stamp.st_atim,stamp.st_mtim};Require(utimensat(AT_FDCWD,file.c_str(),times,0)==0,"fault fixture restore mtime");
        bool rejected=false;try{IntegralProductionSources(source,seeds,threads);}catch(const std::exception&){rejected=true;}
        Require(rejected,"same-size same-inode restored-mtime source corruption accepted");
        Require(activeValidationBatch==&parent&&parent.checked==checked&&parent.reused==reused&&parent.compounds.Bytes()==bytes&&MemoDigest(parent)==unchanged,"failed integral published partial parent memo/counters");
        Atomic(file,original);SameIntegral(serial,IntegralProductionSources(source,seeds,threads));
    }
    Require(activeValidationBatch==nullptr,"parallel integral scope leaked");
    // Every new round starts a fresh graph check, even when file hashes are warm.
    const auto file=victim+"/result/result_integrals.tsv",saved=Read(file);Atomic(file,saved+" ");
    {ValidationBatch parent("new round");bool rejected=false;try{IntegralProductionSources(source,seeds,4);}catch(const std::exception&){rejected=true;}Require(rejected&&parent.compounds.Bytes()==0&&parent.checked==0,"failed fresh round published graph authority");}
    Atomic(file,saved);Require(activeValidationBatch==nullptr,"fresh-round failure leaked scope");
    {ValidationBatch parent("tiny parent memo");parent.compounds.SetLimit(1);SameIntegral(serial,IntegralProductionSources(source,seeds,4));Require(parent.compounds.Bytes()==0,"tiny parent cache exceeded limit");}
    // Existing mixed/all legacy tests cover eight replicas, both projected-area
    // formats, quantity-sensitive caching and integrity failure during fallback.
    std::string legacy=root+"/legacy";char* legacyArgs[]={argv[0],&legacy[0]};Require(fullauto_legacy_scope_fixture_main(2,legacyArgs)==0,"legacy regression failed");
    std::cout<<"PASS: all8 integral scalars/provenance exact; window/depth arithmetic; shared round memo; <=64MiB private memos; serial+parallel all-or-none; same-size restored-mtime and integral corruption; tiny memo; legacy fallback\n";return 0;
}catch(const std::exception&e){std::cerr<<e.what()<<'\n';return 1;}}
