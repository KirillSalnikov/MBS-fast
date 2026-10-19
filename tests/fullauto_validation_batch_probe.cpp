#include "../src/FullAuto.cpp"
using namespace FullAuto;
int main(int argc,char**argv){try{
    Require(argc==3,"OUTPUT BINARY");ArgPP args;RunConfig config;Engine e(args,config);
    e.root=argv[1];e.signature="graph-scope-v1";Directory(e.root);Directory(e.root+"/jobs");e.exe=argv[2];e.binaryHash=FileHash(e.exe);
    e.controls=true;e.sharedMeans=true;e.gpus={""};e.seeds.assign(production,production+8);e.workerThreads=1;e.analyticMeanThreads=1;e.hostBudgetMb=512;
    e.phi=1;e.target=6;e.theta={0,5,170,175,180};e.initial=128;e.orientationBlock=32;e.workerPolicy="persistent";e.resourcePolicy="off";e.workQueuePolicy="on";
    e.base={"--method","po","--backend","cpu","--particle","1",".2",".1","--ri","1.3116","0","--wavelength-um",".532","--symmetry","1","1","--cutoff-profile","off"};
    WriteGrid(e.root+"/requested_theta.csv",e.theta);
    const auto original=e.PrefixStage(128,e.seeds.data(),8,6,e.theta);e.payloads.SetLimit(1);e.forbidNativeWork=true;
    std::vector<size_t> ids={0,1,2,3,4};Budgets counts(1,std::vector<int>(5,128));
    const auto evaluated=EvaluateBudget(e,{6},counts,ids,e.seeds.data());Require(evaluated.levels[0]==e.Data(original),"evaluation changed all16");
    const auto selected=original[0];const auto originalData=Read(selected.data);
    size_t checked=0,hits=0;
    {
        ValidationBatch batch;
        for(int rep=0;rep<4;++rep)for(const auto& j:original){const auto loaded=e.Materialize(e.cache.at(j.identity));Require(loaded.matrix==j.matrix&&loaded.controls==j.controls&&loaded.controlRefinement==j.controlRefinement,"scoped reload changed estimator");}
        Require(batch.reused>0,"no compound graph reused");checked=batch.checked;hits=batch.reused;
        Atomic(selected.data,originalData+"corrupt\n");bool rejected=false;try{e.Materialize(e.cache.at(selected.identity));}catch(const std::exception&){rejected=true;}Atomic(selected.data,originalData);Require(rejected,"direct payload mutation bypassed validation batch");
        std::istringstream rows(Read(selected.folder+"/angular_sources.tsv"));std::string child,digest;Require(bool(rows>>child>>digest),"no child");const auto marker=Read(child+"/complete.txt");Atomic(child+"/complete.txt",marker+"corrupt\n");
        rejected=false;try{e.Materialize(e.cache.at(selected.identity));}catch(const std::exception&){rejected=true;}Atomic(child+"/complete.txt",marker);Require(rejected,"direct source marker mutation bypassed batch");
    }
    Require(activeValidationBatch==nullptr,"scope leaked");
    // A new evaluation verifies the complete graph afresh, including leaves
    // whose parent matrices were already loaded in the previous evaluation.
    Job victim;bool found=false;for(const auto& item:e.cache)if(!item.second.assembled){victim=item.second;found=true;break;}Require(found,"no native leaf");
    const auto saved=Read(victim.data);Atomic(victim.data,saved+"corrupt\n");
    {
        ValidationBatch batch;for(int attempt=0;attempt<2;++attempt){bool rejected=false;
            try{for(const auto& j:original)e.Materialize(e.cache.at(j.identity));}catch(const std::exception&){rejected=true;}
            Require(rejected,"corrupt descendant or failed validation cached as successful");}
    }
    Atomic(victim.data,saved);Require(activeValidationBatch==nullptr,"exception scope leaked");
    {ValidationBatch batch;batch.compounds.SetLimit(1);for(const auto& j:original)Require(e.Materialize(e.cache.at(j.identity)).matrix==j.matrix,"eviction changed result");Require(batch.compounds.Bytes()==0,"graph cache not bounded");}
    std::cout<<"PASS graph validation: all16 and controls exact; direct payload/marker corruption rejected; fresh-evaluation descendant rejection; failed validation not reused; bounded cache; checked="<<checked<<" reused="<<hits<<'\n';
    return 0;
}catch(const std::exception& x){std::cerr<<x.what()<<'\n';return 1;}}
