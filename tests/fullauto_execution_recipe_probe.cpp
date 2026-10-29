// Replay controller planning from a COPY of a real metadata checkpoint.
// Never run this probe in a live calculation directory: journal repair/status
// writes follow the usual controller rules. No optical payloads are required.
#include "../src/FullAuto.cpp"
using namespace FullAuto;

std::string ReferenceRecipe(const Engine& e){
    std::ostringstream s;s<<e.workerThreads<<' '<<e.workerChunk<<'\n';
    for(const auto& gpu:e.gpus)s<<gpu<<'\n';for(const auto& v:e.base)s<<v<<'\n';
    if(e.phiCutoffPilot)s<<"matched-cutoff-phi-pilot-v1\n";
    if(!e.angleBatchPilotRecipe.empty())s<<"angle-batch-pilot "<<e.angleBatchPilotRecipe<<'\n';
    return Hash(s.str());
}
void Recipes(){
    ArgPP args;RunConfig config;Engine e(args,config);
    auto check=[&](){for(int i=0;i<3;++i)Require(e.ExecutionHash()==ReferenceRecipe(e),"execution recipe changed");};
    check();e.base={"--cutoff-profile","off"};check();e.base.back()="balanced";check();
    e.workerThreads=7;check();e.workerChunk=256;check();e.gpus={"0","1"};check();
    std::swap(e.gpus[0],e.gpus[1]);check();e.phiCutoffPilot=true;check();
    e.angleBatchPilotRecipe="short-packet-candidate";check();
    e.phiCutoffPilot=false;check();e.angleBatchPilotRecipe.clear();check();
    e.base.clear();e.gpus.clear();check();
}
int main(int argc,char** argv){try{
    Recipes();if(argc==1){std::cout<<"PASS exact recipe and every input mutation\n";return 0;}
    Require(argc==2,"usage: probe [COPIED metadata directory]");
    ArgPP args;RunConfig config;config.useGpu=true;Engine e(args,config);e.root=argv[1];
    std::istringstream(Read(e.root+"/configuration.txt"))>>e.signature;
    e.theta=ReadGrid(e.root+"/requested_theta.csv");e.gpuCandidates={"0","1","2","3",""};
    e.seeds.assign(production,production+8);e.compactIndex=true;e.forbidNativeWork=true;
    Require(e.RestoreTuning(),"missing tuning");
    AdaptiveState state;state.counts=Budgets(3,std::vector<int>(e.theta.size(),e.minimum));
    state.previous=Matrix(e.theta.size());state.valid.assign(e.theta.size(),false);state.coarse=state.valid;
    Require(RestoreAdaptive(e,state,false),"missing adaptive state");
    const double start=Seconds();e.restoringIndex=true;const auto restored=e.RestoreIndexJournal();e.restoringIndex=false;
    std::cerr<<"restored "<<restored<<" records in "<<Seconds()-start<<" s\n";
    Require(!e.cache.empty(),"no saved jobs");std::set<int> levels;
    for(const auto& item:e.cache)if(e.SameSampling(item.second,item.second.count,item.second.seed,item.second.depth,item.second.offset))levels.insert(item.second.depth);
    Require(levels.size()==3,"replay expects three saved depth levels");
    const std::vector<int> depths(levels.begin(),levels.end());
    std::ostringstream result;result<<std::setprecision(17);
    const double begin=Seconds();size_t queries=0;
    AdditionalWork work(e,e.seeds.data());
    for(size_t id:state.knots)for(size_t level=0;level<3;++level){
        auto counts=state.counts;
        if(counts[level][id]<=e.maxCount/2)counts[level][id]*=2;
        const auto forecast=work.Forecast(BudgetRequests(e,depths,counts,{id}));
        result<<id<<' '<<level<<' '<<forecast.seconds<<' '<<forecast.jobs<<' '<<forecast.profiled<<'\n';++queries;
    }
    const auto total=work.Forecast(BudgetRequests(e,depths,state.counts,state.knots));
    result<<total.seconds<<' '<<total.jobs<<' '<<total.profiled<<'\n';
    std::cout<<"{\"queries\":"<<queries+1<<",\"jobs\":"<<e.cache.size()<<",\"seconds\":"<<Number(Seconds()-begin)
        <<",\"forecast_hash\":"<<Quoted(Hash(result.str()))<<",\"recipe\":"<<Quoted(e.ExecutionHash())<<"}\n";
    return 0;
}catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}}
