// Real optics, frozen unequal orientation budgets, no adaptive acceptance.
#include "../src/FullAuto.cpp"
int main(int argc,const char* argv[]){try{
    ArgPP args;RegisterCliOptions(args);args.Parse(argc,argv);
    RunConfig config=RunConfig::FromCommandLine(args,false,false);
    FullAuto::Engine e(args,config);e.Init(argc,argv);
    const char* worker=getenv("TASK_MBS_WORKER");FullAuto::Require(worker,"TASK_MBS_WORKER required");
    e.exe=worker;e.binaryHash=FullAuto::FileHash(e.exe);e.controls=false;e.sharedMeans=false;e.SetMirror(false);
    e.rotationBeta=e.rotationGamma=1;e.seeds.assign(FullAuto::production,FullAuto::production+8);
    e.phi=args.GetIntValue("nphi");e.depthSharing=e.coherent;e.depthVerified=true;
    FullAuto::Directory(e.root+"/jobs");FullAuto::WriteGrid(e.root+"/requested_theta.csv",e.theta);
    e.WarmCache();
    FullAuto::Samples samples;std::string sources;
    if(e.go)samples=FullAuto::GoSample(e,e.initial,e.seeds.data(),-1,&sources);
    else{
        FullAuto::Budgets counts(3,std::vector<int>(e.theta.size()));std::vector<size_t> ids;
        for(size_t i=0;i<e.theta.size();++i){ids.push_back(i);counts[0][i]=e.initial;counts[1][i]=i%2?e.initial/2:e.initial/4;counts[2][i]=e.initial/4;}
        const auto data=FullAuto::EvaluateBudget(e,{1,2,e.target},counts,ids,e.seeds.data());samples=FullAuto::Combine(data.levels);sources=data.sources;
    }
    FullAuto::Atomic(e.root+"/final_stage_sources.tsv",sources);
    auto save=[&](const std::string& path,const FullAuto::Matrix& data){
        FullAuto::WriteMatrix(path,e.theta,FullAuto::Matrix(data.begin(),data.begin()+e.theta.size()));
        if(e.go){std::ofstream file(path.c_str(),std::ios::app);file<<std::setprecision(17);
            for(size_t i=e.theta.size();i<data.size();++i){file<<(i==e.theta.size()?0:180)<<" 0";for(double v:data[i])file<<' '<<v;file<<'\n';}}
    };
    for(size_t i=0;i<8;++i)save(e.root+"/replica"+FullAuto::Number(i)+".dat",samples[i]);
    save(e.root+"/frozen_matrix.dat",FullAuto::Mean(samples));
    const auto serial=FullAuto::IntegralProductionSources(sources,e.seeds,1),parallel=FullAuto::IntegralProductionSources(sources,e.seeds,4);
    FullAuto::Require(serial.files==parallel.files&&serial.selectedSources==parallel.selectedSources&&serial.directTargetFallback==parallel.directTargetFallback,"parallel integral provenance differs");
    for(size_t i=0;i<8;++i)FullAuto::Require(serial.samples[i].area==parallel.samples[i].area&&serial.samples[i].ext==parallel.samples[i].ext,"parallel integral scalar differs");
    FullAuto::WriteFullAutoIntegrals(e,e.root,samples,sources,true);FullAuto::ValidateIntegralSources(e.root);
    // New native integral files must be covered by their completed-block receipt.
    auto entry=std::find_if(e.cache.begin(),e.cache.end(),[](const std::pair<const std::string,FullAuto::Job>& item){return !item.second.assembled;});
    FullAuto::Require(entry!=e.cache.end(),"missing native work");
    const auto path=entry->second.folder+"/result/result_integrals.tsv",saved=FullAuto::Read(path);
    FullAuto::Atomic(path,saved+" ");bool rejected=false;
    try{std::set<std::string> visited;FullAuto::ValidateStoredJob(entry->second.folder,visited);}catch(const std::exception&){rejected=true;}
    FullAuto::Atomic(path,saved);FullAuto::Require(rejected,"changed native integrals accepted");
    return 0;
}catch(const std::exception& error){std::cerr<<error.what()<<std::endl;return 1;}}
