#define main old_lazy_main
#include "fullauto_lazy_restore_probe.cpp"
#undef main

void SparseRows(){
    Require(Engine::PrefixPlanEstimate(16777216,512,8,6)<48*1024*1024,
            "large prefixes disabled batching by retaining repeated path strings");
    std::vector<double> wide;for(int i=0;i<12012;++i)wide.push_back(double(i)/10);
    for(const auto& request:std::vector<std::vector<double>>{{},{0},{0,42.3,170,9999},wide}){
        std::vector<size_t> reference,actual;std::set<double> wanted(request.begin(),request.end());
        for(size_t r=0;r<wide.size();++r)if(wanted.count(wide[r]))reference.push_back(r);
        Engine::MatchingAngularRows(wide,request,[&](size_t r){actual.push_back(r);});
        Require(actual==reference,"sparse/dense row lookup changed selected coordinates");
    }
}

void PartitionWork(const std::string& root,bool traceBound){
    ArgPP args;RunConfig config;Engine e(args,config);SetupMemory(e,root,"/bin/true",false);
    e.orientationBlock=256;e.initial=8192;e.maxCount=65536;e.local=true;e.depthSharing=false;
    e.theta.clear();std::vector<size_t> ids;for(size_t i=0;i<32;++i){e.theta.push_back(170+i*.25);ids.push_back(i);}
    const std::vector<double> half(e.theta.begin(),e.theta.begin()+16);
    // All current prefixes are complete; half the requested future angles
    // also have saved work. Cost per requested angle-orientation is therefore
    // misleading, although partitioning adds almost no diffraction cost.
    for(int offset=0;offset<16384;offset+=256)for(unsigned seed:e.seeds){
        const auto& grid=offset<8192?e.theta:half;Job j;e.SetSampling(j,256,seed,6,grid,false);
        j.offset=offset;j.identity=Number(offset)+"_"+Number(seed);j.folder=root+"/fake_"+j.identity;
        j.traceSeconds=traceBound?100:.01;j.diffractionSeconds=(traceBound?.01:100)*grid.size()/32.;
        j.seconds=j.traceSeconds+j.diffractionSeconds;e.CacheJob(j);
    }
    Budgets current(1,std::vector<int>(32,8192)),proposed(1,std::vector<int>(32,16384));
    Statistics stats;stats.interval.resize(32);stats.ratioInterval.resize(32);
    for(size_t r=0;r<32;++r)stats.interval[r][0]=.1+r*.01;
    AdditionalWork work(e,e.seeds.data());const auto full=proposed;
    BoundRoundWork(e,work,current,proposed,{6},ids,stats);
    const auto forecast=work.Forecast(BudgetRequests(e,{6},proposed,ids));
    if(traceBound)Require(proposed==full,"time limit multiplied trace-bound work");
    else Require(proposed!=current&&proposed!=full&&forecast.profiled&&forecast.seconds<=1800,
                 "asymmetric cache coverage prevented a bounded useful subround");
    for(size_t id:ids)Require(proposed[0][id]>=current[0][id]&&proposed[0][id]<=full[0][id]
        &&(proposed[0][id]&(proposed[0][id]-1))==0,"partition changed prefix validity");
    std::cout<<"PASS partition trace_bound="<<traceBound<<" selected_seconds="<<forecast.seconds<<'\n';
}

void MixedPhi(const std::string& root,const std::string& binary,bool controls){
    ArgPP args;RunConfig config;Engine old(args,config),batch(args,config);
    SetupMemory(old,root+"/old",binary,controls);SetupMemory(batch,root+"/batch",binary,controls);
    old.theta=batch.theta={0,5,10,15,170,175,180};old.orientationBlock=batch.orientationBlock=32;
    for(auto* e:{&old,&batch}){e->thetaPhi={{0,1}};e->phi=2;WriteGrid(e->root+"/requested_theta.csv",e->theta);
        if(controls){e->weightTheta=e->theta;e->weights.resize(e->theta.size());
            for(size_t r=0;r<e->theta.size();++r)e->weights[r]={{.3+r*.01,.2-r*.005}};}}
    for(int phi:{1,2}){PhiScope scope(old,phi);for(int offset=0;offset<256;offset+=32){
        auto jobs=old.DirectStage(32,old.seeds.data(),8,6,old.theta,true,false,offset);
        for(const auto& j:jobs)batch.CacheJob(j);}}
    old.forbidNativeWork=batch.forbidNativeWork=true;
    old.payloads.Clear();batch.payloads.Clear();old.payloads.SetLimit(48*1024);batch.payloads.SetLimit(48*1024);
    const std::map<int,std::vector<std::vector<double>>> grids={
        {1,{{0,5},{0,10},{0,15}}},{2,{{5,10},{170,175},{175,180}}}};
    std::map<int,std::vector<std::vector<Job>>> references;
    const auto oldReads=old.payloads.misses,newReads=batch.payloads.misses;
    for(const auto& group:grids){PhiScope scope(old,group.first);
        for(const auto& grid:group.second)references[group.first].push_back(old.PrefixStage(256,old.seeds.data(),8,6,grid));}
    batch.ClearPrefixAssemblyPlan();
    for(const auto& group:grids){PhiScope scope(batch,group.first);std::vector<Engine::WorkRequest> requests;
        for(const auto& grid:group.second)requests.push_back({256,0,6,0,grid});batch.AddPrefixAssemblyPlan(requests,batch.seeds.data());}
    Require(batch.prefixAssemblyPlans.size()==6,"mixed phi plan disabled or conflated");
    for(const auto& group:grids){PhiScope scope(batch,group.first);
        for(size_t k=0;k<group.second.size();++k){const auto result=batch.PrefixStage(256,batch.seeds.data(),8,6,group.second[k]);
            for(size_t s=0;s<8;++s){const auto& a=result[s];const auto& b=references[group.first][k][s];
                Require(a.matrix==b.matrix&&a.controls==b.controls&&a.controlSetup==b.controlSetup
                    &&a.controlRefinement==b.controlRefinement,"phi-routed batching changed optical estimator");
                Require(Read(a.folder+"/angular_sources.tsv")==Read(b.folder+"/angular_sources.tsv"),"phi-routed batching changed source order");}
            if(k==0)Require(batch.preparedPrefixes.size()==2,"unlike phi groups were mixed or batching was disabled");}}
    Require(batch.preparedPrefixes.empty()&&batch.preparedPrefixBytes==0,"mixed phi batches retained payloads");
    Require((batch.payloads.misses-newReads)*2<old.payloads.misses-oldReads,"mixed phi batching did not reuse wide payloads");
    std::cout<<"PASS mixed phi controls="<<controls<<": all16/controls/source order exact; payload loads "
        <<old.payloads.misses-oldReads<<" -> "<<batch.payloads.misses-newReads<<'\n';
}

int main(int argc,char** argv){try{
    Require(argc==3,"OUTPUT WORKER_BINARY");SparseRows();
    PartitionWork(std::string(argv[1])+"/diffraction",false);PartitionWork(std::string(argv[1])+"/trace",true);
    for(bool controls:{false,true})MixedPhi(std::string(argv[1])+"/mixed_phi_"+Number(controls),argv[2],controls);
    std::cout<<"PASS: sparse lookup, bounded native-work partition and mixed-phi exact prefix reuse\n";return 0;
}catch(const std::exception& error){std::cerr<<error.what()<<'\n';return 1;}}
