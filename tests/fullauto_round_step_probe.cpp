// Scheduling regressions only. Synthetic timing/coverage fixtures do not
// claim optical convergence or measured GPU speedup.
#include "../src/FullAuto.cpp"
using namespace FullAuto;
namespace {
void Init(Engine& e,const std::string& root){
    e.root=root;e.signature="round-prefix-step-probe";Directory(root);
    e.exe="/bin/true";e.binaryHash=FileHash(e.exe);e.controls=false;
    e.depthSharing=false;e.gpus={""};e.seeds.assign(production,production+8);
    e.local=true;e.incremental=true;e.angularCache=true;e.orientationBlock=256;
    e.phi=1;e.target=6;e.initial=2048;e.maxCount=65536;
    e.theta.clear();for(size_t i=0;i<32;++i)e.theta.push_back(150+i*.1);
}
void Coverage(Engine& e,double seconds){
    for(int offset=0;offset<2048;offset+=256)for(unsigned seed:e.seeds){
        Job job;e.SetSampling(job,256,seed,6,e.theta,false);job.offset=offset;
        job.identity=Number(offset)+"_"+Number(seed);job.folder=e.root+"/fake_"+job.identity;
        job.seconds=seconds;job.traceSeconds=.01;job.diffractionSeconds=seconds-.01;
        e.CacheJob(job);
    }
}
void CheckBounds(const Budgets& before,const Budgets& wanted,const Budgets& after){
    for(size_t l=0;l<before.size();++l)for(size_t r=0;r<before[l].size();++r)
        Require(after[l][r]>=before[l][r]&&after[l][r]<=wanted[l][r]
            &&!(after[l][r]&(after[l][r]-1)),"invalid/shrinking selected prefix");
}
void Run(const std::string& root,double timing,bool pair,bool expectWithin){
    ArgPP args;RunConfig config;Engine e(args,config);Init(e,root);
    e.nonconvexGeometry=true;
    if(timing>0)Coverage(e,timing);
    std::vector<size_t> ids;Statistics stats;stats.interval.resize(e.theta.size());stats.ratioInterval.resize(e.theta.size());
    std::vector<RoundOutputPriority> groups;
    for(size_t r=0;r<e.theta.size();++r){ids.push_back(r);stats.interval[r][0]=.1;}
    for(size_t r=0;r<e.theta.size();r+=pair?2:1)
        groups.push_back({1.-r*.01,pair?std::vector<size_t>{r,r+1}:std::vector<size_t>{r}});
    Budgets current(1,std::vector<int>(32,2048)),wanted(1,std::vector<int>(32,8192)),selected=wanted;
    const auto seeds=e.seeds;AdditionalWork work(e,e.seeds.data());
    BoundRoundWork(e,work,current,selected,{6},ids,stats,{}, {},nullptr,groups);
    CheckBounds(current,wanted,selected);Require(e.seeds==seeds&&e.seeds.size()==8,"replica selection changed");
    for(const auto& group:groups)for(size_t r:group.nodes)
        Require(selected[0][r]==selected[0][group.nodes.front()],"complete paired cohort split");
    const auto receipt=Read(root+"/round_work_limit.json");
    if(!timing){
        Require(selected[0][0]==8192&&selected[0][2]==2048,"unknown-cost complete-cohort fallback changed");
        Require(receipt.find("\"bounded_prefix_step_applied\":false")!=std::string::npos,"unknown cost was treated as measured");
    }else{
        const auto forecast=work.Forecast(BudgetRequests(e,{6},selected,ids));
        Require(selected[0][0]==4096,"oversized highest-priority cohort did not receive one 2x step");
        Require(forecast.profiled&&(forecast.seconds<=1800)==expectWithin,"forecast budget receipt inconsistent");
        Require(receipt.find("\"bounded_prefix_step_applied\":true")!=std::string::npos,"shortened-prefix receipt absent");
        if(expectWithin&&!pair)Require(selected[0][1]==4096&&selected[0][2]==4096,
            "remaining budget did not widen to other whole ranked cohorts");
        if(!expectWithin&&!pair){
            auto singleton=current;singleton[0][0]=4096;
            const auto one=work.Forecast(BudgetRequests(e,{6},singleton,ids));
            Require(selected[0][1]==4096&&selected[0][2]==2048&&forecast.jobs<=one.jobs&&forecast.seconds<=one.seconds,
                "irreducible prefix wasted its second native row or increased forecast cost");
        }
        for(int count:selected[0])Require(count==2048||count==4096,"one-step growth bound lost");
        std::cout<<"known_timing="<<timing<<" paired_cohort="<<pair<<" forecast_seconds="<<forecast.seconds<<'\n';
    }
    // Replay of the same frozen metadata must choose exactly the same plan.
    auto repeated=wanted;BoundRoundWork(e,work,current,repeated,{6},ids,stats,{}, {},nullptr,groups);
    Require(repeated==selected,"selection depends on proposal traversal history");
}
void UnchangedCases(const std::string& root,bool local){
    ArgPP args;RunConfig config;Engine e(args,config);Init(e,root);Coverage(e,840);e.local=local;
    std::vector<size_t> ids;Statistics stats;stats.interval.resize(e.theta.size());stats.ratioInterval.resize(e.theta.size());
    for(size_t r=0;r<e.theta.size();++r){ids.push_back(r);stats.interval[r][0]=.1;}
    Budgets current(1,std::vector<int>(32,2048)),wanted(1,std::vector<int>(32,local?4096:8192)),selected=wanted;
    AdditionalWork work(e,e.seeds.data());
    BoundRoundWork(e,work,current,selected,{6},ids,stats,{}, {},nullptr,{{1.,ids}});
    Require(selected==wanted,local?"already 2x irreducible whole cohort changed":"nonlocal allocation changed");
    Require(Read(root+"/round_work_limit.json").find("\"bounded_prefix_step_applied\":false")!=std::string::npos,
        "inactive prefix-step policy reported as applied");
}
void CompatibleBackfill(const std::string& root){
    ArgPP args;RunConfig config;Engine e(args,config);Init(e,root);
    std::vector<size_t> ids;Statistics stats;stats.interval.resize(32);stats.ratioInterval.resize(32);
    for(size_t r=0;r<32;++r){ids.push_back(r);stats.interval[r][0]=.1;}
    Budgets current(1,std::vector<int>(32,2048));current[0][0]=current[0][31]=65536;current[0][2]=8192;
    auto wanted=current;for(size_t id:{size_t(1),size_t(2),size_t(3),size_t(4),size_t(5)})wanted[0][id]*=4;
    Job timing;e.SetSampling(timing,256,9999,6,e.theta,false);timing.identity="timing";timing.folder=root+"/timing";
    timing.seconds=280;timing.traceSeconds=.01;timing.diffractionSeconds=279.99;e.CacheJob(timing);
    for(int offset=0;offset<65536;offset+=256){std::vector<double> rows;
        for(size_t id:ids)if(current[0][id]>=offset+256)rows.push_back(e.theta[id]);
        for(unsigned seed:e.seeds){Job job;e.SetSampling(job,256,seed,6,rows,false);job.offset=offset;job.assembled=true;
            job.identity=Number(offset)+"_"+Number(seed);job.folder=root+"/"+job.identity;e.CacheJob(job);}
    }
    // The second-ranked cohort occupies a disjoint, much longer prefix.
    // Later whole cohorts can share the first one's exact Sobol interval.
    const std::vector<RoundOutputPriority> groups={{1.,{1}},{.9,{2}},{.8,{3,6}},{.7,{4,7}},{.6,{5,8}}};
    auto selected=wanted;AdditionalWork work(e,e.seeds.data());
    BoundRoundWork(e,work,current,selected,{6},ids,stats,{}, {},nullptr,groups);
    CheckBounds(current,wanted,selected);
    Require(selected[0][1]==4096&&selected[0][2]==8192&&selected[0][3]==4096&&selected[0][4]==4096
        &&selected[0][5]==2048,"costly intervening cohort hid compatible higher-priority backfill");
    for(size_t id:{size_t(6),size_t(7),size_t(8)})Require(selected[0][id]==current[0][id],"backfill invented unwanted stencil growth");
    const auto forecast=work.Forecast(BudgetRequests(e,{6},selected,ids));
    Require(forecast.profiled&&forecast.seconds<=1800,"compatible backfill exceeded exact budget");
    const auto receipt=Read(root+"/round_work_limit.json");
    Require(receipt.find("\"sharing_backfill_cohorts\":2")!=std::string::npos
        &&receipt.find("\"selected_cohorts\":3")!=std::string::npos,"non-prefix admissions omitted from receipt");
    std::cout<<"PASS blocked-rank backfill: worst cohort retained, free padding replacement then three useful rows, "<<forecast.seconds<<" forecast seconds\n";
}

void TraceSharingBatch(const std::string& root){
    ArgPP args;RunConfig config;Engine e(args,config);Init(e,root);
    const size_t rows=192;e.theta.clear();std::vector<size_t> ids;Statistics stats;
    stats.interval.resize(rows);stats.ratioInterval.resize(rows);
    for(size_t i=0;i<rows;++i){ids.push_back(i);e.theta.push_back(i*.5);stats.interval[i][0]=1.+double(rows-i)/rows;}
    for(unsigned seed:e.seeds){
        Job j;e.SetSampling(j,256,seed,6,e.theta,false);j.identity=Number(seed);j.folder=root+"/"+j.identity;
        j.seconds=750;j.traceSeconds=60;j.diffractionSeconds=690;e.CacheJob(j);
    }
    Budgets current(1,std::vector<int>(rows,256)),wanted(1,std::vector<int>(rows,512));
    AdditionalWork work(e,e.seeds.data());
    auto narrow=wanted;BoundRoundWork(e,work,current,narrow,{6},ids,stats);
    e.nonconvexGeometry=true;auto wide=wanted;BoundRoundWork(e,work,current,wide,{6},ids,stats);
    const auto before=work.Forecast(BudgetRequests(e,{6},narrow,ids));
    const auto after=work.Forecast(BudgetRequests(e,{6},wide,ids));
    const auto oldRows=std::count(narrow[0].begin(),narrow[0].end(),512);
    const auto newRows=std::count(wide[0].begin(),wide[0].end(),512);
    std::cout<<"compatible rows "<<oldRows<<" -> "<<newRows<<" forecast "<<before.seconds<<" -> "<<after.seconds
             <<" jobs "<<before.jobs<<" -> "<<after.jobs<<'\n';
    CheckBounds(current,wanted,wide);
    Require(newRows>oldRows&&before.profiled&&after.profiled&&after.jobs==before.jobs&&after.seconds<=1800,
            "same-prefix batch failed to widen within the existing jobs and time budget");
    for(size_t id:ids)Require(wide[0][id]==(id<size_t(newRows)?512:256),"compatible rows lost ranked order");
}

}
int main(int argc,char** argv){try{
    Require(argc==2,"OUTPUT");const std::string root=argv[1];
    Run(root+"/wide",280,false,true);Run(root+"/paired",280,true,true);
    Run(root+"/irreducible",840,false,false);Run(root+"/unknown",0,true,false);
    UnchangedCases(root+"/nonlocal",false);UnchangedCases(root+"/already2x",true);
    CompatibleBackfill(root+"/compatible_backfill");
    TraceSharingBatch(root+"/trace_sharing_batch");
    std::cout<<"PASS: oversized 4x cohort becomes bounded 2x; measured breadth; whole cohorts; eight seeds; no rollback; deterministic; unknown/irreducible costs honest\n";
    return 0;
}catch(const std::exception& error){std::cerr<<error.what()<<'\n';return 1;}}
