#include "../src/FullAuto.cpp"
using namespace FullAuto;
int main(int argc,char** argv){try{
    Require(argc==3,"expected actual front fixture and private output directory");
    size_t rows,levels,left,right;std::ifstream input(argv[1]);
    Require(bool(input>>rows>>levels>>left>>right)&&rows==33&&levels==3,"invalid actual front fixture");
    std::vector<double> theta(rows);std::vector<int> phi(rows);std::vector<size_t> ids;
    Budgets counts(levels,std::vector<int>(rows));
    for(size_t r=0;r<rows;++r){ids.push_back(r);Require(bool(input>>theta[r]>>counts[0][r]>>counts[1][r]>>counts[2][r]>>phi[r]),"short front fixture");}
    const auto original=counts;const std::vector<size_t> stencil={left,right};
    auto build=[&](const Budgets& budget,const std::vector<size_t>& nodes,const std::vector<size_t>& evaluated){
        return FullAutoAngularCohort::Build(ids,evaluated,budget,{0,1,2},nodes,1<<24,[&](size_t id){return phi[id];});};
    const auto actual=build(counts,stencil,ids);
    Require(actual.nodes==std::vector<size_t>({6,7,8,9,10,11,12,13}),"actual front did not spend bounded cohort on seven low-prefix nodes");
    Require(actual.targets==std::vector<std::pair<size_t,int>>({{0,2048}}),"equal-prefix correction levels recruited");
    Require(counts==original,"planner mutated frozen counts");
    Require(build(counts,{13,14},ids).targets.empty(),"proposal recursively propagated across original equal prefixes");
    for(int barrier:{512,4096}){auto b=counts;b[0][8]=barrier;const auto plan=build(b,stencil,ids);
        Require(std::find(plan.nodes.begin(),plan.nodes.end(),8)==plan.nodes.end()&&plan.nodes.back()==7,"planner skipped a low/high-prefix barrier");}
    phi[8]=64;Require(build(counts,stencil,ids).nodes.back()==7,"planner crossed phi boundary");phi[8]=32;
    auto missing=ids;missing.erase(missing.begin()+8);
    Require(build(counts,stencil,missing).nodes.back()==7,"planner crossed unevaluated support");
    missing=ids;missing.erase(missing.begin()+5);
    Require(build(counts,{4,6,8},missing).targets.empty(),"span-two stencil crossed missing interior support");
    auto endpoint=counts;for(int& n:endpoint[0])n=1024;endpoint[0][0]=2048;
    Require(build(endpoint,{0,1},ids).nodes==std::vector<size_t>({0,1,2,3,4,5,6,7}),"left endpoint underflow or incomplete useful halo");
    for(int& n:endpoint[0])n=1024;endpoint[0][32]=2048;
    Require(build(endpoint,{31,32},ids).nodes==std::vector<size_t>({25,26,27,28,29,30,31,32}),"right endpoint overflow");
    auto multi=counts;multi[1][6]=2048;multi[2][6]=2048;
    Require(build(multi,stencil,ids).nodes.size()==8,"eight-node cap multiplied by level count");
    auto large=counts;for(int& n:large[0])n=1073741824;large[0][6]=2147483647;
    auto maxPlan=FullAutoAngularCohort::Build(ids,ids,large,{0},stencil,2147483647,[&](size_t id){return phi[id];});
    Require(maxPlan.nodes.size()==8&&maxPlan.targets[0].second==2147483647,"large integer target overflow");
    Require(FullAutoAngularCohort::Build(ids,ids,counts,{0},stencil,1024,[&](size_t id){return phi[id];}).targets.empty(),"maximum count exceeded");

    ArgPP args;RunConfig config;Engine e(args,config);
    e.root=argv[2];Directory(e.root);e.signature="angular-cohort-probe";e.exe="/bin/true";e.binaryHash=FileHash(e.exe);
    e.theta=theta;e.phi=32;e.eps=.04;e.interpolation=.01;e.go=false;e.coherent=false;e.local=true;e.ratios=true;
    e.maxCount=1<<24;e.depthSharing=false;e.controls=false;e.gpus={""};e.seeds.assign(production,production+8);
    e.incremental=true;e.angularCache=true;e.orientationBlock=256;e.target=8;
    AdaptiveState state;state.counts=counts;state.knots=ids;
    Samples samples(8,Matrix(rows));
    for(size_t s=0;s<8;++s)for(size_t r=0;r<rows;++r){samples[s][r][0]=1;samples[s][r][15]=.3+(r<=left?(s<4?-.01:.01):0);}
    AdaptiveData stage;stage.levels.push_back(samples);stage.levels.push_back(Samples(8,Matrix(rows)));stage.levels.push_back(Samples(8,Matrix(rows)));
    const auto samplesBefore=samples;const auto stats=Estimate(samples);
    const auto decision=CheckAngularSampling(e,state,stage,ids,samples,stats);
    Require(!decision.passed&&decision.haloCohorts>0&&decision.maximumCohortNodes==8,"actual gate did not recruit bounded halo");
    Require(decision.proposed[0][13]==2048&&decision.proposed[0][14]==1024,"integration failed to include useful halo or spread recursively");
    Require(decision.proposed[1]==counts[1]&&decision.proposed[2]==counts[2],"non-dominating levels changed");
    Require(samples==samplesBefore&&state.counts==counts,"acceptance samples or current allocation modified");
    for(const auto& group:decision.priorities)Require(group.nodes.size()<=8,"union priority exceeds node cap");
    // Every newly grown halo node belongs to its full scheduling cohort.
    for(size_t r=0;r<rows;++r)if(decision.proposed[0][r]>counts[0][r]){
        bool found=false;for(const auto& group:decision.priorities)found=found||std::find(group.nodes.begin(),group.nodes.end(),r)!=group.nodes.end();
        Require(found,"new halo growth lost its priority cohort");
    }
    SaveAngularSampling(e,decision,e.root);
    auto restored=SavedAngularPriorities(e,state,e.root);
    Require(restored.size()==decision.priorities.size(),"saved cohort count changed");
    for(size_t k=0;k<restored.size();++k)Require(restored[k].nodes==decision.priorities[k].nodes&&restored[k].error==decision.priorities[k].error,"resume lost expanded cohort");
    std::vector<RoundOutputPriority> untouched={{1,{0,1}}};
    const auto before=untouched;
    Require(!ReadSavedAngularCohorts(e,state,"error_scaled_M11,node_count,expanded,theta_degrees\n",restored.size(),untouched)
        &&untouched.size()==before.size()&&untouched[0].nodes==before[0].nodes,"truncated table published partial priorities");
    const std::string header="error_scaled_M11,node_count,expanded,theta_degrees\n";
    for(const auto& row:std::vector<std::string>{"1,2,0,150,150.625\n","1,3,0,150,150.3125,150.9375\n","1,4,1,150,150.3125,150.9375,151.25\n","1,9,1,150\n"})
        Require(!ReadSavedAngularCohorts(e,state,header+row,1,untouched),"malformed saved cohort accepted");
    Atomic(e.root+"/angular_cohorts.csv",header);
    const auto fallback=SavedAngularPriorities(e,state,e.root);
    Require(fallback.size()==decision.priorities.size(),"truncated cohort table disabled legacy fallback");
    for(const auto& group:fallback)Require(group.nodes.size()<=3,"legacy fallback did not use original stencils");
    Atomic(e.root+"/angular_cohorts.csv",decision.cohortTable);
    // Unknown-cost admission must retain a COMPLETE expanded cohort.
    Budgets selected=decision.proposed;AdditionalWork work(e,e.seeds.data());
    BoundRoundWork(e,work,counts,selected,{4,6,8},ids,stats,ids,{},nullptr,restored);
    bool complete=false;
    for(const auto& group:restored){bool all=true,grew=false;
        for(size_t id:group.nodes)for(size_t level=0;level<counts.size();++level){
            all=all&&selected[level][id]==decision.proposed[level][id];grew=grew||selected[level][id]>counts[level][id];}
        complete=complete||(all&&grew&&group.nodes.size()>3);
    }
    Require(complete,"unknown-cost round limiter split/discarded expanded cohort");
    e.coherent=true;const auto disabled=CheckAngularSampling(e,state,stage,ids,samples,stats);
    Require(disabled.cohortTable==header&&!disabled.enabled&&disabled.priorities.empty(),"disabled mode omitted header or recruited halo");
    std::cout<<"PASS actual round5 front: seven low-prefix nodes in one <=8-node cohort; barriers, frozen counts, integer bounds, same phi, evaluated support, gate integration, durable hint restore/fallback, unknown-cost whole-cohort admission\n";
    return 0;
}catch(const std::exception& ex){std::cerr<<ex.what()<<'\n';return 1;}}
