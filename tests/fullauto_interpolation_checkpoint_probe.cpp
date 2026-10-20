#define main lazy_restore_test_main
#include "fullauto_lazy_restore_probe.cpp"
#undef main
int main(int argc,char** argv){try{
    Require(argc==3,"OUTPUT BINARY");
    for(bool controls:{false,true}){
        ArgPP args;RunConfig config;Engine e(args,config);SetupMemory(e,std::string(argv[1])+"/controls"+Number(controls),argv[2],controls);
        e.workQueuePolicy="off";e.workerPolicy="process";e.orientationBlock=64;e.initial=128;e.minimum=64;e.depthSharing=false;
        e.theta={170,172,174,176,178,179,180};e.weightTheta=e.theta;e.weights.assign(e.theta.size(),{{0.25,0.125}});
        WriteGrid(e.root+"/requested_theta.csv",e.theta);
        e.polarWeights.resize(e.theta.size());for(auto& row:e.polarWeights)for(auto& control:row)for(double& value:control)value=0.125;
        const std::vector<int> depths={4,6};const std::vector<size_t> ids={0,1,2,3,4,5,6};
        AdaptiveState state;state.knots={0,2,4,6};state.counts=Budgets(2,std::vector<int>(e.theta.size(),64));state.counts[0].assign(e.theta.size(),128);
        auto original=CommonPrefixInterpolationImpl(e,state,depths,ids,e.seeds.data());
        const std::string sources=Read(e.root+"/interpolation_production_sources.tsv"),budgets=Read(e.root+"/interpolation_production_budgets.csv");
        const auto groups=FullAutoAuditPlan::Stencils(state.knots,state.counts,ids);Require(!groups.empty(),"fixture has no interpolation groups");
        const auto& group=*groups.begin();const std::vector<size_t> rows(group.second.rows.begin(),group.second.rows.end());
        const auto counts=StencilBudgets(e,state.counts,group);const auto key=InterpolationGroupKey(e,depths,counts,rows,e.seeds.data());
        AdaptiveData saved;Require(RestoreInterpolationGroup(e,key,depths.size(),rows.size(),saved),"missing saved group");
        const auto direct=EvaluateBudget(e,depths,counts,rows,e.seeds.data(),false);
        Require(saved.levels==direct.levels&&saved.costs==direct.costs&&saved.analytic==direct.analytic&&saved.sources==direct.sources,"checkpoint roundtrip changed eight-replica data");
        e.payloads.SetLimit(1);std::vector<double> assemblyTimes,receiptTimes;
        for(int trial=0;trial<5;++trial){
            const double begin=Seconds();const auto assembled=EvaluateBudget(e,depths,counts,rows,e.seeds.data(),false);assemblyTimes.push_back(Seconds()-begin);
            AdaptiveData receipt;const double resumedAt=Seconds();Require(RestoreInterpolationGroup(e,key,depths.size(),rows.size(),receipt),"benchmark receipt absent");receiptTimes.push_back(Seconds()-resumedAt);
            Require(receipt.levels==assembled.levels,"benchmark changed replicas");
        }
        std::sort(assemblyTimes.begin(),assemblyTimes.end());std::sort(receiptTimes.begin(),receiptTimes.end());
        Atomic(e.root+"/group_resume_benchmark.json","{\"trials\":5,\"assembly_median_seconds\":"+Number(assemblyTimes[2])
            +",\"receipt_median_seconds\":"+Number(receiptTimes[2])+",\"speedup\":"+Number(assemblyTimes[2]/receiptTimes[2])
            +",\"native_work_repeated\":false,\"all16_and_replicas_exact\":true,\"scope\":\"small saved CPU prism fixture; controller-only replay\"}");
        // A resumed interpolation must not re-enter prefix assembly or workers.
        e.forbidNativeWork=true;e.payloads.Clear();e.cache.clear();e.samplingIndex.clear();e.folderIndex.clear();e.warmedCacheControls.insert(controls);
        const auto resumed=CommonPrefixInterpolationImpl(e,state,depths,ids,e.seeds.data());
        Require(resumed.errors==original.errors&&resumed.matrixErrors==original.matrixErrors&&resumed.ratioErrors==original.ratioErrors
            &&resumed.bad==original.bad&&resumed.maximum==original.maximum,"resume changed all16/ratio interpolation bounds");
        Require(e.cache.empty(),"resume reassembled native prefixes");
        Require(Read(e.root+"/interpolation_production_sources.tsv")==sources&&Read(e.root+"/interpolation_production_budgets.csv")==budgets,"resume changed provenance");
        ++e.phi;Require(InterpolationGroupKey(e,depths,counts,rows,e.seeds.data())!=key,"phi not frozen");--e.phi;
        auto otherCounts=counts;otherCounts[0][rows[0]]*=2;Require(InterpolationGroupKey(e,depths,otherCounts,rows,e.seeds.data())!=key,"budget not frozen");
        auto otherSeeds=e.seeds;++otherSeeds[0];Require(InterpolationGroupKey(e,depths,counts,rows,otherSeeds.data())!=key,"seed not frozen");
        if(controls){e.weights[0][0]+=0.01;Require(InterpolationGroupKey(e,depths,counts,rows,e.seeds.data())!=key,"control coefficient not frozen");e.weights[0][0]-=0.01;}
        // Both snapshot corruption and descendant corruption must be rejected.
        const auto path=InterpolationGroupPath(e,key),text=Read(path);Atomic(path,text+"corrupt\n");
        bool rejected=false;try{AdaptiveData bad;RestoreInterpolationGroup(e,key,depths.size(),rows.size(),bad);}catch(const std::exception&){rejected=true;}
        Atomic(path,text);Require(rejected,"corrupt group accepted");
        const auto roots=InterpolationGroupRoots(saved.sources);std::string leaf=*roots.begin();
        for(;;){if(!Exists(leaf+"/angular_sources.tsv"))break;
            std::istringstream input(Read(leaf+"/angular_sources.tsv"));std::string digest;Require(bool(input>>leaf>>digest),"invalid fixture descendant");}
        const auto sourceData=Read(leaf+"/result/result.dat");Atomic(leaf+"/result/result.dat",sourceData+"corrupt\n");rejected=false;
        try{ValidationBatch fresh("corruption test",true);AdaptiveData bad;RestoreInterpolationGroup(e,key,depths.size(),rows.size(),bad);}catch(const std::exception&){rejected=true;}
        Atomic(leaf+"/result/result.dat",sourceData);Require(rejected,"corrupt descendant accepted by group receipt");
        std::cout<<"PASS controls="<<controls<<": exact eight replicas/all16/ratios/provenance, no reassembly, frozen inputs, corrupt snapshot and source rejected\n";
    }
    return 0;
}catch(const std::exception& error){std::cerr<<error.what()<<'\n';return 1;}}
