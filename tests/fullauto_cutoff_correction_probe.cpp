#include "../src/FullAuto.cpp"
using namespace FullAuto;
void Setup(Engine& e,const std::string& root,const std::string& binary,bool controls,bool absorbing){
    e.root=root;e.signature="cutoff-correction-regression-v1";Directory(root);Directory(root+"/jobs");e.exe=binary;e.binaryHash=FileHash(binary);
    e.controls=controls;e.depthSharing=true;e.gpus={""};e.seeds.assign(production,production+8);e.workerThreads=1;e.hostBudgetMb=512;
    e.phi=4;e.target=8;e.theta={170,175,180};e.initial=1024;e.orientationBlock=512;e.workerPolicy="process";e.resourcePolicy="off";
    e.base={"--method","po","--backend","cpu","--particle","1","2","1","--ri","1.3116",absorbing?".001":"0","--wavelength-um",".532","--symmetry","1","1","--cutoff-profile","off"};e.selectedCutoff="off";
    WriteGrid(root+"/requested_theta.csv",e.theta);
}
void Close(double a,double b,const char* message){Require(std::fabs(a-b)<=2e-12*std::max(1.,std::max(std::fabs(a),std::fabs(b))),message);}
int main(int argc,char** argv){try{
    Require(argc==3,"OUTPUT BINARY");ArgPP args;RunConfig config;
    for(bool controls:{false,true}){
        const std::string root=std::string(argv[1])+"/controls"+Number(controls);Engine e(args,config);Setup(e,root,argv[2],controls,controls);e.WarmCache();
        {AdaptiveState state;state.round=2;state.counts={{1048576}};Statistics noisy;noisy.maxInterval=1;
            const int initial=e.initial;e.initial=524288;e.statistical=.01;
            Require(WorkingPhiPilotDue(e,state,noisy),"severely noisy rule waits eight initial blocks");
            noisy.maxInterval=.05;Require(!WorkingPhiPilotDue(e,state,noisy),"ordinary rule rechecked too early");
            state.round=1;noisy.maxInterval=1;Require(!WorkingPhiPilotDue(e,state,noisy),"one-round pilot triggered");e.initial=initial;}
        auto old=e.DirectStage(512,e.seeds.data(),8,8,e.theta);const auto oldHash=FileHash(old[0].folder+"/complete.txt");
        e.balancedCorrection=true;auto reused=e.DirectStage(512,e.seeds.data(),8,8,e.theta);Require(!reused[0].assembled&&FileHash(old[0].folder+"/complete.txt")==oldHash,"old off blocks rewritten");
        // Queue batches five native legs per paired block, followed by assembly.
        e.collecting=true;e.DirectStage(512,e.seeds.data(),8,8,e.theta,true,false,512,6);e.collecting=false;
        std::vector<size_t> todo;for(size_t i=0;i<e.pendingJobs.size();++i)todo.push_back(i);e.ExecutePrepared(e.pendingJobs,e.pendingCommands,todo);e.pendingJobs.clear();e.pendingCommands.clear();
        auto corrected=e.DirectStage(512,e.seeds.data(),8,8,e.theta,true,false,512,6);
        for(const auto& j:corrected){Require(j.assembled&&j.pairedDepth==6,"paired correction missing");IntegralReader reader;
            for(int depth:{8,6}){const Job target=depth==8?j:e.LowerView(j,6);Matrix expected(3,Mueller{});IntegralScalar integrals;
                std::istringstream sources(Read(j.folder+"/cutoff_sources.tsv"));int d,sign;std::string folder;
                while(sources>>d>>sign>>folder)if(d==depth){const auto source=e.Materialize(*e.folderIndex.at(folder));const auto v=reader.Resolve(folder,depth);
                    if(depth==8)Require(source.pairedDepth==6,"balanced cutoff duplicated the depth trace");
                    integrals.area+=sign*v.area;integrals.ext+=sign*v.ext;
                    for(size_t r=0;r<3;++r)for(int k=0;k<16;++k)expected[r][k]+=sign*source.matrix[r][k];}
                for(size_t r=0;r<3;++r)for(int k=0;k<16;++k)Close(target.matrix[r][k],expected[r][k],"signed matrix wrong");
                const auto actual=reader.Resolve(target.folder,depth);Close(actual.area,integrals.area,"signed projected area wrong");Close(actual.ext,integrals.ext,"signed extinction wrong");
            }
        }
        // A fresh pilot seed must measure the SAME signed estimator as the
        // corrected production, while ordinary holdouts stay direct-off.
        unsigned pilotSeeds[8];for(int i=0;i<8;++i)pilotSeeds[i]=1758000000u+17*i;
        const auto productionHash=e.ExecutionHash();
        auto directPilot=e.Stage(512,pilotSeeds,8,8,e.theta);
        for(const auto& j:directPilot)Require(!j.assembled,"ordinary pilot unexpectedly corrected");
        std::vector<Job> matched;std::string matchedHash;
        {PhiCutoffPilotScope scope(e);matchedHash=e.ExecutionHash();Require(e.ExecutionHash()!=productionHash,"pilot cache identity not separated");
            matched=e.Stage(512,pilotSeeds,8,8,e.theta);
            for(size_t i=0;i<matched.size();++i){const auto& j=matched[i];Require(j.assembled&&j.folder!=directPilot[i].folder,"matched pilot reused off block");
                Matrix expected(3,Mueller{});std::istringstream sources(Read(j.folder+"/cutoff_sources.tsv"));int d,sign;std::string folder;int terms=0;
                while(sources>>d>>sign>>folder){Require(d==8,"wrong pilot depth");const auto part=e.Materialize(*e.folderIndex.at(folder));++terms;
                    Require(part.seed==pilotSeeds[i]&&part.offset==j.offset,"pilot pairing lost");
                    for(size_t r=0;r<3;++r)for(int k=0;k<16;++k)expected[r][k]+=sign*part.matrix[r][k];}
                Require(terms==3,"pilot correction source count");for(size_t r=0;r<3;++r)for(int k=0;k<16;++k)Close(j.matrix[r][k],expected[r][k],"pilot signed assembly wrong");
            }}
        Require(!e.phiCutoffPilot&&e.ExecutionHash()==productionHash,"pilot scope changed production cache identity");
        // Mixed small direct-off and large corrected production prefixes.
        const Matrix reference=Mean(e.Data(directPilot));const Budgets budgets{{256,1024,256}};const std::vector<size_t> ids{0,1,2};
        const auto forecast=EvaluateWorkingPhi(e,{8},budgets,ids,reference,512,pilotSeeds);
        Samples expected(8,reference);std::istringstream trace(forecast.sources);int level,sign,count;unsigned seed;std::string pilotFolder;size_t assembledRows=0,directRows=0;
        while(trace>>level>>sign>>count>>seed>>pilotFolder){const Job job=e.Materialize(*e.folderIndex.at(pilotFolder));const auto group=e.Data(std::vector<Job>{job});
            const size_t replica=(seed-pilotSeeds[0])/17;Require(replica<8,"unknown pilot scramble");
            for(size_t row=0;row<3;++row){if((row==1)!=(job.executionHash==matchedHash))continue;const auto at=std::lower_bound(job.angles.begin(),job.angles.end(),e.theta[row]);
                if(at==job.angles.end()||*at!=e.theta[row])continue;
                if(job.executionHash==matchedHash)++assembledRows;else ++directRows;
                for(int k=0;k<16;++k)expected[replica][row][k]+=std::sqrt(512./budgets[0][row])*group[0][at-job.angles.begin()][k];}}
        Require(assembledRows==8&&directRows==16,"pilot did not split cutoff estimators by row");
        const auto center=Mean(expected);for(auto& sample:expected)for(size_t row=0;row<3;++row)for(int k=0;k<16;++k)sample[row][k]+=reference[row][k]-center[row][k];
        Close(forecast.error,StatisticalError(e,Estimate(expected)),"pilot projected variance differs from literal source reconstruction");
        Require(!e.phiCutoffPilot&&e.ExecutionHash()==productionHash,"evaluation leaked pilot state");
        // Mixed old off/new corrected prefix, both on the same unpruned target.
        auto prefix=e.PrefixStage(1024,e.seeds.data(),8,8,e.theta);IntegralReader reader;
        for(size_t i=0;i<8;++i){for(size_t r=0;r<3;++r)for(int k=0;k<16;++k)Close(prefix[i].matrix[r][k],.5*(old[i].matrix[r][k]+corrected[i].matrix[r][k]),"mixed checkpoint prefix wrong");
            const auto a=reader.Resolve(old[i].folder,8),b=reader.Resolve(corrected[i].folder,8),v=reader.Resolve(prefix[i].folder,8);Close(v.ext,.5*(a.ext+b.ext),"mixed extinction prefix wrong");}
        const auto before=FileHash(corrected[0].folder+"/complete.txt");Engine restored(args,config);Setup(restored,root,argv[2],controls,controls);restored.balancedCorrection=true;restored.WarmCache();restored.forbidNativeWork=true;
        auto replay=restored.DirectStage(512,e.seeds.data(),8,8,e.theta,true,false,512,6);Require(FileHash(replay[0].folder+"/complete.txt")==before,"completed correction rerun");
        auto mixed=restored.PrefixStage(1024,e.seeds.data(),8,8,e.theta);Require(mixed[0].matrix==prefix[0].matrix,"resume changed mixed estimator");
        auto isolated=e.DirectStage(512,e.seeds.data(),1,8,e.theta,true,false,2048);
        const auto folder=isolated[0].folder;std::remove((folder+"/complete.txt").c_str());
        Engine interrupted(args,config);Setup(interrupted,root,argv[2],controls,controls);interrupted.balancedCorrection=true;interrupted.WarmCache();interrupted.forbidNativeWork=true;
        auto recovered=interrupted.DirectStage(512,e.seeds.data(),1,8,e.theta,true,false,2048);
        Require(recovered[0].matrix==isolated[0].matrix,"interrupted assembly changed or reran children");
        const auto provenance=Read(folder+"/cutoff_sources.tsv");Atomic(folder+"/cutoff_sources.tsv",provenance+"corruption\n");bool rejected=false;
        try{std::set<std::string> seen;ValidateStoredJob(folder,seen);}catch(const std::exception&){rejected=true;}
        Require(rejected,"changed signed estimator was accepted");Atomic(folder+"/cutoff_sources.tsv",provenance);


    }
    std::cout<<"PASS corrected paired matrices, controls, absorption, integrals, mixed prefixes, resume\n";return 0;
}catch(const std::exception& ex){std::cerr<<ex.what()<<'\n';return 1;}}
