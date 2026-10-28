// This probe explicitly tests a native build with prefix-output support.
#define MBS_NATIVE_ORIENTATION_PREFIX_OUTPUT 1
#include "../src/FullAuto.cpp"
using namespace FullAuto;
void Setup(Engine& e,const std::string& root,const std::string& binary,bool controls,bool absorbing,bool incoherent=false){
    e.compactIndex=true;
    e.root=root;e.signature="cutoff-correction-regression-v1";Directory(root);Directory(root+"/jobs");e.exe=binary;e.binaryHash=FileHash(binary);
    e.controls=controls;e.depthSharing=true;e.gpus={""};e.seeds.assign(production,production+8);e.workerThreads=1;e.hostBudgetMb=512;
    e.phi=4;e.target=8;e.theta={170,175,180};e.initial=1024;e.orientationBlock=512;e.workerPolicy="process";e.resourcePolicy="off";
    e.base={"--method","po","--backend","cpu","--particle","1","2","1","--ri","1.3116",absorbing?".001":"0","--wavelength-um",".532","--symmetry","1","1","--cutoff-profile","off"};e.selectedCutoff="off";
    e.coherent=!incoherent;if(incoherent){e.base.push_back("--incoherent");e.base[6]="501.2";e.base[7]="155.8";e.phi=8;e.theta={163.75909223159732,167.13352685050799,170};}
    WriteGrid(root+"/requested_theta.csv",e.theta);
}
// A normal interruption cannot leave a durable journal record without its
// marker: publishing the marker precedes appending/fsyncing the record. Remove
// both in this fixture to simulate stopping before either commit; leaving only
// the durable record is deliberate corruption and is tested separately.
void DropUncommittedJournalJobs(const std::string& root,const std::set<std::string>& folders){
    const auto path=root+"/completed_jobs.journal";if(!Exists(path))return;
    std::istringstream in(Read(path));std::ostringstream kept;std::string header;
    while(std::getline(in,header)){size_t length=0;std::string hash;std::istringstream h(header);Require(bool(h>>length>>hash),"invalid test journal");
        std::string body(length,'\0');in.read(&body[0],length);char newline=0;Require(size_t(in.gcount())==length&&bool(in.get(newline))&&newline=='\n'&&Hash(body)==hash,"damaged test journal");
        const auto name=body.substr(0,body.find('\n'));if(!folders.count(root+"/jobs/"+name))kept<<header<<'\n'<<body<<'\n';
    }
    Atomic(path,kept.str());
}
void Close(double a,double b,const char* message){Require(std::fabs(a-b)<=2e-12*std::max(1.,std::max(std::fabs(a),std::fabs(b))),message);}
int main(int argc,char** argv){try{
    Require(argc==3,"OUTPUT BINARY");ArgPP args;RunConfig config;
    for(int mode:{0,1,2}){const bool controls=mode==1,incoherent=mode==2;
        const std::string root=std::string(argv[1])+"/mode"+Number(mode);Engine e(args,config);Setup(e,root,argv[2],controls,controls,incoherent);e.WarmCache();
        if(incoherent){
            // Equal counts must telescope sample-by-sample to the direct
            // highest depth, including correlations across all eight seeds.
            e.target=12;e.VerifyDepthSharing(e.theta);Require(e.depthSharing,"incoherent sharing pilot rejected");
            const auto levels=EvaluateBudget(e,{8,10,12},Budgets(3,std::vector<int>(3,128)),{0,1,2},e.seeds.data());
            const auto combined=Combine(levels.levels),direct=e.Data(e.DirectStage(128,e.seeds.data(),8,12,e.theta));
            double error=0;for(size_t s=0;s<8;++s){error=std::max(error,Difference(combined[s],direct[s],direct[s]));
                for(size_t r=0;r<3;++r)for(int k=1;k<16;++k)error=std::max(error,std::fabs(combined[s][r][k]/combined[s][r][0]-direct[s][r][k]/direct[s][r][0]));}
            Require(error<1e-7,"incoherent multilevel estimate differs from direct highest depth");
            std::cout<<"incoherent actual column501 equal-count levels versus direct n12: "<<Number(error)<<'\n';e.target=8;
        }
        {AdaptiveState state;state.round=2;state.counts={{1048576}};Statistics noisy;noisy.maxInterval=1;
            const int initial=e.initial;e.initial=524288;e.statistical=.01;
            Require(WorkingPhiPilotDue(e,state,noisy),"severely noisy rule waits eight initial blocks");
            noisy.maxInterval=.05;Require(!WorkingPhiPilotDue(e,state,noisy),"ordinary rule rechecked too early");
            state.round=1;noisy.maxInterval=1;Require(!WorkingPhiPilotDue(e,state,noisy),"one-round pilot triggered");e.initial=initial;}
        auto old=e.DirectStage(512,e.seeds.data(),8,8,e.theta);const auto oldHash=FileHash(old[0].folder+"/complete.txt");
        e.balancedCorrection=true;auto reused=e.DirectStage(512,e.seeds.data(),8,8,e.theta);Require(!reused[0].assembled&&FileHash(old[0].folder+"/complete.txt")==oldHash,"old off blocks rewritten");
        // New transactions need two native legs; the balanced prefix is captured in B_N.
        e.collecting=true;e.DirectStage(512,e.seeds.data(),8,8,e.theta,true,false,512,6);e.collecting=false;
        std::vector<size_t> todo;for(size_t i=0;i<e.pendingJobs.size();++i)todo.push_back(i);e.ExecutePrepared(e.pendingJobs,e.pendingCommands,todo);e.pendingJobs.clear();e.pendingCommands.clear();
        auto corrected=e.DirectStage(512,e.seeds.data(),8,8,e.theta,true,false,512,6);
        for(const auto& j:corrected){Require(j.assembled&&j.pairedDepth==6,"paired correction missing");IntegralReader reader;
            Require(Read(j.folder+"/cutoff_prefix_recipe.txt")=="embedded-prefix-v1\n","new correction did not save its balanced prefix");
            std::istringstream signedRows(Read(j.folder+"/cutoff_sources.tsv"));int sourceDepth,sourceSign;std::string sourceFolder;int prefixViews=0,nativeLegs=0;
            while(signedRows>>sourceDepth>>sourceSign>>sourceFolder)if(sourceDepth==8){const Job source=e.Materialize(*e.folderIndex.at(sourceFolder));
                if(sourceSign==-1){Require(source.assembled&&source.count==64,"balanced prefix was traced twice");++prefixViews;
                    const auto replay=Read(source.folder+"/command.args");Require(replay.find("--cutoff-profile\nbalanced\n")!=std::string::npos&&replay.find("--orientation-prefix-output")==std::string::npos,"prefix standalone replay changed its estimator");}
                else {Require(!source.assembled,"unexpected native correction source");++nativeLegs;}}
            Require(prefixViews==1&&nativeLegs==2,"correction lost its two-leg transaction");
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
        const auto before=FileHash(corrected[0].folder+"/complete.txt");Engine restored(args,config);Setup(restored,root,argv[2],controls,controls,incoherent);restored.balancedCorrection=true;restored.WarmCache();restored.forbidNativeWork=true;
        auto replay=restored.DirectStage(512,e.seeds.data(),8,8,e.theta,true,false,512,6);Require(FileHash(replay[0].folder+"/complete.txt")==before,"completed correction rerun");
        auto mixed=restored.PrefixStage(1024,e.seeds.data(),8,8,e.theta);Require(mixed[0].matrix==prefix[0].matrix,"resume changed mixed estimator");
        auto isolated=e.DirectStage(512,e.seeds.data(),1,8,e.theta,true,false,2048);
        const auto folder=isolated[0].folder;std::remove((folder+"/complete.txt").c_str());
        if(e.compactIndex){Engine missing(args,config);Setup(missing,root,argv[2],controls,controls,incoherent);missing.balancedCorrection=true;missing.WarmCache();missing.forbidNativeWork=true;
            bool rejected=false;try{missing.DirectStage(512,e.seeds.data(),1,8,e.theta,true,false,2048);}catch(const std::exception&){rejected=true;}
            Require(rejected,"missing durably indexed marker was accepted");}
        DropUncommittedJournalJobs(root,{folder});
        Engine interrupted(args,config);Setup(interrupted,root,argv[2],controls,controls,incoherent);interrupted.balancedCorrection=true;interrupted.WarmCache();interrupted.forbidNativeWork=true;
        auto recovered=interrupted.DirectStage(512,e.seeds.data(),1,8,e.theta,true,false,2048);
        Require(recovered[0].matrix==isolated[0].matrix,"interrupted assembly changed or reran children");
        // Rebuild a lost prefix view and parent using only its durable native
        // products. The prefix source is hashed even though no worker ran for it.
        std::istringstream prefixSources(Read(folder+"/cutoff_sources.tsv"));int pd,ps;std::string pf,prefixFolder,balancedFolder;
        while(prefixSources>>pd>>ps>>pf)if(pd==8){if(ps==-1)prefixFolder=pf;else if(e.folderIndex.at(pf)->count==512)balancedFolder=pf;}
        Require(!prefixFolder.empty()&&!balancedFolder.empty(),"missing embedded prefix provenance");
        const auto nativePrefix=balancedFolder+"/result/result_prefix64.dat",prefixOriginal=Read(nativePrefix);
        Atomic(nativePrefix,prefixOriginal+"corruption\n");bool prefixRejected=false;
        try{std::set<std::string> seen;ValidateStoredJob(folder,seen);}catch(const std::exception&){prefixRejected=true;}
        Require(prefixRejected,"changed embedded prefix was accepted");Atomic(nativePrefix,prefixOriginal);
        std::remove((folder+"/complete.txt").c_str());std::remove((prefixFolder+"/complete.txt").c_str());
        DropUncommittedJournalJobs(root,{folder,prefixFolder});
        Require(std::rename((folder+"/result").c_str(),(folder+"/saved_before_prefix_recovery").c_str())==0,"cannot preserve test assembly");
        Engine prefixRecovered(args,config);Setup(prefixRecovered,root,argv[2],controls,controls,incoherent);prefixRecovered.balancedCorrection=true;prefixRecovered.WarmCache();prefixRecovered.forbidNativeWork=true;
        auto rebuilt=prefixRecovered.DirectStage(512,e.seeds.data(),1,8,e.theta,true,false,2048);
        Require(rebuilt[0].matrix==isolated[0].matrix,"lost prefix view caused native replay or changed matrix");
        const auto provenance=Read(folder+"/cutoff_sources.tsv");Atomic(folder+"/cutoff_sources.tsv",provenance+"corruption\n");bool rejected=false;
        try{std::set<std::string> seen;ValidateStoredJob(folder,seen);}catch(const std::exception&){rejected=true;}
        Require(rejected,"changed signed estimator was accepted");Atomic(folder+"/cutoff_sources.tsv",provenance);
        // An interrupted pre-prefix transaction keeps its original five
        // native legs; never silently reinterpret or discard existing work.
        e.collecting=true;e.DirectStage(512,e.seeds.data(),1,8,e.theta,true,false,4096,6);e.collecting=false;
        Require(e.pendingJobs.size()==1,"legacy test queued unexpected jobs");
        const auto legacyIdentity=Hash(e.pendingJobs[0].identity+"\ncutoff-child-v1 512 8 0 1");Directory(root+"/jobs/c"+legacyIdentity);
        e.ExecutePrepared(e.pendingJobs,e.pendingCommands,{0});const Job legacy=e.pendingJobs[0];
        Require(!Exists(legacy.folder+"/cutoff_prefix_recipe.txt"),"interrupted legacy correction changed recipes");
        std::istringstream legacySources(Read(legacy.folder+"/cutoff_sources.tsv"));int ld,ls;std::string lf;int independentBalanced=0;
        while(legacySources>>ld>>ls>>lf){const Job part=e.Materialize(*e.folderIndex.at(lf));if(!part.assembled&&part.pairedDepth==0)++independentBalanced;}
        Require(independentBalanced==4,"legacy independent-depth children were not preserved");e.pendingJobs.clear();e.pendingCommands.clear();


    }
    std::cout<<"PASS coherent/incoherent corrected paired matrices, controls, absorption, integrals, mixed prefixes, resume\n";return 0;
}catch(const std::exception& ex){std::cerr<<ex.what()<<'\n';return 1;}}
