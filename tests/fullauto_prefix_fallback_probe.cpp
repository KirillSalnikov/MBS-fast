#include "../src/FullAuto.cpp"
using namespace FullAuto;

void SetupFallback(Engine& e,const std::string& root,const std::string& binary,bool gpu,bool controls){
    e.root=root;e.signature="prefix-capability-fallback-test-v1";Directory(root);Directory(root+"/jobs");
    e.exe=binary;e.binaryHash=FileHash(binary);e.compactIndex=true;e.controls=controls;e.sharedMeans=controls;
    e.depthSharing=true;e.gpus={gpu?"0":""};e.seeds.assign(production,production+8);e.workerThreads=2;e.hostBudgetMb=1024;
    e.phi=4;e.target=8;e.theta={175.9,176};e.initial=512;e.orientationBlock=512;e.workerPolicy="persistent";e.resourcePolicy="off";
    e.kernel.pipeline=gpu;e.balancedCorrection=true;e.selectedCutoff="off";
    e.base={"--method","po","--backend",gpu?"cuda":"cpu","--particle","1","2","1","--ri","1.3116",controls?".001":"0",
        "--wavelength-um",".532","--symmetry","1","1","--cutoff-profile","off"};
    WriteGrid(root+"/requested_theta.csv",e.theta);
}
void TestRecipeState(Engine& e){
    Job parent;parent.identity="recipe-fixture";parent.count=512;parent.depth=8;
    parent.folder=e.root+"/fixture_parent";Directory(parent.folder);Atomic(parent.folder+"/cutoff_prefix_recipe.txt","embedded-prefix-v1\n");
    const auto sourceId=Hash(parent.identity+"\ncutoff-child-v1 512 8 0 1 prefix-v1 64"),source=e.root+"/jobs/c"+sourceId;
    Directory(source+"/result");Atomic(source+"/result/incomplete.txt","preserve partial source\n");
    Require(!e.RestoreCutoffPrefixRecipe(parent,64),"missing native capability retained incomplete embedded recipe");
    Require(Read(parent.folder+"/cutoff_prefix_recipe.txt")=="three-leg-v1\n","fallback recipe not committed");
    Require(Read(parent.folder+"/cutoff_prefix_recipe.before_fallback.txt")=="embedded-prefix-v1\n","old recipe lost");
    Require(Read(source+"/result/incomplete.txt")=="preserve partial source\n","partial worker data lost");
    const auto receipt=Read(parent.folder+"/cutoff_prefix_fallback.json");
    Require(!e.RestoreCutoffPrefixRecipe(parent,64)&&Read(parent.folder+"/cutoff_prefix_fallback.json")==receipt,"fallback resume changed receipt");

    // A durable, valid embedded native source is readable even when the current
    // executable cannot produce another one. This is a metadata fixture only.
    parent.identity="completed-recipe-fixture";parent.folder=e.root+"/fixture_completed";Directory(parent.folder);
    Atomic(parent.folder+"/cutoff_prefix_recipe.txt","embedded-prefix-v1\n");
    const auto completedId=Hash(parent.identity+"\ncutoff-child-v1 512 8 0 1 prefix-v1 64"),completed=e.root+"/jobs/c"+completedId;
    Directory(completed+"/result");Atomic(completed+"/result/result.dat","validated fixture\n");
    Atomic(completed+"/result/result_prefix64.dat","validated prefix fixture\n");
    Atomic(completed+"/result/result_prefix64_integrals.tsv","validated scalar fixture\n");
    Atomic(completed+"/command.args",e.exe+"\n--orientation-prefix-output\n64\n");
    std::ostringstream manifest;for(const char* suffix:{".dat","_integrals.tsv"}){
        const auto path=completed+"/result/result_prefix64"+suffix;manifest<<path<<' '<<FileHash(path)<<'\n';}
    Atomic(completed+"/integral_outputs.fnv",manifest.str());
    Atomic(completed+"/complete.txt",completedId+" 1 "+FileHash(completed+"/result/result.dat")+" none none "+Hash(manifest.str())+"\n");
    const auto before=FileHash(completed+"/complete.txt");
    Require(e.RestoreCutoffPrefixRecipe(parent,64),"valid completed embedded source was discarded");
    Require(FileHash(completed+"/complete.txt")==before&&Read(parent.folder+"/cutoff_prefix_recipe.txt")=="embedded-prefix-v1\n","completed source rewritten");
    Atomic(completed+"/result/result_prefix64.dat","corrupted prefix\n");bool rejected=false;
    try{e.RestoreCutoffPrefixRecipe(parent,64);}catch(const std::exception&){rejected=true;}
    Require(rejected,"corrupted prefix silently migrated or accepted");
}
int main(int argc,char** argv){try{
    Require(argc==3||argc==4,"OUTPUT BINARY [cpu|cuda]");
    Require(!NativeWorkerCapabilities::OrientationPrefixOutput(),"fallback probe must use default absent native capability");
    const bool gpu=argc==4&&std::string(argv[3])!="cpu";ArgPP args;RunConfig config;config.useGpu=gpu;
    {Engine e(args,config);SetupFallback(e,std::string(argv[1])+"/recipe",argv[2],gpu,false);TestRecipeState(e);
        Job job;job.count=65536;job.inputAngles={175.9,176};job.dispatchChunk=512;
        std::vector<std::string> command={"--sobol-seed","65536","11","--orientation-prefix-output","8192"};const auto original=command;
        Require(!e.VerifyCutoffPrefix(job,command),"unsupported native prefix reached numerical pilot");
        Require(!e.VerifyShortPacketChunk(job,command),"short packet pilot accepted unbounded prefix");
        e.ApplyShortPacketChunk(job,command);Require(command==original&&job.dispatchChunk==0,"short packet pilot rewrote embedded prefix recipe");
    }
    for(bool controls:{false,true})for(int paired:{0,6}){
        const auto root=std::string(argv[1])+"/controls"+Number(controls)+"_paired"+Number(paired);Engine e(args,config);SetupFallback(e,root,argv[2],gpu,controls);
        e.collecting=true;e.DirectStage(512,e.seeds.data(),1,8,e.theta,true,false,0,paired);e.collecting=false;
        Require(e.pendingJobs.size()==1,"unexpected queued estimator count");const auto parentFolder=e.pendingJobs[0].folder;
        Directory(parentFolder);if(paired)Atomic(parentFolder+"/cutoff_prefix_recipe.txt","embedded-prefix-v1\n");
        e.ExecutePrepared(e.pendingJobs,e.pendingCommands,{0});const Job output=e.pendingJobs[0];
        Require(output.assembled&&Exists(parentFolder+"/complete.txt"),"fallback did not finish corrected estimator");
        Require(paired?Read(parentFolder+"/cutoff_prefix_recipe.txt")=="three-leg-v1\n":!Exists(parentFolder+"/cutoff_prefix_recipe.txt"),"fresh or interrupted recipe selected unavailable prefix");
        for(int depth:{8,paired})if(depth){Matrix expected(e.theta.size(),Mueller{});std::vector<std::array<double,5>> expectedControls(e.theta.size(),std::array<double,5>{});size_t terms=0;IntegralReader reader;IntegralScalar scalar;
            std::istringstream sources(Read(parentFolder+"/cutoff_sources.tsv"));int d,sign;std::string folder;
            while(sources>>d>>sign>>folder)if(d==depth){++terms;const auto source=e.Materialize(*e.folderIndex.at(folder));
                Require(Read(folder+"/command.args").find("--orientation-prefix-output")==std::string::npos,"fallback emitted unsupported worker option");
                for(size_t row=0;row<e.theta.size();++row)for(int k=0;k<16;++k)expected[row][k]+=sign*source.matrix[row][k];
                if(controls)for(size_t row=0;row<e.theta.size();++row)for(int k=0;k<5;++k)expectedControls[row][k]+=sign*source.controls[row][k];
                const auto v=reader.Resolve(folder,depth);scalar.ext+=sign*v.ext;scalar.area+=sign*v.area;}
            const auto target=depth==8?output:e.LowerView(output,depth);Require(terms==3&&target.matrix==expected,"fallback changed signed all16 estimator");
            if(controls)Require(target.controls==expectedControls,"fallback changed signed analytic control components");
            const auto actual=reader.Resolve(target.folder,depth);
            Require(std::fabs(actual.ext-scalar.ext)<1e-10*std::max(1.,std::fabs(scalar.ext))&&std::fabs(actual.area-scalar.area)<1e-10*std::max(1.,std::fabs(scalar.area)),"fallback changed signed integral characteristics");
        }
        const auto marker=FileHash(parentFolder+"/complete.txt");for(auto& service:e.services)service.reset();
        Engine resumed(args,config);SetupFallback(resumed,root,argv[2],gpu,controls);resumed.WarmCache();resumed.forbidNativeWork=true;
        const auto replay=resumed.DirectStage(512,resumed.seeds.data(),1,8,resumed.theta,true,false,0,paired);
        Require(replay[0].matrix==output.matrix&&FileHash(parentFolder+"/complete.txt")==marker,"completed fallback reran or changed on resume");
        std::cout<<"PASS three-leg native fallback controls="<<controls<<" paired="<<paired<<", all16/integrals/immutable resume\n";
    }
    std::cout<<"PASS absent capability, durable recipe migration, completed source retention, corruption refusal, prefix/chunk guard\n";return 0;
}catch(const std::exception& error){std::cerr<<error.what()<<'\n';return 1;}}
