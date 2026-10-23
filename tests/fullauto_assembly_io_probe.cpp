#define main original_lazy_main
#include "fullauto_lazy_restore_probe.cpp"
#undef main

void InterruptedFiles(const std::string& root){
    Directory(root+"/result");
    const pid_t child=fork();Require(child>=0,"cannot fork interruption fixture");
    if(!child){AtomicFileBatch files(root);Atomic(root+"/result/result.dat","saved matrix\n");
        Atomic(root+"/angular_sources.tsv","durable parent\n");_exit(86);}
    int status;Require(waitpid(child,&status,0)==child&&WIFEXITED(status)&&WEXITSTATUS(status)==86,"interruption fixture failed");
    Require(!Exists(root+"/complete.txt"),"unfinished files were committed");
    {
        AtomicFileBatch files(root);
        Atomic(root+"/result/result.dat","saved matrix\n");Atomic(root+"/angular_sources.tsv","durable parent\n");
        Require(files.files.size()==2,"matching interrupted files escaped durability fence");
        bool rejected=false;try{Atomic(root+"/complete.txt","premature\n");}catch(const std::exception&){rejected=true;}
        Require(rejected&&!Exists(root+"/complete.txt"),"completion published before data flush");
        files.Commit();Require(!activeAtomicFileBatch,"committed transaction remains active");
        Atomic(root+"/complete.txt","complete\n");
    }
    bool rejected=false;try{AtomicFileBatch files(root);}catch(const std::exception&){rejected=true;}
    Require(rejected&&!activeAtomicFileBatch,"completed job entered a rewrite transaction");
    Require(Read(root+"/result/result.dat")=="saved matrix\n","interruption recovery lost data");
}

void CompactGraphMemo(const std::string& root){
    Directory(root+"/leaf/result");Atomic(root+"/leaf/result/result.dat","leaf\n");
    Atomic(root+"/leaf/complete.txt","leaf 1 "+FileHash(root+"/leaf/result/result.dat")+" none\n");
    std::vector<std::string> roots;
    for(int i=0;i<8;++i){const auto p=root+"/graph"+Number(i);Directory(p+"/result");
        Atomic(p+"/result/result.dat","derived\n");const auto sources=root+"/leaf\t"+FileHash(root+"/leaf/complete.txt")+"\n";
        Atomic(p+"/angular_sources.tsv",sources);
        Atomic(p+"/complete.txt",Number(i)+" 1 "+FileHash(p+"/result/result.dat")+" none "+Hash(sources)+"\n"+std::string(60000,'x'));
        roots.push_back(p);}
    ValidationBatch batch("large compact marker fixture",false,nullptr,8192);
    for(const auto& p:roots){std::set<std::string> seen;ValidateStoredJob(p,seen);}
    const auto before=batch.reused;
    for(const auto& p:roots){std::set<std::string> seen;ValidateStoredJob(p,seen);}
    Require(batch.reused-before==8&&batch.compounds.Bytes()<8192,"sampling grids evicted verified graph memos");
}

void ManyStencils(const std::string& root,const std::string& binary,bool controls){
    ArgPP args;RunConfig config;Engine reference(args,config),batch(args,config);
    SetupMemory(reference,root+"/reference",binary,controls);SetupMemory(batch,root+"/batch",binary,controls);
    reference.theta.clear();for(int i=0;i<48;++i)reference.theta.push_back(170+i*.2);batch.theta=reference.theta;
    reference.orientationBlock=batch.orientationBlock=32;
    for(auto* e:{&reference,&batch}){e->compactIndex=true;WriteGrid(e->root+"/requested_theta.csv",e->theta);
        if(controls){e->weightTheta=e->theta;e->weights.resize(e->theta.size(),{{.4,.3}});}}
    for(int offset=0;offset<256;offset+=32){const auto jobs=reference.DirectStage(32,reference.seeds.data(),8,6,reference.theta,true,false,offset);
        // Importing a source from outside the local root does not append it to
        // this fixture's journal; the real refinement importer maps its parent.
        batch.restoringIndex=true;for(const auto& j:jobs)batch.CacheJob(j);batch.restoringIndex=false;}
    reference.forbidNativeWork=batch.forbidNativeWork=true;
    // Fit one eight-seed block, but not all eight orientation blocks.
    reference.payloads.Clear();batch.payloads.Clear();reference.payloads.SetLimit(160*1024);batch.payloads.SetLimit(160*1024);
    std::vector<Engine::WorkRequest> requests;std::vector<std::vector<Job>> expected;
    for(int i=0;i<24;++i){std::vector<double> x={reference.theta[2*i],reference.theta[2*i+1]};
        requests.push_back({256,0,6,0,x});expected.push_back(reference.PrefixStage(256,reference.seeds.data(),8,6,x));}
    batch.SetPrefixAssemblyPlan(requests,batch.seeds.data());const auto misses=batch.payloads.misses;
    for(size_t i=0;i<requests.size();++i){const auto got=batch.PrefixStage(256,batch.seeds.data(),8,6,requests[i].angles);
        if(!i)Require(batch.preparedPrefixes.size()==23,"narrow batches still capped at eight");
        for(size_t s=0;s<8;++s){const auto& a=expected[i][s];const auto& b=got[s];
            Require(a.matrix==b.matrix&&a.controls==b.controls&&a.controlRefinement==b.controlRefinement&&a.controlSetup==b.controlSetup,"batched output changed estimator");
            Require(Read(a.folder+"/angular_sources.tsv")==Read(b.folder+"/angular_sources.tsv"),"batched output changed sources");
            std::set<std::string> checked;ValidateStoredJob(b.folder,checked);
            IntegralReader reader;const auto ia=reader.Resolve(a.folder,6),ib=reader.Resolve(b.folder,6);
            Require(ia.area==ib.area&&ia.ext==ib.ext,"batched output changed integral characteristics");}
    }
    Require(batch.preparedPrefixes.empty()&&batch.preparedPrefixBytes==0,"batch memory retained after use");
    Require(batch.payloads.misses-misses<=64,"wide blocks reread for narrow groups");
    std::cout<<"PASS controls="<<controls<<": 24 stencils, all16/controls/integrals/source order exact; payload loads="<<batch.payloads.misses-misses<<'\n';
}
int main(int argc,char** argv){try{
    Require(argc==3,"OUTPUT BINARY");const std::string root=argv[1];
    InterruptedFiles(root+"/interrupted");CompactGraphMemo(root+"/graphs");
    for(bool controls:{false,true})ManyStencils(root+"/controls"+Number(controls),argv[2],controls);
    std::cout<<"PASS interrupted derived writes, commit ordering, compact validation memos, bounded wider prefix batching\n";return 0;
}catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}}
