#include "../src/FullAuto.cpp"
using namespace FullAuto;
int main(int argc,char**argv){try{
    Require(argc==2,"existing lazy-regression fixture root");
    for(bool controls:{false,true}){
        ArgPP a;RunConfig c;Engine serial(a,c),parallel(a,c);
        for(auto* e:{&serial,&parallel}){e->root=std::string(argv[1])+(controls?"/controls":"/plain");e->signature="bounded-memory-regression-v1";e->controls=controls;e->forbidNativeWork=true;e->payloads.SetLimit(1);}
        serial.WarmCache(1);parallel.WarmCache(4);
        Require(!serial.cache.empty()&&serial.cache.size()==parallel.cache.size(),"index membership changed");
        Require(serial.samplingIndex.size()==parallel.samplingIndex.size()&&serial.folderIndex.size()==parallel.folderIndex.size(),"lookup membership changed");
        for(const auto& item:serial.cache){const Job& x=item.second;const Job& y=parallel.cache.at(item.first);
            Require(x.inputAngles==y.inputAngles&&x.count==y.count&&x.offset==y.offset&&x.seed==y.seed&&x.depth==y.depth&&x.phi==y.phi&&x.block==y.block&&x.executionHash==y.executionHash&&x.weightHash==y.weightHash&&x.warp==y.warp&&x.hasControls==y.hasControls&&x.assembled==y.assembled&&x.paired==y.paired&&x.pairedDepth==y.pairedDepth&&x.pipeline==y.pipeline&&x.mirrored==y.mirrored,"sampling metadata changed");
            Require(x.controlRefinement==y.controlRefinement&&x.controlSetup==y.controlSetup&&x.seconds==y.seconds&&x.traceSeconds==y.traceSeconds&&x.diffractionSeconds==y.diffractionSeconds&&x.preparationThreads==y.preparationThreads&&x.processingSeconds==y.processingSeconds&&x.preparedWaitSeconds==y.preparedWaitSeconds&&x.producerWaitSeconds==y.producerWaitSeconds&&x.assemblySeconds==y.assemblySeconds,"diagnostics changed");
            const auto first=serial.Materialize(x),second=parallel.Materialize(y);
            Require(first.matrix==second.matrix&&first.controls==second.controls&&first.angles==second.angles,"payload changed");
        }
        Require(parallel.payloads.Bytes()==0,"one-byte cache retained payload");
        const Job victim=serial.cache.begin()->second;const auto path=victim.folder+"/sampling.index";const auto saved=Read(path);Atomic(path,saved+"broken\n");
        bool rejected=false;try{Engine bad(a,c);bad.root=serial.root;bad.signature=serial.signature;bad.controls=controls;bad.WarmCache(4);}catch(const std::exception&){rejected=true;}
        Atomic(path,saved);Require(rejected,"parallel reader accepted changed index");
        std::cout<<"PASS "<<(controls?"controls":"plain")<<": serial/parallel metadata, diagnostics, all16 bitwise equal; corrupt-index exception propagated\n";
    }
    std::cout<<"PASS parallel restore\n";return 0;
}catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}}
