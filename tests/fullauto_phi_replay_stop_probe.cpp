#define main old_lazy_main
#include "fullauto_lazy_restore_probe.cpp"
#undef main

int main(int argc,char** argv){try{
    Require((argc==3||argc==4) && !Exists(argv[1]),"fresh output root, worker binary, optional cpu|cuda backend required");
    const std::string backend=argc==4?argv[3]:"cpu";
    Require(backend=="cpu"||backend=="cuda","unknown backend");
    const double tolerance=backend=="cuda"?1e-5:1e-10;
    ArgPP args;RegisterCliOptions(args);RunConfig config;Engine e(args,config);
    SetupMemory(e,argv[1],argv[2],false);
    e.pilot=16;e.phi=4;e.theta={0.,31.,170.,180.};
    e.base={"--method","po","--backend",backend,"--particle","3","0.2","0.4","--rosette-joined",
            "--ri","1.3194","0","--wavelength-um",".400","--symmetry","2","2","--cutoff-profile","off"};
    e.memoryMode=true;
    e.memoryStore.reset(new FullAutoMemoryStore::Store(e.root,8*1024*1024,e.root));
    e.memoryScope.reset(new FullAutoMemoryStore::Scope(*e.memoryStore));
    std::vector<int> rules={1,2,4};
    std::map<int,Job> reference;
    for(int rule:rules){
        e.phi=rule;
        reference[rule]=e.DirectStage(e.pilot,validation+1,1,4,e.theta,false)[0];
    }
    const auto base=e.base;const auto original=e.ExecutionHash();
    auto before=e.nativeSubmitted;
    std::map<int,std::vector<Job>> partial;
    Require(!e.FillTiledPhi(partial,validation,1,2,rules,4,e.theta,[&]{return e.nativeSubmitted>before;}),
            "a losing replay continued");
    Require(e.nativeSubmitted==before+1 && e.phi==4 && e.base==base && e.ExecutionHash()==original
            && e.kernelTrialNonce.empty(),"stopping changed the native recipe");
    Require(e.memoryStore->List(e.root+"/kernel_trials",true).empty(),"aborted private outputs retained");
    std::map<int,std::vector<Job>> full;
    Require(e.FillTiledPhi(full,validation,1,2,rules,4,e.theta),"unrestricted replay stopped");
    before=e.nativeSubmitted;
    for(int rule:rules){
        e.phi=rule;
        const Job cold=e.DirectStage(e.pilot,validation+1,1,4,e.theta,false)[0];
        Require(cold.matrix==reference.at(rule).matrix,"standalone validation matrix changed");
        const auto& tiled=full.at(rule)[1];
        for(size_t i=0;i<e.theta.size();++i){
            const double scale=cold.matrix[i][0];
            Require(scale>0 && tiled.matrix[i][0]>0,"invalid M11");
            for(size_t j=0;j<16;++j){
                Require(std::fabs(tiled.matrix[i][j]-cold.matrix[i][j])/scale<tolerance,"raw Mueller mismatch");
                Require(std::fabs(tiled.matrix[i][j]/tiled.matrix[i][0]-cold.matrix[i][j]/scale)<tolerance,"ratio mismatch");
            }
        }
    }
    Require(e.nativeSubmitted==before,"cold validation samples were recomputed");
    std::cout<<"PASS stopped after one private request, recovered cold samples, complete replay agrees on all 16 elements\n";
    return 0;
}catch(const std::exception& ex){std::cerr<<ex.what()<<'\n';return 1;}}
