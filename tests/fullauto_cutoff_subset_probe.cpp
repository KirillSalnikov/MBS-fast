#define main fullauto_cutoff_correction_main
#include "fullauto_cutoff_correction_probe.cpp"
#undef main

int main(int argc,char** argv){
    try{
        Require(argc==3,"expected temporary output and native solver");
        ArgPP args;
        RunConfig config;
        Engine engine(args,config);
        Setup(engine,argv[1],argv[2],false,false);
        engine.balancedCorrection=true;
        engine.cutoffPolicy="staged";
        engine.ConfigurePlanning(true);
        const std::vector<double> subset={175,180};
        const auto full=engine.DirectStage(512,engine.seeds.data(),1,8,engine.theta).front();
        const auto part=engine.DirectStage(512,engine.seeds.data(),1,8,subset).front();
        Require(part.angles==subset,"subset output grid changed");
        for(size_t i=0;i<subset.size();++i)
            for(size_t k=0;k<16;++k)
                Close(part.matrix.at(i)[k],full.matrix.at(i+1)[k],"projected Mueller component changed");
        IntegralReader integrals;
        const auto all=integrals.Resolve(full.folder,8),small=integrals.Resolve(part.folder,8);
        Close(small.ext,all.ext,"projected extinction changed");
        Close(small.area,all.area,"projected area changed");
        std::istringstream sources(Read(part.folder+"/cutoff_sources.tsv"));
        int depth,sign;
        std::string folder;
        bool projected=false;
        while(sources>>depth>>sign>>folder){
            const auto& source=*engine.folderIndex.at(folder);
            Require(source.inputAngles==subset,"signed correction points to a different grid");
            if(source.assembled&&Exists(folder+"/angular_sources.tsv"))projected=true;
        }
        Require(projected,"missing source projection");
        Engine restored(args,config);
        Setup(restored,argv[1],argv[2],false,false);
        restored.balancedCorrection=true;
        restored.cutoffPolicy="staged";
        restored.ConfigurePlanning(true);
        restored.WarmCache();
        restored.forbidNativeWork=true;
        const auto warm=restored.DirectStage(512,restored.seeds.data(),1,8,subset).front();
        Require(warm.matrix==part.matrix,"restored projection differs");
        std::cout<<"PASS native direction subset, signed sources, integrals and warm restart\n";
        return 0;
    }catch(const std::exception& ex){std::cerr<<ex.what()<<'\n';return 1;}
}
