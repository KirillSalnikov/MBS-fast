#define main fullauto_cutoff_correction_main
#include "fullauto_cutoff_correction_probe.cpp"
#undef main

static void Configure(Engine& e,const std::string& root,const std::string& binary){
    Setup(e,root,binary,true,false);
    e.balancedCorrection=true;
    e.cutoffPolicy="staged";
    e.ConfigurePlanning(true);
    e.sharedMeans=true;
    e.theta={0,12.5,90,170,175,180};
    e.weightTheta={0,180};
    e.weights={{{.2,-.15}},{{-.1,.6}}};
    WriteGrid(e.MeanReferenceGrid(),e.theta);
}
int main(int argc,char** argv){
    try{
        Require(argc==3||argc==4,"expected output, native solver, and optionally legacy");
        const std::string schema=argc==4?argv[3]:"modern";
        Require(schema=="modern"||schema=="legacy","unknown table schema");
        const bool modern=schema=="modern";
        ArgPP args;RunConfig config;
        Engine wide(args,config),direct(args,config);
        Configure(wide,std::string(argv[1])+"/wide",argv[2]);
        Configure(direct,std::string(argv[1])+"/reference",argv[2]);
        wide.DirectStage(512,wide.seeds.data(),1,8,wide.theta,true,false,0,6).at(0);
        const std::vector<double> subset={0,180};
        wide.forbidNativeWork=modern;
        const auto part=wide.DirectStage(512,wide.seeds.data(),1,8,subset,true,false,0,6).at(0);
        const auto fresh=direct.DirectStage(512,direct.seeds.data(),1,8,subset,true,false,0,6).at(0);
        int projections=0;
        std::istringstream sources(Read(part.folder+"/cutoff_sources.tsv"));
        int depth,sign;std::string folder;
        while(sources>>depth>>sign>>folder){
            if(depth!=8||!Exists(folder+"/angular_sources.tsv"))continue;
            if(Read(folder+"/stdout.log").find("Exact direction subset")==std::string::npos)continue;
            const auto used=wide.Materialize(*wide.folderIndex.at(folder));
            const auto parentFolder=Read(folder+"/angular_sources.tsv");
            const auto parent=wide.Materialize(*wide.folderIndex.at(parentFolder.substr(0,parentFolder.find('\t'))));
            Require(used.inputAngles==subset&&parent.inputAngles==wide.theta,"wrong projection source");
            Require(used.rowRefinement.size()==2&&used.controlRefinement<parent.controlRefinement,"projection did not narrow analytic bound");
            ++projections;
        }
        Require(modern?projections==3:projections==0,"wrong policy for native table schema");
        IntegralReader ir;
        for(int n:{8,6}){
            const Job a=n==8?part:wide.LowerView(part,n),b=n==8?fresh:direct.LowerView(fresh,n);
            Require(a.angles==b.angles&&a.angles==subset,"wrong output grid");
            for(size_t r=0;r<subset.size();++r){
                for(size_t k=0;k<16;++k)Close(a.matrix.at(r)[k],b.matrix.at(r)[k],"matrix differs");
                for(size_t k=0;k<5;++k)Close(a.controls.at(r)[k],b.controls.at(r)[k],"analytic primitive differs");
            }
            Close(a.controlRefinement,b.controlRefinement,"wrong analytic error");
            const auto av=ir.Resolve(a.folder,n),bv=ir.Resolve(b.folder,n);
            Close(av.ext,bv.ext,"extinction differs");Close(av.area,bv.area,"area differs");
        }
        Engine restored(args,config);Configure(restored,wide.root,argv[2]);restored.WarmCache();
        restored.forbidNativeWork=true;
        const auto warm=restored.DirectStage(512,restored.seeds.data(),1,8,subset,true,false,0,6).at(0);
        Require(warm.matrix==part.matrix&&warm.controls==part.controls,"resume changed projection");
        std::cout<<"PASS "<<schema<<" verified native projections="<<projections<<"\n";
    }catch(const std::exception& ex){std::cerr<<ex.what()<<'\n';return 1;}
}
