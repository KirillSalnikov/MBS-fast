// Check prefetch placement when a two-row request has one uncached angle.
#include "../src/FullAuto.cpp"
int main(){
    using namespace FullAuto;
    try{
        ArgPP args;RunConfig config;config.useGpu=true;Engine e(args,config);
        e.nonconvexGeometry=true;e.controls=false;e.coherent=true;e.phi=8;e.seeds.assign(production,production+8);
        setenv("MBS_FULLAUTO_PRODUCTIVE_PADDING","1",1);e.ConfigurePlanning(true);
        Engine::WorkRequest request{512,0,12,0,{0,179.5}};
        auto& demand=e.paddingDemand[{12,8}];demand.counts[179.5]=256;
        demand.angles[256]={100,179.375,179.5,179.625};
        std::vector<std::set<double>> covered(8,std::set<double>{0,179.375});
        auto choose=[&]{return e.PrefixPadding(request,covered,179.5,0,e.seeds.data());};
        Require(choose()==179.625,"compatible uncovered direction missed");
        covered[3].insert(179.5);Require(choose()==0,"partly complete exact request changed");covered[3].erase(179.5);
        covered[0].erase(0);Require(choose()==0,"missing reference row");covered[0].insert(0);
        covered[7].insert(179.625);Require(choose()==0,"partial secondary coverage or distant row selected");covered[7].erase(179.625);
        e.controls=true;Require(choose()==0,"analytic weights require original grid");e.controls=false;
        request.pairedDepth=10;Require(choose()==0,"paired source changed");request.pairedDepth=0;
        Require(e.PrefixPadding(request,covered,179.5,0,training)==0,"pilot row selected");
        e.phi=4;Require(choose()==0,"azimuth mismatch");e.phi=8;
        e.nonconvexGeometry=false;e.ConfigurePlanning(true);Require(choose()==0,"convex schedule changed");
        unsetenv("MBS_FULLAUTO_PRODUCTIVE_PADDING");e.nonconvexGeometry=true;e.ConfigurePlanning(true);
        Require(choose()==0,"option unexpectedly active");
        std::cout<<"PASS fullauto productive singleton padding\n";return 0;
    }catch(const std::exception& ex){std::cerr<<ex.what()<<'\n';return 1;}
}
