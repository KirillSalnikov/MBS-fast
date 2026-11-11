#include "../src/FullAuto.cpp"
using namespace FullAuto;
int main(int argc,char** argv){try{
 Require(argc==2,"fixture directory required");Directory(argv[1]);
 auto policy=[](std::vector<std::string> flags,bool expected){
  flags.insert(flags.begin(),"probe");std::vector<const char*> av;for(const auto& f:flags)av.push_back(f.c_str());
  ArgPP args;RegisterCliOptions(args);args.Parse(int(av.size()),av.data());RunConfig c;Engine e(args,c);
  Require(e.ConfiguredFinalChecks()==expected,"audit opt-in/master priority differs");
  e.finalChecks=expected;Require(e.RequiredProductionPasses()==1,"implicit repeat enabled");
 };
 policy({},false);policy({"--fullauto-final-checks","off"},false);policy({"--fullauto-final-checks","on"},true);
 policy({"--fullauto-final-audit","on"},true);policy({"--fullauto-physics-audit","on"},true);
 policy({"--fullauto-final-checks","off","--fullauto-final-audit","on"},false);
 policy({"--fullauto-final-checks","off","--fullauto-physics-audit","on"},false);
 policy({"--fullauto-final-audit","off"},false);
 const std::string base=argv[1];Atomic(base+"/implicit.conf","theta.max_points = 4097\n");
 Atomic(base+"/explicit.conf","controller.stable_passes = 3\n");
 ArgPP args;RunConfig cfg;cfg.adaptive=LoadAdaptiveConvergenceConfig(base+"/implicit.conf");Engine a(args,cfg);
 Require(a.RequiredProductionPasses()==1,"unrelated config retained implicit 2 passes");
 cfg.adaptive=LoadAdaptiveConvergenceConfig(base+"/explicit.conf");Require(a.RequiredProductionPasses()==3,"explicit config passes ignored");
 ArgPP cli;RegisterCliOptions(cli);const char* av[]={"probe","--convergence-passes","4"};cli.Parse(3,av);cfg.adaptive.stablePasses=4;Engine b(cli,cfg);Require(b.RequiredProductionPasses()==4,"explicit CLI passes ignored");
 std::cout<<"PASS fullauto defaults: opt-in audits, single pass with or without audit, explicit master/child priorities, explicit versus unrelated adaptive config\n";
}catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}}
