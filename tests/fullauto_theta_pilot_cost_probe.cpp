#define main old_lazy_probe_main
#include "fullauto_lazy_restore_probe.cpp"
#undef main
#include <chrono>
int main(int argc,char** argv){try{
 Require(argc==3,"OUTPUT BINARY");const std::string out=argv[1];Directory(out);
 ArgPP args;RunConfig cfg;std::vector<size_t> knots[2];double elapsed[2];
 Engine baseline(args,cfg);SetupMemory(baseline,out+"/n256",argv[2],false);baseline.phi=4;baseline.target=6;baseline.eps=.04;baseline.interpolation=.01;
 baseline.base={"--method","po","--backend","cpu","--particle","2","5","2.5","--ri","1.3116","0","--wavelength-um",".532","--symmetry","1","1","--cutoff-profile","off"};
 baseline.theta.clear();for(int k=0;k<=128;++k)baseline.theta.push_back(k*180./128);baseline.thetaAuditCount=256;
 Engine reduced(args,cfg);SetupMemory(reduced,out+"/n32",argv[2],false);reduced.phi=baseline.phi;reduced.target=baseline.target;reduced.eps=baseline.eps;reduced.interpolation=baseline.interpolation;reduced.theta=baseline.theta;reduced.base=baseline.base;reduced.thetaAuditCount=32;
 for(int k=0;k<2;++k){Engine& e=k?reduced:baseline;std::vector<size_t> initial={0,16,32,48,64,80,96,112,128};const double start=Seconds();knots[k]=e.AuditThetaSupport(initial);elapsed[k]=Seconds()-start;}
 auto data=baseline.Data(baseline.DirectStage(256,angularValidation,8,baseline.target,baseline.theta));const auto ref=Mean(data);
 double error[2]={};
 for(int k=0;k<2;++k){std::vector<double> x;Samples support(8),predicted(8);for(auto row:knots[k]){x.push_back(baseline.theta[row]);for(int seed=0;seed<8;++seed)support[seed].push_back(data[seed][row]);}
  for(int seed=0;seed<8;++seed)predicted[seed]=CubicInterpolate(x,support[seed],baseline.theta);
  const auto delta=Subtract(predicted,data);const auto mean=Mean(delta),variance=Variance(delta);
  for(size_t row=0;row<ref.size();++row)for(int j=0;j<16;++j)error[k]=std::max(error[k],(std::fabs(mean[row][j])+t95*std::sqrt(variance[row][j]/8))/ref[row][0]);
 }
 Require(error[1]<=baseline.interpolation*.5,"reduced pilot failed same reference all16 interpolation budget");
 std::ostringstream receipt;receipt<<"{\"baseline_seconds\":"<<elapsed[0]<<",\"reduced_seconds\":"<<elapsed[1]<<",\"speedup\":"<<elapsed[0]/elapsed[1]<<",\"baseline_knots\":"<<knots[0].size()<<",\"reduced_knots\":"<<knots[1].size()<<",\"reference_bound_baseline\":"<<error[0]<<",\"reference_bound_reduced\":"<<error[1]<<",\"budget\":"<<baseline.interpolation*.5<<"}";
 Atomic(out+"/receipt.json",receipt.str());std::cout<<receipt.str()<<'\n';
}catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}}
