#include "FullAuto.h"
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <stdexcept>
using namespace FullAuto;
int main(int argc,char** argv) {
    if(argc==3 && std::string(argv[1])=="--allocation") {
        std::ifstream manifest(argv[2]);double tolerance;size_t count;
        if(!(manifest>>tolerance>>count))throw std::runtime_error("bad allocation manifest");
        std::vector<Samples> levels;std::vector<double> costs;std::vector<int> counts;
        for(size_t level=0;level<count;++level){int n;double cost;manifest>>n>>cost;counts.push_back(n);costs.push_back(cost);Samples data(8);
            for(int seed=0;seed<8;++seed){std::string path,line;manifest>>path;std::ifstream file(path);std::getline(file,line);
                while(std::getline(file,line)){std::istringstream row(line);double theta,area;Mueller m;if(!(row>>theta>>area))continue;
                    for(double& v:m)if(!(row>>v))throw std::runtime_error("bad allocation matrix");
                    data[seed].push_back(m);}}
            levels.push_back(data);}
        Samples combined=levels[0];for(size_t l=1;l<levels.size();++l)for(int s=0;s<8;++s)for(size_t r=0;r<combined[s].size();++r)
            for(int j=0;j<16;++j)combined[s][r][j]+=levels[l][s][r][j];
        auto budgets=Allocate(levels,costs,counts,Estimate(combined).mean,tolerance);
        for(int n:budgets)std::cout<<n<<' ';
        std::cout<<'\n';return 0;
    }
    if(argc==4){
        std::ifstream input(argv[1]),grid(argv[2]);std::string line;std::getline(input,line);
        std::vector<double> x,knots;Matrix y,selected;
        while(std::getline(input,line)){std::istringstream row(line);double t,area;if(!(row>>t>>area))continue;Mueller m;for(double& v:m)if(!(row>>v))throw std::runtime_error("bad matrix");x.push_back(t);y.push_back(m);}
        double v;while(grid>>v)knots.push_back(v);
        for(double k:knots){size_t i=0;while(i<x.size()&&std::fabs(x[i]-k)>1e-12)++i;if(i==x.size())throw std::runtime_error("knot missing");selected.push_back(y[i]);}
        auto out=CubicInterpolate(knots,selected,x);std::ofstream f(argv[3]);f<<std::setprecision(17);
        for(size_t i=0;i<x.size();++i){f<<x[i];for(double a:out[i])f<<' '<<a;f<<'\n';}return 0;
    }
    // Cubic reproduction on an irregular grid checks boundary conditions;
    // a known sample variance checks the confidence interval convention.
    std::vector<double> x={0,.1,.7,1.9,3.};Matrix y(x.size());
    for(size_t i=0;i<x.size();++i)for(int j=0;j<16;++j)y[i][j]=(j+1)*(2+3*x[i]-x[i]*x[i]+.4*x[i]*x[i]*x[i]);
    std::vector<double> q={0,.05,.4,1.,2.4,3.};auto a=CubicInterpolate(x,y,q);
    for(size_t i=0;i<q.size();++i)for(int j=0;j<16;++j)if(std::fabs(a[i][j]-(j+1)*(2+3*q[i]-q[i]*q[i]+.4*q[i]*q[i]*q[i]))>1e-11)throw std::runtime_error("cubic reproduction failed");
    Samples s(8,Matrix(1));for(int k=0;k<8;++k){s[k][0][0]=10+.1*k;s[k][0][1]=.2*k;}
    auto stats=Estimate(s);double expected=2.3646242510102993*std::sqrt(.06/8)/10.35;
    if(std::fabs(stats.interval[0][0]-expected)>1e-14)throw std::runtime_error("confidence convention mismatch");
    auto delta=s;for(int k=0;k<8;++k)delta[k][0][0]=1+.01*k;
    auto budgets=Allocate({s,delta},{1.,10.},{16,16},stats.mean,.01);
    if(budgets[0]<=budgets[1] || budgets[0]<=0)throw std::runtime_error("allocation failed");
    bool rejected=false;try{CubicInterpolate(x,y,{-1});}catch(const std::exception&){rejected=true;}
    if(!rejected)throw std::runtime_error("extrapolation allowed");
    std::cout<<"fullauto numerical tests passed\n";
}
