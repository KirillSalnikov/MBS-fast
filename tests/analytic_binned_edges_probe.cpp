#include "../src/AnalyticFacetAverage.cpp"
#include <cassert>
#include <iostream>

int main()
{
    using namespace AnalyticFacetAverage;
    struct Patch { std::vector<double> x, y; };
    const std::vector<Patch> patches={{{0, 1, 0.32}, {0, 0.2, 0.73}}};
    const auto& p=patches.front();
    std::vector<double> nodes, weights;
    Gauss(64,nodes,weights);
    for(size_t i=0;i<nodes.size();++i){nodes[i]=(1+nodes[i])/2;weights[i]/=2;}

    const double step=.03;
    BinnedEdgeMeasure measure;
    measure.Build(patches,{0},{2},nodes,weights,step,80.,2.);
    assert(!measure.bins.empty());
    size_t budget=0;
    EdgeDistances direct(p.x,p.y,nodes,budget);
    Basis basis(128,165.,step);
    Moments moments({1.3243,0},128);
    size_t accepted=0;
    for(double q:{.01,.14,1.15,4.9,12.4,39.5,78.7}){
        const auto green=basis.Combine(moments.At(.71),2*q+step);
        const double expected=2*EdgeMean(p.x,p.y,q,nodes,weights,direct,
            [&](double z){return green.Value(z);},[&](double z){return green.Same(z);});
        double actual=0;
        if(measure.Evaluate(green,q,actual)){
            assert(std::fabs(actual-expected)<1.1e-8*std::abs(expected));
            ++accepted;
        }
    }
    assert(accepted>=3);
    BinnedEdgeMeasure oversize;
    oversize.Build(patches,{0},{1},nodes,weights,step,1.e9,2.);
    assert(oversize.bins.empty());
    std::cout<<"PASS: validated "<<accepted<<" binned quadratures and bounded allocation\n";
}
