// Include the implementation to exercise internal radial primitives without
// adding a test-only public API to the physical model.
#include "../src/AnalyticFacetAverage.cpp"
#include <iostream>
int main(){
    AnalyticFacetAverage::Basis basis(320,24000,.015);
    std::cout<<std::setprecision(17);
    int degree;double z;
    while(std::cin>>degree>>z){
        if(degree<0||degree>320||degree%2||z<64)return 2;
        std::vector<double> factors(161,0);factors[degree/2]=1;
        const auto value=basis.Exact(z,factors);
        std::cout<<value[0]<<' '<<value[1]<<'\n';
    }
}
