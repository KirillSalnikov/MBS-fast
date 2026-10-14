#include "FullAuto.h"
#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <stdexcept>
using namespace FullAuto;
int main(int argc,char** argv){
    try{
        if(argc!=2&&(argc!=3||std::string(argv[2])!="--expect-first-correction-only"))
            throw std::runtime_error("expected per-angle level sample file [--expect-first-correction-only]");
        std::ifstream f(argv[1]);size_t rows;double tolerance;
        if(!(f>>rows>>tolerance))throw std::runtime_error("invalid priority input");
        std::cout<<std::setprecision(17)<<"theta,error";
        for(const char* name:{"current","joint","staged"})for(int l=0;l<3;++l)std::cout<<','<<name<<l;
        std::cout<<'\n';
        for(size_t row=0;row<rows;++row){
            double theta;f>>theta;std::vector<int> counts(3);std::vector<double> costs(3);
            std::vector<Samples> levels(3,Samples(8,Matrix(1)));Samples combined(8,Matrix(1));
            for(int l=0;l<3;++l){f>>counts[l]>>costs[l];for(int s=0;s<8;++s)for(int j=0;j<16;++j){
                if(!(f>>levels[l][s][0][j]))throw std::runtime_error("incomplete samples");
                combined[s][0][j]+=levels[l][s][0][j];}}
            const auto stats=Estimate(combined);auto proposed=Allocate(levels,costs,counts,stats.mean,tolerance*.85);
            auto transformed=levels;
            for(auto& level:transformed)for(auto& sample:level)for(int j=1;j<16;++j)
                sample[0][j]-=stats.mean[0][j]/stats.mean[0][0]*sample[0][0];
            auto ratio=Allocate(transformed,costs,counts,stats.mean,tolerance*.85);
            const double error=std::max(stats.maxInterval,stats.maxRatioInterval);
            for(int l=0;l<3;++l)proposed[l]=error>tolerance?std::max(counts[l],std::min(4*counts[l],std::max(proposed[l],ratio[l]))):counts[l];
            const auto saved=levels;auto staged=PrioritizeAllocation(levels,costs,counts,proposed,tolerance);
            if(levels!=saved)throw std::runtime_error("priority changed sampled matrices");
            bool progress=false;for(int l=0;l<3;++l){if(staged[l]<counts[l]||staged[l]>proposed[l])throw std::runtime_error("invalid staged budget");progress=progress||staged[l]>counts[l];}
            if(proposed!=counts&&!progress)throw std::runtime_error("priority prevented all progress");
            if(argc==3&&(staged[0]!=counts[0]||staged[1]!=proposed[1]||staged[1]<=counts[1]||staged[2]!=counts[2]))
                throw std::runtime_error("saved dominant correction was not prioritized");
            if(PrioritizeAllocation(levels,costs,counts,proposed,1e101)!=proposed)throw std::runtime_error("stability work deferred");
            std::cout<<theta<<','<<error;for(const auto& v:{counts,proposed,staged})for(int n:v)std::cout<<','<<n;std::cout<<'\n';
        }
        return 0;
    }catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}
}
