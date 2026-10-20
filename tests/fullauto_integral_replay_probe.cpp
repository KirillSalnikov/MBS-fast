// Read-only replay of a saved production graph. Compare output hashes and use
// strace -c for file-open counts when comparing implementations.
#include "../src/FullAuto.cpp"
using namespace FullAuto;
int main(int argc,char** argv){try{
    Require(argc==3,"source list and repetition count required");
    const auto sources=Read(argv[1]);const int repetitions=std::stoi(argv[2]);
    Require(repetitions>0,"positive repetition count required");
    const std::vector<unsigned> seeds(production,production+8);
    std::string expected;double seconds=0;size_t files=0;
    for(int trial=0;trial<repetitions;++trial){
        const double begun=Seconds();const auto result=IntegralProductionSources(sources,seeds,1);seconds+=Seconds()-begun;
        std::ostringstream serialized;serialized<<std::setprecision(17)<<result.directTargetFallback<<'\n'<<result.selectedSources;
        for(const auto& value:result.samples)serialized<<value.area<<' '<<value.ext<<'\n';
        for(const auto& file:result.files)serialized<<file.first<<' '<<file.second<<'\n';
        const auto hash=Hash(serialized.str());if(trial==0)expected=hash;else Require(hash==expected,"replayed integral changed");
        files=result.files.size();
    }
    std::cout<<std::setprecision(9)<<"{\"result_hash\":"<<Quoted(expected)<<",\"source_files\":"<<files
        <<",\"repetitions\":"<<repetitions<<",\"seconds_per_replay\":"<<seconds/repetitions
        <<",\"checked_cache_opens\":"<<ValidationFiles().opens<<",\"checked_cache_bytes\":"<<ValidationFiles().Bytes()
        <<",\"checked_cache_limit\":"<<ValidationFiles().Limit()<<"}\n";
    return 0;
}catch(const std::exception& error){std::cerr<<error.what()<<'\n';return 1;}}
