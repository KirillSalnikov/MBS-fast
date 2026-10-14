// Replays saved sampling metadata, or exercises incompatible cached recipes.
// Compare the ordered source identities, not just the number of cache hits.
#include "../src/FullAuto.cpp"
int main(int argc,char** argv){
    try{
        using namespace FullAuto;
        ArgPP args;RunConfig config;Engine e(args,config);
        std::vector<std::string> identities;std::map<std::string,std::string> recipes;
        if(argc==2){
            std::ifstream input(argv[1]);Require(bool(input),"cannot read metadata");
            std::string line;
            while(std::getline(input,line)){
                std::istringstream row(line);std::string configuration;Job job;
                Require(bool(row>>job.folder>>configuration>>job.identity>>job.count>>job.offset>>job.seed
                    >>job.depth>>job.phi>>job.block>>job.pipeline>>job.mirrored>>job.paired>>job.assembled
                    >>job.hasControls>>job.pairedDepth>>job.warp>>job.weightHash>>job.executionHash),"invalid metadata");
                recipes[job.identity]=job.executionHash;
                e.base={"audit_original_recipe",job.executionHash};job.executionHash=e.ExecutionHash();
                job.inputAngles={170,180};e.CacheJob(job);identities.push_back(job.identity);
            }
        }else{
            for(int k=0;k<1024;++k){
                Job j;e.base={"--cutoff-profile",k%2?"off":"safe"};
                e.phi=1+(k/2)%4;e.kernel.block=64*(1+(k/8)%2);e.kernel.pipeline=(k/16)%2;
                e.kernel.warp=(k/32)%2?"thread":"auto";e.mirror=(k/64)%2;e.controls=(k/128)%2;
                e.SetSampling(j,64,11,6,{170,180},(k/256)%2);j.offset=(k/512)*64;
                j.identity=Hash(Number(k));j.folder="fixture/"+j.identity;j.assembled=k%3==0;
                recipes[j.identity]=k%2?"off":"safe";e.CacheJob(j);identities.push_back(j.identity);
            }
        }
        Require(!identities.empty(),"empty sampling metadata");
        const size_t queries=std::min(size_t(1024),identities.size());
        double scanSeconds=0,indexSeconds=0;size_t scanned=0,candidates=0,matches=0;
        for(size_t q=0;q<queries;++q){
            const auto& j=e.cache.at(identities[q*(identities.size()-1)/std::max(size_t(1),queries-1)]);
            e.base={argc==2?"audit_original_recipe":"--cutoff-profile",recipes.at(j.identity)};
            e.phi=j.phi;e.kernel.block=j.block;e.kernel.pipeline=j.pipeline;e.kernel.warp=j.warp;
            e.mirror=j.mirrored;e.controls=j.hasControls;
            std::vector<std::string> before,after;double begin=Seconds();
            for(const auto& pair:e.cache){++scanned;if(e.SameSampling(pair.second,j.count,j.seed,j.depth,j.offset))before.push_back(pair.first);}
            scanSeconds+=Seconds()-begin;begin=Seconds();
            for(const auto& pair:e.SamplingCandidates(j.count,j.seed,j.depth,j.offset)){
                ++candidates;if(e.SameSampling(*pair.second,j.count,j.seed,j.depth,j.offset))after.push_back(pair.first);}
            indexSeconds+=Seconds()-begin;Require(before==after,"indexed source selection differs");matches+=after.size();
            Require(e.folderIndex.at(j.folder)==&j,"folder lookup differs");
        }
        Require(e.SamplingCandidates(19,UINT32_MAX,999,7).empty(),"absent prefix matched");
        Job replacement=e.cache.begin()->second;replacement.seconds=123;e.CacheJob(replacement);
        Require(e.folderIndex.at(replacement.folder)->seconds==123,"updated entry not visible");
        replacement.offset++;bool rejected=false;try{e.CacheJob(replacement);}catch(const std::exception&){rejected=true;}
        Require(rejected,"same identity with changed prefix was accepted");
        std::cout<<std::setprecision(17)<<"{\"jobs\":"<<e.cache.size()<<",\"queries\":"<<queries
            <<",\"linear_candidates\":"<<scanned<<",\"indexed_candidates\":"<<candidates
            <<",\"matches\":"<<matches<<",\"linear_seconds\":"<<scanSeconds<<",\"indexed_seconds\":"<<indexSeconds
            <<",\"ordered_sources_identical\":true,\"changed_identity_rejected\":true}\n";
        return 0;
    }catch(const std::exception& error){std::cerr<<error.what()<<'\n';return 1;}
}
