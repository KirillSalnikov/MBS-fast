// Metadata-only replay. No optical work, payload reads, or production changes.
// Queries use each sampled job's grid, not an observed trace of Stage calls.
#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <stdexcept>
#include <string>
#include <tuple>
#include <vector>
#include <sys/resource.h>

namespace {
using Clock=std::chrono::steady_clock;
double Now(){return std::chrono::duration<double>(Clock::now().time_since_epoch()).count();}
void Check(bool yes,const char* message){if(!yes)throw std::runtime_error(message);}
struct Job {
    std::string source,configuration,identity,warp,weightHash,executionHash;
    int count=0,offset=0,depth=0,phi=0,block=0,pipeline=0,mirrored=0,paired=0,assembled=0,hasControls=0,pairedDepth=0;
    unsigned seed=0;std::vector<double> angles;
};
using Namespace=std::pair<std::string,std::string>;
using Coarse=std::tuple<Namespace,int,unsigned,int,int>;
using Complete=std::tuple<Coarse,int,int,int,int,int,std::string,std::string>;
Coarse Bucket(const Job& j){return {{j.source,j.configuration},j.count,j.seed,j.depth,j.offset};}
Complete Key(const Job& j){return {Bucket(j),j.phi,j.block,j.pipeline,j.mirrored,j.hasControls,j.warp,j.executionHash};}
bool Same(const Job& j,const Job& q){
    return j.offset==q.offset&&j.count==q.count&&j.seed==q.seed&&j.depth==q.depth&&j.phi==q.phi&&j.block==q.block
        &&j.pipeline==q.pipeline&&j.warp==q.warp&&j.mirrored==q.mirrored&&!j.paired&&j.hasControls==q.hasControls&&j.executionHash==q.executionHash;
}
using Row=std::pair<const Job*,size_t>;
struct Selection {std::map<double,Row> rows;const Job* exact=nullptr;};
struct Counters {uint64_t queries=0,bypassed=0,candidates=0,eligible=0,rowMatches=0,ownerLookups=0,exact=0;};
struct Owned {std::map<double,Row> first;std::vector<const Job*> exactCandidates;};
struct Replay {
    std::vector<Job> requests;
    std::map<std::pair<Namespace,std::string>,Job> cache;
    std::map<Coarse,std::map<std::string,const Job*>> buckets;
    std::map<Complete,Owned> index;
    size_t parsedAngles=0,duplicates=0;double loadSeconds=0,indexBuildSeconds=0;
    void Load(const std::string& path){const double started=Now();std::ifstream file(path);Check(bool(file),"cannot read replay");
        std::string magic;size_t count=0;file>>magic>>count;Check(magic=="MBS_ANGULAR_METADATA_V1"&&count<=1000000,"invalid replay header");requests.reserve(count);
        for(size_t k=0;k<count;++k){Job j;size_t rows=0;
            Check(bool(file>>j.source>>j.configuration>>j.identity>>j.count>>j.offset>>j.seed>>j.depth>>j.phi>>j.block
                >>j.pipeline>>j.mirrored>>j.paired>>j.assembled>>j.hasControls>>j.pairedDepth>>j.warp>>j.weightHash>>j.executionHash>>rows),"truncated metadata");
            Check(rows>=2&&rows<=1000000&&parsedAngles+rows<=10000000,"invalid sampled angle count");j.angles.resize(rows);parsedAngles+=rows;
            for(auto& angle:j.angles)Check(bool(file>>angle)&&std::isfinite(angle),"invalid sampled angle");
            Check(std::is_sorted(j.angles.begin(),j.angles.end()),"sample grid is not sorted");
            requests.push_back(j);const auto key=std::make_pair(Namespace{j.source,j.configuration},j.identity);duplicates+=cache.count(key);cache[key]=std::move(j);
        }
        std::string extra;Check(!(file>>extra),"unexpected replay tail");
        for(const auto& entry:cache){const Job& j=entry.second;buckets[Bucket(j)][j.identity]=&j;}
        loadSeconds=Now()-started;
    }
    Selection Baseline(const Job& q,Counters& counter)const{Selection result;++counter.queries;if(q.paired){++counter.bypassed;return result;}
        const auto found=buckets.find(Bucket(q));if(found==buckets.end())return result;
        for(const auto& item:found->second){const Job& j=*item.second;++counter.candidates;if(!Same(j,q))continue;++counter.eligible;
            if(j.angles==q.angles&&j.weightHash==q.weightHash)result.exact=&j;
            auto matched=[&](size_t row){++counter.rowMatches;result.rows.insert({j.angles[row],{&j,row}});};
            // Exact copy of the old MatchingAngularRows branch selection.
            if(q.angles.size()<j.angles.size()/8){for(double angle:q.angles){const auto at=std::lower_bound(j.angles.begin(),j.angles.end(),angle);if(at!=j.angles.end()&&*at==angle)matched(size_t(at-j.angles.begin()));}}
            else {size_t wanted=0;for(size_t row=0;row<j.angles.size()&&wanted<q.angles.size();++row){
                while(wanted<q.angles.size()&&q.angles[wanted]<j.angles[row])++wanted;
                if(wanted<q.angles.size()&&q.angles[wanted]==j.angles[row])matched(row);}}
        }
        counter.exact+=result.exact!=nullptr;return result;
    }
    void BuildIndex(){const double started=Now();for(const auto& item:cache){const Job& j=item.second;if(j.paired)continue;
        auto& owned=index[Key(j)];owned.exactCandidates.push_back(&j);for(size_t row=0;row<j.angles.size();++row)owned.first.insert({j.angles[row],{&j,row}});}
        indexBuildSeconds=Now()-started;
    }
    Selection Indexed(const Job& q,Counters& counter)const{Selection result;++counter.queries;if(q.paired){++counter.bypassed;return result;}
        const auto found=index.find(Key(q));if(found==index.end())return result;const auto& owned=found->second;
        for(auto at=owned.exactCandidates.rbegin();at!=owned.exactCandidates.rend();++at){++counter.candidates;if((*at)->angles==q.angles&&(*at)->weightHash==q.weightHash){result.exact=*at;break;}}
        for(double angle:q.angles){++counter.ownerLookups;const auto row=owned.first.find(angle);if(row!=owned.first.end()){result.rows.insert(*row);++counter.rowMatches;}}
        counter.exact+=result.exact!=nullptr;return result;
    }
};
uint64_t HashText(uint64_t value,const std::string& text){for(unsigned char byte:text){value^=byte;value*=1099511628211ULL;}return value;}
uint64_t Digest(const Selection& s){uint64_t value=14695981039346656037ULL;value=HashText(value,s.exact?s.exact->identity:"none");
    value^=s.rows.size();value*=1099511628211ULL;
    // Constant-work anti-elision checksum. Full row-by-row equivalence is
    // checked separately, outside the timed region.
    if(!s.rows.empty())for(const auto* row:{&s.rows.begin()->second,&s.rows.rbegin()->second}){value=HashText(value,row->first->identity);value^=row->second;value*=1099511628211ULL;}
    return value;}
void Equal(const Selection& a,const Selection& b){Check(a.exact==b.exact&&a.rows.size()==b.rows.size(),"index changed exact source or row count");
    auto x=a.rows.begin(),y=b.rows.begin();for(;x!=a.rows.end();++x,++y)Check(x->first==y->first&&x->second==y->second,"index changed row precedence");}
struct Timing {double once=0,repeated=0;size_t repetitions=0;uint64_t digest=0;Counters counters;};
volatile uint64_t sink=0;
Timing Measure(const Replay& replay,bool index,size_t maximum,double cap){Timing out;
    auto pass=[&](Counters& counters){uint64_t value=14695981039346656037ULL;for(const auto& q:replay.requests){const auto selected=index?replay.Indexed(q,counters):replay.Baseline(q,counters);value^=Digest(selected);value*=1099511628211ULL;}sink=value;return value;};
    const double first=Now();out.digest=pass(out.counters);out.once=Now()-first;
    const double started=Now();double last=out.once;
    while(out.repetitions<maximum&&Now()-started+last<=cap){Counters unused;const double begin=Now();const auto digest=pass(unused);last=Now()-begin;Check(digest==out.digest,"replay changed between repetitions");++out.repetitions;}
    out.repeated=Now()-started;return out;
}
void PrintTiming(const Timing& t){std::cout<<"{\"once_seconds\":"<<t.once<<",\"repetitions\":"<<t.repetitions<<",\"repeat_seconds\":"<<t.repeated
    <<",\"repeat_mean_seconds\":"<<(t.repetitions?t.repeated/t.repetitions:0)<<",\"digest\":\""<<std::hex<<t.digest<<std::dec<<"\",\"queries\":"<<t.counters.queries
    <<",\"paired_mirror_bypassed\":"<<t.counters.bypassed<<",\"candidate_comparisons\":"<<t.counters.candidates<<",\"eligible_candidates\":"<<t.counters.eligible
    <<",\"row_matches\":"<<t.counters.rowMatches<<",\"owner_lookups\":"<<t.counters.ownerLookups<<",\"exact_queries\":"<<t.counters.exact<<"}";}
}
int main(int argc,char** argv){try{Check(argc==2,"usage: angular_benchmark replay.txt");Replay replay;replay.Load(argv[1]);
    // Baseline first, with the input already in memory; no index build cost
    // or file parsing can contaminate the reported scan duration.
    const auto baseline=Measure(replay,false,100,5);replay.BuildIndex();
    Counters before,after;for(const auto& request:replay.requests)Equal(replay.Baseline(request,before),replay.Indexed(request,after));
    const auto indexed=Measure(replay,true,100,5);Check(baseline.digest==indexed.digest,"aggregate selection mismatch");
    size_t largest=0,indexRows=0;for(const auto& b:replay.buckets)largest=std::max(largest,b.second.size());for(const auto& entry:replay.index)indexRows+=entry.second.first.size();
    rusage usage{};getrusage(RUSAGE_SELF,&usage);std::cout<<std::setprecision(17)<<"{\"records\":"<<replay.requests.size()<<",\"unique_jobs\":"<<replay.cache.size()
        <<",\"metadata_replacements\":"<<replay.duplicates<<",\"input_angle_values\":"<<replay.parsedAngles<<",\"coarse_buckets\":"<<replay.buckets.size()
        <<",\"largest_coarse_bucket\":"<<largest<<",\"complete_keys\":"<<replay.index.size()<<",\"index_row_owners\":"<<indexRows<<",\"load_and_bucket_build_seconds\":"<<replay.loadSeconds
        <<",\"simulated_index_build_seconds\":"<<replay.indexBuildSeconds<<",\"peak_rss_kib\":"<<usage.ru_maxrss<<",\"repeat_cap_seconds_per_method\":5,\"baseline\":";
    PrintTiming(baseline);std::cout<<",\"simulated_unbounded_index\":";PrintTiming(indexed);
    std::cout<<",\"selection_equal\":true,\"scope\":\"Incomplete head/tail metadata sample; each sampled job grid is a replay query, not a Stage call trace. Timings exclude payload verification/native work/disk. Simulated index is unbounded and includes neither eviction nor mutation cost. Do not extrapolate to the complete campaign.\"}\n";return 0;
}catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}}
