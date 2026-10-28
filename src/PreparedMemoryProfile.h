#pragma once
#include <algorithm>
#include <cstddef>
#include <map>
#include <mutex>
#include <string>

// A memory estimate, never a cached optical result. A persistent worker learns
// across seeds/offsets of the same physics recipe. Live RAM checks still bound
// every chunk, including the first one. Retain at most 32 small records.
namespace PreparedMemoryProfile {
struct Estimate {
    long long safeBytes=0;
    size_t observations=0;
    int largestChunk=0;
    unsigned long long touched=0;
    int Start(int maximum,long long budget)const{
        if(observations<128||safeBytes<=0||largestChunk<32)return std::min(16,maximum);
        // Another factor of two protects new orientations; never start larger
        // than an actually measured chunk or 256 orientations.
        const long long memory=budget/std::max(1LL,safeBytes)/2;
        return int(std::max(1LL,std::min(memory,(long long)std::min(maximum,std::min(256,largestChunk)))));
    }
};
struct Store {std::mutex mutex;std::map<std::string,Estimate> values;unsigned long long clock=0;};
inline Store& Cache(){static Store s;return s;}
inline Estimate Get(const std::string& key){
    if(key.empty())return {};
    auto& s=Cache();std::lock_guard<std::mutex> lock(s.mutex);
    auto it=s.values.find(key);if(it==s.values.end())return {};
    it->second.touched=++s.clock;return it->second;
}
inline void Observe(const std::string& key,long long bytes,size_t count){
    if(key.empty()||bytes<=0||!count)return;
    auto& s=Cache();std::lock_guard<std::mutex> lock(s.mutex);
    if(!s.values.count(key)&&s.values.size()>=32){
        auto oldest=s.values.begin();for(auto it=s.values.begin();it!=s.values.end();++it)
            if(it->second.touched<oldest->second.touched)oldest=it;
        s.values.erase(oldest);
    }
    auto& value=s.values[key];value.safeBytes=std::max(value.safeBytes,bytes);
    value.observations=std::min(size_t(1000000),value.observations+count);
    value.largestChunk=std::max(value.largestChunk,int(count));value.touched=++s.clock;
}
}
