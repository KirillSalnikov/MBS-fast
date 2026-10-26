#pragma once
#include "BoundedCache.h"
#include <cstdint>
#include <cstring>
#include <limits>
#include <vector>

// Dependency graphs repeat the same long output directory in every key.
// Intern that directory without hashing it: keys remain collision-free and
// a bounded memo can retain substantially more already verified graphs.
// Peek and traversal are read-only while independent roots are validated.
class ValidatedGraphCache {
    BoundedCache<std::string> values;
    std::map<std::string,uint32_t> directoryIds;
    std::vector<std::string> directories;
    size_t maximum,directoryBytes=0;
    static size_t Split(const std::string& path){const auto p=path.rfind('/');return p==std::string::npos?0:p+1;}
    static std::string Key(uint32_t id,const std::string& path,size_t split){
        std::string result(reinterpret_cast<const char*>(&id),sizeof(id));result.append(path,split,std::string::npos);return result;
    }
    bool Lookup(const std::string& path,std::string& key)const{
        const auto split=Split(path);const auto found=directoryIds.find(path.substr(0,split));
        if(found==directoryIds.end())return false;
        key=Key(found->second,path,split);
        return true;
    }
public:
    explicit ValidatedGraphCache(size_t limit=256*1024*1024):values(limit),maximum(limit){}
    size_t Bytes()const{return directoryBytes+values.Bytes();}
    size_t Limit()const{return maximum;}
    size_t Evictions()const{return values.evictions;}
    void Clear(){values.Clear();directoryIds.clear();directories.clear();directoryBytes=0;values.SetLimit(maximum);}
    void SetLimit(size_t limit){maximum=limit;if(directoryBytes>limit)Clear();else values.SetLimit(limit-directoryBytes);}
    std::shared_ptr<const std::string> Get(const std::string& path){std::string key;return Lookup(path,key)?values.Get(key):std::shared_ptr<const std::string>{};}
    std::shared_ptr<const std::string> Peek(const std::string& path)const{std::string key;return Lookup(path,key)?values.Peek(key):std::shared_ptr<const std::string>{};}
    void Put(const std::string& path,std::shared_ptr<const std::string> value,size_t bytes){
        const auto split=Split(path);const auto directory=path.substr(0,split);auto found=directoryIds.find(directory);
        if(found==directoryIds.end()){
            // Bound interned names too, including vector and tree overhead.
            const size_t cost=2*directory.size()+256;
            if(cost>maximum/8||directoryBytes>maximum/8-cost||directories.size()>=std::numeric_limits<uint32_t>::max())return;
            const uint32_t id=uint32_t(directories.size());directories.push_back(directory);
            found=directoryIds.emplace(directory,id).first;directoryBytes+=cost;values.SetLimit(maximum-directoryBytes);
        }
        values.Put(Key(found->second,path,split),std::move(value),bytes);
    }
    template<class Visitor> void LeastRecentFirst(Visitor visitor)const{
        values.LeastRecentFirst([&](const std::string& key,std::shared_ptr<const std::string> value){
            uint32_t id=0;std::memcpy(&id,key.data(),sizeof(id));visitor(directories.at(id)+key.substr(sizeof(id)),value);
        });
    }
};
