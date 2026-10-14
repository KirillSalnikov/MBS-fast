#pragma once
#include "handler/HandlerPO.h"
#include <list>
#include <map>
#include <mutex>
#include <memory>
// Process-local, bounded cache of angle-independent CPU beam preparation.
// Keys contain the complete physics/sampling recipe and exact orientation.
namespace PreparedReuse {
struct Value {
    PreparedOrientation upper,lower;double incoming=0,outgoing=0;
    static size_t PreparedBytes(const PreparedOrientation& p){
        size_t n=p.beams.capacity()*sizeof(PreparedBeam);
        for(const auto& b:p.beams){n+=b.absorptionPaths.capacity()*sizeof(double);
            if(b.edgeData.HasOverflow())n+=sizeof(BeamEdgeData);if(b.fallback)n+=sizeof(PreparedBeamFallback);}
        return n;
    }
    size_t Bytes()const{return sizeof(Value)+PreparedBytes(upper)+PreparedBytes(lower);}
};
struct Store {
    std::mutex mutex;size_t limit=0,bytes=0;unsigned long long hits=0,misses=0;
    std::list<std::string> lru;
    struct Entry {std::shared_ptr<const Value> value;std::list<std::string>::iterator position;size_t bytes;};
    std::map<std::string,Entry> entries;
};
inline Store& Cache(){static Store value;return value;}
inline void Configure(size_t limit){auto& c=Cache();std::lock_guard<std::mutex> lock(c.mutex);c.limit=limit;c.hits=c.misses=0;
    while(c.bytes>limit&&!c.lru.empty()){auto i=c.entries.find(c.lru.back());c.bytes-=i->second.bytes;c.entries.erase(i);c.lru.pop_back();}}
inline std::shared_ptr<const Value> Get(const std::string& key){auto& c=Cache();std::lock_guard<std::mutex> lock(c.mutex);
    auto i=c.entries.find(key);if(i==c.entries.end()){++c.misses;return {};}
    ++c.hits;c.lru.splice(c.lru.begin(),c.lru,i->second.position);return i->second.value;}
inline void Put(const std::string& key,Value value){auto& c=Cache();const size_t bytes=value.Bytes()+2*key.size()+256;
    std::lock_guard<std::mutex> lock(c.mutex);if(bytes>c.limit||c.entries.count(key))return;
    while(c.bytes+bytes>c.limit&&!c.lru.empty()){auto i=c.entries.find(c.lru.back());c.bytes-=i->second.bytes;c.entries.erase(i);c.lru.pop_back();}
    auto saved=std::make_shared<Value>(std::move(value));c.lru.push_front(key);c.entries.emplace(key,Store::Entry{saved,c.lru.begin(),bytes});c.bytes+=bytes;}
}
