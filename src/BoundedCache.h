#pragma once
#include <list>
#include <map>
#include <memory>
#include <string>

// Eviction changes residency only. A caller holding a value keeps it alive
// until its current operation ends; disk remains the authoritative checkpoint.
template<class T> class BoundedCache {
    struct Entry {
        std::shared_ptr<const T> value;
        size_t bytes;
        typename std::list<std::string>::iterator position;
    };
    std::map<std::string,Entry> entries;
    std::list<std::string> recent;
    size_t used=0,maximum;
public:
    size_t hits=0,misses=0,evictions=0;
    explicit BoundedCache(size_t bytes):maximum(bytes){}
    size_t Bytes()const{return used;}
    size_t Limit()const{return maximum;}
    void Clear(){entries.clear();recent.clear();used=0;}
    void SetLimit(size_t bytes){maximum=bytes;Trim();}
    void Erase(const std::string& key){
        auto i=entries.find(key);if(i==entries.end())return;
        used-=i->second.bytes;recent.erase(i->second.position);entries.erase(i);
    }
    void Trim(){while(used>maximum&&!recent.empty()){Erase(recent.back());++evictions;}}
    std::shared_ptr<const T> Get(const std::string& key){
        auto i=entries.find(key);if(i==entries.end()){++misses;return {};}
        recent.splice(recent.begin(),recent,i->second.position);++hits;return i->second.value;
    }
    // Read-only lookup/traversal for a frozen cache borrowed by validation
    // workers. The owner must not mutate it until those workers have joined.
    std::shared_ptr<const T> Peek(const std::string& key)const{
        const auto i=entries.find(key);return i==entries.end()?std::shared_ptr<const T>{}:i->second.value;
    }
    template<class Visitor> void LeastRecentFirst(Visitor visitor)const{
        for(auto i=recent.rbegin();i!=recent.rend();++i)visitor(*i,entries.at(*i).value);
    }
    void Put(const std::string& key,std::shared_ptr<const T> value,size_t bytes){
        Erase(key);
        // Account for two owned keys and the containers as well as the value.
        const size_t overhead=2*key.size()+sizeof(Entry)+192;
        if(bytes>maximum||overhead>maximum-bytes)return;
        bytes+=overhead;
        while(used>maximum-bytes&&!recent.empty()){Erase(recent.back());++evictions;}
        recent.push_front(key);entries.emplace(key,Entry{value,bytes,recent.begin()});used+=bytes;
    }
};
