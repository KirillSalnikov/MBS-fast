#pragma once
#include "BoundedCache.h"
#include <cstdint>
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <sys/stat.h>

// Process-local, bounded cache. Every hit is guarded by device/inode, size,
// nanosecond mtime AND ctime; a restart never trusts old cached digests.
// Small validation manifests share their text and digest. Large matrices keep
// only a digest. Compact inode keys avoid retaining two long checkpoint paths.
class CheckedFileCache {
    struct Entry {struct stat stamp;std::string hash,text;bool hasText=false;};
    BoundedCache<Entry> entries;
    static void Check(bool ok,const std::string& message){if(!ok)throw std::runtime_error("fullauto: "+message);}
    static bool Same(const struct stat& a,const struct stat& b){
        return a.st_dev==b.st_dev&&a.st_ino==b.st_ino&&a.st_size==b.st_size
            &&a.st_mtim.tv_sec==b.st_mtim.tv_sec&&a.st_mtim.tv_nsec==b.st_mtim.tv_nsec
            &&a.st_ctim.tv_sec==b.st_ctim.tv_sec&&a.st_ctim.tv_nsec==b.st_ctim.tv_nsec;
    }
    static std::string Key(const struct stat& stamp){
        std::string key(reinterpret_cast<const char*>(&stamp.st_dev),sizeof(stamp.st_dev));
        key.append(reinterpret_cast<const char*>(&stamp.st_ino),sizeof(stamp.st_ino));return key;
    }
    static std::string Hex(uint64_t hash){std::ostringstream out;out<<std::hex<<hash;return out.str();}
    std::shared_ptr<const Entry> Get(const std::string& path,bool text){
        struct stat before,after;Check(stat(path.c_str(),&before)==0,"cannot stat "+path);
        const auto key=Key(before);const auto found=entries.Get(key);
        if(found&&Same(before,found->stamp)&&(!text||found->hasText)){
            if(text)++textHits;else ++hashHits;return found;
        }
        std::ifstream input(path.c_str(),std::ios::binary);Check(bool(input),"cannot read "+path);++opens;
        // A source marker is normally hashed immediately before recursively
        // validating that source. Retain its bytes during the hash read too,
        // so the following Text() needs no second open/read of the same file.
        const std::string marker="/complete.txt";
        const bool markerText=before.st_size<=65536&&path.size()>=marker.size()&&path.compare(path.size()-marker.size(),marker.size(),marker)==0;
        const bool keepText=text||markerText;
        auto value=std::make_shared<Entry>();value->stamp=before;value->hasText=keepText;
        uint64_t hash=14695981039346656037ULL;char buffer[65536];
        while(input){input.read(buffer,sizeof(buffer));const auto size=input.gcount();bytesRead+=size;
            for(std::streamsize i=0;i<size;++i){hash^=static_cast<unsigned char>(buffer[i]);hash*=1099511628211ULL;}
            if(keepText)value->text.append(buffer,size);
        }
        Check(input.eof(),"cannot read "+path);
        Check(stat(path.c_str(),&after)==0&&Same(before,after),"file changed while reading "+path);
        value->hash=Hex(hash);
        // A large manifest can be returned, but must not displace thousands of
        // small files from the cache. Keep only its digest in that case.
        std::shared_ptr<Entry> retained=value;
        if(text&&value->text.size()>65536){retained=std::make_shared<Entry>();retained->stamp=before;retained->hash=value->hash;}
        entries.Put(key,retained,sizeof(Entry)+retained->hash.capacity()+retained->text.capacity());
        return value;
    }
public:
    size_t opens=0,hashHits=0,textHits=0;uint64_t bytesRead=0;
    explicit CheckedFileCache(size_t limit=32*1024*1024):entries(limit){}
    size_t Bytes()const{return entries.Bytes();}
    size_t Limit()const{return entries.Limit();}
    std::string Hash(const std::string& path){return Get(path,false)->hash;}
    std::string Text(const std::string& path){return Get(path,true)->text;}
    static std::string Read(const std::string& path){
        std::ifstream input(path.c_str(),std::ios::binary);Check(bool(input),"cannot read "+path);
        std::string result;char buffer[65536];
        while(input){input.read(buffer,sizeof(buffer));result.append(buffer,input.gcount());}
        Check(input.eof(),"cannot read "+path);return result;
    }
};
