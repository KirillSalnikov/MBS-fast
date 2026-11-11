#pragma once
// Controller working data, independent of durable checkpoints. Small files
// stay in a bounded RAM cache; evicted bytes go to one anonymous temporary
// file. Only Snapshot() performs a durability barrier. The caller owns returned
// Read() strings, which are temporary working buffers outside this cache budget.
#include <algorithm>
#include <cerrno>
#include <cstdint>
#include <cstring>
#include <exception>
#include <fcntl.h>
#include <limits>
#include <list>
#include <map>
#include <mutex>
#include <stdexcept>
#include <string>
#include <sys/stat.h>
#include <unistd.h>
#include <utility>
#include <vector>
#include "FullAutoMemoryCatalog.h"

#ifndef MBS_MEMORY_STORE_EVENT
#define MBS_MEMORY_STORE_EVENT(name,path) ((void)0)
#endif
#ifndef MBS_MEMORY_STORE_SPILL_OPEN
#define MBS_MEMORY_STORE_SPILL_OPEN(fd) ((void)0)
#endif

namespace FullAutoMemoryStore {
inline void Require(bool okay,const std::string& message){if(!okay)throw std::runtime_error("fullauto memory store: "+message);}
inline uint64_t HashBytes(const char* data,size_t length,uint64_t value=UINT64_C(14695981039346656037)){
    for(size_t i=0;i<length;++i){value^=static_cast<unsigned char>(data[i]);value*=UINT64_C(1099511628211);}return value;
}
inline uint64_t Hash(const std::string& bytes){return HashBytes(bytes.data(),bytes.size());}
inline std::string Canonical(const std::string& input){
    Require(!input.empty()&&input.size()<=4096&&input[0]=='/'&&input.find('\0')==std::string::npos,"expected an absolute path of at most 4096 bytes");
    std::string out;size_t first=1;
    while(first<input.size()){
        const size_t last=input.find('/',first);const auto part=input.substr(first,last==std::string::npos?last:last-first);
        Require(part!="..","parent traversal is not allowed");if(!part.empty()&&part!=".")out+="/"+part;
        if(last==std::string::npos)break;
        first=last+1;
    }
    return out.empty()?"/":out;
}
inline void WriteAll(int fd,const char* data,size_t length){
    while(length){const ssize_t n=::write(fd,data,length);if(n<0&&errno==EINTR)continue;
        Require(n>0,"write failed: "+std::string(std::strerror(errno)));data+=n;length-=size_t(n);}
}
inline void ReadAt(int fd,uint64_t offset,char* data,size_t length){
    Require(offset<=uint64_t(std::numeric_limits<off_t>::max()),"file offset overflow");
    while(length){const ssize_t n=::pread(fd,data,length,off_t(offset));if(n<0&&errno==EINTR)continue;
        Require(n>0,"truncated or unreadable spill/checkpoint");data+=n;length-=size_t(n);offset+=uint64_t(n);}
}
inline void Sync(int fd){int rc;do{rc=::fsync(fd);}while(rc<0&&errno==EINTR);Require(rc==0,"checkpoint fsync failed");}
struct Descriptor {
    int fd;
    explicit Descriptor(int value=-1):fd(value){}
    ~Descriptor(){if(fd>=0)::close(fd);}
    Descriptor(const Descriptor&)=delete;Descriptor& operator=(const Descriptor&)=delete;
    int Release(){const int result=fd;fd=-1;return result;}
};
struct Statistics {
    size_t residentBytes=0,metadataBytes=0,limitBytes=0,metadataLimitBytes=0,files=0,directories=0;
    uint64_t spillBytes=0,liveSpillBytes=0,spillReadBytes=0,snapshotReadBytes=0,spillWriteBytes=0,evictions=0,snapshots=0,epoch=0;
    bool catalogActive=false,healthy=true;uint64_t catalogBytes=0;size_t catalogCacheBytes=0;
};
class Store {
    struct Entry {
        bool directory=false,resident=false,spilled=false;
        uint64_t size=0,digest=UINT64_C(14695981039346656037),offset=0,generation=0;
        size_t metadata=0;
        std::string bytes;
        std::list<std::string>::iterator recent;
    };
    std::string root,spillDirectory;
    size_t limit,metadataLimit,residentBytes=0,metadataBytes=0;
    std::map<std::string,Entry> entries;
    std::list<std::string> recent;
    std::unique_ptr<FullAutoMemoryCatalog::Catalog> catalog;
    size_t fileCount=0,directoryCount=0;
    bool poisoned=false,streamCallbackActive=false;
    int spill=-1;
    uint64_t spillEnd=0,liveSpill=0,spillReads=0,snapshotReads=0,spillWrites=0,evictions=0,snapshots=0,epoch=0;
    mutable std::recursive_mutex mutex;
    struct Mutation {
        Store& store;bool alreadyUnwinding,allBackends;
        explicit Mutation(Store& owner,bool all=false):store(owner),alreadyUnwinding(std::uncaught_exception()),allBackends(all){}
        ~Mutation(){if((store.catalog||allBackends)&&!alreadyUnwinding&&std::uncaught_exception())store.poisoned=true;}
    };
    struct StreamCallback {
        Store& store;
        explicit StreamCallback(Store& owner):store(owner){Require(!store.streamCallbackActive,"stream callbacks cannot reenter this store");store.streamCallbackActive=true;}
        ~StreamCallback(){store.streamCallbackActive=false;}
    };
    struct StagedResident {
        Store& store;size_t charge;bool transferred=false;
        StagedResident(Store& owner,size_t bytes):store(owner),charge(bytes){store.residentBytes+=charge;}
        ~StagedResident(){if(!transferred)store.residentBytes-=charge;}
    };
    void RequireHealthy()const{Require(!streamCallbackActive,"stream callbacks cannot reenter this store");Require(!poisoned&&(!catalog||catalog->Healthy()),"working catalog mutation failed; cannot use or checkpoint this state; restore the previous valid checkpoint");}
    static size_t Charge(const std::string& path){
        // Conservatively cover map/list nodes, Entry, owned path strings and
        // allocator overhead. Unlike payload-only caches, millions of empty
        // file names cannot bypass the RAM limit.
        Require(path.size()<(std::numeric_limits<size_t>::max()-512)/3,"metadata size overflow");
        return 512+3*path.size();
    }
    static size_t ResidentCharge(const std::string& value){return value.capacity()+1;}
    std::string Key(const std::string& path)const{
        const auto key=Canonical(path);Require(key==root||(key.size()>root.size()&&key.compare(0,root.size(),root)==0&&(root=="/"||key[root.size()]=='/')),"path outside owned root");return key;
    }
    static std::string Parent(const std::string& path){const auto slash=path.find_last_of('/');return slash==0?"/":path.substr(0,slash);}
    std::string Relative(const std::string& path)const{return path==root?"":path.substr(root=="/"?0:root.size());}
    std::string Absolute(const std::string& relative)const{return relative.empty()?root:(root=="/"?relative:root+relative);}
    uint64_t Next(){Require(epoch<UINT64_MAX,"generation overflow");return ++epoch;}
    size_t CatalogCache()const{return std::min<size_t>(limit/8,1024*1024);}
    size_t CatalogReservation()const{return CatalogCache()+256*1024;}
    static FullAutoMemoryCatalog::Row Row(const Entry& e){FullAutoMemoryCatalog::Row row;row.directory=e.directory;row.size=e.size;row.digest=e.digest;
        // The catalog may index a payload still owned by a resident Entry.
        // Such bytes must be spilled before that cache entry can be dropped.
        row.offset=!e.directory&&!e.spilled?UINT64_MAX:e.offset;row.generation=e.generation;return row;}
    static bool Spilled(const FullAutoMemoryCatalog::Row& row){return !row.directory&&row.offset!=UINT64_MAX;}
    bool Lookup(const std::string& path,FullAutoMemoryCatalog::Row& row)const{
        const auto found=entries.find(path);if(found!=entries.end()){row=Row(found->second);return true;}return catalog&&catalog->Get(path,row);
    }
    int NewSpill()const{
        std::string name=spillDirectory+"/mbs-fullauto-working-XXXXXX";std::vector<char> text(name.begin(),name.end());text.push_back('\0');
        const int fd=::mkstemp(text.data());Require(fd>=0,"cannot create temporary spill in "+spillDirectory);
        if(::unlink(text.data())!=0){::close(fd);throw std::runtime_error("fullauto memory store: cannot unlink temporary spill");}
        ::fcntl(fd,F_SETFD,FD_CLOEXEC);MBS_MEMORY_STORE_SPILL_OPEN(fd);return fd;
    }
    void OpenSpill(){if(spill<0)spill=NewSpill();}
    void SaveSpill(Entry& e,const std::string& path=std::string()){
        if(e.spilled)return;
        Require(e.resident&&!e.directory,"invalid spill source");OpenSpill();
        Require(e.size<=uint64_t(std::numeric_limits<off_t>::max())-spillEnd,"spill size overflow");
        Require(::lseek(spill,off_t(spillEnd),SEEK_SET)==off_t(spillEnd),"cannot seek spill");
        WriteAll(spill,e.bytes.data(),e.bytes.size());e.offset=spillEnd;e.spilled=true;spillEnd+=e.size;liveSpill+=e.size;spillWrites+=e.size;
        if(catalog){Require(!path.empty(),"missing resident catalog spill path");catalog->Put(path,Row(e));}
    }
    void EvictOne(){
        Require(!recent.empty(),"metadata exceeds working-memory limit");const auto path=recent.back();auto& e=entries.at(path);
        SaveSpill(e,path);residentBytes-=ResidentCharge(e.bytes);std::string().swap(e.bytes);e.resident=false;recent.pop_back();++evictions;
        if(catalog){metadataBytes-=e.metadata;entries.erase(path);}
    }
    void Room(size_t extra){
        if(catalog){Require(extra<=limit-CatalogReservation(),"working value exceeds memory cache budget");while(metadataBytes+residentBytes>limit-extra)EvictOne();return;}
        Require(extra<=limit&&metadataBytes<=limit-extra,"metadata exceeds working-memory limit; increase --fullauto-memory-mib");
        while(residentBytes>limit-metadataBytes-extra)EvictOne();
    }
    void Touch(Entry& e){if(e.resident)recent.splice(recent.begin(),recent,e.recent);}
    void Drop(std::map<std::string,Entry>::iterator found){
        auto& e=found->second;if(e.resident){residentBytes-=ResidentCharge(e.bytes);recent.erase(e.recent);}
        if(e.spilled)liveSpill-=e.size;
        if(e.directory)--directoryCount;else --fileCount;metadataBytes-=e.metadata;entries.erase(found);
    }
    void DropResident(const std::string& path){
        const auto found=entries.find(path);if(found==entries.end())return;const auto& e=found->second;
        if(e.resident){residentBytes-=ResidentCharge(e.bytes);recent.erase(e.recent);}metadataBytes-=e.metadata;entries.erase(found);
    }
    void CacheCatalog(const std::string& key,const FullAutoMemoryCatalog::Row& row,const std::string& bytes){
        const size_t charge=Charge(key),dataCharge=ResidentCharge(bytes);
        if(charge>limit-CatalogReservation()||dataCharge>limit-CatalogReservation()-charge)return;
        DropResident(key);while(metadataBytes+residentBytes>limit-charge-dataCharge)EvictOne();
        Entry e;e.size=row.size;e.digest=row.digest;e.offset=row.offset;e.generation=row.generation;e.spilled=true;e.resident=true;e.bytes=bytes;e.metadata=charge;
        recent.push_front(key);e.recent=recent.begin();metadataBytes+=charge;residentBytes+=ResidentCharge(e.bytes);entries.emplace(key,std::move(e));
    }
    void MigrateCatalog(){
        if(catalog)return;
        Require(limit>=1024*1024,"metadata exceeds working-memory limit; use at least 1 MiB for metadata spill");
        std::unique_ptr<FullAutoMemoryCatalog::Catalog> staged(new FullAutoMemoryCatalog::Catalog(spillDirectory,CatalogCache()));
        for(auto& item:entries){if(!item.second.directory)SaveSpill(item.second);staged->Put(item.first,Row(item.second));}
        entries.clear();recent.clear();residentBytes=0;metadataBytes=CatalogReservation();catalog.swap(staged);
    }
    std::vector<std::string> MissingParents(const std::string& path)const{
        std::vector<std::string> result;std::string current=Parent(path);
        while(current!=root){
            FullAutoMemoryCatalog::Row row;if(Lookup(current,row)){Require(row.directory,"parent is a file");break;}
            result.push_back(current);current=Parent(current);
        }
        std::reverse(result.begin(),result.end());return result;
    }
    void AddDirectory(const std::string& path){Entry e;e.directory=true;e.metadata=Charge(path);e.generation=Next();++directoryCount;
        if(catalog){catalog->Put(path,Row(e));return;}metadataBytes+=e.metadata;entries.emplace(path,std::move(e));}
    void Preflight(const std::string& key,const std::vector<std::string>& parents){
        if(catalog)return;
        size_t required=entries.count(key)?0:Charge(key);
        for(const auto& p:parents){Require(required<=limit&&Charge(p)<=limit-required,"metadata exceeds working-memory limit");required+=Charge(p);}
        if(limit>=1024*1024&&required<=limit&&(metadataBytes>limit-required||metadataBytes+required>metadataLimit)){MigrateCatalog();return;}
        Room(required);
    }
    void DirectoryCatalog(const std::string& key){
        FullAutoMemoryCatalog::Row row;if(Lookup(key,row)){Require(row.directory,"directory path is a file");return;}
        const auto parents=MissingParents(key);for(const auto& path:parents)AddDirectory(path);AddDirectory(key);
    }
    void WriteCatalog(const std::string& key,const std::string& value){
        FullAutoMemoryCatalog::Row old;const bool exists=Lookup(key,old);Require(!exists||!old.directory,"cannot replace directory with a file");
        Require(!exists||!Spilled(old)||(old.offset<=spillEnd&&old.size<=spillEnd-old.offset&&old.size<=liveSpill),"invalid overwritten catalog range");
        const auto parents=MissingParents(key);for(const auto& path:parents)AddDirectory(path);OpenSpill();
        Require(value.size()<=uint64_t(std::numeric_limits<off_t>::max())-spillEnd,"spill size overflow");
        Require(::lseek(spill,off_t(spillEnd),SEEK_SET)==off_t(spillEnd),"cannot seek spill");WriteAll(spill,value.data(),value.size());
        FullAutoMemoryCatalog::Row row;row.size=value.size();row.digest=Hash(value);row.offset=spillEnd;row.generation=Next();spillEnd+=value.size();spillWrites+=value.size();
        catalog->Put(key,row);liveSpill+=row.size;if(exists){if(Spilled(old))liveSpill-=old.size;}else ++fileCount;
        DropResident(key);CacheCatalog(key,row,value);Compact();
    }
    void RemoveCatalog(const std::string& key,bool recursive){
        FullAutoMemoryCatalog::Row row;if(!Lookup(key,row))return;auto children=List(key,true);Require(recursive||children.empty(),"directory is not empty");children.push_back(key);
        for(const auto& path:children){FullAutoMemoryCatalog::Row item;Require(Lookup(path,item),"catalog entry disappeared");DropResident(path);catalog->Erase(path);
            if(item.directory)--directoryCount;else{if(Spilled(item)){Require(item.size<=liveSpill,"invalid removed catalog range");liveSpill-=item.size;}--fileCount;}}
        Next();Compact();
    }
    void Compact(){
        // Overwritten/deleted data must not grow the temporary file forever.
        if(spill<0||spillEnd<UINT64_C(64)*1024*1024||liveSpill>spillEnd/2)return;
        if(catalog){
            while(!recent.empty())EvictOne();
            Descriptor replacement(NewSpill());uint64_t end=0;char buffer[65536];
            std::unique_ptr<FullAutoMemoryCatalog::Catalog> next(new FullAutoMemoryCatalog::Catalog(spillDirectory,CatalogCache()));
            catalog->All([&](const std::string& path,FullAutoMemoryCatalog::Row row){
                if(!row.directory){Require(row.offset<=spillEnd&&row.size<=spillEnd-row.offset,"invalid compacted catalog range");uint64_t done=0,digest=UINT64_C(14695981039346656037);const auto oldOffset=row.offset;row.offset=end;
                    while(done<row.size){const auto n=size_t(std::min<uint64_t>(sizeof(buffer),row.size-done));ReadAt(spill,oldOffset+done,buffer,n);digest=HashBytes(buffer,n,digest);WriteAll(replacement.fd,buffer,n);done+=n;spillReads+=n;spillWrites+=n;}
                    Require(digest==row.digest,"spill checksum changed during compaction");end+=row.size;}
                next->Put(path,row);
            });
            ::close(spill);spill=replacement.Release();spillEnd=end;catalog.swap(next);return;
        }
        Descriptor replacement(NewSpill());uint64_t end=0;std::map<std::string,uint64_t> offsets;char buffer[65536];
        for(const auto& item:entries){const auto& e=item.second;if(!e.spilled)continue;offsets[item.first]=end;uint64_t done=0,digest=UINT64_C(14695981039346656037);
            while(done<e.size){const size_t n=size_t(std::min<uint64_t>(sizeof(buffer),e.size-done));ReadAt(spill,e.offset+done,buffer,n);digest=HashBytes(buffer,n,digest);WriteAll(replacement.fd,buffer,n);done+=n;spillReads+=n;spillWrites+=n;}
            Require(digest==e.digest,"spill checksum changed during compaction");end+=e.size;
        }
        ::close(spill);spill=replacement.Release();spillEnd=end;for(const auto& item:offsets)entries.at(item.first).offset=item.second;
    }
    struct Writer {
        int fd;uint64_t digest=UINT64_C(14695981039346656037);
        char buffer[65536];size_t used=0;
        explicit Writer(int value):fd(value){}
        void Flush(){if(used){WriteAll(fd,buffer,used);used=0;}}
        void Raw(const char* data,size_t size){
            while(size){if(!used&&size>=sizeof(buffer)){WriteAll(fd,data,size);return;}
                const size_t n=std::min(size,sizeof(buffer)-used);std::memcpy(buffer+used,data,n);used+=n;data+=n;size-=n;if(used==sizeof(buffer))Flush();}
        }
        void Bytes(const char* data,size_t size){digest=HashBytes(data,size,digest);Raw(data,size);}
        void Number(uint64_t value){char bytes[8];for(unsigned i=0;i<8;++i)bytes[i]=char((value>>(8*i))&255);Bytes(bytes,8);}
        void Text(const std::string& value){Number(value.size());Bytes(value.data(),value.size());}
        void Finish(){const auto saved=digest;char bytes[8];for(unsigned i=0;i<8;++i)bytes[i]=char((saved>>(8*i))&255);Raw(bytes,8);Flush();}
    };
    struct Reader {
        int fd;uint64_t offset=0,length,digest=UINT64_C(14695981039346656037);
        char buffer[65536];size_t cursor=0,available=0;
        Reader(int value,uint64_t bytes):fd(value),length(bytes){}
        void Bytes(char* data,size_t size){Require(offset<=length&&size<=length-offset,"truncated checkpoint");char* const original=data;const size_t total=size;
            while(size){
                if(cursor==available){cursor=available=0;
                    if(size>=sizeof(buffer)){ReadAt(fd,offset,data,size);offset+=size;break;}
                    available=size_t(std::min<uint64_t>(sizeof(buffer),length-offset));ReadAt(fd,offset,buffer,available);}
                const size_t n=std::min(size,available-cursor);std::memcpy(data,buffer+cursor,n);cursor+=n;offset+=n;data+=n;size-=n;
            }
            digest=HashBytes(original,total,digest);}
        uint64_t Number(){unsigned char bytes[8];Bytes(reinterpret_cast<char*>(bytes),8);uint64_t value=0;for(unsigned i=0;i<8;++i)value|=uint64_t(bytes[i])<<(8*i);return value;}
        std::string Text(size_t maximum){const auto size=Number();Require(size<=maximum&&offset<=length&&size<=length-offset,"checkpoint string bound exceeded");std::string value(size,'\0');if(size)Bytes(&value[0],size);return value;}
        void Finish(){Require(offset<=length&&length-offset==8,"checkpoint has trailing or missing bytes");const auto expected=digest;const auto actual=Number();Require(actual==expected,"checkpoint checksum mismatch");}
    };
public:
    explicit Store(const std::string& ownedRoot,size_t limitBytes,const std::string& temporaryDirectory="/tmp",size_t metadataBudget=0)
        :root(Canonical(ownedRoot)),spillDirectory(Canonical(temporaryDirectory)),limit(limitBytes),metadataLimit(metadataBudget?std::min(limitBytes,metadataBudget):limitBytes){
        Require(limit>=Charge(root),"working-memory limit too small");AddDirectory(root);
    }
    ~Store(){if(spill>=0)::close(spill);}
    Store(const Store&)=delete;Store& operator=(const Store&)=delete;
    const std::string& Root()const{return root;}
    size_t Limit()const{return limit;}
    bool Healthy()const{std::lock_guard<std::recursive_mutex> guard(mutex);return !poisoned&&(!catalog||catalog->Healthy());}
    size_t MemoryBytes()const{std::lock_guard<std::recursive_mutex> guard(mutex);return metadataBytes+residentBytes;}
    uint64_t Epoch()const{std::lock_guard<std::recursive_mutex> guard(mutex);return epoch;}
    bool Owns(const std::string& path)const{try{Key(path);return true;}catch(const std::exception&){return false;}}
    bool Exists(const std::string& path)const{std::lock_guard<std::recursive_mutex> guard(mutex);RequireHealthy();FullAutoMemoryCatalog::Row row;return Lookup(Key(path),row);}
    bool IsDirectory(const std::string& path)const{std::lock_guard<std::recursive_mutex> guard(mutex);RequireHealthy();FullAutoMemoryCatalog::Row row;return Lookup(Key(path),row)&&row.directory;}
    uint64_t Generation(const std::string& path)const{std::lock_guard<std::recursive_mutex> guard(mutex);RequireHealthy();FullAutoMemoryCatalog::Row row;Require(Lookup(Key(path),row),"missing generation");return row.generation;}
    uint64_t Digest(const std::string& path)const{std::lock_guard<std::recursive_mutex> guard(mutex);RequireHealthy();FullAutoMemoryCatalog::Row row;Require(Lookup(Key(path),row)&&!row.directory,"missing data digest");return row.digest;}
    uint64_t Size(const std::string& path)const{std::lock_guard<std::recursive_mutex> guard(mutex);RequireHealthy();FullAutoMemoryCatalog::Row row;Require(Lookup(Key(path),row)&&!row.directory,"missing data size");return row.size;}
    void Directory(const std::string& path){
        std::lock_guard<std::recursive_mutex> guard(mutex);RequireHealthy();Mutation mutation(*this);const auto key=Key(path);if(catalog){DirectoryCatalog(key);return;}const auto found=entries.find(key);
        if(found!=entries.end()){Require(found->second.directory,"directory path is a file");return;}
        const auto parents=MissingParents(key);Preflight(key,parents);if(catalog){DirectoryCatalog(key);return;}for(const auto& p:parents)AddDirectory(p);AddDirectory(key);
    }
    void Write(const std::string& path,const std::string& value){
        std::lock_guard<std::recursive_mutex> guard(mutex);RequireHealthy();Mutation mutation(*this);const auto key=Key(path);if(catalog){WriteCatalog(key,value);return;}const auto old=entries.find(key);
        Require(old==entries.end()||!old->second.directory,"cannot replace directory with a file");
        const auto parents=MissingParents(key);Preflight(key,parents);if(catalog){WriteCatalog(key,value);return;}
        // Prepare bytes/spill before replacing a valid value. A failed disk
        // allocation never turns a previously valid matrix into partial data.
        Entry replacement;replacement.size=value.size();replacement.digest=Hash(value);replacement.metadata=Charge(key);
        const size_t additional=old==entries.end()?replacement.metadata:0;size_t parentBytes=0;for(const auto& p:parents)parentBytes+=Charge(p);
        const size_t available=limit-metadataBytes-additional-parentBytes;
        if(value.size()+1<=available){replacement.bytes=value;replacement.resident=ResidentCharge(replacement.bytes)<=available;}
        if(!replacement.resident){std::string().swap(replacement.bytes);OpenSpill();Require(value.size()<=uint64_t(std::numeric_limits<off_t>::max())-spillEnd,"spill size overflow");
            Require(::lseek(spill,off_t(spillEnd),SEEK_SET)==off_t(spillEnd),"cannot seek spill");WriteAll(spill,value.data(),value.size());
            replacement.spilled=true;replacement.offset=spillEnd;spillEnd+=value.size();spillWrites+=value.size();}
        if(old!=entries.end())Drop(old);
        for(const auto& p:parents)AddDirectory(p);
        replacement.generation=Next();metadataBytes+=replacement.metadata;
        if(replacement.spilled)liveSpill+=replacement.size;
        if(replacement.resident){Room(ResidentCharge(replacement.bytes));recent.push_front(key);replacement.recent=recent.begin();residentBytes+=ResidentCharge(replacement.bytes);}
        entries.emplace(key,std::move(replacement));++fileCount;Compact();
    }
    // Reader(dst, capacity) returns 1..capacity bytes, or zero on EOF. Exactly
    // size bytes belong to this frame: never probe the following frame for EOF.
    // The callback must not reenter this Store. Resident staging is charged to
    // the cache and moved into its Entry, never copied as a second full value.
    // Only the fixed spill transfer buffer is outside the retained cache.
    template<class Reader> void WriteStream(const std::string& path,uint64_t size,uint64_t expectedDigest,Reader reader){
        std::lock_guard<std::recursive_mutex> guard(mutex);RequireHealthy();const auto key=Key(path);
        FullAutoMemoryCatalog::Row old;const bool existed=Lookup(key,old);Require(!existed||!old.directory,"cannot replace directory with a file");
        const auto parents=MissingParents(key);
        Require(uint64_t(parents.size())<UINT64_MAX-epoch,"generation overflow");
        // Admission and metadata spill happen before consuming the frame. They
        // may move old bytes, but do not publish a new logical file/generation.
        Preflight(key,parents);Require(!existed||Lookup(key,old),"overwritten stream entry disappeared");
        size_t metadataRequired=catalog?Charge(key):(existed?0:Charge(key));
        if(!catalog)for(const auto& parent:parents)metadataRequired+=Charge(parent);
        const size_t fixedMetadata=catalog?CatalogReservation():metadataBytes;
        if(metadataRequired<=limit-fixedMetadata&&size<limit-fixedMetadata-metadataRequired){
            Room(metadataRequired+size_t(size)+1);
            std::string value(size_t(size),'\0');const size_t charge=ResidentCharge(value);
            if(charge<=limit-fixedMetadata-metadataRequired){
                Room(metadataRequired+charge);StagedResident staging(*this,charge);uint64_t done=0,digest=UINT64_C(14695981039346656037);
                {StreamCallback callback(*this);while(done<size){const auto capacity=size_t(std::min<uint64_t>(65536,size-done));const size_t n=reader(&value[size_t(done)],capacity);
                    Require(n>0&&n<=capacity,"truncated stream or invalid reader count");digest=HashBytes(value.data()+size_t(done),n,digest);done+=n;}}
                Require(digest==expectedDigest,"stream checksum mismatch");MBS_MEMORY_STORE_EVENT("before_stream_publish",key);
                // Room() may have evicted the previous resident value and
                // assigned it a spill offset. Refresh its authoritative row.
                Require(!existed||Lookup(key,old),"overwritten stream entry disappeared");
                Mutation mutation(*this,true);for(const auto& parent:parents)AddDirectory(parent);
                Entry replacement;replacement.size=size;replacement.digest=digest;replacement.generation=Next();replacement.resident=true;
                replacement.metadata=Charge(key);replacement.bytes.swap(value);
                if(catalog){catalog->Put(key,Row(replacement));if(existed&&Spilled(old))liveSpill-=old.size;DropResident(key);if(!existed)++fileCount;}
                else{const auto found=entries.find(key);if(found!=entries.end())Drop(found);++fileCount;}
                recent.push_front(key);replacement.recent=recent.begin();metadataBytes+=replacement.metadata;
                entries.emplace(key,std::move(replacement));staging.transferred=true;Compact();return;
            }
        }
        Require(!existed||Lookup(key,old),"overwritten stream entry disappeared");
        if(catalog&&existed&&Spilled(old))Require(old.offset<=spillEnd&&old.size<=spillEnd-old.offset&&old.size<=liveSpill,"invalid overwritten catalog range");
        OpenSpill();Require(size<=uint64_t(std::numeric_limits<off_t>::max())-spillEnd,"spill size overflow");
        const uint64_t start=spillEnd;Require(::lseek(spill,off_t(start),SEEK_SET)==off_t(start),"cannot seek streamed spill");
        uint64_t done=0,digest=UINT64_C(14695981039346656037);char buffer[65536];
        try{
            {StreamCallback callback(*this);while(done<size){const auto capacity=size_t(std::min<uint64_t>(sizeof(buffer),size-done));const size_t n=reader(buffer,capacity);
                Require(n>0&&n<=capacity,"truncated stream or invalid reader count");digest=HashBytes(buffer,n,digest);WriteAll(spill,buffer,n);done+=n;spillWrites+=n;}}
            Require(digest==expectedDigest,"stream checksum mismatch");MBS_MEMORY_STORE_EVENT("before_stream_publish",key);
        }catch(...){
            // Nothing references this append yet. Restore its physical length
            // where possible; even a failed truncate cannot expose its bytes,
            // and the next append overwrites from the unchanged spillEnd.
            int rc;do{rc=::ftruncate(spill,off_t(start));}while(rc<0&&errno==EINTR);
            throw;
        }
        // From here a failed metadata mutation is fail-closed, including the
        // map backend: no checkpoint may publish an incomplete installation.
        Mutation mutation(*this,true);spillEnd=start+size;
        for(const auto& parent:parents)AddDirectory(parent);
        FullAutoMemoryCatalog::Row row;row.size=size;row.digest=digest;row.offset=start;row.generation=Next();
        if(catalog){catalog->Put(key,row);liveSpill+=size;if(existed){if(Spilled(old))liveSpill-=old.size;}else ++fileCount;DropResident(key);}
        else{auto found=entries.find(key);if(found!=entries.end())Drop(found);Entry replacement;replacement.size=size;replacement.digest=digest;replacement.offset=start;
            replacement.generation=row.generation;replacement.spilled=true;replacement.metadata=Charge(key);metadataBytes+=replacement.metadata;
            entries.emplace(key,std::move(replacement));liveSpill+=size;++fileCount;}
        Compact();
    }
    // Writer(data, size) consumes the entire chunk or throws. A receiver must
    // stage its value until the advertised Digest() is verified, since a bad
    // spill checksum is reported after the last streamed chunk. No promotion
    // to the resident cache or whole-value allocation is performed.
    template<class Writer> void ReadStream(const std::string& path,Writer writer){
        std::lock_guard<std::recursive_mutex> guard(mutex);RequireHealthy();const auto key=Key(path);
        FullAutoMemoryCatalog::Row row;Require(Lookup(key,row)&&!row.directory,"missing file "+key);auto found=entries.find(key);
        const bool resident=found!=entries.end()&&found->second.resident;
        if(resident){Require(found->second.bytes.size()==row.size,"invalid resident payload size");Touch(found->second);}
        else Require(spill>=0&&row.offset<=spillEnd&&row.size<=spillEnd-row.offset,"invalid streamed spill range");
        StreamCallback callback(*this);uint64_t done=0,digest=UINT64_C(14695981039346656037);char buffer[65536];
        while(done<row.size){const auto n=size_t(std::min<uint64_t>(sizeof(buffer),row.size-done));const char* data=nullptr;
            if(resident)data=found->second.bytes.data()+size_t(done);
            else{ReadAt(spill,row.offset+done,buffer,n);spillReads+=n;data=buffer;}
            digest=HashBytes(data,n,digest);writer(data,n);done+=n;
        }
        Require(digest==row.digest,"spill checksum mismatch during streaming read");
    }
    std::string Read(const std::string& path){
        std::lock_guard<std::recursive_mutex> guard(mutex);RequireHealthy();const auto key=Key(path);auto found=entries.find(key);
        if(catalog&&found==entries.end()){
            FullAutoMemoryCatalog::Row row;Require(catalog->Get(key,row)&&!row.directory,"missing file "+key);Require(row.size<=std::numeric_limits<size_t>::max()&&row.offset<=spillEnd&&row.size<=spillEnd-row.offset,"invalid catalog payload range");
            std::string value(size_t(row.size),'\0');if(row.size)ReadAt(spill,row.offset,&value[0],value.size());spillReads+=row.size;Require(Hash(value)==row.digest,"spill checksum mismatch");CacheCatalog(key,row,value);return value;
        }
        Require(found!=entries.end(),"missing file "+key);auto& e=found->second;
        Require(!e.directory,"cannot read directory");if(e.resident){Touch(e);return e.bytes;}
        Require(e.spilled&&e.size<=std::numeric_limits<size_t>::max(),"invalid spill record");std::string value(size_t(e.size),'\0');
        if(e.size)ReadAt(spill,e.offset,&value[0],value.size());
        spillReads+=e.size;Require(Hash(value)==e.digest,"spill checksum mismatch");
        if(ResidentCharge(value)<=limit-metadataBytes){Room(ResidentCharge(value));e.bytes=value;e.resident=true;recent.push_front(key);e.recent=recent.begin();residentBytes+=ResidentCharge(e.bytes);}
        return value;
    }
    void Remove(const std::string& path,bool recursive=false){
        std::lock_guard<std::recursive_mutex> guard(mutex);RequireHealthy();Mutation mutation(*this);const auto key=Key(path);Require(key!=root,"cannot remove owned root");if(catalog){RemoveCatalog(key,recursive);return;}auto found=entries.find(key);if(found==entries.end())return;
        auto children=List(key,true);Require(recursive||children.empty(),"directory is not empty");for(const auto& child:children)Drop(entries.find(child));Drop(entries.find(key));Next();Compact();
    }
    void Link(const std::string& from,const std::string& to){
        std::lock_guard<std::recursive_mutex> guard(mutex);Require(!Exists(to),"link destination exists");Write(to,Read(from));
    }
    void Rename(const std::string& from,const std::string& to){
        std::lock_guard<std::recursive_mutex> guard(mutex);RequireHealthy();Mutation mutation(*this);const auto a=Key(from),b=Key(to);if(a==b)return;Require(a!=root&&b!=root,"cannot rename owned root");
        if(catalog){FullAutoMemoryCatalog::Row source;Require(Lookup(a,source),"rename source missing");if(!source.directory){Write(b,Read(a));Remove(a);return;}
            Require(!Exists(b)&&b.compare(0,a.size()+1,a+"/")!=0,"invalid directory rename destination");const auto parents=MissingParents(b);for(const auto& path:parents)AddDirectory(path);
            auto paths=List(a,true);paths.insert(paths.begin(),a);for(const auto& path:paths){const auto cached=entries.find(path);if(cached!=entries.end()&&!cached->second.directory)SaveSpill(cached->second,path);
                FullAutoMemoryCatalog::Row row;Require(Lookup(path,row),"rename source missing");const auto destination=b+path.substr(a.size());Require(!Exists(destination),"rename destination exists");
                row.generation=Next();catalog->Put(destination,row);catalog->Erase(path);DropResident(path);}return;}
        const auto found=entries.find(a);Require(found!=entries.end(),"rename source missing");
        if(!found->second.directory){Write(b,Read(a));Remove(a);return;}
        Require(!Exists(b)&&b.compare(0,a.size()+1,a+"/")!=0,"invalid directory rename destination");
        auto paths=List(a,true);paths.insert(paths.begin(),a);const auto parents=MissingParents(b);size_t oldCharge=0,newCharge=0;
        for(const auto& path:paths){oldCharge+=entries.at(path).metadata;const auto destination=b+path.substr(a.size());Require(!Exists(destination),"directory rename destination exists");
            Require(Charge(destination)<=limit-newCharge,"renamed metadata exceeds working-memory limit");newCharge+=Charge(destination);}
        for(const auto& p:parents){Require(Charge(p)<=limit-newCharge,"renamed metadata exceeds working-memory limit");newCharge+=Charge(p);}
        if(newCharge>oldCharge)Room(newCharge-oldCharge);
        for(const auto& p:parents)AddDirectory(p);
        for(const auto& path:paths){auto old=entries.find(path);const auto destination=b+path.substr(a.size());Entry moved=std::move(old->second);entries.erase(old);
            metadataBytes-=moved.metadata;moved.metadata=Charge(destination);metadataBytes+=moved.metadata;moved.generation=Next();if(moved.resident)*moved.recent=destination;entries.emplace(destination,std::move(moved));}
    }
    std::vector<std::string> List(const std::string& path,bool recursive=true)const{
        std::lock_guard<std::recursive_mutex> guard(mutex);RequireHealthy();const auto key=Key(path);const auto prefix=key=="/"?key:key+"/";std::vector<std::string> result;
        if(catalog){catalog->Children(key,recursive,[&](const std::string& name,const FullAutoMemoryCatalog::Row&){result.push_back(name);});return result;}
        for(auto it=entries.lower_bound(prefix);it!=entries.end()&&it->first.compare(0,prefix.size(),prefix)==0;++it)
            if(it->first!=key&&(recursive||it->first.find('/',prefix.size())==std::string::npos))result.push_back(it->first);
        return result;
    }
    Statistics Stats()const{
        std::lock_guard<std::recursive_mutex> guard(mutex);Statistics result;result.residentBytes=residentBytes;result.metadataBytes=metadataBytes;result.limitBytes=limit;result.metadataLimitBytes=metadataLimit;
        result.spillBytes=spillEnd;result.liveSpillBytes=liveSpill;result.spillReadBytes=spillReads;result.snapshotReadBytes=snapshotReads;result.spillWriteBytes=spillWrites;result.evictions=evictions;result.snapshots=snapshots;result.epoch=epoch;
        result.files=fileCount;result.directories=directoryCount;result.catalogActive=bool(catalog);result.healthy=Healthy();result.catalogBytes=catalog?catalog->DiskBytes():0;result.catalogCacheBytes=catalog?CatalogCache():0;
        return result;
    }
    void Snapshot(const std::string& path,const std::string& metadata){
        std::lock_guard<std::recursive_mutex> guard(mutex);RequireHealthy();Require(Canonical(path)==path,"checkpoint path must be canonical and absolute");Require(metadata.size()<=16*1024*1024,"checkpoint metadata too large");
        std::string name=path+".tmp.XXXXXX";std::vector<char> temp(name.begin(),name.end());temp.push_back('\0');Descriptor fd(::mkstemp(temp.data()));Require(fd.fd>=0,"cannot create checkpoint temporary file");
        ::fcntl(fd.fd,F_SETFD,FD_CLOEXEC);
        try{
            Writer writer(fd.fd);writer.Text("MBS_FULLAUTO_MEMORY_STORE_V1");writer.Text(root);writer.Text(metadata);writer.Number(epoch);writer.Number(fileCount+directoryCount);char buffer[65536];
            auto save=[&](const std::string& path,const FullAutoMemoryCatalog::Row& row,const Entry* resident){
                Require(Key(path)==path&&row.generation>0&&row.generation<=epoch,"invalid working catalog record");
                writer.Text(Relative(path));writer.Number(row.directory?1:0);writer.Number(row.generation);writer.Number(row.size);writer.Number(row.digest);
                if(row.directory){Require(row.size==0&&row.digest==Hash(""),"invalid working directory data");return;}
                if(resident&&resident->resident){writer.Bytes(resident->bytes.data(),resident->bytes.size());return;}
                Require(row.offset<=spillEnd&&row.size<=spillEnd-row.offset,"invalid spill range in working catalog");
                uint64_t done=0,digest=UINT64_C(14695981039346656037);while(done<row.size){const auto n=size_t(std::min<uint64_t>(sizeof(buffer),row.size-done));ReadAt(spill,row.offset+done,buffer,n);digest=HashBytes(buffer,n,digest);writer.Bytes(buffer,n);done+=n;spillReads+=n;snapshotReads+=n;}
                Require(digest==row.digest,"spill checksum mismatch while saving checkpoint");
            };
            if(catalog)catalog->All([&](const std::string& path,const FullAutoMemoryCatalog::Row& row){const auto cached=entries.find(path);save(path,row,cached==entries.end()?nullptr:&cached->second);});
            else for(const auto& item:entries)save(item.first,Row(item.second),&item.second);
            writer.Finish();MBS_MEMORY_STORE_EVENT("before_snapshot_fsync",path);Sync(fd.fd);MBS_MEMORY_STORE_EVENT("before_snapshot_publish",path);Require(::rename(temp.data(),path.c_str())==0,"cannot publish checkpoint");
            Descriptor directory(::open(Parent(Canonical(path)).c_str(),O_RDONLY|O_DIRECTORY));Require(directory.fd>=0,"cannot open checkpoint directory");Sync(directory.fd);++snapshots;
        }catch(...){::unlink(temp.data());throw;}
    }
    std::string Restore(const std::string& path){
        std::lock_guard<std::recursive_mutex> guard(mutex);Require(!streamCallbackActive,"stream callbacks cannot reenter this store");Descriptor fd(::open(path.c_str(),O_RDONLY|O_CLOEXEC|O_NOFOLLOW));Require(fd.fd>=0,"cannot open checkpoint");struct stat st{};
        Require(::fstat(fd.fd,&st)==0&&S_ISREG(st.st_mode)&&st.st_size>=8,"invalid checkpoint file");Reader reader(fd.fd,uint64_t(st.st_size));
        Require(reader.Text(64)=="MBS_FULLAUTO_MEMORY_STORE_V1","unsupported checkpoint schema");Require(reader.Text(4096)==root,"checkpoint root differs");const auto metadata=reader.Text(16*1024*1024);const auto savedEpoch=reader.Number(),count=reader.Number();
        Require(count>=1&&count<=uint64_t(st.st_size)/40,"invalid checkpoint record count");
        Store staged(root,limit,spillDirectory,metadataLimit);staged.entries.clear();staged.metadataBytes=0;staged.epoch=0;staged.directoryCount=staged.fileCount=0;std::string previous;char buffer[65536];
        for(uint64_t i=0;i<count;++i){const auto relative=reader.Text(4096);Require(relative.empty()||relative[0]=='/',"invalid checkpoint relative path");const auto absolute=staged.Absolute(relative);const auto key=staged.Key(absolute);Require(key==absolute&&(i==0?key==root:key>previous),"checkpoint path order or spelling changed");previous=key;
            const auto directory=reader.Number();Entry e;e.generation=reader.Number();e.size=reader.Number();e.digest=reader.Number();e.directory=directory==1;e.metadata=Charge(key);
            Require(directory<=1&&e.generation>0&&e.generation<=savedEpoch,"invalid checkpoint generation");Require(e.size<=uint64_t(st.st_size)-std::min(reader.offset,uint64_t(st.st_size)),"invalid checkpoint payload length");
            if(i==0)Require(e.directory,"checkpoint root is not a directory");else{FullAutoMemoryCatalog::Row parent;Require(staged.Lookup(Parent(key),parent)&&parent.directory,"checkpoint parent missing");}
            staged.Preflight(key,{});if(!staged.catalog)staged.metadataBytes+=e.metadata;
            if(e.directory){Require(e.size==0&&e.digest==Hash(""),"invalid checkpoint directory data");++staged.directoryCount;}
            else if(!staged.catalog&&e.size+1<=staged.limit-staged.metadataBytes){e.bytes.resize(size_t(e.size));if(e.size)reader.Bytes(&e.bytes[0],e.bytes.size());Require(Hash(e.bytes)==e.digest,"checkpoint member checksum mismatch");e.resident=true;
                staged.Room(ResidentCharge(e.bytes));staged.recent.push_front(key);e.recent=staged.recent.begin();staged.residentBytes+=ResidentCharge(e.bytes);}
            else{staged.OpenSpill();Require(::lseek(staged.spill,off_t(staged.spillEnd),SEEK_SET)==off_t(staged.spillEnd),"cannot seek restored spill");e.offset=staged.spillEnd;e.spilled=true;uint64_t done=0,digest=UINT64_C(14695981039346656037);
                while(done<e.size){const auto n=size_t(std::min<uint64_t>(sizeof(buffer),e.size-done));reader.Bytes(buffer,n);digest=HashBytes(buffer,n,digest);WriteAll(staged.spill,buffer,n);done+=n;}
                Require(digest==e.digest,"checkpoint member checksum mismatch");staged.spillEnd+=e.size;staged.liveSpill+=e.size;staged.spillWrites+=e.size;
            }
            if(!e.directory)++staged.fileCount;
            if(staged.catalog)staged.catalog->Put(key,Row(e));else staged.entries.emplace(key,std::move(e));
        }
        reader.Finish();staged.epoch=savedEpoch;
        // Swapping both containers preserves Entry's list iterators. A failed
        // restore never changes the existing live computation.
        entries.swap(staged.entries);recent.swap(staged.recent);catalog.swap(staged.catalog);std::swap(spill,staged.spill);std::swap(spillEnd,staged.spillEnd);std::swap(liveSpill,staged.liveSpill);
        std::swap(residentBytes,staged.residentBytes);std::swap(metadataBytes,staged.metadataBytes);std::swap(epoch,staged.epoch);spillReads=0;snapshotReads=0;spillWrites=staged.spillWrites;evictions=0;
        std::swap(fileCount,staged.fileCount);std::swap(directoryCount,staged.directoryCount);
        std::swap(poisoned,staged.poisoned);
        return metadata;
    }
};
inline Store*& Active(){static Store* value=nullptr;return value;}
struct Scope {
    Store* previous;
    explicit Scope(Store& value):previous(Active()){Active()=&value;}
    ~Scope(){Active()=previous;}
    Scope(const Scope&)=delete;Scope& operator=(const Scope&)=delete;
};
}
