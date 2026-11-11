#include <string>
#include <stdexcept>
#include <unistd.h>
static std::string failEvent,crashEvent;
static std::string failCatalogEvent,failCatalogPath;
static int lastSpillFd=-1;
static void MemoryStoreEvent(const char* name,const std::string&){if(crashEvent==name)_exit(91);if(failEvent==name)throw std::runtime_error("injected checkpoint failure");}
static void MemoryCatalogEvent(const char* name,const std::string& path){if(failCatalogEvent==name&&failCatalogPath==path)throw std::runtime_error("injected catalog mutation failure");}
#define MBS_MEMORY_STORE_EVENT(name,path) MemoryStoreEvent(name,path)
#define MBS_MEMORY_STORE_SPILL_OPEN(fd) (lastSpillFd=(fd))
#define MBS_MEMORY_CATALOG_EVENT(name,path) MemoryCatalogEvent(name,path)
#include "../src/FullAutoMemoryStore.h"
#include <cassert>
#include <dirent.h>
#include <fstream>
#include <iostream>
#include <thread>
#include <chrono>
#include <sys/wait.h>
using FullAutoMemoryStore::Store;
template<class F> static bool Throws(F f){try{f();return false;}catch(const std::exception&){return true;}}
static std::string Bytes(const std::string& path){std::ifstream stream(path,std::ios::binary);return std::string(std::istreambuf_iterator<char>(stream),std::istreambuf_iterator<char>());}
static void Bytes(const std::string& path,const std::string& value){std::ofstream stream(path,std::ios::binary);stream.write(value.data(),value.size());assert(stream.good());}
static size_t CountFiles(const std::string& path){DIR* dir=opendir(path.c_str());assert(dir);size_t n=0;while(auto* entry=readdir(dir))if(std::string(entry->d_name)!="."&&std::string(entry->d_name)!="..")++n;closedir(dir);return n;}
static size_t RssBytes(){std::ifstream stream("/proc/self/status");std::string line;while(std::getline(stream,line))if(line.compare(0,6,"VmRSS:")==0)return size_t(std::stoull(line.substr(6)))*1024;return 0;}
static char Pattern(uint64_t offset){return char((offset*17+offset/251)%251);}
static uint64_t PatternHash(uint64_t size){char block[65536];uint64_t digest=UINT64_C(14695981039346656037),done=0;
    while(done<size){const auto n=size_t(std::min<uint64_t>(sizeof(block),size-done));for(size_t i=0;i<n;++i)block[i]=Pattern(done+i);digest=FullAutoMemoryStore::HashBytes(block,n,digest);done+=n;}return digest;}
static void StreamChecks(const std::string& root,const std::string& disk,const std::string& snapshot,bool catalog){
    Store streamed(root,8*1024*1024,disk);
    if(catalog){for(int i=0;i<20000;++i)streamed.Write(root+"/metadata"+std::to_string(i),"x");assert(streamed.Stats().catalogActive);}
    int workingSpill=catalog?lastSpillFd:-1;
    auto patternWrite=[&](const std::string& path,uint64_t bytes,uint64_t expected){uint64_t cursor=0;
        streamed.WriteStream(path,bytes,expected,[&](char* data,size_t capacity)->size_t{
            for(size_t i=0;i<capacity;++i)data[i]=Pattern(cursor+i);
            cursor+=capacity;assert(streamed.MemoryBytes()<=streamed.Limit());return capacity;});};
    auto patternRead=[&](const std::string& path,uint64_t bytes){uint64_t cursor=0;
        streamed.ReadStream(path,[&](const char* data,size_t n){for(size_t i=0;i<n;++i)assert(data[i]==Pattern(cursor+i));cursor+=n;});assert(cursor==bytes);};
    const uint64_t residentSize=1024*1024,residentHash=PatternHash(residentSize);const auto writesBefore=streamed.Stats().spillWriteBytes;
    patternWrite(root+"/resident/data",residentSize,residentHash);
    assert(streamed.Stats().spillWriteBytes==writesBefore&&streamed.Stats().residentBytes>=residentSize);
    patternRead(root+"/resident/data",residentSize);const auto residentGeneration=streamed.Generation(root+"/resident/data");
    assert(Throws([&]{patternWrite(root+"/resident/data",residentSize,residentHash^1);}));
    assert(streamed.Healthy()&&streamed.Generation(root+"/resident/data")==residentGeneration&&streamed.Stats().spillWriteBytes==writesBefore);
    patternRead(root+"/resident/data",residentSize);
    // Catalog snapshots include resident-only rows without first forcing them
    // into spill. Restore still authenticates the same portable V1 payload.
    streamed.Snapshot(snapshot,"resident stream");assert(streamed.Stats().spillWriteBytes==writesBefore);
    {Store check(root,8*1024*1024,disk);assert(check.Restore(snapshot)=="resident stream"&&check.Digest(root+"/resident/data")==residentHash);}
    const uint64_t part=512*1024,partHash=PatternHash(part);
    for(int i=0;i<22;++i)patternWrite(root+"/resident/f"+std::to_string(i),part,partHash);
    if(workingSpill<0)workingSpill=lastSpillFd;
    assert(streamed.Stats().spillWriteBytes>writesBefore&&streamed.MemoryBytes()<=streamed.Limit());
    for(int i=0;i<22;++i)patternRead(root+"/resident/f"+std::to_string((i*7)%22),part);
    streamed.Rename(root+"/resident",root+"/moved_resident");patternRead(root+"/moved_resident/f21",part);
    patternWrite(root+"/overwrite",residentSize,residentHash);streamed.Write(root+"/overwrite","ordinary replacement");assert(streamed.Read(root+"/overwrite")=="ordinary replacement");streamed.Remove(root+"/overwrite");
    patternWrite(root+"/removed",residentSize,residentHash);streamed.Remove(root+"/removed");assert(!streamed.Exists(root+"/removed"));
    streamed.Remove(root+"/moved_resident",true);
    const std::string path=root+"/large/data";const uint64_t size=UINT64_C(80)*1024*1024+17,digest=PatternHash(size);
    streamed.Write(path,"old value");const auto generation=streamed.Generation(path),epoch=streamed.Epoch();const auto baselineRss=RssBytes();uint64_t position=0;size_t largestChunk=0;
    auto reader=[&](char* data,size_t capacity)->size_t{assert(capacity<=65536);largestChunk=std::max(largestChunk,capacity);const auto n=std::min(capacity,size_t(37013));
        for(size_t i=0;i<n;++i)data[i]=Pattern(position+i);
        position+=n;assert(streamed.MemoryBytes()<=streamed.Limit());return n;};
    streamed.WriteStream(path,size,digest,reader);assert(position==size&&largestChunk==65536&&streamed.Generation(path)>generation&&streamed.Epoch()>epoch);
    assert(streamed.Size(path)==size&&streamed.Digest(path)==digest&&streamed.MemoryBytes()<=streamed.Limit());
#if !defined(__SANITIZE_ADDRESS__)
    assert(RssBytes()<=baselineRss+8*1024*1024); // payload stays out of process RAM
#endif
    auto verify=[&](Store& source){uint64_t seen=0;source.ReadStream(path,[&](const char* data,size_t n){assert(n<=65536);for(size_t i=0;i<n;++i)assert(data[i]==Pattern(seen+i));seen+=n;});assert(seen==size);};
    verify(streamed);const auto savedGeneration=streamed.Generation(path),savedEpoch=streamed.Epoch(),savedBytes=streamed.Stats().spillBytes;
    // A broken frame, bad checksum, caller exception or publication hook
    // leaves the previous generation usable, even with a SQLite catalog.
    for(int failure=0;failure<5;++failure){position=0;if(failure==4)failEvent="before_stream_publish";
        assert(Throws([&]{streamed.WriteStream(path,size,digest^(failure==1?1:0),[&](char* data,size_t capacity)->size_t{
            if(failure==0&&position>100000)return 0;
            if(failure==2&&position>100000)throw std::runtime_error("reader failure");
            if(failure==3)return capacity+1;
            return reader(data,capacity);
        });}));failEvent.clear();assert(streamed.Healthy()&&streamed.Generation(path)==savedGeneration&&streamed.Epoch()==savedEpoch);
        assert(streamed.Size(path)==size&&streamed.Digest(path)==digest&&streamed.Stats().spillBytes==savedBytes);
        struct stat st{};assert(fstat(workingSpill,&st)==0&&uint64_t(st.st_size)==savedBytes);
    }
    assert(Throws([&]{streamed.ReadStream(path,[&](const char*,size_t){throw std::runtime_error("writer failure");});}));assert(streamed.Healthy());
    assert(Throws([&]{streamed.WriteStream(path,100000,digest,[&](char*,size_t)->size_t{streamed.Write(root+"/forbidden","x");return 0;});}));
    assert(streamed.Healthy()&&!streamed.Exists(root+"/forbidden")&&streamed.Generation(path)==savedGeneration);
    assert(Throws([&]{streamed.WriteStream(path,UINT64_MAX,0,[](char*,size_t)->size_t{assert(false);return 0;});}));assert(streamed.Healthy());
    assert(Throws([&]{streamed.Size(root);}));assert(Throws([&]{streamed.WriteStream(root,0,FullAutoMemoryStore::Hash(""),[](char*,size_t)->size_t{assert(false);return 0;});}));assert(streamed.Healthy());
    // Small and empty frames preserve the existing resident Write fast path.
    const std::string small("small\0binary",12);position=0;
    streamed.WriteStream(root+"/small",small.size(),FullAutoMemoryStore::Hash(small),[&](char* data,size_t capacity)->size_t{const size_t n=std::min(size_t(2),capacity);std::memcpy(data,small.data()+position,n);position+=n;return n;});assert(streamed.Read(root+"/small")==small);
    streamed.WriteStream(root+"/empty",0,FullAutoMemoryStore::Hash(""),[](char*,size_t)->size_t{assert(false);return 0;});assert(streamed.Size(root+"/empty")==0);
    assert(Throws([&]{streamed.WriteStream(root+"/absent/new",10,0,[](char*,size_t)->size_t{return 0;});}));assert(streamed.Healthy()&&!streamed.Exists(root+"/absent"));
    streamed.Snapshot(snapshot,"large streamed frame");Store restored(root,8*1024*1024,disk);assert(restored.Restore(snapshot)=="large streamed frame");verify(restored);assert(restored.Size(path)==size&&restored.Digest(path)==digest);
    // Detect actual spill corruption during bounded reads and snapshotting.
    const auto offset=restored.Stats().spillBytes;position=0;const uint64_t corruptionSize=9*1024*1024;
    restored.WriteStream(root+"/corrupt",corruptionSize,PatternHash(corruptionSize),reader);
    const int restoredSpill=lastSpillFd;assert(restoredSpill>=0);const char wrong=char(Pattern(17)^1);assert(pwrite(restoredSpill,&wrong,1,off_t(offset+17))==1);
    assert(Throws([&]{restored.ReadStream(root+"/corrupt",[](const char*,size_t){});}));
    struct stat saved{},after{};assert(stat(snapshot.c_str(),&saved)==0);assert(Throws([&]{restored.Snapshot(snapshot,"bad streamed spill");}));assert(stat(snapshot.c_str(),&after)==0&&saved.st_ino==after.st_ino&&saved.st_size==after.st_size);
    if(catalog){position=0;failCatalogEvent="before_put";failCatalogPath=root+"/sql_failure";
        assert(Throws([&]{streamed.WriteStream(root+"/sql_failure",100000,PatternHash(100000),reader);}));failCatalogEvent.clear();failCatalogPath.clear();
        assert(!streamed.Healthy()&&Throws([&]{streamed.Snapshot(snapshot,"partial publication");}));assert(stat(snapshot.c_str(),&after)==0&&saved.st_ino==after.st_ino&&saved.st_size==after.st_size);}
    std::cout<<"bounded streaming: catalog="<<catalog<<" bytes="<<size<<" budget="<<streamed.Limit()<<" max_chunk="<<largestChunk<<" checksum/failure/restore PASS\n";
}
static int CatalogScaleOnly(int records){
    assert(records>=100000&&records<=10000000);char temp[]="/tmp/mbs-memory-million-test-XXXXXX";assert(mkdtemp(temp));const std::string disk=temp,root=disk+"/logical";const size_t baseline=RssBytes();const auto start=std::chrono::steady_clock::now();
    {
        Store store(root,8*1024*1024,disk);auto path=[&](int i){return root+"/jobs/j"+std::to_string(i/16)+"/file"+std::to_string(i%16);};
        auto value=[](int i){return "matrix "+std::to_string(i)+" seed "+std::to_string(i%8)+" count 4096\n";};
        for(int i=0;i<records;++i)store.Write(path(i),value(i));
        for(int i=0;i<1000;++i){const int index=int((uint64_t(i)*7919)%uint64_t(records));assert(store.Read(path(index))==value(index));assert(store.Digest(path(index))==FullAutoMemoryStore::Hash(value(index)));}
        const auto stats=store.Stats();const auto rss=RssBytes();assert(stats.catalogActive&&stats.files==size_t(records)&&store.MemoryBytes()<=store.Limit()&&CountFiles(disk)==0);
#if !defined(__SANITIZE_ADDRESS__)
        assert(rss<=baseline+32*1024*1024);
#endif
        const auto elapsed=std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count();
        std::cout<<"{\"records\":"<<records<<",\"accounted_bytes\":"<<store.MemoryBytes()<<",\"limit_bytes\":"<<store.Limit()<<",\"rss_delta_bytes\":"<<int64_t(rss)-int64_t(baseline)<<",\"catalog_bytes\":"<<stats.catalogBytes<<",\"spill_bytes\":"<<stats.spillBytes<<",\"elapsed_seconds\":"<<elapsed<<"}\n";
    }
    assert(CountFiles(disk)==0);assert(rmdir(disk.c_str())==0);return 0;
}
int main(int argc,char** argv){
    if(argc==2)return CatalogScaleOnly(std::stoi(argv[1]));
    char temp[]="/tmp/mbs-memory-store-test-XXXXXX";assert(mkdtemp(temp));const std::string disk=temp,root=disk+"/logical",snapshot=disk+"/checkpoint.bin";
    assert(FullAutoMemoryStore::Hash("hello")==UINT64_C(0xa430d84680aabd0b));
    Store store(root,24576,disk);assert(store.Exists(root)&&store.IsDirectory(root));
    assert(store.Owns(root+"/a")&&!store.Owns(root+"x/a")&&!store.Owns(root+"/../a"));
    store.Directory(root+"/jobs/a/result");store.Write(root+"/jobs/a/result/data","hello");
    assert(store.Read(root+"/jobs/a/result/data")=="hello");assert(store.Digest(root+"/jobs/a/result/data")==FullAutoMemoryStore::Hash("hello"));
    auto generation=store.Generation(root+"/jobs/a/result/data"),epoch=store.Epoch();store.Write(root+"/jobs/a/result/data","world");
    assert(store.Generation(root+"/jobs/a/result/data")>generation&&store.Epoch()>epoch);
    assert(store.Read(root+"/jobs/a/result/data")=="world");
    store.Link(root+"/jobs/a/result/data",root+"/jobs/a/result/copy");store.Write(root+"/jobs/a/result/data","new");
    assert(store.Read(root+"/jobs/a/result/copy")=="world"); // links are immutable value aliases, never writable POSIX aliases
    store.Rename(root+"/jobs/a/result/copy",root+"/jobs/a/result/moved");assert(!store.Exists(root+"/jobs/a/result/copy"));
    store.Rename(root+"/jobs/a",root+"/jobs/longer-name");assert(!store.Exists(root+"/jobs/a"));
    assert(store.Read(root+"/jobs/longer-name/result/moved")=="world");
    assert(Throws([&]{store.Remove(root+"/jobs");}));
    assert(Throws([&]{store.Directory(root+"/jobs/longer-name/result/moved/child");}));
    const auto shallow=store.List(root+"/jobs",false);assert(shallow.size()==1&&shallow[0]==root+"/jobs/longer-name");
    assert(CountFiles(disk)==0); // virtual directories/files never touch disk
    std::string binary(8192,'q');binary[20]='\0';store.Write(root+"/large",binary);
    for(int i=0;i<5;++i)store.Write(root+"/blob"+std::to_string(i),std::string(7000,char('a'+i)));
    assert(store.Stats().spillWriteBytes>0&&store.Stats().evictions>0);assert(store.MemoryBytes()<=store.Limit());
    assert(store.Read(root+"/large")==binary);assert(store.MemoryBytes()<=store.Limit());assert(CountFiles(disk)==0); // spill is anonymous
    const auto names=store.List(root);const auto savedEpoch=store.Epoch();std::map<std::string,std::string> expected;
    for(const auto& name:names)if(!store.IsDirectory(name))expected[name]=store.Read(name);
    store.Snapshot(snapshot,"configuration and adaptive-state metadata");const auto good=Bytes(snapshot);assert(!good.empty()&&CountFiles(disk)==1);
    Store restored(root,24576,disk);assert(restored.Restore(snapshot)=="configuration and adaptive-state metadata");assert(restored.Epoch()==savedEpoch);
    assert(restored.List(root)==names);for(const auto& item:expected){assert(restored.Read(item.first)==item.second);assert(restored.Digest(item.first)==FullAutoMemoryStore::Hash(item.second));assert(restored.Generation(item.first)==store.Generation(item.first));}
    assert(restored.MemoryBytes()<=restored.Limit());
    // Malformed, truncated and trailing-byte snapshots cannot replace any part
    // of a running state. The old valid state survives each rejection.
    restored.Write(root+"/sentinel","keep me");
    for(const auto& bad:std::vector<std::string>{good.substr(0,good.size()-1),good+"extra",std::string(good.size(),'x')}){
        Bytes(snapshot,bad);assert(Throws([&]{restored.Restore(snapshot);}));assert(restored.Read(root+"/sentinel")=="keep me");
    }
    auto corrupted=good;corrupted[corrupted.size()-25]^=1;Bytes(snapshot,corrupted);assert(Throws([&]{restored.Restore(snapshot);}));assert(restored.Read(root+"/sentinel")=="keep me");Bytes(snapshot,good);
    // A member checksum is checked independently of the outer generation hash.
    corrupted=good;const auto member=corrupted.find("world");assert(member!=std::string::npos);corrupted[member]^=1;
    const auto forgedOuter=FullAutoMemoryStore::HashBytes(corrupted.data(),corrupted.size()-8);for(unsigned i=0;i<8;++i)corrupted[corrupted.size()-8+i]=char((forgedOuter>>(8*i))&255);
    Bytes(snapshot,corrupted);assert(Throws([&]{restored.Restore(snapshot);}));assert(restored.Read(root+"/sentinel")=="keep me");Bytes(snapshot,good);
    Store wrongRoot(root+"-other",24576,disk);assert(Throws([&]{wrongRoot.Restore(snapshot);}));
    Store small(root,2048,disk);assert(Throws([&]{small.Restore(snapshot);}));assert(small.List(root).empty());
    // Failed checkpoint publication leaves the previous complete generation.
    for(const auto& event:{"before_snapshot_fsync","before_snapshot_publish"}){
        failEvent=event;assert(Throws([&]{restored.Snapshot(snapshot,"replacement");}));failEvent.clear();assert(Bytes(snapshot)==good&&CountFiles(disk)==1);
    }
    crashEvent="before_snapshot_publish";const pid_t child=fork();assert(child>=0);if(child==0){restored.Snapshot(snapshot,"crashed generation");_exit(92);}
    int status=0;assert(waitpid(child,&status,0)==child&&WIFEXITED(status)&&WEXITSTATUS(status)==91);crashEvent.clear();assert(Bytes(snapshot)==good);
    {DIR* dir=opendir(disk.c_str());assert(dir);while(auto* entry=readdir(dir)){const std::string name=entry->d_name;if(name.find("checkpoint.bin.tmp.")==0)assert(unlink((disk+"/"+name).c_str())==0);}closedir(dir);}
    // Existing stale temporary file, such as after SIGKILL, cannot be mistaken
    // for a committed generation or be overwritten by the next save.
    Bytes(snapshot+".tmp.interrupted","incomplete");restored.Snapshot(snapshot,"new generation");assert(Bytes(snapshot+".tmp.interrupted")=="incomplete");
    Store newest(root,24576,disk);assert(newest.Restore(snapshot)=="new generation"&&newest.Read(root+"/sentinel")=="keep me");unlink((snapshot+".tmp.interrupted").c_str());
    // Metadata remains bounded even when every payload is empty.
    Store metadata(root,4096,disk);int added=0;while(!Throws([&]{metadata.Write(root+"/m"+std::to_string(added),"");}))++added;
    assert(added>0&&added<20&&metadata.MemoryBytes()<=metadata.Limit());for(int i=0;i<added;++i)assert(metadata.Exists(root+"/m"+std::to_string(i)));
    assert(!metadata.Exists(root+"/m"+std::to_string(added)));
    // A failed oversized directory tree is not partially installed.
    const auto before=metadata.List(root);assert(Throws([&]{metadata.Write(root+"/a/b/c/d/e/f/g/value","x");}));assert(metadata.List(root)==before);
    // Concurrent read-only validation safely shares LRU/spill state.
    std::thread readerA([&]{for(int i=0;i<50;++i)assert(restored.Read(root+"/blob0")==std::string(7000,'a'));});
    std::thread readerB([&]{for(int i=0;i<50;++i)assert(restored.Read(root+"/blob4")==std::string(7000,'e'));});readerA.join();readerB.join();assert(restored.MemoryBytes()<=restored.Limit());
    Store slash("/",4096,disk);slash.Write("/a","root test");slash.Snapshot(snapshot,"slash root");Store slashRestored("/",4096,disk);assert(slashRestored.Restore(snapshot)=="slash root"&&slashRestored.Read("/a")=="root test");
    // Disk corruption is detected when evicted bytes are materialized or
    // checkpointed; a cached digest itself does not silently reread the disk.
    {Store corruptSpill(root,2048,disk);corruptSpill.Write(root+"/huge",std::string(4000,'z'));const auto recorded=corruptSpill.Digest(root+"/huge");assert(lastSpillFd>=0);
        assert(pwrite(lastSpillFd,"x",1,0)==1);assert(corruptSpill.Digest(root+"/huge")==recorded);
        assert(Throws([&]{corruptSpill.Read(root+"/huge");}));const auto previous=Bytes(snapshot);
        assert(Throws([&]{corruptSpill.Snapshot(snapshot,"corrupt spill");}));assert(Bytes(snapshot)==previous);}
    // Repeated mutable state updates cannot leave an ever-growing anonymous
    // log of dead values. Compaction preserves other live entries and hashes.
    {Store compact(root,4096,disk);compact.Write(root+"/stable",std::string(5000,'s'));
        for(int i=0;i<68;++i)compact.Write(root+"/changing",std::string(1024*1024,char('a'+i%20)));
        assert(compact.Stats().spillBytes<64*1024*1024&&compact.Stats().liveSpillBytes==1024*1024+5000);
        assert(compact.Read(root+"/changing")==std::string(1024*1024,char('a'+67%20)));
        assert(compact.Read(root+"/stable")==std::string(5000,'s')&&compact.MemoryBytes()<=compact.Limit());}
    // Catalog spill removes the previous O(number-of-files) RAM floor. More
    // than one hundred thousand names stay within an 8 MiB working budget.
    assert(Throws([&]{FullAutoMemoryCatalog::Api absent("/definitely-missing-mbs-sqlite-library.so");}));
    {
        const size_t rssBefore=RssBytes(),limit=8*1024*1024;Store scale(root,limit,disk);const int records=120000;
        auto path=[&](int i){return root+"/many/j"+std::to_string(i/12)+"/file"+std::to_string(i%12);};
        auto value=[](int i){return "mueller="+std::to_string(i)+";seed="+std::to_string(i%8)+";count=4096";};
        for(int i=0;i<records;++i)scale.Write(path(i),value(i));
        const auto stats=scale.Stats();const size_t rssAfter=RssBytes();assert(stats.catalogActive&&stats.catalogBytes>0&&stats.catalogCacheBytes<=1024*1024);
        assert(stats.files==size_t(records)&&stats.directories==10002&&scale.MemoryBytes()<=limit);
        assert(CountFiles(disk)==1); // existing checkpoint only; both working spills are anonymous
        assert(scale.List(root+"/many",false).size()==10000&&scale.List(root+"/many/j13",false).size()==12);
        for(int i=0;i<1000;++i){const int index=(i*7919)%records;assert(scale.Read(path(index))==value(index));assert(scale.Digest(path(index))==FullAutoMemoryStore::Hash(value(index)));}
#if !defined(__SANITIZE_ADDRESS__)
        assert(rssAfter<=rssBefore+32*1024*1024); // includes loaded SQLite code and allocator retention
#endif
        scale.Snapshot(snapshot,"large bounded catalog");Store scaleRestored(root,limit,disk);assert(scaleRestored.Restore(snapshot)=="large bounded catalog");assert(scaleRestored.Stats().catalogActive&&scaleRestored.Stats().files==size_t(records));
        for(int i=0;i<1000;++i){const int index=(i*7907)%records;assert(scaleRestored.Read(path(index))==value(index));assert(scaleRestored.Generation(path(index))==scale.Generation(path(index)));}
        const auto generationBefore=scale.Generation(path(9));scale.Write(path(9),"new generation");assert(scale.Read(path(9))=="new generation"&&scale.Generation(path(9))>generationBefore);
        scale.Rename(root+"/many/j0",root+"/many/moved");assert(!scale.Exists(root+"/many/j0")&&scale.Read(root+"/many/moved/file9")=="new generation");
        scale.Remove(root+"/many/moved",true);assert(scale.Stats().files==size_t(records-12)&&scale.Stats().directories==10001);
        for(int i=0;i<68;++i)scale.Write(root+"/changing",std::string(1024*1024,char('a'+i%20)));
        assert(scale.Stats().spillBytes<64*1024*1024&&scale.Read(root+"/changing")==std::string(1024*1024,char('a'+67%20)));
        assert(scale.Read(path(111111))==value(111111)&&scale.MemoryBytes()<=limit&&scale.Stats().files==size_t(records-11));
        scale.Snapshot(snapshot,"catalog after compaction");assert(scaleRestored.Restore(snapshot)=="catalog after compaction");
        assert(scaleRestored.Read(path(111111))==value(111111)&&!scaleRestored.Exists(root+"/many/j0"));
        auto damagedLarge=Bytes(snapshot);damagedLarge.back()^=1;Bytes(snapshot,damagedLarge);
        assert(Throws([&]{scaleRestored.Restore(snapshot);}));assert(scaleRestored.Read(path(111111))==value(111111));
        std::cout<<"catalog scale: records="<<records<<" budget="<<limit<<" accounted="<<stats.residentBytes+stats.metadataBytes<<" rss_delta="<<int64_t(rssAfter)-int64_t(rssBefore)<<" catalog_bytes="<<stats.catalogBytes<<"\n";
    }
    // Both temporary files are anonymous even if the controller is killed
    // before its destructor can close the SQLite connection.
    {int ready[2];assert(pipe(ready)==0);const pid_t killed=fork();assert(killed>=0);
        if(killed==0){close(ready[0]);Store disposable(root+"/killed",1024*1024,disk);for(int i=0;i<5000;++i)disposable.Write(root+"/killed/f"+std::to_string(i),"working value");
            assert(disposable.Stats().catalogActive);assert(write(ready[1],"r",1)==1);for(;;)pause();}
        close(ready[1]);char signal=0;assert(read(ready[0],&signal,1)==1&&signal=='r');close(ready[0]);assert(kill(killed,SIGKILL)==0);int status=0;assert(waitpid(killed,&status,0)==killed&&WIFSIGNALED(status));assert(CountFiles(disk)==1);}
    // Journal-OFF working mutations are fail-closed: an interrupted multirow
    // rename may not publish a recovery snapshot of a partial operation.
    {Store failed(root,1024*1024,disk);for(int i=0;i<5000;++i)failed.Write(root+"/f"+std::to_string(i),"value "+std::to_string(i));assert(failed.Stats().catalogActive);
        failed.Write(root+"/literal_%/quote'/file","literal");assert(failed.List(root+"/literal_%",true).size()==2);
        failed.Snapshot(snapshot,"before catalog failure");const auto previous=Bytes(snapshot);
        failCatalogEvent="before_erase";failCatalogPath=root+"/f1";assert(Throws([&]{failed.Rename(root+"/f1",root+"/renamed");}));failCatalogEvent.clear();failCatalogPath.clear();
        assert(!failed.Healthy()&&!failed.Stats().healthy);assert(Throws([&]{failed.Read(root+"/f2");}));assert(Throws([&]{failed.Snapshot(snapshot,"must not publish");}));assert(Bytes(snapshot)==previous);
        assert(failed.Restore(snapshot)=="before catalog failure"&&failed.Healthy());assert(failed.Read(root+"/f1")=="value 1"&&!failed.Exists(root+"/renamed"));}
    {Store slashCatalog("/",1024*1024,disk);for(int i=0;i<2500;++i)slashCatalog.Write("/f"+std::to_string(i),"x");assert(slashCatalog.Stats().catalogActive);
        assert(slashCatalog.List("/",false).size()==2500&&slashCatalog.List("/",true).size()==2500);}
    StreamChecks(root+"/stream_map",disk,snapshot,false);
    StreamChecks(root+"/stream_catalog",disk,snapshot,true);
    assert(FullAutoMemoryStore::Active()==nullptr);{FullAutoMemoryStore::Scope first(store);assert(FullAutoMemoryStore::Active()==&store);{FullAutoMemoryStore::Scope second(restored);assert(FullAutoMemoryStore::Active()==&restored);}assert(FullAutoMemoryStore::Active()==&store);}assert(FullAutoMemoryStore::Active()==nullptr);
    store.Remove(root+"/jobs",true);assert(!store.Exists(root+"/jobs"));assert(Throws([&]{store.Remove(root,true);}));
    unlink(snapshot.c_str());assert(CountFiles(disk)==0);rmdir(disk.c_str());
    std::cout<<"fullauto memory store: RAM bounds, anonymous spill, generations, concurrent reads and atomic checked snapshots passed\n";
}
