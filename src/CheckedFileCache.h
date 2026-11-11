#pragma once
#include "BoundedCache.h"
#include "PackedJobStore.h"
#include "ReadOnlyFile.h"
#include <cstdint>
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <sys/stat.h>
#include <time.h>

// Process-local, bounded cache. Every hit is guarded by device/inode, size,
// nanosecond mtime AND ctime; a restart never trusts old cached digests.
// Small validation manifests share their text and digest. Large matrices keep
// only a digest. Loose keys also bind the logical path: packing can delete a
// file and the filesystem can immediately reuse its inode for another job.
class CheckedFileCache {
    struct Entry {struct stat stamp;std::string hash,text;bool hasText=false;};
    struct LargeEntry {std::shared_ptr<const Entry> value;std::string prefix;size_t decodedSize=0;};
    BoundedCache<Entry> entries;
    // Large source manifests are revisited when an already verified graph is
    // reused. Keep them in a separate, bounded part of the SAME budget: they
    // must neither be reread every time nor evict all small leaf digests.
    BoundedCache<LargeEntry> largeTexts;
    size_t maximum;
    static size_t LargeBudget(size_t limit){return limit>=1024*1024?limit/4:0;}
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
    // Some ext4/WSL kernels assign several writes the same coarse timestamp.
    // A same-path rewrite can then retain all of Same()'s fields. Never use a
    // loose-file hit inside that current timestamp interval. One second also
    // covers filesystems exposing whole-second timestamps; future stamps and
    // a failed clock read are conservative misses. Packed immutable members
    // keep their separate verified generation/offset identity below.
    static bool Mature(const struct stat& stamp){
        struct timespec now{};
        if(clock_gettime(CLOCK_REALTIME,&now)!=0||now.tv_sec<=0
                ||now.tv_nsec<0||now.tv_nsec>=1000000000)return false;
        auto oldEnough=[&](const struct timespec& value){
            return value.tv_nsec>=0&&value.tv_nsec<1000000000
                &&(value.tv_sec<now.tv_sec-1
                    ||(value.tv_sec==now.tv_sec-1&&value.tv_nsec<=now.tv_nsec));
        };
        return oldEnough(stamp.st_mtim)&&oldEnough(stamp.st_ctim);
    }
    static std::string Hex(uint64_t hash){std::ostringstream out;out<<std::hex<<hash;return out.str();}
    std::shared_ptr<const Entry> Find(const std::string& key,bool text){
        const auto found=entries.Get(key);
        if(found&&(!text||found->hasText))return found;
        const auto large=largeTexts.Get(key);if(!large)return found;
        if(!text||large->prefix.empty())return large->value;
        auto decoded=std::make_shared<Entry>();decoded->stamp=large->value->stamp;decoded->hash=large->value->hash;decoded->hasText=true;
        decoded->text.reserve(large->decodedSize);const auto& body=large->value->text;
        // A final line may consist solely of the factored prefix, with no LF:
        // then its encoded suffix is empty but the original length retains it.
        for(size_t at=0;at<body.size()||decoded->text.size()<large->decodedSize;){const auto newline=body.find('\n',at);const size_t end=newline==std::string::npos?body.size():newline+1;
            decoded->text+=large->prefix;decoded->text.append(body,at,end-at);at=end;}
        Check(decoded->text.size()==large->decodedSize,"invalid resident manifest length");return decoded;
    }
    static std::shared_ptr<LargeEntry> Large(const std::shared_ptr<Entry>& value){
        auto result=std::make_shared<LargeEntry>();result->value=value;result->decodedSize=value->text.size();
        // Source lists repeat the same long absolute directory on every line.
        // Factor that prefix in RAM only. The remainder includes its original
        // newline bytes, so CRLF, NUL and an unterminated final line survive.
        const auto& text=value->text;const auto first=text.find('\n');if(first==std::string::npos)return result;
        size_t common=first,lines=1;
        for(size_t at=first+1;at<text.size();){++lines;const auto newline=text.find('\n',at);const size_t end=newline==std::string::npos?text.size():newline;
            common=std::min(common,end-at);size_t k=0;while(k<common&&text[k]==text[at+k])++k;common=k;
            if(common<16)return result;at=end+(newline!=std::string::npos);}
        if(lines<2||common<16||common*(lines-1)<text.size()/4)return result;
        auto encoded=std::make_shared<Entry>();encoded->stamp=value->stamp;encoded->hash=value->hash;
        encoded->text.reserve(text.size()-common*lines);result->prefix=text.substr(0,common);
        for(size_t at=0;at<text.size();){const auto newline=text.find('\n',at);const size_t end=newline==std::string::npos?text.size():newline+1;
            encoded->text.append(text,at+common,end-at-common);at=end;}
        result->value=encoded;return result;
    }
    void Retain(const std::string& key,const std::shared_ptr<Entry>& value){
        largeTexts.Erase(key);entries.SetLimit(maximum-largeTexts.Bytes());std::shared_ptr<Entry> retained=value;
        if(value->text.size()>65536){
            retained=std::make_shared<Entry>();retained->stamp=value->stamp;retained->hash=value->hash;
            // Bound both one decoded manifest and its resident representation.
            if(value->text.size()<=largeTexts.Limit()/2){const auto large=Large(value);
                const size_t bytes=sizeof(LargeEntry)+large->prefix.capacity()+sizeof(Entry)
                    +large->value->hash.capacity()+large->value->text.capacity();
                const size_t charge=BoundedCache<LargeEntry>::RetainedBytes(key,bytes);
                if(bytes<=largeTexts.Limit()/2&&charge<=largeTexts.Limit()){
                    // Reserve before retaining the new entry; then return any
                    // slack after large-LRU eviction to the small-file region.
                    // A small-only workload keeps the original entire budget.
                    entries.SetLimit(maximum-std::min(largeTexts.Limit(),largeTexts.Bytes()+charge));
                    largeTexts.Put(key,large,bytes);entries.SetLimit(maximum-largeTexts.Bytes());}}
        }
        entries.Put(key,retained,sizeof(Entry)+retained->hash.capacity()+retained->text.capacity());
    }
    std::shared_ptr<const Entry> Get(const std::string& path,bool text){
        struct stat before,after;
        if(stat(path.c_str(),&before)!=0){
            const int error=errno;
            Check(error==ENOENT||error==ENOTDIR,"cannot stat "+path);
            PackedJobStore::Member member;
            Check(PackedJobStore::Locate(path,member)&&!member.directory,"cannot stat "+path);
            // Our own verified appends preserve the member generation. Any
            // unexplained archive stamp change invalidates the whole generation.
            const auto key=Key(member.stamp)+":"+std::to_string(member.generation)+":"+std::to_string(member.offset);
            const auto found=Find(key,text);
            if(found&&(!text||found->hasText)){if(text)++textHits;else ++hashHits;return found;}
            const auto bytes=PackedJobStore::Read(member);++opens;bytesRead+=bytes.size();
            auto value=std::make_shared<Entry>();value->stamp=member.stamp;value->hash=member.hash;
            const bool marker=path.size()>=13&&path.compare(path.size()-13,13,"/complete.txt")==0;
            value->hasText=text||marker;if(value->hasText)value->text=bytes;
            Retain(key,value);return value;
        }
        // The full path is charged by BoundedCache with the key's capacity.
        // Keep raw and packed key namespaces distinct after the inode bytes.
        const auto key=Key(before)+std::string(1,'\0')+path;
        const bool mature=Mature(before);
        const auto found=mature?Find(key,text):std::shared_ptr<const Entry>();
        if(!mature){
            // A young read is not reusable later: the file can change again
            // in the same timestamp tick before becoming mature. Discard any
            // earlier value too, including a retained large manifest.
            entries.Erase(key);largeTexts.Erase(key);
            entries.SetLimit(maximum-largeTexts.Bytes());
        }
        if(found&&Same(before,found->stamp)&&(!text||found->hasText)){
            if(text)++textHits;else ++hashHits;return found;
        }
        CheckpointRead::File input(path.c_str());Check(input.fd>=0,"cannot read "+path);++opens;
        // A source marker is normally hashed immediately before recursively
        // validating that source. Retain its bytes during the hash read too,
        // so the following Text() needs no second open/read of the same file.
        const std::string marker="/complete.txt";
        const bool markerText=before.st_size<=65536&&path.size()>=marker.size()&&path.compare(path.size()-marker.size(),marker.size(),marker)==0;
        const bool keepText=text||markerText;
        auto value=std::make_shared<Entry>();value->stamp=before;value->hasText=keepText;
        uint64_t hash=14695981039346656037ULL;char buffer[65536];
        auto consume=[&](size_t size){bytesRead+=size;
            for(size_t i=0;i<size;++i){hash^=static_cast<unsigned char>(buffer[i]);hash*=1099511628211ULL;}
            if(keepText)value->text.append(buffer,size);
        };
        // Checkpoint files have a known length. Avoid a second empty read for
        // every small file; the unchanged post-read stamp still rejects size,
        // inode, mtime or ctime changes. Zero-length virtual files (e.g. procfs)
        // and non-regular streams retain the read-to-EOF path.
        if(S_ISREG(before.st_mode)&&before.st_size>0){
            uint64_t remaining=uint64_t(before.st_size);
            while(remaining){const auto wanted=size_t(std::min<uint64_t>(remaining,sizeof(buffer)));
                const auto size=CheckpointRead::Read(input.fd,buffer,wanted);
                Check(size>=0,"cannot read "+path);consume(size_t(size));
                // Some virtual regular files advertise a larger extent than
                // their contents. A real truncation still changes the stamp.
                if(size==0)break;
                remaining-=uint64_t(size);
            }
        }else{
            for(;;){const auto size=CheckpointRead::Read(input.fd,buffer,sizeof(buffer));
                Check(size>=0,"cannot read "+path);if(size==0)break;consume(size_t(size));}
        }
        Check(stat(path.c_str(),&after)==0&&Same(before,after),"file changed while reading "+path);
        value->hash=Hex(hash);
        if(mature)Retain(key,value);
        return value;
    }
public:
    size_t opens=0,hashHits=0,textHits=0;uint64_t bytesRead=0;
    explicit CheckedFileCache(size_t limit=32*1024*1024):entries(limit),largeTexts(LargeBudget(limit)),maximum(limit){}
    size_t Bytes()const{return entries.Bytes()+largeTexts.Bytes();}
    size_t Limit()const{return maximum;}
    std::string Hash(const std::string& path){return Get(path,false)->hash;}
    std::string Text(const std::string& path){return Get(path,true)->text;}
    // Verify and return the same bytes with one freshness check. A cached
    // digest alone must never authenticate a separately reopened manifest.
    std::string VerifiedText(const std::string& path,const std::string& expected){
        const auto value=Get(path,true);
        Check(value->hash==expected,"cached manifest changed "+path);
        return value->text;
    }
    static std::string Read(const std::string& path){
        CheckpointRead::File input(path.c_str());
        if(input.fd<0){const int error=errno;
            // A pack may replace an absent loose member, not hide an I/O or
            // permission failure while opening an existing source.
            Check(error==ENOENT||error==ENOTDIR,"cannot read "+path);
            PackedJobStore::Member member;Check(PackedJobStore::Locate(path,member)&&!member.directory,"cannot read "+path);return PackedJobStore::Read(member);}
        std::string result;char buffer[65536];
        for(;;){const auto size=CheckpointRead::Read(input.fd,buffer,sizeof(buffer));
            Check(size>=0,"cannot read "+path);if(size==0)break;result.append(buffer,size_t(size));}
        return result;
    }
};
