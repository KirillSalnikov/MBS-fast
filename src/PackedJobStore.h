#pragma once
#include "BoundedCache.h"
#include <algorithm>
#include <cerrno>
#include <cstdint>
#include <cstring>
#include <dirent.h>
#include <fcntl.h>
#include <iomanip>
#include <mutex>
#include <sstream>
#include <stdexcept>
#include <sys/stat.h>
#include <unistd.h>
#include <vector>

// Lossless backing store for immutable fullauto job trees. Logical paths and
// their original hashes do not change. A task uses at most 256 pack files;
// completed jobs no longer consume a directory and many inodes apiece.
// The controller lock is the writer lock. Readers never repair storage.
namespace PackedJobStore {
inline void Check(bool good,const std::string& message){if(!good)throw std::runtime_error("fullauto pack: "+message);}
inline uint64_t Digest(const std::string& bytes){uint64_t h=14695981039346656037ULL;for(unsigned char c:bytes){h^=c;h*=1099511628211ULL;}return h;}
inline std::string Hex(uint64_t value){std::ostringstream s;s<<std::hex<<value;return s.str();}
inline uint64_t HexNumber(const std::string& text){Check(!text.empty()&&text.size()<=16,"invalid stored member hash");uint64_t value=0;for(char c:text){Check((c>='0'&&c<='9')||(c>='a'&&c<='f'),"invalid stored member hash");value=16*value+uint64_t(c<='9'?c-'0':c-'a'+10);}return value;}
inline bool Same(const struct stat& a,const struct stat& b){return a.st_dev==b.st_dev&&a.st_ino==b.st_ino&&a.st_size==b.st_size&&a.st_mtim.tv_sec==b.st_mtim.tv_sec&&a.st_mtim.tv_nsec==b.st_mtim.tv_nsec&&a.st_ctim.tv_sec==b.st_ctim.tv_sec&&a.st_ctim.tv_nsec==b.st_ctim.tv_nsec;}
struct Descriptor {int fd;explicit Descriptor(int value):fd(value){}~Descriptor(){if(fd>=0)close(fd);}Descriptor(const Descriptor&)=delete;};
inline std::string Bytes(int fd,uint64_t offset,size_t size){
    std::string result(size,'\0');size_t at=0;
    while(at<size){const auto n=pread(fd,&result[at],size-at,off_t(offset+at));if(n<0&&errno==EINTR)continue;Check(n>0,"short archive read");at+=size_t(n);}return result;
}
inline void Write(int fd,const std::string& bytes){size_t at=0;while(at<bytes.size()){const auto n=write(fd,bytes.data()+at,bytes.size()-at);if(n<0&&errno==EINTR)continue;Check(n>0,"archive write failed");at+=size_t(n);}}
inline void Flush(int fd){int result;do{result=fsync(fd);}while(result<0&&errno==EINTR);Check(result==0,"archive fsync failed");}
inline void FlushDirectory(const std::string& path){Descriptor fd(open(path.c_str(),O_RDONLY|O_DIRECTORY));Check(fd.fd>=0,"cannot open directory "+path);Flush(fd.fd);}
inline bool JobName(const std::string& name){if(name.size()<2||name.size()>17||name[0]<'a'||name[0]>'z'||(name.size()>2&&name[1]=='0'))return false;for(size_t k=1;k<name.size();++k)if(!((name[k]>='0'&&name[k]<='9')||(name[k]>='a'&&name[k]<='f')))return false;return true;}
inline bool Relative(const std::string& path){if(path.empty()||path[0]=='/'||path.back()=='/')return false;size_t at=0;while(at<path.size()){const auto end=path.find('/',at);const auto part=path.substr(at,end==std::string::npos?end:end-at);if(part.empty()||part=="."||part==".."||part.find('\0')!=std::string::npos)return false;if(end==std::string::npos)break;at=end+1;}return true;}
inline bool Split(const std::string& path,std::string& root,std::string& name,std::string& member){
    const auto begin=path.rfind("/jobs/");if(begin==std::string::npos)return false;
    const auto end=path.find('/',begin+6);root=path.substr(0,begin);name=path.substr(begin+6,end==std::string::npos?end:end-begin-6);member=end==std::string::npos?"":path.substr(end+1);
    return !root.empty()&&JobName(name)&&(member.empty()||Relative(member));
}
inline std::string CanonicalRoot(const std::string& root){
    // A frozen parent may still be referenced through its old symlink after a
    // disk migration. Both paths must share one append generation and index.
    static std::mutex mutex;static std::map<std::string,std::string> roots;std::lock_guard<std::mutex> lock(mutex);
    const auto found=roots.find(root);if(found!=roots.end())return found->second;
    char* resolved=realpath(root.c_str(),nullptr);if(!resolved)return root;const std::string actual=resolved;free(resolved);
    if(roots.size()>=64)roots.clear();roots[root]=actual;return actual;
}
inline std::string Shard(const std::string& root,const std::string& name){std::ostringstream out;out<<CanonicalRoot(root)<<"/job_packs/"<<std::hex<<std::setw(2)<<std::setfill('0')<<(Digest(name)&255)<<".pack";return out.str();}
inline std::string Encode(const std::string& value){const char* digits="0123456789abcdef";std::string out;out.reserve(2*value.size());for(unsigned char c:value){out+=digits[c>>4];out+=digits[c&15];}return out;}
inline std::string Decode(const std::string& value){Check(value.size()%2==0,"bad member name encoding");std::string out;for(size_t k=0;k<value.size();k+=2){unsigned v=0;for(size_t j=k;j<k+2;++j){const auto c=value[j];Check((c>='0'&&c<='9')||(c>='a'&&c<='f'),"bad member name encoding");v=16*v+(c<='9'?c-'0':c-'a'+10);}out+=char(v);}Check(Relative(out),"unsafe member path");return out;}
const size_t HeaderBytes=256,TrailerBytes=64;
const size_t FrontierSlotBytes=4096,DataBegin=2*FrontierSlotBytes;
const uint64_t MaximumRecordBytes=64ULL*1024*1024;
struct Oversize : std::runtime_error {explicit Oversize(const std::string& message):std::runtime_error("fullauto pack: "+message){}};
inline std::string Fixed(const std::string& line,size_t width){Check(line.size()<width,"record header overflow");return line+std::string(width-line.size()-1,' ')+'\n';}
struct Frontier {uint64_t sequence=0,end=0;};
inline std::string FrontierText(uint64_t sequence,uint64_t end){const auto text=std::to_string(sequence)+" "+std::to_string(end);return Fixed("MBS_JOB_FRONTIER_V1 "+text+" "+Hex(Digest(text)),FrontierSlotBytes);}
inline Frontier ReadFrontier(int fd,uint64_t size){
    Check(size>=DataBegin,"truncated archive frontier");Frontier best;
    for(unsigned slot=0;slot<2;++slot){const auto text=Bytes(fd,slot*FrontierSlotBytes,FrontierSlotBytes);std::istringstream row(text);std::string magic,hash,extra;Frontier candidate;
        const bool valid=bool(row>>magic>>candidate.sequence>>candidate.end>>hash)&&!(row>>extra)&&magic=="MBS_JOB_FRONTIER_V1"&&candidate.sequence>0&&candidate.end>=DataBegin&&hash==Hex(Digest(std::to_string(candidate.sequence)+" "+std::to_string(candidate.end)));
        if(valid){Check(candidate.end<=size,"committed archive extent was truncated");Check(candidate.sequence!=best.sequence||candidate.end==best.end,"conflicting archive frontiers");if(candidate.sequence>best.sequence)best=candidate;}
    }
    Check(best.sequence>0,"both archive frontier copies are corrupt");return best;
}
inline bool MirrorFrontier(int fd,uint64_t sequence,uint64_t end){
    const auto text=FrontierText(sequence,end);bool changed=false;
    for(unsigned slot=0;slot<2;++slot)if(Bytes(fd,slot*FrontierSlotBytes,FrontierSlotBytes)!=text){
        Check(lseek(fd,off_t(slot*FrontierSlotBytes),SEEK_SET)==off_t(slot*FrontierSlotBytes),"cannot seek mirrored frontier");Write(fd,text);Flush(fd);changed=true;
    }
    Check(Bytes(fd,0,FrontierSlotBytes)==text&&Bytes(fd,FrontierSlotBytes,FrontierSlotBytes)==text,"mirrored frontier readback differs");return changed;
}
struct Record {uint64_t offset=0,metadata=0,payload=0;std::string metadataHash,payloadHash;};
// Each generated name is one letter plus a canonical 64-bit hexadecimal hash.
// Store only its identity and record offset, not 2 million string/map nodes.
// At 2.3 million jobs the 256 shard indices occupy about 64 MiB in total.
class CompactIndex {
    struct Slot {uint64_t identity=0,location=0;};
    std::vector<Slot> slots;size_t used=0;
    static uint64_t Identity(const std::string& name){uint64_t value=0;for(size_t k=1;k<name.size();++k)value=16*value+uint64_t(name[k]<='9'?name[k]-'0':name[k]-'a'+10);return value;}
    static uint64_t Mix(uint64_t value,unsigned prefix){value^=uint64_t(prefix)*0x9e3779b97f4a7c15ULL;value^=value>>30;value*=0xbf58476d1ce4e5b9ULL;value^=value>>27;value*=0x94d049bb133111ebULL;return value^(value>>31);}
    void Add(uint64_t identity,uint64_t location){size_t at=size_t(Mix(identity,unsigned(location&255)))&(slots.size()-1);while(slots[at].location)at=(at+1)&(slots.size()-1);slots[at].identity=identity;slots[at].location=location;++used;}
public:
    size_t size()const{return used;}
    size_t Bytes()const{return sizeof(CompactIndex)+slots.capacity()*sizeof(Slot);}
    bool Find(const std::string& name,uint64_t& offset)const{if(slots.empty())return false;const auto id=Identity(name);size_t at=size_t(Mix(id,unsigned(name[0])))&(slots.size()-1);while(slots[at].location){const auto& slot=slots[at];if(slot.identity==id&&(slot.location&255)==unsigned(name[0])){offset=slot.location>>8;return true;}at=(at+1)&(slots.size()-1);}return false;}
    size_t count(const std::string& name)const{uint64_t offset;return Find(name,offset)?1:0;}
    bool Insert(const std::string& name,uint64_t offset){Check(JobName(name)&&offset<(1ULL<<56),"invalid compact index identity");if(count(name))return false;
        if(slots.empty()||(used+1)*10>slots.size()*7){std::vector<Slot> old;old.swap(slots);slots.resize(old.empty()?16:old.size()*2);used=0;for(const auto& slot:old)if(slot.location)Add(slot.identity,slot.location);}
        Add(Identity(name),(offset<<8)|unsigned(name[0]));return true;}
    template<class Consumer> void Each(Consumer consume)const{for(const auto& slot:slots)if(slot.location)consume(std::string(1,char(slot.location&255))+Hex(slot.identity),slot.location>>8);}
};
struct Index {struct stat stamp{};uint64_t validEnd=DataBegin,generation=0,frontierSequence=0;CompactIndex jobs;};
struct Member {std::string archive,hash;uint64_t offset=0,size=0,generation=0;bool directory=false;struct stat stamp{};};
// Cache the extent only. Repeating a long archive path and a struct stat in
// every file entry made a few thousand ordinary jobs evict each other.
// Archive identity/freshness belongs to the shard index, once per lookup.
struct MemberExtent {uint64_t offset=0,size=0,hash=0;bool directory=false;};
struct MemberIndex {std::map<std::string,MemberExtent> entries;};
struct Store {std::mutex mutex;uint64_t generation=0;BoundedCache<Index> indices{128*1024*1024};BoundedCache<MemberIndex> members{8*1024*1024};};
inline Store& Global(){static Store store;return store;}
inline Record ParseHeader(const std::string& header,std::string& name){
    std::istringstream in(header);std::string magic,extra;Record record;
    Check(bool(in>>magic>>name>>record.metadata>>record.payload>>record.metadataHash>>record.payloadHash)&&!(in>>extra)&&magic=="MBS_JOB_PACK_V1"&&JobName(name),"invalid record header");
    Check(record.metadata>0&&record.metadata<=16*1024*1024&&record.payload<=MaximumRecordBytes&&record.metadata+record.payload<=MaximumRecordBytes,"invalid record length");return record;
}
inline std::shared_ptr<const Index> LoadIndex(const std::string& path){
    // Caller holds Global().mutex. Only a few shard indices reside in RAM.
    struct stat stamp{};if(stat(path.c_str(),&stamp)!=0){Check(errno==ENOENT,"cannot stat archive "+path);return {};}
    auto found=Global().indices.Get(path);if(found&&Same(found->stamp,stamp))return found;
    Descriptor fd(open(path.c_str(),O_RDONLY));Check(fd.fd>=0,"cannot open archive "+path);
    auto index=std::make_shared<Index>();index->stamp=stamp;index->generation=++Global().generation;
    const auto frontier=ReadFrontier(fd.fd,uint64_t(stamp.st_size));index->frontierSequence=frontier.sequence;
    uint64_t at=DataBegin;const auto size=frontier.end;
    while(at<size){
        Check(size-at>=HeaderBytes,"truncated committed record header");
        const auto header=Bytes(fd.fd,at,HeaderBytes);std::string name;auto record=ParseHeader(header,name);record.offset=at;
        const uint64_t total=HeaderBytes+record.metadata+record.payload+TrailerBytes;
        Check(size-at>=total,"truncated committed record body");
        const auto trailer=Bytes(fd.fd,at+total-TrailerBytes,TrailerBytes);
        // The commit authenticates the header, including the metadata hash.
        // Metadata and payload remain lazy and are checked on first use.
        Check(trailer==Fixed("MBS_JOB_COMMIT_V1 "+Hex(Digest(header)),TrailerBytes),"committed record trailer changed");
        Check(index->jobs.Insert(name,at),"duplicate packed job");at+=total;
    }
    index->validEnd=at;struct stat after{};Check(fstat(fd.fd,&after)==0&&Same(stamp,after),"archive changed during index read");
    Global().indices.Put(path,index,sizeof(Index)+index->jobs.Bytes());return index;
}
inline Record ReadRecord(int fd,const std::string& name,uint64_t offset){std::string saved;auto record=ParseHeader(Bytes(fd,offset,HeaderBytes),saved);Check(name==saved,"packed index points to another job");record.offset=offset;return record;}
inline std::map<std::string,Member> Members(int fd,const std::string& archive,const Record& record,const struct stat& stamp){
    const auto metadata=Bytes(fd,record.offset+HeaderBytes,size_t(record.metadata));Check(Hex(Digest(metadata))==record.metadataHash,"member index changed");
    std::istringstream rows(metadata);std::string line;std::map<std::string,Member> members;uint64_t next=0;
    while(std::getline(rows,line)){std::istringstream row(line);std::string type,name,extra;Member value;value.archive=archive;value.stamp=stamp;
        Check(bool(row>>type>>name),"invalid member metadata");name=Decode(name);
        if(type=="D"){value.directory=true;Check(!(row>>extra),"invalid directory metadata");}
        else{uint64_t offset=0;Check(type=="F"&&bool(row>>offset>>value.size>>value.hash)&&!(row>>extra)&&offset==next&&value.size<=record.payload-next,"invalid member extent");value.offset=record.offset+HeaderBytes+record.metadata+offset;next+=value.size;}
        Check(members.emplace(name,value).second,"duplicate member");
    }
    Check(next==record.payload&&members.count("complete.txt")&&!members.at("complete.txt").directory,"packed job has no completion marker");return members;
}
inline bool Locate(const std::string& path,Member& member){
    std::string root,name,part;if(!Split(path,root,name,part))return false;
    const auto archive=Shard(root,name);std::lock_guard<std::mutex> lock(Global().mutex);auto index=LoadIndex(archive);if(!index)return false;
    uint64_t offset=0;if(!index->jobs.Find(name,offset))return false;
    if(part.empty()){member.archive=archive;member.directory=true;member.stamp=index->stamp;member.generation=index->generation;return true;}
    std::string key(reinterpret_cast<const char*>(&index->stamp.st_dev),sizeof(index->stamp.st_dev));
    key.append(reinterpret_cast<const char*>(&index->stamp.st_ino),sizeof(index->stamp.st_ino));
    key.append(reinterpret_cast<const char*>(&index->generation),sizeof(index->generation));key.append(reinterpret_cast<const char*>(&offset),sizeof(offset));
    auto members=Global().members.Get(key);
    if(!members){Descriptor fd(open(archive.c_str(),O_RDONLY));Check(fd.fd>=0,"cannot open archive");auto loaded=std::make_shared<MemberIndex>();const auto decoded=Members(fd.fd,archive,ReadRecord(fd.fd,name,offset),index->stamp);
        size_t size=sizeof(MemberIndex);for(const auto& item:decoded){MemberExtent extent;extent.offset=item.second.offset;extent.size=item.second.size;extent.directory=item.second.directory;if(!extent.directory)extent.hash=HexNumber(item.second.hash);loaded->entries.emplace(item.first,extent);size+=sizeof(MemberExtent)+128+item.first.size();}Global().members.Put(key,loaded,size);members=loaded;}
    const auto item=members->entries.find(part);if(item==members->entries.end())return false;
    member.archive=archive;member.offset=item->second.offset;member.size=item->second.size;member.directory=item->second.directory;member.hash=member.directory?"":Hex(item->second.hash);member.stamp=index->stamp;member.generation=index->generation;return true;
}
inline bool Exists(const std::string& path){Member member;return Locate(path,member);}
inline std::string Read(const Member& member){
    Check(!member.directory,"cannot read packed directory");Descriptor fd(open(member.archive.c_str(),O_RDONLY));Check(fd.fd>=0,"cannot read packed member");
    struct stat before{},after{};Check(fstat(fd.fd,&before)==0&&Same(before,member.stamp),"archive changed before member read");
    const auto bytes=Bytes(fd.fd,member.offset,size_t(member.size));Check(Hex(Digest(bytes))==member.hash,"packed member checksum mismatch");
    Check(fstat(fd.fd,&after)==0&&Same(before,after),"archive changed during member read");return bytes;
}
inline std::vector<std::string> Jobs(const std::string& root){
    std::vector<std::string> result;std::lock_guard<std::mutex> lock(Global().mutex);
    for(unsigned k=0;k<256;++k){std::ostringstream path;path<<CanonicalRoot(root)<<"/job_packs/"<<std::hex<<std::setw(2)<<std::setfill('0')<<k<<".pack";
        auto index=LoadIndex(path.str());if(index)index->jobs.Each([&](const std::string& name,uint64_t){result.push_back(name);});}
    return result;
}
struct Source {std::string name,bytes;bool directory=false;struct stat stamp{};};
inline void Collect(const std::string& folder,const std::string& prefix,std::vector<Source>& sources,uint64_t& bytes){
    DIR* directory=opendir((folder+(prefix.empty()?"":"/"+prefix)).c_str());Check(directory,"cannot list job tree");std::vector<std::string> names;
    while(auto* entry=readdir(directory))if(std::string(entry->d_name)!="."&&std::string(entry->d_name)!="..")names.push_back(entry->d_name);closedir(directory);std::sort(names.begin(),names.end());
    for(const auto& name:names){Source source;source.name=prefix.empty()?name:prefix+"/"+name;Check(Relative(source.name),"unsafe job member");const auto path=folder+"/"+source.name;
        Check(lstat(path.c_str(),&source.stamp)==0,"job member disappeared");source.directory=S_ISDIR(source.stamp.st_mode);
        Check(source.directory||S_ISREG(source.stamp.st_mode),"symlink or special file in completed job");
        if(source.directory){sources.push_back(source);Collect(folder,source.name,sources,bytes);}
        else{Check(source.stamp.st_size>=0,"negative source size");if(uint64_t(source.stamp.st_size)>MaximumRecordBytes-bytes)throw Oversize("job exceeds 64 MiB pack record bound; original retained");Descriptor fd(open(path.c_str(),O_RDONLY|O_NOFOLLOW));Check(fd.fd>=0,"cannot open source member");source.bytes=Bytes(fd.fd,0,size_t(source.stamp.st_size));bytes+=source.bytes.size();struct stat after{};Check(fstat(fd.fd,&after)==0&&Same(source.stamp,after),"source changed while packing");sources.push_back(std::move(source));}
    }
}
inline size_t Pack(const std::string& folder){
    std::string root,name,part;Check(Split(folder,root,name,part)&&part.empty(),"not a generated job path");
    struct stat folderStamp{};if(lstat(folder.c_str(),&folderStamp)!=0){Check(errno==ENOENT,"cannot inspect job directory");return 0;}Check(S_ISDIR(folderStamp.st_mode),"job path is not a real directory");
    const auto archive=Shard(root,name);const auto packRoot=CanonicalRoot(root)+"/job_packs";
    std::vector<Source> sources;uint64_t bytes=0;Collect(folder,"",sources,bytes);
    std::lock_guard<std::mutex> lock(Global().mutex);auto index=LoadIndex(archive);
    bool already=index&&index->jobs.count(name);bool complete=false;for(const auto& source:sources)if(source.name=="complete.txt"&&!source.directory)complete=true;
    Check(complete||already,"cannot pack incomplete job");
    if(!already){
        Check(mkdir(packRoot.c_str(),0755)==0||errno==EEXIST,"cannot create archive directory");
        const bool newShard=!index||index->jobs.size()==0;
        if(newShard)FlushDirectory(root);
        if(!index){
            // Publish a valid empty archive atomically. Failed initialization
            // leaves only this temporary file and cannot hide completed data.
            const auto temporary=archive+".init.tmp";Descriptor init(open(temporary.c_str(),O_WRONLY|O_CREAT|O_TRUNC,0644));Check(init.fd>=0,"cannot initialize archive");
            Write(init.fd,FrontierText(1,DataBegin));Write(init.fd,std::string(FrontierSlotBytes,'\0'));Flush(init.fd);
            Check(rename(temporary.c_str(),archive.c_str())==0,"cannot publish empty archive");FlushDirectory(packRoot);index=LoadIndex(archive);
        }
        std::ostringstream metadata;std::string payload;payload.reserve(size_t(bytes));
        for(const auto& source:sources){if(source.directory)metadata<<"D "<<Encode(source.name)<<'\n';else{metadata<<"F "<<Encode(source.name)<<' '<<payload.size()<<' '<<source.bytes.size()<<' '<<Hex(Digest(source.bytes))<<'\n';payload+=source.bytes;}}
        const auto meta=metadata.str();if(meta.size()+payload.size()>MaximumRecordBytes)throw Oversize("packed record exceeds 64 MiB bound; original retained");
        const auto header=Fixed("MBS_JOB_PACK_V1 "+name+" "+std::to_string(meta.size())+" "+std::to_string(payload.size())+" "+Hex(Digest(meta))+" "+Hex(Digest(payload)),HeaderBytes);
        Descriptor fd(open(archive.c_str(),O_RDWR,0644));Check(fd.fd>=0,"cannot open archive for append");const uint64_t end=index->validEnd;
        // Repair a damaged redundant copy before selecting a slot to overwrite.
        // Otherwise the next slot could happen to be the ONLY surviving old
        // frontier, and a torn write would destroy both copies at once.
        MirrorFrontier(fd.fd,index->frontierSequence,index->validEnd);
        Check(ftruncate(fd.fd,off_t(end))==0&&lseek(fd.fd,off_t(end),SEEK_SET)==off_t(end),"cannot repair uncommitted archive tail");
        Write(fd.fd,header);Write(fd.fd,meta);Write(fd.fd,payload);Write(fd.fd,Fixed("MBS_JOB_COMMIT_V1 "+Hex(Digest(header)),TrailerBytes));Flush(fd.fd);
        const uint64_t published=end+HeaderBytes+meta.size()+payload.size()+TrailerBytes;
        Check(Bytes(fd.fd,end,HeaderBytes)==header&&Bytes(fd.fd,end+HeaderBytes,meta.size())==meta&&Bytes(fd.fd,end+HeaderBytes+meta.size(),payload.size())==payload&&Bytes(fd.fd,published-TrailerBytes,TrailerBytes)==Fixed("MBS_JOB_COMMIT_V1 "+Hex(Digest(header)),TrailerBytes),"archive append readback differs");
        const uint64_t sequence=index->frontierSequence+1;Check(sequence>index->frontierSequence,"archive frontier sequence overflow");
        // Alternate two independent filesystem blocks. A torn update keeps the
        // previous frontier; bytes beyond it are uncommitted, regardless of
        // file length or partially persisted headers/trailers.
        const auto watermark=FrontierText(sequence,published);const auto slot=(sequence-1)%2;
        Check(lseek(fd.fd,off_t(slot*FrontierSlotBytes),SEEK_SET)==off_t(slot*FrontierSlotBytes),"cannot seek archive frontier");Write(fd.fd,watermark);Flush(fd.fd);
        const auto verifiedFrontier=ReadFrontier(fd.fd,published);Check(verifiedFrontier.sequence==sequence&&verifiedFrontier.end==published,"archive frontier readback differs");
        // Before losing the loose copy, BOTH independently flushed slots must
        // describe the latest prefix. Later corruption of one slot must not
        // roll back a completed job and make a future append truncate it.
        MirrorFrontier(fd.fd,sequence,published);
        // This writer owns the controller lock and just verified the exact
        // append. Keep the generation of older members; external stamp changes
        // still rebuild the index and invalidate every cached hash.
        auto updated=index?std::const_pointer_cast<Index>(index):std::make_shared<Index>();
        if(!index)updated->generation=++Global().generation;
        updated->jobs.Insert(name,end);updated->validEnd=published;updated->frontierSequence=sequence;
        Check(fstat(fd.fd,&updated->stamp)==0&&uint64_t(updated->stamp.st_size)==updated->validEnd,"archive append length changed");
        Descriptor verify(open(archive.c_str(),O_RDONLY));Check(verify.fd>=0&&Bytes(verify.fd,end,HeaderBytes)==header&&Bytes(verify.fd,updated->validEnd-TrailerBytes,TrailerBytes)==Fixed("MBS_JOB_COMMIT_V1 "+Hex(Digest(header)),TrailerBytes),"archive commit readback differs");
        index=updated;Global().indices.Put(archive,index,sizeof(Index)+index->jobs.Bytes());
    }
    // Read back every stored byte before removing any source. Interrupted
    // removals are idempotent: remaining files are compared to the committed
    // archive. No active/incomplete source is ever considered disposable.
    Descriptor fd(open(archive.c_str(),O_RDWR));Check(fd.fd>=0,"cannot verify archive");uint64_t offset=0;Check(index->jobs.Find(name,offset),"packed job missing");const auto members=Members(fd.fd,archive,ReadRecord(fd.fd,name,offset),index->stamp);
    for(const auto& source:sources){const auto found=members.find(source.name);Check(found!=members.end()&&found->second.directory==source.directory,"source tree changed after archive commit");
        if(!source.directory){const auto saved=Bytes(fd.fd,found->second.offset,size_t(found->second.size));Check(saved==source.bytes&&Hex(Digest(saved))==found->second.hash,"archive readback differs");}
        struct stat now{};Check(lstat((folder+"/"+source.name).c_str(),&now)==0&&Same(now,source.stamp),"source changed before archive cleanup");}
    struct stat archiveAfter{};Check(fstat(fd.fd,&archiveAfter)==0&&Same(archiveAfter,index->stamp),"archive changed before source cleanup");
    if(MirrorFrontier(fd.fd,index->frontierSequence,index->validEnd)){
        auto updated=std::const_pointer_cast<Index>(index);Check(fstat(fd.fd,&updated->stamp)==0,"cannot stat mirrored archive");
        Global().indices.Put(archive,updated,sizeof(Index)+updated->jobs.Bytes());
    }
    size_t removed=0;for(auto at=sources.rbegin();at!=sources.rend();++at){const auto path=folder+"/"+at->name;struct stat current{};Check(lstat(path.c_str(),&current)==0,"source disappeared during cleanup");Check((at->directory?rmdir(path.c_str()):unlink(path.c_str()))==0,"cannot remove verified archived member");if(at->directory||current.st_nlink==1)++removed;}
    Check(rmdir(folder.c_str())==0,"cannot remove archived job directory");FlushDirectory(root+"/jobs");return removed+1;
}
// Explicit export for inspection or a legacy binary. Archives are retained.
// Never overwrite an existing different file or traverse an in-job symlink.
inline size_t Extract(const std::string& folder){
    std::string root,name,part;Check(Split(folder,root,name,part)&&part.empty(),"not a generated job path");
    const auto archive=Shard(root,name);std::lock_guard<std::mutex> lock(Global().mutex);auto index=LoadIndex(archive);uint64_t offset=0;Check(index&&index->jobs.Find(name,offset),"packed job absent");
    Descriptor fd(open(archive.c_str(),O_RDONLY));Check(fd.fd>=0,"cannot export archive");const auto members=Members(fd.fd,archive,ReadRecord(fd.fd,name,offset),index->stamp);
    auto directory=[](const std::string& path){struct stat st{};if(mkdir(path.c_str(),0755)!=0)Check(errno==EEXIST&&lstat(path.c_str(),&st)==0&&S_ISDIR(st.st_mode),"unsafe existing export directory");};
    directory(root+"/jobs");directory(folder);std::vector<std::string> directories{folder};size_t created=0;
    for(const auto& item:members){const auto path=folder+"/"+item.first;if(item.second.directory){directory(path);directories.push_back(path);continue;}
        const auto bytes=Bytes(fd.fd,item.second.offset,size_t(item.second.size));Check(Hex(Digest(bytes))==item.second.hash,"corrupt archived export member");struct stat old{};
        if(lstat(path.c_str(),&old)==0){Check(S_ISREG(old.st_mode),"unsafe existing export member");Descriptor previous(open(path.c_str(),O_RDONLY|O_NOFOLLOW));Check(previous.fd>=0&&uint64_t(old.st_size)==item.second.size&&Bytes(previous.fd,0,size_t(old.st_size))==bytes,"existing export member differs");continue;}
        Check(errno==ENOENT,"cannot inspect export path");const auto temporary=path+".unpack.tmp";Descriptor output(open(temporary.c_str(),O_WRONLY|O_CREAT|O_TRUNC|O_NOFOLLOW,0644));Check(output.fd>=0,"cannot create exported member");Write(output.fd,bytes);Flush(output.fd);Check(rename(temporary.c_str(),path.c_str())==0,"cannot publish exported member");++created;
    }
    for(auto at=directories.rbegin();at!=directories.rend();++at)FlushDirectory(*at);FlushDirectory(root+"/jobs");return created;
}
} // namespace PackedJobStore
