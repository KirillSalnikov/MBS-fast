#pragma once
// Private controller/worker transport. Length-delimited arguments preserve
// spaces and newlines without shell interpretation; payloads are bounded.
#include <sys/socket.h>
#include <unistd.h>
#include <cerrno>
#include <cstdint>
#include <cstring>
#include <deque>
#include <algorithm>
#include <functional>
#include <stdexcept>
#include <string>
#include <set>
#include <vector>
#include "GpuChunkRetry.h"
namespace WorkerProtocol {
// Retry the identical request once in a fresh process. Never accept its partial output.
const int MemoryPressure = 75;
const int StagingCapacity = 77;
// A wide GPU packet refused its device path before unsafe host fallback.
// A controller may replay identical sampling once with a smaller packet cap.
const int GpuChunkRetry = ::GpuChunkRetry::StatusCode;
inline bool Transfer(int fd,void* buffer,size_t bytes,bool sendData) {
    char* p=static_cast<char*>(buffer);
    while(bytes){ssize_t n=sendData?send(fd,p,bytes,MSG_NOSIGNAL):recv(fd,p,bytes,0);
        if(n<0&&errno==EINTR)continue;
        if(n<=0)return false;
        p+=n;
        bytes-=size_t(n);
        }
    return true;
}
inline bool Send(int fd,const std::vector<std::string>& words){
    uint32_t count=uint32_t(words.size());if(count>4096)return false;
    if(!Transfer(fd,&count,sizeof(count),true))return false;
    for(const auto& s:words){uint32_t n=uint32_t(s.size());if(s.size()>1048576)return false;
        if(!Transfer(fd,&n,sizeof(n),true)||!Transfer(fd,const_cast<char*>(s.data()),n,true))return false;}
    return true;
}
inline bool Receive(int fd,std::vector<std::string>& words){
    words.clear();uint32_t count=0;if(!Transfer(fd,&count,sizeof(count),false)||count>4096)return false;
    size_t total=0;for(uint32_t i=0;i<count;++i){uint32_t n=0;
        if(!Transfer(fd,&n,sizeof(n),false)||n>1048576||(total+=n)>16777216)return false;
        std::string s(n,'\0');if(n&&!Transfer(fd,&s[0],n,false))return false;
        if(s.find('\0')!=std::string::npos)return false;
        words.push_back(std::move(s));
        }
    return true;
}
inline bool Reply(int fd,int code){int32_t value=code;return Transfer(fd,&value,sizeof(value),true);}
inline bool Result(int fd,int& code){int32_t value=0;if(!Transfer(fd,&value,sizeof(value),false))return false;code=value;return true;}

// The optional memory request is identified inside the existing argument
// envelope. Legacy workers/controllers continue to use the status-only reply.
// Every returned member is relative to the controller output root. ReceiveFiles
// publishes no partial bundle on EOF, a malformed member, or a checksum failure.
struct File { std::string path,bytes; };
const uint64_t MaximumBundleBytes=64ULL*1024*1024;
const uint32_t MaximumBundleFiles=4096;
inline const char* MemoryRequestTag(){return "MBS_WORKER_MEMORY_REQUEST_V1";}
inline uint64_t FileDigest(const std::string& bytes){
    uint64_t hash=14695981039346656037ULL;
    for(unsigned char byte:bytes){hash^=byte;hash*=1099511628211ULL;}return hash;
}
inline uint64_t MemberDigest(const File& file){
    uint64_t hash=FileDigest(file.path);hash*=1099511628211ULL; // embedded NUL separates path and payload
    for(unsigned char byte:file.bytes){hash^=byte;hash*=1099511628211ULL;}return hash;
}
inline bool RelativeFilePath(const std::string& path){
    if(path.empty()||path.size()>4096||path.front()=='/'||path.back()=='/'||path.find('\0')!=std::string::npos)return false;
    size_t begin=0;while(begin<path.size()){
        const auto end=path.find('/',begin);const auto name=path.substr(begin,end==std::string::npos?end:end-begin);
        if(name.empty()||name=="."||name=="..")return false;
        if(end==std::string::npos)break;
        begin=end+1;
    }return true;
}
inline bool CheckFiles(const std::vector<File>& files,uint64_t& total){
    total=0;if(files.size()>MaximumBundleFiles)return false;std::set<std::string> paths;
    for(const auto& f:files){if(!RelativeFilePath(f.path)||!paths.insert(f.path).second||f.bytes.size()>MaximumBundleBytes)return false;
        const uint64_t add=uint64_t(f.path.size())+uint64_t(f.bytes.size());
        if(add>MaximumBundleBytes-total)return false;
        total+=add;}
    return true;
}
inline bool SendFiles(int fd,const std::vector<File>& files){
    uint64_t total=0;if(!CheckFiles(files,total))return false;
    uint32_t magic=0x4d425342,version=1,count=uint32_t(files.size());
    if(!Transfer(fd,&magic,sizeof(magic),true)||!Transfer(fd,&version,sizeof(version),true)
        ||!Transfer(fd,&count,sizeof(count),true)||!Transfer(fd,&total,sizeof(total),true))return false;
    for(const auto& f:files){uint32_t pathBytes=uint32_t(f.path.size());uint64_t bytes=f.bytes.size(),digest=MemberDigest(f);
        if(!Transfer(fd,&pathBytes,sizeof(pathBytes),true)||!Transfer(fd,&bytes,sizeof(bytes),true)
            ||!Transfer(fd,&digest,sizeof(digest),true)||!Transfer(fd,const_cast<char*>(f.path.data()),pathBytes,true)
            ||(bytes&&!Transfer(fd,const_cast<char*>(f.bytes.data()),size_t(bytes),true)))return false;}
    return true;
}
inline bool ReceiveFiles(int fd,std::vector<File>& files){
    files.clear();uint32_t magic=0,version=0,count=0;uint64_t expected=0,total=0;
    if(!Transfer(fd,&magic,sizeof(magic),false)||!Transfer(fd,&version,sizeof(version),false)
        ||!Transfer(fd,&count,sizeof(count),false)||!Transfer(fd,&expected,sizeof(expected),false)
        ||magic!=0x4d425342||version!=1||count>MaximumBundleFiles||expected>MaximumBundleBytes)return false;
    std::vector<File> result;std::set<std::string> paths;result.reserve(count);
    for(uint32_t i=0;i<count;++i){uint32_t pathBytes=0;uint64_t bytes=0,digest=0;
        if(!Transfer(fd,&pathBytes,sizeof(pathBytes),false)||!Transfer(fd,&bytes,sizeof(bytes),false)
            ||!Transfer(fd,&digest,sizeof(digest),false)||pathBytes==0||pathBytes>4096||bytes>MaximumBundleBytes
            ||uint64_t(pathBytes)+bytes>expected-total)return false;
        File item;item.path.resize(pathBytes);item.bytes.resize(size_t(bytes));
        if(!Transfer(fd,&item.path[0],pathBytes,false)||!RelativeFilePath(item.path)||!paths.insert(item.path).second
            ||(bytes&&!Transfer(fd,&item.bytes[0],size_t(bytes),false))||MemberDigest(item)!=digest)return false;
        total+=uint64_t(pathBytes)+bytes;result.push_back(std::move(item));
    }
    if(total!=expected)return false;
    files.swap(result);return true;
}
inline bool SendMemoryRequest(int fd,const std::string& logicalRoot,const std::vector<std::string>& ordinaryRequest,const std::vector<File>& inputs){
    uint64_t bytes=0;
    if(logicalRoot.size()<2||logicalRoot.front()!='/'||logicalRoot.back()=='/'||ordinaryRequest.size()<4||!CheckFiles(inputs,bytes))return false;
    std::vector<std::string> words{MemoryRequestTag(),logicalRoot};words.insert(words.end(),ordinaryRequest.begin(),ordinaryRequest.end());
    return Send(fd,words)&&SendFiles(fd,inputs);
}
inline bool ReplyMemory(int fd,int code,const std::vector<File>& files){return Reply(fd,code)&&SendFiles(fd,files);}
inline bool ResultMemory(int fd,int& code,std::vector<File>& files){return Result(fd,code)&&ReceiveFiles(fd,files);}

// V2 transports file bodies in bounded pieces. Metadata is bounded separately;
// a legitimate large angular grid no longer needs a bundle-sized RAM string.
// Consumers must publish a file only after its expected size/hash is verified.
// A complete job is published by the controller only after the whole frame.
const uint64_t MaximumStreamBytes=8ULL*1024*1024*1024;
// Return zero when the remaining frame space cannot accommodate a larger retry.
inline uint64_t LargerOutputBudget(uint64_t previous,uint64_t inputBytes){
    if(inputBytes>=MaximumStreamBytes||!previous)return 0;
    const uint64_t ceiling=MaximumStreamBytes-inputBytes;
    if(previous>=ceiling)return 0;
    return previous>ceiling/4?ceiling:4*previous;
}
const size_t StreamChunkBytes=64*1024;
struct StreamFile {std::string path;uint64_t size=0,digest=14695981039346656037ULL;
    StreamFile()=default;StreamFile(const std::string& p,uint64_t n,uint64_t h):path(p),size(n),digest(h){}};
using StreamReader=std::function<size_t(char*,size_t)>;
using StreamWriter=std::function<void(const char*,size_t)>;
inline const char* MemoryStreamRequestTag(){return "MBS_WORKER_MEMORY_REQUEST_V2";}
inline uint64_t HashBlock(uint64_t hash,const char* bytes,size_t count){
    for(size_t i=0;i<count;++i){hash^=static_cast<unsigned char>(bytes[i]);hash*=1099511628211ULL;}return hash;
}
inline uint64_t StreamDescriptorDigest(const StreamFile& f){
    uint64_t hash=FileDigest(f.path);hash*=1099511628211ULL;
    for(uint64_t value:{f.size,f.digest})for(unsigned i=0;i<8;++i){hash^=(value>>(8*i))&255;hash*=1099511628211ULL;}
    return hash;
}
inline bool CheckStreamFiles(const std::vector<StreamFile>& files,uint64_t& total){
    total=0;if(files.size()>MaximumBundleFiles)return false;std::set<std::string> paths;
    for(const auto& f:files){if(!RelativeFilePath(f.path)||!paths.insert(f.path).second||f.size>MaximumStreamBytes)return false;
        const uint64_t add=uint64_t(f.path.size())+f.size;if(add>MaximumStreamBytes-total)return false;total+=add;}
    return true;
}
template<class Producer> bool SendStreamFiles(int fd,const std::vector<StreamFile>& files,Producer producer){
    uint64_t total=0;if(!CheckStreamFiles(files,total))return false;
    uint32_t magic=0x4d425353,version=2,count=uint32_t(files.size());
    if(!Transfer(fd,&magic,4,true)||!Transfer(fd,&version,4,true)||!Transfer(fd,&count,4,true)||!Transfer(fd,&total,8,true))return false;
    for(const auto& file:files){uint32_t pathBytes=uint32_t(file.path.size());uint64_t size=file.size,digest=file.digest,descriptor=StreamDescriptorDigest(file);
        if(!Transfer(fd,&pathBytes,4,true)||!Transfer(fd,&size,8,true)||!Transfer(fd,&digest,8,true)||!Transfer(fd,&descriptor,8,true)
                ||!Transfer(fd,const_cast<char*>(file.path.data()),pathBytes,true))return false;
        uint64_t sent=0,hash=14695981039346656037ULL;
        const StreamWriter writer=[&](const char* bytes,size_t length){
            if(uint64_t(length)>file.size-sent)throw std::runtime_error("stream producer exceeded member size");
            while(length){const size_t n=std::min(length,StreamChunkBytes);
                if(!Transfer(fd,const_cast<char*>(bytes),n,true))throw std::runtime_error("stream member send failed");
                hash=HashBlock(hash,bytes,n);sent+=n;bytes+=n;length-=n;}
        };
        producer(file,writer);
        if(sent!=file.size||hash!=file.digest)return false;
    }return true;
}
template<class Consumer,class Header> bool ReceiveStreamFiles(int fd,Consumer consumer,Header header){
    uint32_t magic=0,version=0,count=0;uint64_t expected=0,total=0;
    if(!Transfer(fd,&magic,4,false)||!Transfer(fd,&version,4,false)||!Transfer(fd,&count,4,false)||!Transfer(fd,&expected,8,false)
            ||magic!=0x4d425353||version!=2||count>MaximumBundleFiles||expected>MaximumStreamBytes)return false;
    header(expected,count);std::set<std::string> paths;
    for(uint32_t i=0;i<count;++i){uint32_t pathBytes=0;StreamFile file;uint64_t descriptor=0;
        if(!Transfer(fd,&pathBytes,4,false)||!Transfer(fd,&file.size,8,false)||!Transfer(fd,&file.digest,8,false)||!Transfer(fd,&descriptor,8,false)
                ||pathBytes==0||pathBytes>4096||file.size>MaximumStreamBytes||uint64_t(pathBytes)+file.size>expected-total)return false;
        file.path.resize(pathBytes);
        if(!Transfer(fd,&file.path[0],pathBytes,false)||!RelativeFilePath(file.path)||!paths.insert(file.path).second
                ||descriptor!=StreamDescriptorDigest(file))return false;
        uint64_t received=0,hash=14695981039346656037ULL;
        const StreamReader reader=[&](char* bytes,size_t requested)->size_t{
            const auto n=size_t(std::min<uint64_t>(file.size-received,std::min(requested,StreamChunkBytes)));
            if(!n)return 0;
            if(!Transfer(fd,bytes,n,false))throw std::runtime_error("truncated stream member");
            hash=HashBlock(hash,bytes,n);received+=n;return n;
        };
        consumer(file,reader);
        if(received!=file.size||hash!=file.digest)return false;
        total+=uint64_t(pathBytes)+file.size;
    }return total==expected;
}
template<class Consumer> bool ReceiveStreamFiles(int fd,Consumer consumer){
    return ReceiveStreamFiles(fd,consumer,[](uint64_t,uint32_t){});
}
template<class Producer> bool SendMemoryStreamRequest(int fd,const std::string& logicalRoot,uint64_t outputBytes,
        const std::vector<std::string>& request,const std::vector<StreamFile>& inputs,Producer producer){
    uint64_t bytes=0;if(logicalRoot.size()<2||logicalRoot.front()!='/'||logicalRoot.back()=='/'||request.size()<4
            ||!outputBytes||outputBytes>MaximumStreamBytes||!CheckStreamFiles(inputs,bytes)||outputBytes>MaximumStreamBytes-bytes)return false;
    std::vector<std::string> words{MemoryStreamRequestTag(),logicalRoot,std::to_string(outputBytes)};words.insert(words.end(),request.begin(),request.end());
    return Send(fd,words)&&SendStreamFiles(fd,inputs,producer);
}
template<class Producer> bool ReplyMemoryStream(int fd,int code,const std::vector<StreamFile>& files,Producer producer){
    return Reply(fd,code)&&SendStreamFiles(fd,files,producer);
}
template<class Consumer> bool ResultMemoryStream(int fd,int& code,Consumer consumer){return Result(fd,code)&&ReceiveStreamFiles(fd,consumer);}
inline bool ReplyMemoryStreamSmall(int fd,int code,const std::vector<File>& files){
    uint64_t total=0;if(!CheckFiles(files,total))return false;std::vector<StreamFile> descriptions;
    for(const auto& file:files)descriptions.emplace_back(file.path,file.bytes.size(),FileDigest(file.bytes));
    return ReplyMemoryStream(fd,code,descriptions,[&](const StreamFile& member,const StreamWriter& writer){
        for(const auto& file:files)if(file.path==member.path){writer(file.bytes.data(),file.bytes.size());return;}
    });
}

struct CompletionDispatchStats {
    size_t readyBatches=0,readyBatchPeak=0,pendingPeak=0;
};
// A ready GPU must not wait behind another GPU's checkpoint/callback. First
// consume the complete ready set (reap() also refills its worker), then commit
// one result. Recheck between commits so a long CPU callback cannot postpone
// every subsequent dispatch. Bound the queue to < 2*slots completed requests:
// no unbounded backlog of uncommitted outputs or extra concurrent physics jobs.
// pollReady(block) returns slot IDs; reap(slot,index) returns false for a
// recovered request that was relaunched and has not completed successfully.
template<class PollReady,class Reap,class Commit>
CompletionDispatchStats DispatchCompletions(size_t total,size_t slots,PollReady pollReady,Reap reap,Commit commit) {
    CompletionDispatchStats stats;std::deque<size_t> pending;size_t committed=0;
    while(committed<total){
        if(pending.size()<slots){
            const auto ready=pollReady(pending.empty());
            if(!ready.empty()){++stats.readyBatches;if(ready.size()>stats.readyBatchPeak)stats.readyBatchPeak=ready.size();}
            for(size_t slot:ready){size_t index=0;if(reap(slot,index))pending.push_back(index);}
            if(pending.size()>stats.pendingPeak)stats.pendingPeak=pending.size();
        }
        if(!pending.empty()){const size_t index=pending.front();commit(index);pending.pop_front();++committed;}
    }
    return stats;
}
}
