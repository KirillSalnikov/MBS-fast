#pragma once
// A request-local RAM filesystem adapter for the existing ostream-based native
// solver. Bounded, checksummed chunks cross the worker socket. Small requests
// use tmpfs; large declared outputs use disposable disk staging without fsync.
#include "WorkerProtocol.h"
#include "HostMemoryReserve.h"
#include <sys/stat.h>
#include <sys/statfs.h>
#include <sys/statvfs.h>
#include <sys/resource.h>
#include <signal.h>
#include <fcntl.h>
#include <dirent.h>
#include <algorithm>
#include <map>
#include <stdexcept>
#include <limits>
#include <fstream>
#include <ostream>

namespace WorkerMemoryStage {
inline void Require(bool okay,const std::string& message){if(!okay)throw std::runtime_error("memory worker: "+message);}
struct CapacityExceeded : std::runtime_error {
    explicit CapacityExceeded(const char* reason):std::runtime_error(std::string("memory worker: ")+reason){}
};
inline void RemoveTree(const std::string& path) noexcept {
    struct stat st{};if(lstat(path.c_str(),&st)!=0)return;
    if(S_ISDIR(st.st_mode)){
        DIR* directory=opendir(path.c_str());if(directory){while(auto* entry=readdir(directory)){
            const std::string name=entry->d_name;if(name!="."&&name!="..")RemoveTree(path+"/"+name);}
            closedir(directory);}
        rmdir(path.c_str());
    }else unlink(path.c_str());
}
// A killed solver cannot run C++ destructors. Its owning controller may call
// this only after reaping that worker, before reusing its service slot.
inline bool CleanupForPid(pid_t pid) noexcept {
    if(pid<=0)return false;
    const std::string prefix="mbs-worker-memory-"+std::to_string(pid)+"-";
    DIR* directory=opendir("/dev/shm");if(!directory)return false;bool okay=true;
    while(auto* entry=readdir(directory)){
        const std::string name=entry->d_name;if(name.compare(0,prefix.size(),prefix)!=0)continue;
        const std::string path="/dev/shm/"+name;struct stat st{};
        if(lstat(path.c_str(),&st)!=0||!S_ISDIR(st.st_mode)||st.st_uid!=geteuid()||(st.st_mode&0777)!=0700){okay=false;continue;}
        RemoveTree(path);if(lstat(path.c_str(),&st)==0||errno!=ENOENT)okay=false;
    }
    closedir(directory);return okay;
}
// A controller killed with its workers cannot execute CleanupForPid(). A later
// memory controller can reclaim only unmistakable, owned private trees whose
// PID no longer exists. A live/reused PID, EPERM, symlink or unexpected mode is
// deliberately preserved. The directory argument isolates the regression test;
// production always uses the fixed /dev/shm wrapper below.
inline bool CleanupDeadWorkersAt(const std::string& base,const std::string& prefix="mbs-worker-memory-") noexcept {
    DIR* directory=opendir(base.c_str());if(!directory)return false;bool okay=true;
    try{
        for(;;){
            errno=0;auto* entry=readdir(directory);if(!entry){if(errno)okay=false;break;}
            const std::string name=entry->d_name;
            if(name.compare(0,prefix.size(),prefix)!=0)continue;
            const size_t dash=name.find('-',prefix.size());
            if(dash==std::string::npos||dash==prefix.size()||name.size()-dash-1!=6)continue;
            bool valid=true;uint64_t number=0;
            for(size_t i=prefix.size();i<dash;++i){
                const unsigned char c=name[i];
                if(c<'0'||c>'9'||number>(uint64_t(std::numeric_limits<pid_t>::max())-uint64_t(c-'0'))/10){valid=false;break;}
                number=number*10+uint64_t(c-'0');
            }
            for(size_t i=dash+1;i<name.size();++i){const unsigned char c=name[i];
                if(!((c>='0'&&c<='9')||(c>='a'&&c<='z')||(c>='A'&&c<='Z')))valid=false;}
            if(!valid||number==0)continue;
            const auto pid=pid_t(number);errno=0;
            if(kill(pid,0)==0||errno!=ESRCH)continue;
            const std::string path=base+"/"+name;struct stat before{},after{};
            if(lstat(path.c_str(),&before)!=0){if(errno!=ENOENT)okay=false;continue;}
            if(!S_ISDIR(before.st_mode)||before.st_uid!=geteuid()||(before.st_mode&07777)!=0700)continue;
            // Recheck both identity and liveness immediately before removal.
            errno=0;if(kill(pid,0)==0||errno!=ESRCH)continue;
            if(lstat(path.c_str(),&after)!=0){if(errno!=ENOENT)okay=false;continue;}
            if(before.st_dev!=after.st_dev||before.st_ino!=after.st_ino
                ||!S_ISDIR(after.st_mode)||after.st_uid!=geteuid()||(after.st_mode&07777)!=0700)continue;
            RemoveTree(path);
            if(lstat(path.c_str(),&after)==0||errno!=ENOENT)okay=false;
        }
    }catch(...){okay=false;}
    closedir(directory);return okay;
}
inline bool CleanupDeadWorkers() noexcept {return CleanupDeadWorkersAt("/dev/shm");}
inline bool CleanupDeadWorkers(const std::string& logicalRoot) noexcept {
    const bool memory=CleanupDeadWorkers();return CleanupDeadWorkersAt(logicalRoot,".mbs-worker-memory-")&&memory;
}
inline bool CleanupForPid(pid_t pid,const std::string& logicalRoot) noexcept {
    // Restrict the disk scan to the reaped worker's own exact PID prefix.
    const bool memory=CleanupForPid(pid);if(pid<=0)return false;
    const auto prefix=".mbs-worker-memory-"+std::to_string(pid)+"-";DIR* dir=opendir(logicalRoot.c_str());if(!dir)return false;
    bool okay=true;while(auto* entry=readdir(dir)){const std::string name=entry->d_name;
        if(name.compare(0,prefix.size(),prefix)!=0)continue;
        const auto path=logicalRoot+"/"+name;struct stat st{};
        if(lstat(path.c_str(),&st)!=0||!S_ISDIR(st.st_mode)||st.st_uid!=geteuid()||(st.st_mode&07777)!=0700){okay=false;continue;}
        RemoveTree(path);if(lstat(path.c_str(),&st)==0||errno!=ENOENT)okay=false;}
    closedir(dir);return okay&&memory;
}
inline volatile sig_atomic_t& FileLimitHit(){static volatile sig_atomic_t value=0;return value;}
inline void FileLimitSignal(int){FileLimitHit()=1;}
class Request {
    std::string logicalRoot,root;
    std::map<std::string,std::pair<uint64_t,uint64_t>> originals;
    rlimit savedLimit{};bool limitChanged=false;
    bool streaming=false,signalChanged=false;struct sigaction savedSignal{};
    uint64_t stagingBudget=WorkerProtocol::MaximumBundleBytes,outputBudget=WorkerProtocol::MaximumBundleBytes,inputWritten=0;
    std::string backend="tmpfs";
    static uint64_t Available(const std::string& path){struct statvfs space{};
        Require(statvfs(path.c_str(),&space)==0,"cannot inspect staging capacity");
        Require(!space.f_frsize||uint64_t(space.f_bavail)<=UINT64_MAX/uint64_t(space.f_frsize),"staging capacity overflow");
        return uint64_t(space.f_bavail)*uint64_t(space.f_frsize);}
    void Parents(const std::string& relative){
        size_t slash=0;unsigned depth=0;
        while((slash=relative.find('/',slash))!=std::string::npos){
            Require(++depth<=64,"path nesting exceeds bound");const auto directory=root+"/"+relative.substr(0,slash);
            if(mkdir(directory.c_str(),0700)!=0){struct stat st{};Require(errno==EEXIST&&lstat(directory.c_str(),&st)==0&&S_ISDIR(st.st_mode),"unsafe staging directory");}
            ++slash;
        }
    }
    void Write(const WorkerProtocol::File& file){
        Parents(file.path);const auto path=root+"/"+file.path;
        const int fd=open(path.c_str(),O_WRONLY|O_CREAT|O_EXCL|O_NOFOLLOW|O_CLOEXEC,0600);
        Require(fd>=0,"cannot create staging input "+file.path);
        size_t offset=0;while(offset<file.bytes.size()){
            const ssize_t count=write(fd,file.bytes.data()+offset,file.bytes.size()-offset);
            if(count<0&&errno==EINTR)continue;
            if(count<=0){close(fd);throw std::runtime_error("memory worker: cannot write staging input");}offset+=size_t(count);
        }
        Require(close(fd)==0,"cannot close staging input");
        originals.emplace(file.path,std::make_pair(uint64_t(file.bytes.size()),WorkerProtocol::FileDigest(file.bytes)));
    }
    void CollectDirectory(const std::string& relative,std::vector<WorkerProtocol::File>& result,uint64_t& total,size_t& visited,unsigned depth)const{
        Require(depth<=64,"output nesting exceeds bound");const std::string folder=relative.empty()?root:root+"/"+relative;
        DIR* directory=opendir(folder.c_str());Require(directory,"cannot enumerate staging output");
        std::vector<std::string> names;
        errno=0;while(auto* entry=readdir(directory)){
            const std::string name=entry->d_name;if(name=="."||name=="..")continue;
            if(++visited>WorkerProtocol::MaximumBundleFiles*2){closedir(directory);throw std::runtime_error("memory worker: too many staging entries");}
            names.push_back(name);
        }
        const int error=errno;closedir(directory);Require(error==0,"cannot read staging directory");std::sort(names.begin(),names.end());
        for(const auto& name:names){const auto member=relative.empty()?name:relative+"/"+name;const auto path=root+"/"+member;
            struct stat before{},after{};Require(lstat(path.c_str(),&before)==0,"staging output disappeared");
            if(S_ISDIR(before.st_mode)){CollectDirectory(member,result,total,visited,depth+1);continue;}
            Require(S_ISREG(before.st_mode)&&before.st_nlink==1&&before.st_size>=0,"unsafe staging output "+member);
            const uint64_t bytes=uint64_t(before.st_size);
            Require(bytes<=WorkerProtocol::MaximumBundleBytes,"output member exceeds 64 MiB bound");
            const int fd=open(path.c_str(),O_RDONLY|O_NOFOLLOW|O_CLOEXEC);Require(fd>=0,"cannot read staging output");
            WorkerProtocol::File file;file.path=member;file.bytes.resize(size_t(bytes));size_t offset=0;
            while(offset<file.bytes.size()){const ssize_t count=read(fd,&file.bytes[offset],file.bytes.size()-offset);
                if(count<0&&errno==EINTR)continue;
                if(count<=0){close(fd);throw std::runtime_error("memory worker: truncated staging output");}offset+=size_t(count);}
            const bool stable=fstat(fd,&after)==0&&before.st_dev==after.st_dev&&before.st_ino==after.st_ino&&before.st_size==after.st_size
                &&before.st_mtim.tv_sec==after.st_mtim.tv_sec&&before.st_mtim.tv_nsec==after.st_mtim.tv_nsec;
            close(fd);Require(stable,"staging output changed while reading");
            const auto original=originals.find(member);
            if(original!=originals.end()&&original->second.first==bytes&&original->second.second==WorkerProtocol::FileDigest(file.bytes))continue;
            Require(uint64_t(member.size())+bytes<=WorkerProtocol::MaximumBundleBytes-total,"output bundle exceeds 64 MiB bound");
            total+=uint64_t(member.size())+bytes;result.push_back(std::move(file));
            Require(result.size()<=WorkerProtocol::MaximumBundleFiles,"too many output members");
        }
    }
    static bool Stable(const struct stat& a,const struct stat& b){return a.st_dev==b.st_dev&&a.st_ino==b.st_ino&&a.st_size==b.st_size
        &&a.st_mtim.tv_sec==b.st_mtim.tv_sec&&a.st_mtim.tv_nsec==b.st_mtim.tv_nsec;}
    void CollectStreamDirectory(const std::string& relative,std::vector<WorkerProtocol::StreamFile>& result,
            uint64_t& total,uint64_t& staged,size_t& visited,unsigned depth)const{
        Require(depth<=64,"output nesting exceeds bound");const auto folder=relative.empty()?root:root+"/"+relative;
        DIR* directory=opendir(folder.c_str());Require(directory,"cannot enumerate staging output");std::vector<std::string> names;
        errno=0;while(auto* entry=readdir(directory)){const std::string name=entry->d_name;if(name=="."||name=="..")continue;
            if(++visited>WorkerProtocol::MaximumBundleFiles*2){closedir(directory);throw std::runtime_error("memory worker: too many staging entries");}names.push_back(name);}
        const int error=errno;closedir(directory);Require(!error,"cannot read staging directory");std::sort(names.begin(),names.end());
        for(const auto& name:names){const auto member=relative.empty()?name:relative+"/"+name,path=root+"/"+member;
            struct stat before{};Require(lstat(path.c_str(),&before)==0,"staging output disappeared");
            if(S_ISDIR(before.st_mode)){CollectStreamDirectory(member,result,total,staged,visited,depth+1);continue;}
            Require(S_ISREG(before.st_mode)&&before.st_nlink==1&&before.st_size>=0,"unsafe staging output "+member);
            const uint64_t bytes=uint64_t(before.st_size);if(bytes>stagingBudget-staged)throw CapacityExceeded("native outputs exceeded declared staging budget");
            staged+=bytes;
            uint64_t digest=14695981039346656037ULL;StreamFile(member,[&](const char* data,size_t n){digest=WorkerProtocol::HashBlock(digest,data,n);});
            const auto old=originals.find(member);if(old!=originals.end()&&old->second.first==bytes&&old->second.second==digest)continue;
            if(uint64_t(member.size())+bytes>outputBudget-total)throw CapacityExceeded("output exceeds declared stream budget");
            total+=uint64_t(member.size())+bytes;
            result.emplace_back(member,bytes,digest);Require(result.size()<=WorkerProtocol::MaximumBundleFiles,"too many output members");
        }
    }
public:
    Request(const std::string& logical,const std::vector<WorkerProtocol::File>& inputs):logicalRoot(logical){
        Require(logical.size()>1&&logical.front()=='/'&&logical.back()!='/'&&logical.find('\0')==std::string::npos,"invalid logical output root");
        uint64_t inputBytes=0;Require(WorkerProtocol::CheckFiles(inputs,inputBytes),"invalid input bundle");
        struct statfs fs{};struct statvfs space{};
        Require(statfs("/dev/shm",&fs)==0&&static_cast<unsigned long>(fs.f_type)==0x01021994UL,"/dev/shm must be a RAM filesystem");
        Require(statvfs("/dev/shm",&space)==0,"cannot inspect RAM filesystem capacity");
        const uint64_t available=uint64_t(space.f_bavail)*uint64_t(space.f_frsize);
        Require(available>=inputBytes+4ULL*1024*1024,"insufficient free RAM filesystem space");
        std::string pattern="/dev/shm/mbs-worker-memory-"+std::to_string(getpid())+"-XXXXXX";
        std::vector<char> buffer(pattern.begin(),pattern.end());buffer.push_back('\0');
        char* created=mkdtemp(buffer.data());Require(created,"cannot create private RAM staging directory");root=created;
        try{
            // A corrupt output size cannot consume the whole RAM filesystem in
            // a single file. The complete bundle also has a 64 MiB sum bound.
            Require(getrlimit(RLIMIT_FSIZE,&savedLimit)==0,"cannot query file-size bound");
            rlimit limited=savedLimit;const rlim_t bound=rlim_t(std::min<uint64_t>(WorkerProtocol::MaximumBundleBytes,available-inputBytes));
            limited.rlim_cur=std::min(savedLimit.rlim_cur,bound);
            Require(setrlimit(RLIMIT_FSIZE,&limited)==0,"cannot bound staging file size");limitChanged=true;
            for(const auto& file:inputs)Write(file);
        }catch(...){if(limitChanged)setrlimit(RLIMIT_FSIZE,&savedLimit);RemoveTree(root);root.clear();throw;}
    }
    // V2 uses RAM only for small declared jobs. Large outputs use disposable
    // files beside the controller's anonymous spill, without durability flushes.
    Request(const std::string& logical,uint64_t inputBytes,uint64_t expectedOutput,bool streamMode):logicalRoot(logical),streaming(streamMode){
        Require(streaming&&logical.size()>1&&logical.front()=='/'&&WorkerProtocol::RelativeFilePath(logical.substr(1)),"invalid streaming root");
        Require(expectedOutput>0&&expectedOutput<=WorkerProtocol::MaximumStreamBytes&&inputBytes<=WorkerProtocol::MaximumStreamBytes-expectedOutput,"invalid staging estimate");
        // A shared input cache may be rewritten and returned in full. Its
        // previous size is already in the input allowance; do not count only
        // newly produced bytes when bounding the response.
        stagingBudget=inputBytes+expectedOutput;outputBudget=stagingBudget;
        const auto memory=HostMemoryReserve::Read();const uint64_t reserve=uint64_t(std::max(0LL,HostMemoryReserve::ReserveKb(memory.totalKb)))*1024;
        const uint64_t ram=uint64_t(std::max(0LL,memory.availableKb))*1024;
        struct statfs fs{};bool ramFs=statfs("/dev/shm",&fs)==0&&static_cast<unsigned long>(fs.f_type)==0x01021994UL;
        uint64_t available=0;if(ramFs)available=Available("/dev/shm");
        const uint64_t margin=64ULL*1024*1024;
        const bool useRam=ramFs&&stagingBudget<=WorkerProtocol::MaximumBundleBytes
            &&available>=stagingBudget+margin&&ram>=reserve+stagingBudget+margin;
        backend=useRam?"tmpfs":"disk";const auto base=useRam?std::string("/dev/shm"):logical;
        Require(Available(base)>=stagingBudget+margin,"insufficient free staging space; requires declared inputs/outputs plus 64 MiB reserve");
        const auto pattern=base+(useRam?"/mbs-worker-memory-":"/.mbs-worker-memory-")+std::to_string(getpid())+"-XXXXXX";
        std::vector<char> name(pattern.begin(),pattern.end());name.push_back('\0');char* created=mkdtemp(name.data());Require(created,"cannot create streaming staging directory");root=created;
        try{
            Require(getrlimit(RLIMIT_FSIZE,&savedLimit)==0,"cannot query file-size bound");rlimit limit=savedLimit;
            limit.rlim_cur=std::min(savedLimit.rlim_cur,rlim_t(stagingBudget));Require(setrlimit(RLIMIT_FSIZE,&limit)==0,"cannot bound staging file size");limitChanged=true;
            struct sigaction action{};action.sa_handler=FileLimitSignal;sigemptyset(&action.sa_mask);FileLimitHit()=0;
            Require(sigaction(SIGXFSZ,&action,&savedSignal)==0,"cannot guard staging file-size signal");signalChanged=true;
        }catch(...){if(limitChanged)setrlimit(RLIMIT_FSIZE,&savedLimit);RemoveTree(root);root.clear();throw;}
    }
    Request(const Request&)=delete;Request& operator=(const Request&)=delete;
    ~Request(){if(signalChanged)sigaction(SIGXFSZ,&savedSignal,nullptr);if(limitChanged)setrlimit(RLIMIT_FSIZE,&savedLimit);if(!root.empty())RemoveTree(root);}
    const std::string& Root()const{return root;}
    const std::string& Backend()const{return backend;}
    uint64_t Budget()const{return stagingBudget;}
    void WriteStream(const WorkerProtocol::StreamFile& file,const WorkerProtocol::StreamReader& reader){
        Require(streaming&&WorkerProtocol::RelativeFilePath(file.path)&&!originals.count(file.path),"invalid/duplicate streaming input");
        Require(file.size<=stagingBudget-inputWritten,"input exceeds staging budget");Parents(file.path);const auto path=root+"/"+file.path;
        int fd=open(path.c_str(),O_WRONLY|O_CREAT|O_EXCL|O_NOFOLLOW|O_CLOEXEC,0600);Require(fd>=0,"cannot create streaming input");
        try{uint64_t count=0,digest=14695981039346656037ULL;char buffer[WorkerProtocol::StreamChunkBytes];
            while(count<file.size){const size_t want=size_t(std::min<uint64_t>(sizeof(buffer),file.size-count)),n=reader(buffer,want);
                Require(n>0&&n<=want,"truncated streaming input");digest=WorkerProtocol::HashBlock(digest,buffer,n);size_t offset=0;
                while(offset<n){const auto done=write(fd,buffer+offset,n-offset);if(done<0&&errno==EINTR)continue;Require(done>0,"cannot write streaming input");offset+=size_t(done);}count+=n;}
            Require(digest==file.digest,"streaming input checksum changed");const int closed=close(fd);fd=-1;Require(closed==0,"cannot close streaming input");
            originals.emplace(file.path,std::make_pair(file.size,digest));inputWritten+=file.size;
        }catch(...){if(fd>=0)close(fd);unlink(path.c_str());throw;}
    }
    void StreamFile(const std::string& member,const WorkerProtocol::StreamWriter& writer)const{
        Require(WorkerProtocol::RelativeFilePath(member),"invalid streaming output path");const auto path=root+"/"+member;struct stat before{},after{};
        Require(lstat(path.c_str(),&before)==0&&S_ISREG(before.st_mode)&&before.st_nlink==1&&before.st_size>=0&&uint64_t(before.st_size)<=stagingBudget,"invalid streaming output");
        int fd=open(path.c_str(),O_RDONLY|O_NOFOLLOW|O_CLOEXEC);Require(fd>=0,"cannot open streaming output");
        try{Require(fstat(fd,&after)==0&&Stable(before,after),"staging output identity changed");char buffer[WorkerProtocol::StreamChunkBytes];uint64_t offset=0;
            while(offset<uint64_t(before.st_size)){const size_t want=size_t(std::min<uint64_t>(sizeof(buffer),uint64_t(before.st_size)-offset));
                const ssize_t n=read(fd,buffer,want);if(n<0&&errno==EINTR)continue;Require(n>0,"truncated streaming output");writer(buffer,size_t(n));offset+=uint64_t(n);}
            Require(fstat(fd,&after)==0&&Stable(before,after),"staging output changed during read");close(fd);
        }catch(...){close(fd);throw;}
    }
    std::vector<WorkerProtocol::StreamFile> CollectStream()const{
        Require(streaming,"stream collection requires streaming mode");
        if(FileLimitHit()){
            Require(savedLimit.rlim_cur>stagingBudget,"native output reached inherited file-size limit");
            throw CapacityExceeded("native output reached staging file-size limit");
        }
        std::vector<WorkerProtocol::StreamFile> result;uint64_t total=0,staged=0;size_t visited=0;CollectStreamDirectory("",result,total,staged,visited,0);return result;
    }
    std::string InputIdentityPath(const std::string& physical)const{
        if(physical.compare(0,root.size()+1,root+"/")!=0)return physical;
        const auto relative=physical.substr(root.size()+1);
        Require(WorkerProtocol::RelativeFilePath(relative)&&originals.count(relative),"prepared particle identity requires a declared staging input");
        return logicalRoot+"/"+relative;
    }
    bool Contains(const std::string& path)const{return path==logicalRoot||(path.size()>logicalRoot.size()&&path.compare(0,logicalRoot.size()+1,logicalRoot+"/")==0);}
    std::string Rewrite(const std::string& path)const{
        if(!Contains(path))return path;
        Require(path==logicalRoot||WorkerProtocol::RelativeFilePath(path.substr(logicalRoot.size()+1)),"unsafe logical input/output path");
        return root+path.substr(logicalRoot.size());
    }
    void PreparePath(const std::string& path){
        if(!Contains(path)||path==logicalRoot)return;
        const auto relative=path.substr(logicalRoot.size()+1);Require(WorkerProtocol::RelativeFilePath(relative),"unsafe logical input/output path");Parents(relative);
    }
    void PrepareLog(const std::string& logical){
        Require(Contains(logical)&&logical!=logicalRoot,"stdout path is outside logical root");const auto relative=logical.substr(logicalRoot.size()+1);
        Require(WorkerProtocol::RelativeFilePath(relative),"unsafe stdout path");Parents(relative);
    }
    std::vector<WorkerProtocol::File> Collect()const{
        std::vector<WorkerProtocol::File> result;uint64_t total=0;size_t visited=0;CollectDirectory("",result,total,visited,0);return result;
    }
};
// Only the prepared-particle identity uses this mapping. The native solver
// continues opening the physical tmpfs files, and geometry/visibility bytes
// remain part of the cache key. The scope spans exactly one RunCalculation().
inline const Request*& IdentityRequest(){static thread_local const Request* value=nullptr;return value;}
struct IdentityScope {
    const Request* previous;
    explicit IdentityScope(const Request* request):previous(IdentityRequest()){IdentityRequest()=request;}
    ~IdentityScope(){IdentityRequest()=previous;}
    IdentityScope(const IdentityScope&)=delete;IdentityScope& operator=(const IdentityScope&)=delete;
};
inline std::string PreparedInputPath(const std::string& physical){const auto* request=IdentityRequest();return request?request->InputIdentityPath(physical):physical;}
inline void AppendPreparedParticleIdentity(std::ostream& key,const std::string& physical){
    const auto identity=PreparedInputPath(physical);key<<identity.size()<<':'<<identity;
    std::ifstream geometry(physical,std::ios::binary);const std::string bytes((std::istreambuf_iterator<char>(geometry)),std::istreambuf_iterator<char>());
    key<<bytes.size()<<':'<<bytes;
    std::ifstream sidecar(physical+".visibility",std::ios::binary);const std::string visibility((std::istreambuf_iterator<char>(sidecar)),std::istreambuf_iterator<char>());
    key<<visibility.size()<<':'<<visibility;
}
}
