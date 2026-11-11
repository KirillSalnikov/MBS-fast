#include "../src/WorkerProtocol.h"
#include "../src/WorkerMemoryStage.h"
#include "../src/FullAutoMemoryStore.h"
#include <sys/wait.h>
#include <cassert>
#include <iostream>
#include <fstream>

using WorkerProtocol::StreamFile;
using WorkerProtocol::StreamReader;
using WorkerProtocol::StreamWriter;
static const uint64_t MiB=1024*1024;
template<class F> bool Throws(F action){try{action();return false;}catch(const std::exception&){return true;}}
static bool Exists(const std::string& path){struct stat st{};return lstat(path.c_str(),&st)==0;}
static uint64_t RepeatHash(uint64_t size,char value){char bytes[65536];std::memset(bytes,value,sizeof(bytes));uint64_t hash=14695981039346656037ULL;
    while(size){const auto n=size_t(std::min<uint64_t>(size,sizeof(bytes)));hash=WorkerProtocol::HashBlock(hash,bytes,n);size-=n;}return hash;}
static void Repeat(uint64_t size,char value,const StreamWriter& writer){char bytes[65536];std::memset(bytes,value,sizeof(bytes));while(size){const auto n=size_t(std::min<uint64_t>(size,sizeof(bytes)));writer(bytes,n);size-=n;}}
template<class Writer,class Reader> void Channel(Writer writer,Reader reader){int pair[2];assert(socketpair(AF_UNIX,SOCK_STREAM,0,pair)==0);pid_t pid=fork();assert(pid>=0);
    if(!pid){close(pair[0]);writer(pair[1]);close(pair[1]);_exit(0);}close(pair[1]);reader(pair[0]);close(pair[0]);int status=0;assert(waitpid(pid,&status,0)==pid&&WIFEXITED(status)&&WEXITSTATUS(status)==0);}
static void Frame(int fd,const StreamFile& file,const std::string& bytes,bool corruptDescriptor=false){
    uint32_t magic=0x4d425353,version=2,count=1,pathBytes=file.path.size();uint64_t total=file.path.size()+file.size,size=file.size,digest=file.digest,descriptor=WorkerProtocol::StreamDescriptorDigest(file)+(corruptDescriptor?1:0);
    for(const auto& piece:std::vector<std::pair<void*,size_t>>{{&magic,4},{&version,4},{&count,4},{&total,8},{&pathBytes,4},{&size,8},{&digest,8},{&descriptor,8}})assert(WorkerProtocol::Transfer(fd,piece.first,piece.second,true));
    assert(WorkerProtocol::Transfer(fd,const_cast<char*>(file.path.data()),file.path.size(),true));
    if(!bytes.empty())WorkerProtocol::Transfer(fd,const_cast<char*>(bytes.data()),bytes.size(),true);
}
int main(){
    char temporary[]="/tmp/mbs-worker-stream-test-XXXXXX";assert(mkdtemp(temporary));const std::string root=temporary;
    const uint64_t inputSize=80*MiB,outputSize=140*MiB;
    const StreamFile input("grids/large",inputSize,RepeatHash(inputSize,'I'));
    const auto expectedOutputHash=RepeatHash(outputSize,'O');
    // Real request/reply stream: 80 MiB input -> native filesystem adapter ->
    // 140 MiB output -> 8 MiB controller cache. Neither endpoint builds a body
    // string or a vector bundle; the output stage must choose disk fallback.
    FullAutoMemoryStore::Store store(root,8*MiB,root);int sockets[2];assert(socketpair(AF_UNIX,SOCK_STREAM,0,sockets)==0);const auto child=fork();assert(child>=0);
    if(!child){close(sockets[0]);std::vector<std::string> words;assert(WorkerProtocol::Receive(sockets[1],words)&&words[0]==WorkerProtocol::MemoryStreamRequestTag());
        std::unique_ptr<WorkerMemoryStage::Request> stage;
        assert(WorkerProtocol::ReceiveStreamFiles(sockets[1],[&](const StreamFile& file,const StreamReader& reader){stage->WriteStream(file,reader);},
            [&](uint64_t total,uint32_t count){assert(count==1);stage.reset(new WorkerMemoryStage::Request(root,total,std::stoull(words[2]),true));}));
        assert(stage->Backend()=="disk"&&stage->CollectStream().empty());
        const auto path=stage->Root()+"/result.dat";const int fd=open(path.c_str(),O_CREAT|O_EXCL|O_WRONLY,0600);assert(fd>=0);
        Repeat(outputSize,'O',[&](const char* bytes,size_t n){size_t offset=0;while(offset<n){const auto done=write(fd,bytes+offset,n-offset);assert(done>0);offset+=size_t(done);}});assert(close(fd)==0);
        const auto files=stage->CollectStream();assert(files.size()==1&&files[0].size==outputSize&&files[0].digest==expectedOutputHash);
        assert(WorkerProtocol::ReplyMemoryStream(sockets[1],0,files,[&](const StreamFile& file,const StreamWriter& writer){stage->StreamFile(file.path,writer);}));
        const auto tree=stage->Root();stage.reset();assert(!Exists(tree));close(sockets[1]);_exit(0);
    }
    close(sockets[1]);assert(WorkerProtocol::SendMemoryStreamRequest(sockets[0],root,160*MiB,{root+"/stdout.log","64","auto","--output",root},
        {input},[&](const StreamFile& file,const StreamWriter& writer){Repeat(file.size,'I',writer);}));
    int status=-1;assert(WorkerProtocol::ResultMemoryStream(sockets[0],status,[&](const StreamFile& file,const StreamReader& reader){store.WriteStream(root+"/"+file.path,file.size,file.digest,reader);})&&status==0);
    close(sockets[0]);int exitStatus=0;assert(waitpid(child,&exitStatus,0)==child&&WIFEXITED(exitStatus)&&WEXITSTATUS(exitStatus)==0);
    assert(store.Size(root+"/result.dat")==outputSize&&store.Digest(root+"/result.dat")==expectedOutputHash&&store.MemoryBytes()<=8*MiB);
    uint64_t readSize=0;store.ReadStream(root+"/result.dat",[&](const char* bytes,size_t n){assert(n<=65536);for(size_t i=0;i<n;++i)assert(bytes[i]=='O');readSize+=n;});assert(readSize==outputSize);
    // Every failed member preserves a previously published Store value.
    const auto target=root+"/small";store.Write(target,"previous");const auto generation=store.Generation(target);
    const std::string payload="new complete bytes";const StreamFile small("small",payload.size(),WorkerProtocol::FileDigest(payload));
    for(int mode=0;mode<3;++mode)Channel([&](int fd){Frame(fd,small,mode==0?payload.substr(0,3):mode==1?std::string(payload.size(),'x'):payload,mode==2);},[&](int fd){
        bool accepted=false;try{accepted=WorkerProtocol::ReceiveStreamFiles(fd,[&](const StreamFile& file,const StreamReader& reader){store.WriteStream(root+"/"+file.path,file.size,file.digest,reader);});}catch(const std::exception&){}
        assert(!accepted&&store.Read(target)=="previous"&&store.Generation(target)==generation);});
    Channel([&](int fd){uint32_t magic=0x4d425353,version=2,count=1;uint64_t total=WorkerProtocol::MaximumStreamBytes+1;
        assert(WorkerProtocol::Transfer(fd,&magic,4,true)&&WorkerProtocol::Transfer(fd,&version,4,true)&&WorkerProtocol::Transfer(fd,&count,4,true)&&WorkerProtocol::Transfer(fd,&total,8,true));},
        [&](int fd){bool headerCalled=false;assert(!WorkerProtocol::ReceiveStreamFiles(fd,[](const StreamFile&,const StreamReader&){assert(false);},[&](uint64_t,uint32_t){headerCalled=true;}));assert(!headerCalled);});
    uint64_t total=0;assert(!WorkerProtocol::CheckStreamFiles({StreamFile("../escape",1,0)},total));
    assert(!WorkerProtocol::CheckStreamFiles({small,small},total));
    // A large shared mean cache can be rewritten in place. The response
    // allowance includes declared inputs, not only newly estimated output.
    {WorkerMemoryStage::Request stage(root,inputSize+32,MiB,true);uint64_t received=0;
        stage.WriteStream(StreamFile("means/cache",inputSize,RepeatHash(inputSize,'I')),[&](char* dst,size_t wanted){
            const auto n=size_t(std::min<uint64_t>(wanted,inputSize-received));std::memset(dst,'I',n);received+=n;return n;});
        const int fd=open((stage.Root()+"/means/cache").c_str(),O_WRONLY|O_TRUNC);assert(fd>=0);
        Repeat(inputSize,'R',[&](const char* bytes,size_t n){assert(write(fd,bytes,n)==ssize_t(n));});close(fd);
        const auto changed=stage.CollectStream();assert(changed.size()==1&&changed[0].size==inputSize&&changed[0].digest==RepeatHash(inputSize,'R'));
    }
    // Small jobs retain tmpfs when space is available; output-size failure is
    // reported as an exception rather than killing the service with SIGXFSZ.
    {WorkerMemoryStage::Request stage(root,0,1024,true);const auto path=stage.Root()+"/oversize";int fd=open(path.c_str(),O_CREAT|O_WRONLY,0600);assert(fd>=0);
        char bytes[2048]{};const auto first=write(fd,bytes,sizeof(bytes));assert(first==1024);assert(write(fd,bytes,1)==-1&&errno==EFBIG);close(fd);assert(Throws([&]{stage.CollectStream();}));}
    // A rejected input is removed before its path can become an input identity.
    {WorkerMemoryStage::Request stage(root,32,1024,true);const StreamFile wrong("input",3,0);size_t offset=0;
        assert(Throws([&]{stage.WriteStream(wrong,[&](char* dst,size_t n){const auto count=std::min(n,size_t(3)-offset);std::memset(dst,'x',count);offset+=count;return count;});}));assert(!Exists(stage.Root()+"/input"));}
    int paths[2];assert(pipe(paths)==0);const auto orphan=fork();assert(orphan>=0);
    if(!orphan){close(paths[0]);WorkerMemoryStage::Request stage(root,0,80*MiB,true);assert(stage.Backend()=="disk");
        const auto path=stage.Root();assert(write(paths[1],path.data(),path.size())==ssize_t(path.size()));_exit(0);}
    close(paths[1]);char pathBytes[4096];const auto pathLength=read(paths[0],pathBytes,sizeof(pathBytes));close(paths[0]);assert(pathLength>0);
    const std::string orphanPath(pathBytes,size_t(pathLength));assert(waitpid(orphan,&exitStatus,0)==orphan&&WIFEXITED(exitStatus)&&WEXITSTATUS(exitStatus)==0&&Exists(orphanPath));
    assert(WorkerMemoryStage::CleanupForPid(orphan,root)&&!Exists(orphanPath));
    assert(WorkerMemoryStage::CleanupDeadWorkersAt(root,".mbs-worker-memory-"));assert(rmdir(root.c_str())==0);
    std::cout<<"PASS stream V2: 80MiB input, 140MiB output, bounded8MiB Store, disk fallback, checksum/descriptor/EOF rejection, exact frame limits, SIGXFSZ recovery, normal/orphan cleanup; V1 unaffected\n";
}
