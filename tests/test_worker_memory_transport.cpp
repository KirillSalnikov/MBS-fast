#include "../src/WorkerProtocol.h"
#include "../src/WorkerMemoryStage.h"
#include <sys/wait.h>
#include <cassert>
#include <fstream>
#include <iostream>
#include <functional>

static void Channel(const std::function<void(int)>& writer,const std::function<void(int)>& reader){
    int pair[2];assert(socketpair(AF_UNIX,SOCK_STREAM,0,pair)==0);const pid_t pid=fork();assert(pid>=0);
    if(pid==0){close(pair[0]);writer(pair[1]);close(pair[1]);_exit(0);}
    close(pair[1]);reader(pair[0]);close(pair[0]);int status=0;assert(waitpid(pid,&status,0)==pid&&WIFEXITED(status)&&WEXITSTATUS(status)==0);
}
static bool Exists(const std::string& path){struct stat st{};return lstat(path.c_str(),&st)==0;}
static void Put(const std::string& path,const std::string& bytes){std::ofstream output(path.c_str(),std::ios::binary);assert(output);output.write(bytes.data(),bytes.size());output.close();assert(output);}
static void Frame(int fd,const std::string& path,const std::string& bytes,uint64_t hash,bool truncate=false){
    uint32_t magic=0x4d425342,version=1,count=1,pathBytes=path.size();uint64_t total=path.size()+bytes.size(),size=bytes.size();
    assert(WorkerProtocol::Transfer(fd,&magic,4,true)&&WorkerProtocol::Transfer(fd,&version,4,true)
        &&WorkerProtocol::Transfer(fd,&count,4,true)&&WorkerProtocol::Transfer(fd,&total,8,true)
        &&WorkerProtocol::Transfer(fd,&pathBytes,4,true)&&WorkerProtocol::Transfer(fd,&size,8,true)
        &&WorkerProtocol::Transfer(fd,&hash,8,true)&&WorkerProtocol::Transfer(fd,const_cast<char*>(path.data()),path.size(),true));
    assert(WorkerProtocol::Transfer(fd,const_cast<char*>(bytes.data()),truncate?bytes.size()/2:bytes.size(),true));
}
int main(){
    using WorkerProtocol::File;
    const std::vector<File> files={{"jobs/job/result/result.dat",std::string("123\0\n456",8)},{"jobs/job/stdout.log","ready\n"},{"empty",""}};
    Channel([&](int fd){assert(WorkerProtocol::ReplyMemory(fd,75,files));},[&](int fd){std::vector<File> result;int status=-1;
        assert(WorkerProtocol::ResultMemory(fd,status,result)&&status==75&&result.size()==files.size());
        for(size_t i=0;i<files.size();++i)assert(result[i].path==files[i].path&&result[i].bytes==files[i].bytes);});
    Channel([&](int fd){assert(WorkerProtocol::SendMemoryRequest(fd,"/logical/root",{"/logical/root/jobs/j/stdout.log","64","auto","--output","/logical/root/jobs/j/result"},files));},[&](int fd){
        std::vector<std::string> words;std::vector<File> result;assert(WorkerProtocol::Receive(fd,words)&&words[0]==WorkerProtocol::MemoryRequestTag()&&words[1]=="/logical/root");
        assert(WorkerProtocol::ReceiveFiles(fd,result)&&result.size()==files.size());});
    Channel([&](int fd){assert(WorkerProtocol::Send(fd,{"legacy","64","auto","--close"})&&WorkerProtocol::Reply(fd,0));},[&](int fd){
        std::vector<std::string> words;int status=-1;assert(WorkerProtocol::Receive(fd,words)&&words[0]=="legacy"&&WorkerProtocol::Result(fd,status)&&status==0);});
    uint64_t bytes=0;for(const std::string bad:{"../evil","/absolute","a/../evil","a//b","a/./b","a/"})assert(!WorkerProtocol::CheckFiles({{bad,""}},bytes));
    assert(!WorkerProtocol::CheckFiles({{"a",""},{"a","duplicate"}},bytes));
    File valid={"valid/result","some payload"};
    for(int kind=0;kind<3;++kind)Channel([&](int fd){Frame(fd,kind==2?"other/result":valid.path,valid.bytes,WorkerProtocol::MemberDigest(valid)+(kind==0?1:0),kind==1);},[&](int fd){
        std::vector<File> result={{"old","must disappear"}};assert(!WorkerProtocol::ReceiveFiles(fd,result)&&result.empty());});
    rlimit original{},after{};assert(getrlimit(RLIMIT_FSIZE,&original)==0);std::string stagePath;
    {
        WorkerMemoryStage::Request stage("/logical/root",{{"grids/theta.txt","1\n2\n"},{"means/shared.cache","before"}});stagePath=stage.Root();
        assert(stage.Rewrite("/logical/root/grids/theta.txt")==stagePath+"/grids/theta.txt");
        assert(stage.Rewrite("/logical/root-other/file")=="/logical/root-other/file");
        bool unsafe=false;try{stage.Rewrite("/logical/root/../escape");}catch(const std::runtime_error&){unsafe=true;}assert(unsafe);
        stage.PreparePath("/logical/root/derived/nested/cache");assert(Exists(stagePath+"/derived/nested"));
        stage.PrepareLog("/logical/root/jobs/j/stdout.log");Put(stagePath+"/jobs/j/stdout.log","done");Put(stagePath+"/means/shared.cache","after");
        assert(mkdir((stagePath+"/jobs/j/result").c_str(),0700)==0);Put(stagePath+"/jobs/j/result/result.dat","150 1 2 3\n");
        const auto result=stage.Collect();assert(result.size()==3);
        for(const auto& file:result)assert(file.path!="grids/theta.txt");
    }
    assert(!Exists(stagePath)&&getrlimit(RLIMIT_FSIZE,&after)==0&&original.rlim_cur==after.rlim_cur&&original.rlim_max==after.rlim_max);
    bool rejected=false;
    try{WorkerMemoryStage::Request stage("/logical/root",{});stagePath=stage.Root();assert(symlink("/etc/passwd",(stagePath+"/unsafe").c_str())==0);stage.Collect();}
    catch(const std::runtime_error&){rejected=true;}assert(rejected&&!Exists(stagePath));
    rejected=false;
    try{WorkerMemoryStage::Request stage("/logical/root",{});stagePath=stage.Root();
        for(const std::string name:{"first","second"}){int fd=open((stagePath+"/"+name).c_str(),O_WRONLY|O_CREAT,0600);assert(fd>=0);
            assert(ftruncate(fd,WorkerProtocol::MaximumBundleBytes/2+1)==0);close(fd);}stage.Collect();}
    catch(const std::runtime_error&){rejected=true;}assert(rejected&&!Exists(stagePath));
    // A process killed before stack unwinding leaves only its own private tree;
    // the controller's post-wait cleanup can remove it without following links.
    int pipefd[2];assert(pipe(pipefd)==0);const pid_t child=fork();assert(child>=0);
    if(child==0){close(pipefd[0]);WorkerMemoryStage::Request stage("/logical/root",{{"input","bytes"}});
        const auto root=stage.Root();assert(write(pipefd[1],root.data(),root.size())==ssize_t(root.size()));_exit(0);}
    close(pipefd[1]);char buffer[512];const auto count=read(pipefd[0],buffer,sizeof(buffer));close(pipefd[0]);assert(count>0);stagePath.assign(buffer,size_t(count));int status=0;
    assert(waitpid(child,&status,0)==child&&WIFEXITED(status)&&WEXITSTATUS(status)==0&&Exists(stagePath));
    assert(WorkerMemoryStage::CleanupForPid(child)&&!Exists(stagePath));
    // Isolate startup cleanup from actual /dev/shm worker trees. The same UID
    // may have independent live calculations running during these tests.
    char cleanupPattern[]="/dev/shm/mbs-cleanup-test-XXXXXX";char* made=mkdtemp(cleanupPattern);assert(made);
    const std::string cleanupRoot=made;
    const auto dead=cleanupRoot+"/mbs-worker-memory-"+std::to_string(child)+"-aB1234";
    const auto live=cleanupRoot+"/mbs-worker-memory-"+std::to_string(getpid())+"-cD5678";
    const auto wrongMode=cleanupRoot+"/mbs-worker-memory-"+std::to_string(child)+"-eF9012";
    const auto link=cleanupRoot+"/mbs-worker-memory-"+std::to_string(child)+"-gH3456";
    const auto invalid=cleanupRoot+"/mbs-worker-memory-99999999999999999999999999999-abc123";
    assert(mkdir(dead.c_str(),0700)==0&&mkdir(live.c_str(),0700)==0&&mkdir(wrongMode.c_str(),0755)==0&&mkdir(invalid.c_str(),0700)==0);
    assert(chmod(wrongMode.c_str(),0755)==0); // Preserve this guard even under umask 077.
    Put(dead+"/unfinished","discard only this dead worker");Put(live+"/keep","active PID");Put(wrongMode+"/keep","unexpected mode");
    assert(symlink(live.c_str(),link.c_str())==0);
    assert(symlink(live.c_str(),(dead+"/do-not-follow").c_str())==0);
    assert(WorkerMemoryStage::CleanupDeadWorkersAt(cleanupRoot));
    assert(!Exists(dead)&&Exists(live+"/keep")&&Exists(wrongMode+"/keep")&&Exists(link)&&Exists(invalid));
    assert(WorkerMemoryStage::CleanupDeadWorkersAt(cleanupRoot));
    WorkerMemoryStage::RemoveTree(cleanupRoot);assert(!Exists(cleanupRoot));
    std::cout<<"worker memory transport: binary/legacy roundtrip, bounded authenticated input, EOF/corruption, tmpfs delta/cleanup, oversize, killed-worker cleanup, startup dead/live/symlink/mode guards PASS\n";
}
