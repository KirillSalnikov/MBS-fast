#include "../src/GpuChunkRetry.h"
#include "../src/WorkerProtocol.h"
#include <cassert>
#include <thread>
#include <iostream>
using Words=std::vector<std::string>;
Words Physics(Words words){for(size_t i=0;i+1<words.size();)if(words[i]=="--orientation-chunk")words.erase(words.begin()+i,words.begin()+i+2);else ++i;return words;}
int main(){
    const Words physical={"--particle-file","shape with spaces.dat","--sobol-seed","8192","17","--sobol-offset","1024","--max-reflections","22","--paired-depth","20","--orientation-prefix-output","4096","--theta-grid-file","grid.csv","--phi-points","32","--incoherent"};
    for(int explicitCap:{0,1,64,128,512,2147483647}){auto words=physical;if(explicitCap)words.insert(words.end(),{"--orientation-chunk",std::to_string(explicitCap)});
        const auto before=Physics(words);const int cap=GpuChunkRetry::LimitCommand(words,0);assert(cap==std::min(128,explicitCap?explicitCap:128));assert(Physics(words)==before);
        assert(GpuChunkRetry::LimitCommand(words,1024)==cap&&Physics(words)==before);}
    auto duplicate=physical;duplicate.insert(duplicate.end(),{"--orientation-chunk","1024","--orientation-chunk","32"});
    assert(GpuChunkRetry::LimitCommand(duplicate,64)==32);assert(Physics(duplicate)==physical);
    for(const std::string invalid:{"0","-1","2147483648","bad",""}){auto words=physical;words.insert(words.end(),{"--orientation-chunk",invalid});const auto old=words;
        bool rejected=false;try{GpuChunkRetry::LimitCommand(words,0);}catch(const std::runtime_error&){rejected=true;}assert(rejected&&words==old);}
    assert(GpuChunkRetry::ExitStatus(false)==1&&GpuChunkRetry::ExitStatus(true)==76);
    assert(WorkerProtocol::GpuChunkRetry==76&&WorkerProtocol::GpuChunkRetry!=WorkerProtocol::MemoryPressure);
    try{throw GpuChunkRetry::Required(2048,"test refusal");}catch(const GpuChunkRetry::Required& e){assert(std::string(e.what()).find("2048")!=std::string::npos);}
    for(int version:{0,1,2}){int fd[2];assert(socketpair(AF_UNIX,SOCK_STREAM,0,fd)==0);const std::vector<WorkerProtocol::File> files{{"job/stdout.log","refused before host preparation\n"}};
        std::thread writer([&]{bool ok=version==0?WorkerProtocol::Reply(fd[0],76):version==1?WorkerProtocol::ReplyMemory(fd[0],76,files):WorkerProtocol::ReplyMemoryStreamSmall(fd[0],76,files);assert(ok);close(fd[0]);});
        int status=-1;std::string diagnostic;bool ok=false;
        if(!version)ok=WorkerProtocol::Result(fd[1],status);
        else if(version==1){std::vector<WorkerProtocol::File> received;ok=WorkerProtocol::ResultMemory(fd[1],status,received);assert(received.size()==1);diagnostic=received[0].bytes;}
        else ok=WorkerProtocol::ResultMemoryStream(fd[1],status,[&](const WorkerProtocol::StreamFile&,const WorkerProtocol::StreamReader& read){char block[16];for(size_t n;(n=read(block,sizeof(block)));)diagnostic.append(block,n);});
        assert(ok&&status==76&&(!version||diagnostic==files[0].bytes));close(fd[1]);writer.join();}
    {int fd[2];assert(socketpair(AF_UNIX,SOCK_STREAM,0,fd)==0);assert(WorkerProtocol::Reply(fd[0],76));close(fd[0]);int status=-1;
        assert(!WorkerProtocol::ResultMemoryStream(fd[1],status,[](const WorkerProtocol::StreamFile&,const WorkerProtocol::StreamReader&){}));assert(status==76);close(fd[1]);}
    std::cout<<"PASS chunk limits, scientific argv, typed classification, three wire formats and truncated status76 frame\n";
}
