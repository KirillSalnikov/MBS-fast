#include "../src/WorkerProtocol.h"
#include "../src/WorkerMemoryStage.h"
#include <sys/wait.h>
#include <signal.h>
#include <cassert>
#include <cmath>
#include <fstream>
#include <iostream>
#include <sstream>
#include <stdexcept>
static void Require(bool okay,const std::string& message){if(!okay)throw std::runtime_error(message);}
static std::string Read(const std::string& path){std::ifstream file(path,std::ios::binary);Require(bool(file),"cannot read "+path);return std::string(std::istreambuf_iterator<char>(file),std::istreambuf_iterator<char>());}
static void Write(const std::string& path,const std::string& value){std::ofstream file(path,std::ios::binary);Require(bool(file),"cannot write "+path);file.write(value.data(),value.size());file.close();Require(bool(file),"write failed "+path);}
struct Worker {
    int fd=-1;pid_t pid=-1;
    explicit Worker(const std::string& executable){int pair[2];Require(socketpair(AF_UNIX,SOCK_STREAM,0,pair)==0,"socketpair failed");pid=fork();Require(pid>=0,"fork failed");
        if(pid==0){close(pair[0]);const auto socket=std::to_string(pair[1]);execl(executable.c_str(),executable.c_str(),"--fullauto-worker-server",socket.c_str(),static_cast<char*>(nullptr));_exit(127);}
        close(pair[1]);fd=pair[0];timeval timeout{60,0};setsockopt(fd,SOL_SOCKET,SO_RCVTIMEO,&timeout,sizeof(timeout));}
    ~Worker(){if(fd>=0)close(fd);if(pid>0){kill(pid,SIGTERM);int status;while(waitpid(pid,&status,0)<0&&errno==EINTR){}WorkerMemoryStage::CleanupForPid(pid);}}
    void Finish(){Require(WorkerProtocol::Send(fd,{}),"shutdown send failed");close(fd);fd=-1;int status=0;Require(waitpid(pid,&status,0)==pid&&WIFEXITED(status)&&WEXITSTATUS(status)==0,"worker exited abnormally");pid=-1;}
};
static std::vector<std::vector<double>> Matrix(const std::string& text){std::istringstream input(text);std::string line;std::vector<std::vector<double>> rows;
    while(std::getline(input,line)){std::istringstream row(line);double number=0;if(!(row>>number))continue;std::vector<double> values{number};while(row>>number){Require(std::isfinite(number),"nonfinite matrix value");values.push_back(number);}Require(values.size()==18,"unexpected Mueller row width");rows.push_back(values);}Require(rows.size()==3,"unexpected theta count");return rows;}
static std::pair<unsigned long long,unsigned long long> Reuse(const std::string& log){const std::string prefix="Prepared trace reuse: hits=";const auto first=log.rfind(prefix);Require(first!=std::string::npos,"missing prepared-cache counters\n"+log);const auto hits=std::stoull(log.substr(first+prefix.size()));const auto second=log.find(" misses=",first);Require(second!=std::string::npos,"missing miss count");return {hits,std::stoull(log.substr(second+8))};}
int main(int argc,char** argv){try{
    Require(argc==4,"usage: probe BINARY CUBE_PARTICLE NEW_SCRATCH_DIRECTORY");alarm(120);const std::string binary=argv[1],geometry=Read(argv[2]),disk=argv[3],logical="/mbs-prepared-identity-fixture";
    Require(mkdir(disk.c_str(),0700)==0,"scratch directory must be new");const std::string grid="150\n155\n160\n",visibility="6\n1 1\n1 1\n1 1\n1 1\n1 1\n1 1\n",changedVisibility="6\n0 1\n1 1\n1 1\n1 1\n1 1\n1 1\n";
    Write(disk+"/particle.dat",geometry);Write(disk+"/particle.dat.visibility",visibility);Write(disk+"/theta.txt",grid);Worker worker(binary);
    auto command=[&](const std::string& root,const std::string& label,bool incoherent){std::vector<std::string> words{root+"/jobs/"+label+"/stdout.log","64","auto","--method","po","--backend","cpu","--particle-file",root+"/particle.dat","--resize-dmax-um","5","--wavelength-um","0.532","--refractive-index","1.3116","0","--max-reflections","4","--sobol-seed","32","11","--haar-alpha","--theta-grid-file",root+"/theta.txt","--phi-points","4","--threads","1","--trace-cache-mib","16","--cutoff-profile","off","--output",root+"/jobs/"+label+"/result","--close"};if(incoherent)words.push_back("--incoherent");return words;};
    for(bool incoherent:{false,true}){
        const std::string mode=incoherent?"incoherent":"coherent";std::vector<std::vector<double>> reference;
        for(int index=0;index<5;++index){const auto label=mode+std::to_string(index);const auto data=index==2?"# changed byte identity, same geometry\n"+geometry:geometry;const auto side=index==3?changedVisibility:visibility;
            Require(WorkerProtocol::SendMemoryRequest(worker.fd,logical,command(logical,label,incoherent),{{"particle.dat",data},{"particle.dat.visibility",side},{"theta.txt",grid}}),"memory request failed");
            std::vector<WorkerProtocol::File> files;int status=-1;Require(WorkerProtocol::ResultMemory(worker.fd,status,files),"memory reply failed");std::string log,matrix;
            for(const auto& file:files){if(file.path=="jobs/"+label+"/stdout.log")log=file.bytes;if(file.path=="jobs/"+label+"/result/result.dat")matrix=file.bytes;}
            Write(disk+"/"+label+".log",log);Require(status==0,"native request failed: "+std::to_string(status)+"\n"+log);Require(!matrix.empty(),"matrix missing in bundle");Write(disk+"/"+label+".dat",matrix);
            const auto counts=Reuse(log);const bool hit=index==1||index==4;Require(hit?(counts.first==32&&counts.second==0):(counts.first==0&&counts.second==32),"unexpected cache counts for "+label+": "+std::to_string(counts.first)+"/"+std::to_string(counts.second));
            const auto rows=Matrix(matrix);if(index==0)reference=rows;else if(index!=3)Require(rows==reference,"same-physics memory requests changed all16 matrix elements");
            std::cout<<label<<" hits="<<counts.first<<" misses="<<counts.second<<" rows="<<rows.size()<<'\n';
        }
        // Compare to the ordinary file transport with the same actual particle
        // and sidecar; no path rewriting is involved in this request.
        const auto label=mode+"legacy";Require(mkdir((disk+"/jobs").c_str(),0700)==0||errno==EEXIST,"legacy jobs mkdir failed");Require(mkdir((disk+"/jobs/"+label).c_str(),0700)==0,"legacy job mkdir failed");
        Require(WorkerProtocol::Send(worker.fd,command(disk,label,incoherent)),"legacy request failed");int status=-1;Require(WorkerProtocol::Result(worker.fd,status)&&status==0,"legacy request failed");
        Require(Matrix(Read(disk+"/jobs/"+label+"/result/result.dat"))==reference,"memory transport changed native physics relative to legacy file transport");
        std::cout<<mode<<" legacy comparison: all16 exact PASS\n";
    }
    worker.Finish();std::cout<<"native prepared identity: stable cache hits across tmpfs requests, geometry/sidecar invalidation and all16 legacy equivalence PASS\n";return 0;
}catch(const std::exception& error){std::cerr<<error.what()<<'\n';return 1;}}
