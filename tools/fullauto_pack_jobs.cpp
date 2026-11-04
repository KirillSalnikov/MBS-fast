#include "../src/PackedJobStore.h"
#include <fstream>
#include <iostream>
#include <sys/file.h>
// Build: c++ -std=c++11 -O2 -pthread tools/fullauto_pack_jobs.cpp -o fullauto_pack_jobs
// Run only after stopping the queue runner as well as its native controller.
int main(int argc,char** argv){try{
    using namespace PackedJobStore;
    Check(argc==2||argc==3,"usage: fullauto_pack_jobs /absolute/output [maximum_jobs | --unpack]");
    const std::string root=argv[1];Check(!root.empty()&&root[0]=='/'&&root.back()!='/',"use an absolute output without trailing slash");
    struct stat info{};Check(stat((root+"/configuration.txt").c_str(),&info)==0,"not a fullauto output (configuration.txt absent)");
    Descriptor lock(open((root+"/controller.lock").c_str(),O_WRONLY|O_CREAT|O_CLOEXEC,0644));Check(lock.fd>=0&&flock(lock.fd,LOCK_EX|LOCK_NB)==0,"controller is active; stop the queue before migration");
    if(argc==3&&std::string(argv[2])=="--unpack"){size_t jobs=0,files=0;for(const auto& name:Jobs(root)){files+=Extract(root+"/jobs/"+name);if(++jobs%256==0)std::cout<<"{\"unpacked_jobs\":"<<jobs<<",\"created_files\":"<<files<<",\"finished\":false}"<<std::endl;}
        std::cout<<"{\"unpacked_jobs\":"<<jobs<<",\"created_files\":"<<files<<",\"finished\":true}"<<std::endl;return 0;}
    size_t maximum=0;if(argc==3){size_t end=0;maximum=std::stoull(argv[2],&end);Check(end==std::strlen(argv[2])&&maximum>0,"invalid maximum_jobs");}
    // Persist the reader requirement before removing any ordinary directory.
    const auto policy=root+"/job_pack_format.txt",temporary=policy+".tmp";
    if(stat(policy.c_str(),&info)==0){Descriptor existing(open(policy.c_str(),O_RDONLY));Check(existing.fd>=0&&Bytes(existing.fd,0,size_t(info.st_size))=="sharded_jobs_v1\n","unsupported existing pack policy");}
    Descriptor marker(open(temporary.c_str(),O_WRONLY|O_CREAT|O_TRUNC,0644));Check(marker.fd>=0,"cannot publish reader requirement");Write(marker.fd,"sharded_jobs_v1\n");Flush(marker.fd);Check(rename(temporary.c_str(),policy.c_str())==0,"cannot publish reader requirement");FlushDirectory(root);
    DIR* directory=opendir((root+"/jobs").c_str());Check(directory,"cannot list jobs");size_t jobs=0,inodes=0,skipped=0,oversize=0;
    while(auto* entry=readdir(directory)){
        const std::string name=entry->d_name;if(!JobName(name))continue;
        const auto folder=root+"/jobs/"+name;
        if(stat((folder+"/complete.txt").c_str(),&info)!=0&&!Exists(folder+"/complete.txt")){++skipped;continue;}
        try{inodes+=Pack(folder);++jobs;}catch(const Oversize&){++oversize;continue;}
        if(jobs%256==0)std::cout<<"{\"packed_jobs\":"<<jobs<<",\"released_inodes\":"<<inodes<<",\"incomplete_skipped\":"<<skipped<<",\"finished\":false}"<<std::endl;
        if(maximum&&jobs>=maximum)break;
    }
    closedir(directory);std::cout<<"{\"packed_jobs\":"<<jobs<<",\"released_inodes\":"<<inodes<<",\"incomplete_skipped\":"<<skipped<<",\"oversize_retained\":"<<oversize<<",\"finished\":true}"<<std::endl;return 0;
}catch(const std::exception& error){std::cerr<<error.what()<<'\n';return 1;}}
