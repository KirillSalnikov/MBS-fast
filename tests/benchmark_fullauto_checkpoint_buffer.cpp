#include "../src/FullAutoMemoryStore.h"
#include <chrono>
#include <iostream>
using namespace FullAutoMemoryStore;
int main(int argc,char**argv){if(argc!=2)return 2;std::string path=argv[1];Store s("/checkpoint-buffer-benchmark",64*1024*1024);for(int i=0;i<30000;++i)s.Write("/checkpoint-buffer-benchmark/job"+std::to_string(i),std::string(128,char(i%251)));
 auto start=std::chrono::steady_clock::now();s.Snapshot(path,"fixture");auto written=std::chrono::steady_clock::now();Store r("/checkpoint-buffer-benchmark",64*1024*1024);Require(r.Restore(path)=="fixture","metadata");auto read=std::chrono::steady_clock::now();for(int i=0;i<30000;++i)Require(r.Read("/checkpoint-buffer-benchmark/job"+std::to_string(i))==std::string(128,char(i%251)),"payload");std::cout<<"{\"write_seconds\":"<<std::chrono::duration<double>(written-start).count()<<",\"restore_seconds\":"<<std::chrono::duration<double>(read-written).count()<<"}\n";}
