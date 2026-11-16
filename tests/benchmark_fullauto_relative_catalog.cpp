#include "FullAutoMemoryStore.h"
#include <cassert>
#include <iostream>
#include <iomanip>
#include <chrono>
#include <sstream>
#include <unistd.h>
static std::string Name(int i){std::ostringstream s;s<<"/jobs/"<<std::hex<<std::setw(16)<<std::setfill('0')<<(uint64_t(i)*11400714819323198485ULL)<<"/result.dat";return s.str();}
int main(int argc,char** argv){
 assert(argc==3||argc==4);const int count=argc==4?std::stoi(argv[3]):20000;assert(count>1&&count<=1000000);const std::string root="/home/worker/campaigns/bullet_fullauto4_20261109/gpu1_gpu0_parallel_20261112/results/coherent_followup_column_L446.7_D147.1_w0.532_theta150_170_restart_20261116";
 FullAutoMemoryStore::Store store(root,1024*1024,"/tmp",65536);
 const auto start=std::chrono::steady_clock::now();
 if(std::string(argv[1])=="write"){
  for(int i=0;i<count;++i)store.Write(root+Name(i),"payload:"+std::to_string(i));
  assert(store.Stats().catalogActive);
  store.Snapshot(argv[2],"relative catalog compatibility");
 }else assert(store.Restore(argv[2])=="relative catalog compatibility");
 for(int i=0;i<count;++i)assert(store.Read(root+Name(i))=="payload:"+std::to_string(i));
 assert(store.List(root+"/jobs",false).size()==size_t(count));
 assert(store.List(root+"/jobs",true).size()==size_t(count)*2);
 store.Rename(root+Name(0),root+"/renamed");assert(store.Read(root+"/renamed")=="payload:0");
 store.Remove(root+Name(1));assert(!store.Exists(root+Name(1)));
 auto st=store.Stats();std::cout<<argv[1]<<" catalog_bytes="<<st.catalogBytes<<" seconds="<<std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count()<<" verified="<<count<<"\n";
}
