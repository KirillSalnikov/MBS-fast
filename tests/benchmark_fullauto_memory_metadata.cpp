#include "../src/FullAutoMemoryStore.h"
#include <chrono>
#include <iostream>
using namespace FullAutoMemoryStore;
double now(){return std::chrono::duration<double>(std::chrono::steady_clock::now().time_since_epoch()).count();}
int main(int argc,char** argv){try{
 Require(argc==2,"temporary directory required");
 for(bool early:{false,true}){
  const std::string root=std::string(argv[1])+(early?"/early":"/old");const size_t limit=32*1024*1024;
  Store s(root,limit,argv[1],early?limit/4:0);double start=now();
  for(int i=0;i<14000;++i)s.Write(root+"/jobs/"+std::string(80,'a')+std::to_string(i)+"/metadata", "unchanged record");
  std::string value(32768,'x');for(int i=0;i<640;++i)s.Write(root+"/hot/"+std::to_string(i),value);
  auto initial=s.Stats();const double ready=now();
  for(int repeat=0;repeat<12;++repeat)for(int i=0;i<640;++i)Require(s.Read(root+"/hot/"+std::to_string(i))==value,"payload changed");
  auto final=s.Stats();const double read=now();
  s.Snapshot(std::string(argv[1])+(early?"/early.bin":"/old.bin"),"test");auto snap=s.Stats();
  Store restored(root,limit,argv[1],early?limit/4:0);Require(restored.Restore(std::string(argv[1])+(early?"/early.bin":"/old.bin"))=="test","restore failed");
  Require(restored.Read(root+"/hot/639")==value,"restore payload changed");
  Require(final.residentBytes+final.metadataBytes<=limit,"RAM budget exceeded");
  Require(snap.spillReadBytes-final.spillReadBytes==snap.snapshotReadBytes-final.snapshotReadBytes,"checkpoint read counter mismatch");
  std::cout<<"{\"early\":"<<early<<",\"setup_seconds\":"<<ready-start<<",\"read_seconds\":"<<read-ready<<",\"read_bytes\":"<<final.spillReadBytes-initial.spillReadBytes<<",\"metadata_bytes\":"<<final.metadataBytes<<",\"resident_bytes\":"<<final.residentBytes<<",\"catalog\":"<<final.catalogActive<<",\"snapshot_read_bytes\":"<<snap.snapshotReadBytes<<"}\n";
 }
}catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}}
