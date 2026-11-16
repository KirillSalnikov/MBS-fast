#include "../src/FullAutoMemoryStore.h"
#include <cassert>
#include <iostream>
using namespace FullAutoMemoryStore;
int main(){char folder[]="/tmp/mbs-child-files-XXXXXX";assert(mkdtemp(folder));for(bool spill:{false,true}){
 const std::string root=std::string(folder)+(spill?"/sql":"/map");Store s(root,8*1024*1024,folder,spill?1024:0);
 for(int i=0;i<200;++i)s.Write(root+"/receipt"+std::to_string(i),"data");
 s.Write(root+"/nested/hidden","hidden");s.Directory(root+"/empty");
 auto verify=[&]{std::vector<std::pair<std::string,uint64_t>> expected;for(const auto& p:s.List(root,false))if(!s.IsDirectory(p))expected.emplace_back(p,s.Generation(p));assert(s.ChildFiles(root)==expected);assert(s.ChildFiles(root+"/empty").empty());assert(s.ChildFiles(root+"/nested").size()==1);};
 verify();assert(s.Stats().catalogActive==spill);s.Write(root+"/receipt3","updated");s.Remove(root+"/receipt4");s.Rename(root+"/receipt5",root+"/renamed");verify();
 const std::string cp=std::string(folder)+"/snapshot";s.Snapshot(cp,"unchanged");Store restored(root,8*1024*1024,folder,spill?1024:0);assert(restored.Restore(cp)=="unchanged");assert(restored.ChildFiles(root)==s.ChildFiles(root));unlink(cp.c_str());
 }assert(rmdir(folder)==0);std::cout<<"map/catalog publication listings, mutations and checkpoint roundtrip PASS\n";}
