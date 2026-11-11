#include "../src/FullAutoValidationRootPlan.h"
#include <cassert>
#include <iostream>
#include <algorithm>
int main(){
 const std::vector<std::string> roots={"/r/a","/r/b","/r/c","/r/d"};
 std::vector<std::string> text={"","/r/a\thash-a\n","/r/b\thash-b\n",""};size_t reads=0;
 auto reader=[&](size_t i,size_t budget,std::string& out){++reads;if(text[i].size()>budget)return false;out=text[i];return true;};
 auto p=PlanValidationRoots(roots,reader);assert(p.Ordered()&&p.first==std::vector<size_t>({2,3})&&p.second==std::vector<size_t>({0,1})&&p.edges==2&&reads==4);
 // Duplicate edges must not inflate indegrees; unknown roots never become jobs.
 text[2]+="/r/b\thash-b\n/unknown\thash-other\n";p=PlanValidationRoots(roots,reader);assert(p.Ordered()&&p.edges==2);
 text[0]="/r/c\thash-c\n";p=PlanValidationRoots(roots,reader);assert(!p.Ordered()&&p.first==std::vector<size_t>({0,1,2,3}));
 text[0]="/r/a\thash-a\n";assert(!PlanValidationRoots(roots,reader).Ordered());
 text[0]="/r/a\t\n";assert(!PlanValidationRoots(roots,reader).Ordered());
 text[0]="bad row\n";assert(!PlanValidationRoots(roots,reader).Ordered());
 text[0]="";assert(!PlanValidationRoots(roots,reader,512,2).Ordered());
 p=PlanValidationRoots(roots,reader,512,8*1024*1024,1);
 assert(p.Ordered()&&p.incomplete&&p.edges==1&&p.first==std::vector<size_t>({2,3,1})&&p.second==std::vector<size_t>({0}));
 // Keep fully discovered links when the next file cannot fit. Unread roots
 // are scheduled first, but every root still appears once for real validation.
 p=PlanValidationRoots(roots,reader,512,text[1].size());
 assert(p.Ordered()&&p.incomplete&&p.manifestBytes==text[1].size()&&p.edges==1);
 assert(p.first==std::vector<size_t>({2,3,1})&&p.second==std::vector<size_t>({0}));
 auto unavailable=[&](size_t i,size_t budget,std::string& out){return i==2?false:reader(i,budget,out);};
 p=PlanValidationRoots(roots,unavailable);assert(p.incomplete&&p.Ordered()&&p.second==std::vector<size_t>({0}));
 auto oversized=[&](size_t i,size_t budget,std::string& out){if(i==2){out.assign(budget+1,'x');return true;}return reader(i,budget,out);};
 p=PlanValidationRoots(roots,oversized,512,64);assert(p.incomplete&&p.Ordered()&&p.manifestBytes<=64);
 // Complete and truncated DAGs preserve the complete requested-root set at
 // all bounds. No source mentioned only in a manifest becomes an extra job.
 for(size_t budget=0;budget<100;++budget)for(size_t edges=0;edges<5;++edges){
  p=PlanValidationRoots(roots,reader,512,budget,edges);auto ids=p.first;ids.insert(ids.end(),p.second.begin(),p.second.end());std::sort(ids.begin(),ids.end());
  assert(ids==std::vector<size_t>({0,1,2,3})&&p.manifestBytes<=budget&&p.edges<=edges);
 }
 reads=0;assert(!PlanValidationRoots(roots,reader,3).Ordered()&&reads==0);
 auto duplicate=roots;duplicate.push_back(roots[0]);assert(!PlanValidationRoots(duplicate,reader).Ordered());
 text.assign(4,"");p=PlanValidationRoots(roots,reader);assert(!p.Ordered()&&p.first.size()==4);
 std::cout<<"PASS direct root planning: chains, branches, duplicates, cycles, invalid manifests, bounded partial DAGs/unavailable/oversized reads, byte/root/edge bounds; every root retained\n";
}
