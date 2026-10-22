#include "../src/ValidatedGraphCache.h"
#include <cassert>
#include <iostream>
#include <set>
int main(){
    const size_t budget=2*1024*1024;
    ValidatedGraphCache compact(budget);BoundedCache<std::string> legacy(budget);
    const std::string directory="/home/dima/campaigns/list1_fullauto_2pct_20261014_8ae9e85/dima_gpu0/results/prism_L446.7_D147.1_theta170_180/refine_theta170_176_eps0p05/jobs/";
    auto marker=std::make_shared<std::string>("0123456789abcdef 1 abcdef0123456789 none 123456789abcdef\n");
    for(int i=0;i<4000;++i){const auto path=directory+"p"+std::to_string(i);
        compact.Put(path,marker,sizeof(std::string)+marker->capacity());legacy.Put(path,marker,sizeof(std::string)+marker->capacity());}
    size_t oldHits=0,newHits=0;
    for(int i=0;i<4000;++i){const auto path=directory+"p"+std::to_string(i);oldHits+=bool(legacy.Peek(path));newHits+=bool(compact.Peek(path));}
    assert(newHits==4000&&newHits>oldHits&&compact.Bytes()<=budget);
    // Equal job names in different parent outputs must never alias.
    const auto other=std::make_shared<std::string>("other receipt");compact.Put("/parent/jobs/p17",other,64);
    assert(*compact.Peek(directory+"p17")==*marker&&*compact.Peek("/parent/jobs/p17")==*other);
    std::set<std::string> visited;compact.LeastRecentFirst([&](const std::string& path,std::shared_ptr<const std::string> value){
        assert(visited.insert(path).second);assert(value==compact.Peek(path));});assert(visited.size()==4001);
    // Limit includes directory interning and behaves correctly after shrinking.
    compact.SetLimit(1);assert(compact.Bytes()==0);compact.Put("/new/entry",marker,128);assert(compact.Bytes()==0);
    compact.SetLimit(8192);for(int i=0;i<1000;++i){compact.Put("/directory"+std::to_string(i)+"/x",marker,128);assert(compact.Bytes()<=8192);}
    compact.Clear();assert(compact.Bytes()==0&&!compact.Peek("/parent/jobs/p17"));
    std::cout<<"PASS graph cache: identical paths, separate parents, bounded directory pool, shrink and clear; retained "<<newHits<<" versus "<<oldHits<<" with the same memory budget\n";
}
