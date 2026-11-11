#include "../src/FullAutoThetaPilotPlan.h"
#include <cassert>
#include <iostream>
int main(){using FullAutoThetaPilotPlan::Count;
 assert(Count(256,12839,64,false)==32);
 assert(Count(256,12839,64,true)==256);
 assert(Count(16,12839,64,false)==16);
 assert(Count(256,65,4,false)==256);
 assert(Count(0,12839,64,false)==0);
 assert(Count(256,SIZE_MAX,2147483647,false)==32);
 for(size_t n:{1,65,1025,12839})for(int p:{1,4,64,1024})for(int requested:{0,16,32,64,128,256}){
  const int count=Count(requested,n,p,false);assert(count<=requested&&count>=std::min(requested,32));
  assert(!count||(count&(count-1))==0);assert(Count(requested,n,p,true)==requested);
 }
 std::cout<<"PASS implicit theta pilot work bound, explicit priority, small grids, overflow\n";
}
