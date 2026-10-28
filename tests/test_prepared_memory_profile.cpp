#include "PreparedMemoryProfile.h"
#include <cassert>
#include <iostream>
int main(){
    using namespace PreparedMemoryProfile;
    assert(Get("fresh").Start(512,1LL<<30)==16);
    Observe("droxtal-depth14",4096,64);
    assert(Get("droxtal-depth14").Start(512,1LL<<30)==16);
    Observe("droxtal-depth14",8192,256);
    const auto e=Get("droxtal-depth14");
    assert(e.Start(4096,1LL<<30)==256);
    assert(e.Start(64,1LL<<30)==64); // explicit chunk still wins
    assert(e.Start(512,8192*2*7)==7); // changed free RAM shrinks warm start
    assert(e.Start(512,1)==1);
    assert(Get("different-depth").Start(512,1LL<<30)==16);
    Observe("droxtal-depth14",1024,512);
    assert(Get("droxtal-depth14").safeBytes==8192); // never forget heavy beams
    for(int i=0;i<40;++i)Observe("bounded-"+std::to_string(i),16384,128);
    assert(Cache().values.size()<=32);
    assert(Get("droxtal-depth14").observations==0);
    std::cout << "Prepared memory profile: PASS\n";
}
