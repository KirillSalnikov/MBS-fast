#include "math/Sobol.h"
#include <cassert>
#include <iostream>
int main(){
    for(unsigned seed:{0u,11u,32768u})for(bool scramble:{false,true}){
        Sobol3D sequential(seed,scramble);
        for(unsigned i=0;i<8192;++i){
            double a,b,c,x,y,z;sequential.next(a,b,c);
            if(i<16 || i%127==0 || i==8191){
                Sobol3D block(seed,scramble);block.seek(i);block.next(x,y,z);
                assert(a==x && b==y && c==z);
            }
        }
        Sobol3D large(seed,scramble);large.seek(134217727u);
        double a,b,c,x,y,z;large.next(a,b,c);large.seek(134217728u);large.next(x,y,z);
        Sobol3D repeat(seed,scramble);repeat.seek(134217727u);repeat.next(x,y,z);
        assert(a==x && b==y && c==z);
    }
    std::cout<<"Sobol block offsets preserve exact prefix samples\n";
}
