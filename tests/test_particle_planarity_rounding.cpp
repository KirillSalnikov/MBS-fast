#include "particle/Particle.h"
#include <iostream>
#include <exception>
int main(int argc,char**argv){if(argc!=2)return 2;try{Particle p;p.SetFromFile(argv[1]);std::cout<<p.nFacets<<" "<<p.MaximalDimention()<<"\n";return 0;}catch(const std::exception& e){std::cerr<<e.what()<<"\n";return 1;}}
