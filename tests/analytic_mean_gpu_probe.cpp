#include "AnalyticFacetAverage.h"
#include <cmath>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <chrono>
using AnalyticBackscatter::Vec;
using AnalyticBackscatter::Face;
int main(int argc,char**argv){try{
    if(argc!=8)return 2;
    const bool gpu=std::string(argv[1])=="cuda";
    const double length=std::atof(argv[2]),diameter=std::atof(argv[3]),lo=std::atof(argv[4]),hi=std::atof(argv[5]);
    const int count=std::atoi(argv[6]),threads=std::atoi(argv[7]);const double pi=std::acos(-1.);
    std::vector<Face> faces;Face top,bottom;top.normal=Vec(0,0,1);bottom.normal=Vec(0,0,-1);
    for(int i=0;i<6;++i){const double a=i*pi/3;top.vertices.push_back(Vec(diameter/2*std::cos(a),diameter/2*std::sin(a),length/2));}
    for(int i=5;i>=0;--i){Vec v=top.vertices[i];v.z=-v.z;bottom.vertices.push_back(v);}
    faces.push_back(top);faces.push_back(bottom);
    for(int i=0;i<6;++i){const int j=(i+1)%6;Face side;const double a=(i+.5)*pi/3;side.normal=Vec(std::cos(a),std::sin(a),0);
        const Vec p=top.vertices[i],q=top.vertices[j];side.vertices={q,p,Vec(p.x,p.y,-p.z),Vec(q.x,q.y,-q.z)};faces.push_back(side);}
    std::vector<double>theta;for(int i=0;i<count;++i)theta.push_back((lo+(hi-lo)*i/(count-1))*pi/180);
    const auto start=std::chrono::steady_clock::now();
    AnalyticFacetAverage::Control c(faces,{1.3116,0},.532,theta,"circular","","",threads,gpu);
    std::cerr<<"BENCHMARK_SECONDS="<<std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count()<<'\n';
    std::cout<<std::setprecision(17);
    for(int i=0;i<count;++i){auto m=c.Mean(i);std::cout<<theta[i]<<' '<<m.reflection<<' '<<m.shadow<<' '<<c.RefinementError()<<'\n';}
}catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}}
