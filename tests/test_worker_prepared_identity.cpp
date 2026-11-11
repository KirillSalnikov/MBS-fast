#include "../src/WorkerMemoryStage.h"
#include <cassert>
#include <sstream>
#include <iostream>
static std::string Identity(const std::string& path){std::ostringstream out;WorkerMemoryStage::AppendPreparedParticleIdentity(out,path);return out.str();}
static void Put(const std::string& path,const std::string& bytes){std::ofstream file(path,std::ios::binary);file.write(bytes.data(),bytes.size());file.close();assert(file);}
int main(){
    const std::string logical="/logical/output/particle.dat",geometry="geometry bytes\n",visibility="6\n1 1\n1 1\n1 1\n1 1\n1 1\n1 1\n";
    std::string first,second,physicalFirst,physicalSecond;
    auto recipe=[&](const std::string& root,const std::string& data,const std::string& side,std::string* physical){
        WorkerMemoryStage::Request stage(root,{{"particle.dat",data},{"particle.dat.visibility",side}});
        const auto actual=stage.Rewrite(root+"/particle.dat");if(physical)*physical=actual;WorkerMemoryStage::IdentityScope scope(&stage);
        assert(WorkerMemoryStage::PreparedInputPath(actual)==root+"/particle.dat");
        assert(WorkerMemoryStage::PreparedInputPath(stage.Root()+"-other/particle.dat")==stage.Root()+"-other/particle.dat");
        bool rejected=false;try{WorkerMemoryStage::PreparedInputPath(stage.Root()+"/not-an-input");}catch(const std::runtime_error&){rejected=true;}assert(rejected);
        const auto identity=Identity(actual);assert(stage.Collect().empty());return identity;
    };
    first=recipe("/logical/output",geometry,visibility,&physicalFirst);second=recipe("/logical/output",geometry,visibility,&physicalSecond);
    assert(physicalFirst!=physicalSecond&&first==second&&WorkerMemoryStage::IdentityRequest()==nullptr);
    assert(first!=recipe("/logical/output",geometry+"changed",visibility,nullptr));
    assert(first!=recipe("/logical/output",geometry,"6\n0 1\n1 1\n1 1\n1 1\n1 1\n1 1\n",nullptr));
    assert(first!=recipe("/different/output",geometry,visibility,nullptr));
    {WorkerMemoryStage::Request stage("/logical/output",{{"particle.dat",geometry}});WorkerMemoryStage::IdentityScope scope(&stage);assert(Identity(stage.Rewrite(logical))!=first);}
    // The actual input bytes remain authoritative even if modified within a
    // request; normal transport never changes them after declaring the input.
    {WorkerMemoryStage::Request stage("/logical/output",{{"particle.dat",geometry},{"particle.dat.visibility",visibility}});WorkerMemoryStage::IdentityScope scope(&stage);
        const auto actual=stage.Rewrite(logical);assert(Identity(actual)==first);Put(actual,geometry+"changed");assert(Identity(actual)!=first);}
    char temp[]="/tmp/mbs-prepared-identity-XXXXXX";const int fd=mkstemp(temp);assert(fd>=0);close(fd);const std::string legacy=temp;Put(legacy,geometry);Put(legacy+".visibility",visibility);
    std::ostringstream expected;expected<<legacy.size()<<':'<<legacy<<geometry.size()<<':'<<geometry<<visibility.size()<<':'<<visibility;
    assert(Identity(legacy)==expected.str());
    {WorkerMemoryStage::Request stage("/logical/output",{});WorkerMemoryStage::IdentityScope scope(&stage);assert(Identity(legacy)==expected.str());
        {WorkerMemoryStage::IdentityScope native(nullptr);assert(WorkerMemoryStage::IdentityRequest()==nullptr);}assert(WorkerMemoryStage::IdentityRequest()==&stage);}
    unlink(legacy.c_str());unlink((legacy+".visibility").c_str());assert(WorkerMemoryStage::IdentityRequest()==nullptr);
    std::cout<<"prepared particle identity: stable across staging requests, geometry/visibility changes invalidate, declared-input boundary, unchanged legacy format PASS\n";
}
