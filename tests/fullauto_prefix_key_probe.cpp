#include "../src/FullAuto.cpp"
using namespace FullAuto;
int main(){
 ArgPP args;RunConfig config;Engine e(args,config);
 e.binaryHash="binary";e.eps=.04;e.workerChunk=256;e.seeds={11,23};
 Job j;j.count=1024;j.depth=18;j.phi=64;j.executionHash="execution";j.inputAngles={150,155,170};
 const auto key=e.CutoffPrefixValidationKey(j);
 for(int n:{1024,2048,32768,1048576}){j.count=n;Require(e.CutoffPrefixValidationKey(j)==key,"identical capped pilot was invalidated by production count");}
 j.count=512;Require(e.CutoffPrefixValidationKey(j)!=key,"different pilot count reused");j.count=1024;
 j.inputAngles[1]=156;Require(e.CutoffPrefixValidationKey(j)!=key,"different interior angle reused");j.inputAngles[1]=155;
 j.depth++;Require(e.CutoffPrefixValidationKey(j)!=key,"different depth reused");j.depth--;
 e.eps=.02;Require(e.CutoffPrefixValidationKey(j)!=key,"tighter accuracy reused");e.eps=.04;
 e.seeds.push_back(8803);Require(e.CutoffPrefixValidationKey(j)!=key,"changed independent seeds reused");e.seeds.pop_back();
 j.pairedDepth=16;Require(e.CutoffPrefixValidationKey(j)!=key,"changed paired depth reused");j.pairedDepth=0;
 j.phi=32;Require(e.CutoffPrefixValidationKey(j)!=key,"changed phi reused");j.phi=64;
 j.executionHash="other";Require(e.CutoffPrefixValidationKey(j)!=key,"changed physics or execution reused");j.executionHash="execution";
 e.binaryHash="other";Require(e.CutoffPrefixValidationKey(j)!=key,"changed binary reused");e.binaryHash="binary";
 e.weights={{{.2,.3}}};Require(e.CutoffPrefixValidationKey(j)!=key,"changed control weights reused");
 std::cout<<"PASS capped sample reuse; invalidation for count, interior angle, depth, accuracy, seeds and controls\n";
}
