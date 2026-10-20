#include "../src/FullAuto.cpp"
using namespace FullAuto;
int main(int argc,char**argv){try{
 Require(argc==4,"OUTPUT BINARY NATIVE_COMMAND_FILE");ArgPP args;RunConfig config;config.useGpu=true;Engine e(args,config);
 e.root=argv[1];Directory(e.root);Directory(e.root+"/jobs");e.exe=argv[2];e.binaryHash=FileHash(e.exe);e.signature="short-packet-integration-v4";
 std::istringstream input(Read(argv[3]));std::string line;std::getline(input,line);std::vector<std::string> command;while(std::getline(input,line))command.push_back(line);
 auto value=[&](const std::string& key){const auto i=std::find(command.begin(),command.end(),key);Require(i!=command.end()&&std::next(i)!=command.end(),"missing native argument "+key);return *std::next(i);};
 e.base=command;e.phi=std::stoi(value("--phi-points"));e.workerThreads=std::stoi(value("--threads"));e.target=std::stoi(value("--max-reflections"));
 e.controls=std::find(command.begin(),command.end(),"--analytic-facet-average")!=command.end();e.sharedMeans=true;e.gpus={"0"};e.hostBudgetMb=8192;
 e.workerPolicy="persistent";e.resourcePolicy="auto";e.kernel.pipeline=true;e.eps=.02;e.ratios=true;
 std::istringstream grid(Read(value("--theta-grid-file")));double theta;while(grid>>theta)e.theta.push_back(theta);Require(e.theta.size()>1,"invalid grid");
 Job job;e.SetSampling(job,4096,71,e.target,e.theta,false);job.weightHash=Hash("");
 if(std::find(command.begin(),command.end(),"--paired-depth")!=command.end())job.pairedDepth=std::stoi(value("--paired-depth"));
#ifdef MBS_GPU_FP32
 Require(e.ShortPacketChunk(job)==512,"float cap not selected");
#else
 Require(e.ShortPacketChunk(job)==(job.inputAngles.size()<=8?512:1024),"double packet cap incorrect");
#endif
 const bool accepted=e.VerifyShortPacketChunk(job,command);const auto history=Read(e.root+"/fullauto_history.jsonl");
 Require(e.VerifyShortPacketChunk(job,command)==accepted,"cached decision changed");Require(Read(e.root+"/fullauto_history.jsonl")==history,"cached decision launched work");
 Require(history.find("short_packet_validation")!=std::string::npos,"validation not logged");
 std::cout<<"PASS native chunk validation, all16 and ratios, two independent pilots, cached decision; accepted="<<accepted<<" cap="<<e.ShortPacketCandidate(job.inputAngles.size())<<'\n';return 0;
}catch(const std::exception& error){std::cerr<<error.what()<<'\n';return 1;}}
