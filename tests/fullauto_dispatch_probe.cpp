#include "../src/FullAuto.cpp"
using namespace FullAuto;
int main(int argc,char** argv){try{
    Require(argc==4,"OUTPUT BINARY SOURCE_JOB");ArgPP args;RegisterCliOptions(args);RunConfig config;config.useGpu=true;Engine e(args,config);
    e.root=argv[1];Directory(e.root);Directory(e.root+"/jobs");e.exe=argv[2];e.binaryHash=FileHash(e.exe);e.signature="dispatch-regression-v1";e.gpus={"0"};e.hostBudgetMb=512;e.workerPolicy="persistent";
    Job j;std::istringstream index(Read(std::string(argv[3])+"/sampling.index"));std::string sig;
    Require(bool(index>>sig>>j.identity>>j.count>>j.offset>>j.seed>>j.depth>>j.phi>>j.block>>j.pipeline>>j.mirrored>>j.paired>>j.assembled>>j.hasControls>>j.pairedDepth>>j.warp>>j.weightHash>>j.executionHash),"bad source metadata");
    double angle;while(index>>angle)j.inputAngles.push_back(angle);Require(index.eof()&&j.inputAngles.size()>=2,"bad source angles");
    std::istringstream commandText(Read(std::string(argv[3])+"/command.args"));std::string line;std::getline(commandText,line);std::vector<std::string> command;while(std::getline(commandText,line))command.push_back(line);
    e.controls=j.hasControls;e.phi=j.phi;e.mirror=j.mirrored;e.kernel.block=j.block;e.kernel.pipeline=j.pipeline;e.kernel.warp=j.warp;
    for(size_t i=0;i+1<command.size();++i){if(command[i]=="--threads")e.workerThreads=std::stoi(command[i+1]);if(command[i]=="--orientation-chunk")e.workerChunk=std::stoi(command[i+1]);if(command[i]=="--po-prepare-threads")j.preparationThreads=std::stoi(command[i+1]);}
    e.shortPacketVerifying=false;e.VerifyShortPacketChunk(j,command);
    std::cout<<"PASS actual native recipe packet validation; inspect receipt for performance acceptance\n";return 0;
}catch(const std::exception& ex){std::cerr<<ex.what()<<'\n';return 1;}}
