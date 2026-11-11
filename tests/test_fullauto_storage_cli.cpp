#include "ArgPP.h"
#include "CliOptions.h"
#include "RunConfig.h"
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

namespace {
void Check(const std::vector<std::string>& extra,const std::string& failure="",
           bool fullauto=true,const std::string& method="po") {
    std::vector<std::string> words={"mbs_po_test","--method",method,"--backend","cpu",
        "--particle","1",".2",".1","--ri","1.31","0","--wavelength-um","1"};
    if(fullauto)words.insert(words.end(),{"--fullauto",".02"});
    else words.insert(words.end(),{"--euler-grid","4","4"});
    words.insert(words.end(),extra.begin(),extra.end());
    std::vector<const char*> argv;for(const auto& word:words)argv.push_back(word.c_str());
    std::string message;
    try {ArgPP args;RegisterCliOptions(args);args.Parse(int(argv.size()),argv.data());
        (void)RunConfig::FromCommandLine(args,false,false);}
    catch(const std::exception& error){message=error.what();}
    if(failure.empty()&&!message.empty())throw std::runtime_error(message);
    if(!failure.empty()&&(message.find(failure)==std::string::npos||message.find("Fix:")==std::string::npos))
        throw std::runtime_error("Expected actionable failure containing "+failure+"; got "+message);
}
}
int main(){
    try {
        for(const std::string method:{"po","go"}){
            Check({},"",true,method);
            Check({"--fullauto-storage","memory","--fullauto-memory-mib","512"},"",true,method);
            Check({"--checkpoint"},"",true,method);
            Check({"--fullauto-layout","compact"},"",true,method);
            Check({"--checkpoint","--fullauto-checkpoint-seconds","0"},"",true,method);
            Check({"--fullauto-checkpoint-seconds","600"},"",true,method);
            Check({"--fullauto-storage","legacy","--checkpoint"},"",true,method);
        }
        Check({"--fullauto-layout","bad"},"expects compact or legacy");
        Check({"--fullauto-layout","compact"},"requires --fullauto",false);
        Check({"--fullauto-storage","memry"},"expects memory or legacy");
        Check({"--fullauto-memory-mib","0"},"positive integer");
        Check({"--fullauto-memory-mib","-1"},"positive integer");
        Check({"--fullauto-memory-mib","1.5"},"integer");
        Check({"--fullauto-checkpoint-seconds","-1"},"must not be negative");
        Check({"--fullauto-checkpoint-seconds","0.5"},"integer");
        Check({"--fullauto-storage","legacy","--fullauto-memory-mib","512"},"unused");
        Check({"--fullauto-storage","legacy","--fullauto-checkpoint-seconds","0"},"unused");
        Check({"--fullauto-storage","memory"},"requires --fullauto",false);
        Check({"--fullauto-checkpoint-seconds","0"},"requires --fullauto",false);
        Check({"--fullauto-memory-mib","512"},"requires --fullauto",false);
        Check({"--checkpoint"},"",false);
        Check({"--checkpoint"},"not implemented",false,"go");
        for(const std::string key:{"fullauto_storage","fullauto_checkpoint_seconds","fullauto_memory_mib","checkpoint"})
            if(std::string(FullAutoOptionPolicy(key))!="controller")throw std::runtime_error("not a controller setting: "+key);
        std::cout<<"PASS: fullauto memory/legacy/checkpoint CLI validation and standalone compatibility\n";
    }catch(const std::exception& error){std::cerr<<"FAIL: "<<error.what()<<'\n';return 1;}
}
