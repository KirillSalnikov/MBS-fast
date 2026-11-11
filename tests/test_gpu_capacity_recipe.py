#!/usr/bin/env python3
"""Compile the actual PreparedRecipe function and parser; no CUDA/optical run."""
from pathlib import Path
import subprocess
import sys

root = Path(__file__).resolve().parents[1]
out = Path(sys.argv[1]).resolve()
out.mkdir(parents=True, exist_ok=True)
main = (root / 'src/main.cpp').read_text()
start = main.index('static std::string PreparedRecipe(')
end = main.index('\nint RunCalculation(', start)
function = main[start:end]
probe = r'''
#include "ArgPP.h"
#include "CliOptions.h"
#include "WorkerMemoryStage.h"
#include "RuntimeInfo.h"
#include <set>
#include <cassert>
#include <iostream>
#include <algorithm>
// The recipe's environment suffix is independently tested by RuntimeInfo;
// hold it fixed here to isolate geometry/material/sampling key semantics.
RuntimeEnvironmentSnapshot QueryRuntimeEnvironmentSnapshot(){return {};}
std::string FormatRuntimeEnvironmentReport(const RuntimeEnvironmentSnapshot&){return "fixed environment";}
''' + function + r'''
static std::string Recipe(const std::vector<std::string>& words){
    ArgPP parser;RegisterCliOptions(parser);std::vector<const char*> argv;
    for(const auto& word:words)argv.push_back(word.c_str());
    parser.Parse(int(argv.size()),argv.data());
    return PreparedRecipe(parser,true);
}
static void Set(std::vector<std::string>& args,const std::string& flag,size_t index,const std::string& value){
    auto at=std::find(args.begin(),args.end(),flag);assert(at!=args.end()&&size_t(args.end()-at)>index);at[index]=value;
}
int main(int argc,char** argv){
    assert(argc==2);const std::string file=std::string(argv[1])+"/particle.dat";
    std::vector<std::string> base={"worker","--method","po","--particle","1","501.2","155.8",
        "--ri","1.3116","0","--wavelength-um","0.532","--max-reflections","22","--paired-depth","20",
        "--cutoff-profile","off","--trace-cutoff-jones","0.001","--po-trace-backend","resident",
        "--sobol-seed","512","17","--sobol-offset","0","--phi-points","32",
        "--theta-grid-file","first-grid","--threads","8","--orientation-chunk","64","--output","first"};
    const auto reference=Recipe(base);
    for(const auto& flag: {"--ri","--wavelength-um","--max-reflections","--paired-depth","--trace-cutoff-jones"}){
        auto changed=base;Set(changed,flag,1,"9");assert(Recipe(changed)!=reference);
    }
    auto changed=base;Set(changed,"--particle",2,"500");assert(Recipe(changed)!=reference);
    changed=base;Set(changed,"--particle",3,"156");assert(Recipe(changed)!=reference);
    changed=base;Set(changed,"--particle",1,"10");assert(Recipe(changed)!=reference);
    changed=base;Set(changed,"--ri",2,"0.001");assert(Recipe(changed)!=reference);
    changed=base;Set(changed,"--cutoff-profile",1,"balanced");assert(Recipe(changed)!=reference);
    changed=base;changed.push_back("--no-shadow-beam");assert(Recipe(changed)!=reference);
    changed=base;changed.push_back("--incoherent");assert(Recipe(changed)!=reference);
    changed=base;
    Set(changed,"--sobol-seed",1,"1024");Set(changed,"--sobol-seed",2,"999");Set(changed,"--sobol-offset",1,"4096");
    Set(changed,"--phi-points",1,"64");Set(changed,"--theta-grid-file",1,"second-grid");
    Set(changed,"--threads",1,"4");Set(changed,"--orientation-chunk",1,"115");Set(changed,"--output",1,"second");
    assert(Recipe(changed)==reference);
    std::ofstream(file)<<"geometry A";
    auto loaded=base;loaded.insert(loaded.end(),{"--particle-file",file});const auto first=Recipe(loaded);
    std::ofstream(file)<<"geometry B";assert(Recipe(loaded)!=first);
    std::ofstream(file)<<"geometry A";assert(Recipe(loaded)==first);
    std::ofstream(file+".visibility")<<"visibility A";const auto visible=Recipe(loaded);assert(visible!=first);
    std::ofstream(file+".visibility")<<"visibility B";assert(Recipe(loaded)!=visible);
    std::cout<<"PASS actual PreparedRecipe: geometry, material, wavelength, depth, cutoffs, shadow, coherence and PF/visibility invalidate; seeds/count/offset/angles/output/threads/chunk reuse\n";
}
'''
source = out / 'recipe_probe.cpp'
source.write_text(probe)
binary = out / 'recipe_probe'
subprocess.run(['g++', '-std=c++11', '-O2', '-Wall', '-Wextra', '-I' + str(root / 'src'),
                str(source), str(root / 'src/CliOptions.cpp'), '-o', str(binary)], check=True)
subprocess.run([str(binary), str(out)], check=True)
