#!/usr/bin/env python3
"""Compile the actual consumer scheduling fragment with a checked producer mock.

No GPU or optical calculation. An active resident+facet packet must finish its
consumer work before the next GPU request, including the exception path.
"""
from pathlib import Path
import argparse
import subprocess
import tempfile

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument('--sanitize', action='store_true')
args = parser.parse_args()
root = Path(__file__).resolve().parents[1]
source = (root/'src/scattering/TracerPOTotal.cpp').read_text()
start = source.index('                const bool more=nextStart<myEnd;')
end = source.index('                if(recyclePrepared){', start)
fragment = source[start:end]
pre = r'''
#include <algorithm>
#include <cassert>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>
struct Block { bool residentFull; };
struct Chunks { int Next(int remaining) { return std::min(4, remaining); } };
static void test(bool orderedResidentFacet, bool resident, bool moreExpected, bool fail) {
    int nextStart = 4, myEnd = moreExpected ? 12 : 4;
    int prefixBoundary = 6;
    Chunks memoryChunks;
    std::unique_ptr<Block> block(new Block{resident});
    std::vector<std::string> events;
    bool nextGpuStarted = false;
    const auto prefixBoundedChunk = [&](int start, int size) {
        return start < prefixBoundary ? std::min(size, prefixBoundary-start) : size;
    };
    const auto request = [&](int start, int count) {
        assert(start == 4 && count == 2); events.push_back("request");
        nextGpuStarted = true; // next packet can succeed on GPU after current fallback
    };
    const auto processBlock = [&](Block&) {
        if(orderedResidentFacet) assert(!nextGpuStarted);
        events.push_back("controls");
        if(fail) throw std::runtime_error("injected analytic failure");
        events.push_back("commit");
    };
    try {
'''
post = r'''
    } catch(const std::runtime_error&) { assert(fail); }
    const bool ordered = orderedResidentFacet && moreExpected;
    if(ordered) {
        if(fail) assert(events == std::vector<std::string>{"controls"} && nextStart == 4);
        else assert(events == std::vector<std::string>({"controls", "commit", "request"}) && nextStart == 6);
    } else if(moreExpected) {
        assert(events.front() == "request" && nextStart == 6);
    } else {
        assert(events.front() == "controls" && nextStart == 4);
    }
}
int main() {
    for(bool mode : {false,true}) for(bool resident : {false,true})
        for(bool more : {false,true}) for(bool fail : {false,true}) test(mode,resident,more,fail);
}
'''
with tempfile.TemporaryDirectory(prefix='mbs-resident-facet-schedule-') as tmp:
    cpp = Path(tmp)/'probe.cpp'; binary = Path(tmp)/'probe'
    cpp.write_text(pre+fragment+post)
    command = ['g++', '-std=c++11', '-O1', '-g', '-Wall', '-Wextra', '-Werror', str(cpp), '-o', str(binary)]
    if args.sanitize:
        command += ['-fsanitize=address,undefined', '-fno-omit-frame-pointer']
    subprocess.run(command, check=True)
    subprocess.run([str(binary)], check=True)
print('PASS actual scheduler: ordered resident controls, exception no-next-launch, fallback-to-GPU transition, legacy CPU order, prefix boundary')
