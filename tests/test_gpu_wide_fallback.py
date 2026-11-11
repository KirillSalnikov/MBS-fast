#!/usr/bin/env python3
"""Check the actual wide-packet guard before any CPU beam preparation.

The probe compiles the production guard with a mocked CPU allocation. Source
placement checks cover all current resident input/packing/trace failure exits.
No CUDA device or optical calculation is used.
"""
from pathlib import Path
import argparse
import subprocess
import tempfile

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument('--sanitize', action='store_true')
args = parser.parse_args()
root = Path(__file__).resolve().parents[1]
source = (root / 'src/scattering/TracerPOTotal.cpp').read_text()
start = source.index('    auto prepareBlock = ')
cpu = source.index('        // std::async workers do not inherit', start)
before_cpu = source[start:cpu]
guard_start = before_cpu.index('        if(wideResident && allowResident')
guard = before_cpu[guard_start:]

# These GPU failures all continue to this single guard, before construction of
# local CPU tracing/handler state and before checking CPU prepared-beam caches.
for prerequisite in ('convex->PrepareGpuInput(', 'EncodeGpuBeam(',
                     'const bool success=packed&&GpuConvexQueue::Trace(',
                     'block.residentFull.reset();block.residentNS.reset();'):
    assert before_cpu.index(prerequisite) < guard_start, prerequisite
assert before_cpu.count('return block;') == 1
success = before_cpu.index('                if(success){')
assert success < before_cpu.index('return block;') < guard_start
assert before_cpu.rfind('#endif') < guard_start
assert source.index('            Scattering *localScatter = ', cpu) > cpu
assert source.index('            HandlerPO localHandler(', cpu) > cpu
assert source.index('PreparedReuse::Get(', cpu) > cpu
assert 'residentEligible && GpuConvexQueue::Capacity::WideChunksEnabled()' in source
assert '&& !EnvEnabled("MBS_GPU_TRACE_HYBRID");' in source
assert 'm_scattering->GetMaxReflections(),wideResident)' in source
assert 'maximumChunk=std::min(maximumChunk,64);' in source

pre = r'''
#include "GpuChunkRetry.h"
#include <cassert>
static void test(bool wideResident, bool allowResident, int thisChunk) {
    bool preparedCpu = false, retried = false;
    try {
'''
post = r'''
        preparedCpu = true; // expensive CPU beam allocation begins here
    } catch(const GpuChunkRetry::Required&) { retried = true; }
    const bool expected = wideResident && allowResident && thisChunk > 128;
    assert(retried == expected);
    assert(preparedCpu == !expected);
}
int main() {
    for(bool enabled : {false, true}) for(bool resident : {false, true})
        for(int count : {1, 64, 128, 129, 256, 512}) test(enabled, resident, count);
}
'''
with tempfile.TemporaryDirectory(prefix='mbs-wide-fallback-') as tmp:
    cpp = Path(tmp) / 'probe.cpp'
    binary = Path(tmp) / 'probe'
    cpp.write_text(pre + guard + post)
    command = ['g++', '-std=c++11', '-O1', '-g', '-Wall', '-Wextra', '-Werror',
               '-I' + str(root / 'src'), str(cpp), '-o', str(binary)]
    if args.sanitize:
        command += ['-fsanitize=address,undefined', '-fno-omit-frame-pointer']
    subprocess.run(command, check=True)
    subprocess.run([str(binary)], check=True)
print('PASS actual wide fallback guard: input/packing/trace exits, no CPU beam allocation, narrow/global unsupported unchanged')
