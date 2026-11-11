#!/usr/bin/env python3
"""Exercise actual trace ownership classes with a checked CUDA mock; no GPU."""
from pathlib import Path
import argparse
import subprocess
import tempfile

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument('--sanitize', action='store_true')
parser.add_argument('--output', type=Path)
args = parser.parse_args()
root = Path(__file__).resolve().parents[1]
temporary = tempfile.TemporaryDirectory(prefix='mbs-trace-pool-') if args.output is None else None
out = args.output or Path(temporary.name)
out.mkdir(parents=True, exist_ok=True)
text=(root/'src/cuda/GpuConvexQueue.cu').read_text();start=text.index('struct Buffers {');end=text.index('#include "GpuTraceWorkspace.cuh"',start)
preamble=r'''
#include "cuda/GpuTraceCapacity.h"
#include <cassert>
#include <cstdio>
#include <cstdlib>
#include <map>
#include <memory>
#include <atomic>
#include <stdexcept>
using cudaError_t=int;static const int cudaSuccess=0;using cudaStream_t=void*;using cudaEvent_t=void*;
struct Allocation {size_t bytes;int device;};
static std::map<void*,Allocation> live;static size_t sequence=0,allocations=0,frees=0;
static int currentDevice=0,failSetTarget=-1;static bool failFree=false,failSync=false,failGet=false,failQuery=false,failEvents=false;
static void* failOne=nullptr;
static int cudaGetDevice(int* device){if(failGet)return 1;*device=currentDevice;return 0;}
static int cudaSetDevice(int device){if(device==failSetTarget)return 1;currentDevice=device;return 0;}
static int cudaMemGetInfo(size_t* free,size_t* total){if(failQuery)return 1;*free=size_t(8)*1024*1024*1024;*total=*free*2;return 0;}
template<class T>static int cudaMalloc(T**pointer,size_t bytes){*pointer=reinterpret_cast<T*>((++sequence)*64);live[*pointer]={bytes,currentDevice};++allocations;return 0;}
static int cudaFree(void* pointer){if(!pointer)return 0;assert(live.count(pointer)&&live[pointer].device==currentDevice);if(failFree||pointer==failOne)return 1;live.erase(pointer);++frees;return 0;}
static int cudaStreamSynchronize(cudaStream_t){return failSync?1:0;}
static int cudaEventDestroy(cudaEvent_t){return failEvents?1:0;}
static int cudaStreamDestroy(cudaStream_t){return 0;}
namespace GpuConvexQueue {struct QueueItem {Beam beam;int orientation;};
'''
body=text[start:end]+'\n#include "cuda/GpuTraceWorkspace.cuh"\n}\n'
tests=r'''
using namespace GpuConvexQueue;
static void allocate(Buffers& b){assert(b.Allocate(&b.inputs,128,Buffers::Inputs)==0);assert(b.Allocate(&b.output,256,Buffers::OutputData)==0);assert(b.Allocate(&b.queues[0],512,Buffers::Queue0)==0);}
static void resetQuarantine(){failFree=false;failOne=nullptr;failSetTarget=-1;failGet=false;failSync=false;failEvents=false;
 if(SavedWorkspace()){assert(SavedWorkspace()->ReleaseAllOwned());}
 SavedWorkspace().reset();WorkspaceQuarantined().store(false);assert(live.empty());}
int main(){
 {Buffers b(0);allocate(b);auto* original=b.output;failFree=true;
  assert(!b.ReleaseAllOwned()&&b.Bytes()==896&&b.output==original&&live.size()==3);
  assert(!b.ReleaseAllOwned()&&b.Bytes()==896&&live.size()==3);
  failFree=false;assert(b.ReleaseAllOwned()&&b.Bytes()==0&&live.empty());}
 {Buffers b(0);allocate(b);failOne=b.output;
  assert(!b.ReleaseAllOwned()&&b.Bytes()==256&&live.size()==1&&b.output==failOne);
  assert(!b.ReleaseAllOwned()&&b.Bytes()==256&&live.size()==1);
  failOne=nullptr;assert(b.ReleaseAllOwned()&&live.empty());}
 {currentDevice=1;Buffers b(1);allocate(b);currentDevice=0;failSetTarget=0;
  assert(!b.ReleaseAllOwned()&&live.empty()&&currentDevice==1);
  failSetTarget=-1;currentDevice=0;assert(b.ReleaseAllOwned()&&currentDevice==0);}
 {Buffers b(0);allocate(b);b.stream=reinterpret_cast<void*>(1);failSync=true;
  assert(!b.ReleaseAllOwned()&&live.size()==3&&b.Bytes()==896);
  failSync=false;assert(b.ReleaseAllOwned()&&live.empty()&&!b.stream);}
 {Buffers b(0);allocate(b);auto* original=b.inputs;failFree=true;
  assert(b.Allocate(&b.inputs,1024,Buffers::Inputs)!=0&&b.inputs==original&&b.sizes[Buffers::Inputs]==128&&!b.healthy);
  assert(b.DiscardMismatch(b.inputs,Buffers::Inputs,1024)!=0&&b.inputs==original&&b.sizes[Buffers::Inputs]==128);
  failFree=false;assert(b.ReleaseAllOwned());}
 {Buffers b(0);b.start=reinterpret_cast<void*>(1);b.stop=reinterpret_cast<void*>(2);b.stream=reinterpret_cast<void*>(3);failEvents=true;
  assert(!b.ReleaseAllOwned()&&b.start&&b.stop&&b.stream);failEvents=false;assert(b.ReleaseAllOwned()&&!b.start&&!b.stop&&!b.stream);}
 setenv("MBS_GPU_TRACE_BUFFER_POOL","1",1);unsetenv("MBS_GPU_TRACE_CONVEX_WORKSPACE_MIB");
 {WorkspaceLease lease(0);allocate(*lease.buffer);failFree=true;}
 assert(WorkspaceQuarantined().load()&&SavedWorkspace()&&SavedWorkspace()->Bytes()==896&&live.size()==3);
 size_t before=allocations;bool threw=false;try{WorkspaceLease lease(0);}catch(const std::runtime_error&){threw=true;}
 assert(threw&&allocations==before&&live.size()==3);
 // Persistent failure cannot discard the quarantined owner's pointers.
 assert(!RetireWorkspace(SavedWorkspace())&&SavedWorkspace()&&live.size()==3);resetQuarantine();
 {WorkspaceLease lease(0);allocate(*lease.buffer);lease.success=true;}
 assert(SavedWorkspace()&&!WorkspaceQuarantined().load());failFree=true;setenv("MBS_GPU_TRACE_CONVEX_WORKSPACE_MIB","0",1);threw=false;
 try{WorkspaceLease lease(0);}catch(const std::runtime_error&){threw=true;}
 assert(threw&&WorkspaceQuarantined().load()&&SavedWorkspace()&&SavedWorkspace()->Bytes()==896);resetQuarantine();
 // Failure to restore the previous device is quarantined even after all bytes freed.
 unsetenv("MBS_GPU_TRACE_CONVEX_WORKSPACE_MIB");currentDevice=1;
 {WorkspaceLease lease(1);allocate(*lease.buffer);lease.success=true;}
 currentDevice=0;failSetTarget=0;threw=false;
 try{WorkspaceLease lease(0);}catch(const std::runtime_error&){threw=true;}
 assert(threw&&WorkspaceQuarantined().load()&&SavedWorkspace()&&live.empty());resetQuarantine();
 assert(allocations==frees);std::puts("PASS actual Buffers + WorkspaceLease: repeated/partial free failures, sync/event/device-restore failure, quarantine ownership and no next allocation");
}
'''
source=out/'actual_buffers_probe.cpp';source.write_text(preamble+body+tests);binary=out/('actual_buffers_san'if args.sanitize else'actual_buffers')
cmd=['g++','-std=c++11','-O1','-g','-Wall','-Wextra','-Werror','-I'+str(root/'src'),str(source),'-o',str(binary)]
if args.sanitize:cmd+=['-fsanitize=address,undefined','-fno-omit-frame-pointer']
subprocess.run(cmd,check=True);subprocess.run([str(binary)],check=True)
