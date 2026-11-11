#!/usr/bin/env python3
"""Compile actual resident host ownership classes against a checked CUDA mock.

No GPU is accessed. This exercises cross-thread cache lifetime, opt-out,
device changes, pressure, oversize, failed allocation/free/sync and exceptions.
"""
import argparse
from pathlib import Path
import subprocess
import tempfile

PREAMBLE = r'''
#include <algorithm>
#include <atomic>
#include <cassert>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <map>
#include <memory>
#include <mutex>
#include <stdexcept>
#include <thread>
#include <vector>
#include "cuda/GpuResidentPoolPolicy.h"
using cudaError_t=int;
static const int cudaSuccess=0;
using GpuReal=double;
struct GpuBeam {double values[32];};
struct ResidentSummary {double values[8];};
struct Allocation {size_t bytes;int device;};
static std::map<void*,Allocation> live;
static size_t allocated=0,freed=0,sequence=0,externalBytes=0,memsets=0;
static bool failAllocate=false,failFree=false,failSync=false,failQuery=false;
static size_t getDeviceCalls=0;
static thread_local int currentDevice=0;
static const size_t MiB=1024*1024;
static int cudaGetDevice(int* device){++getDeviceCalls;*device=currentDevice;return 0;}
static int cudaSetDevice(int device){currentDevice=device;return 0;}
static const char* cudaGetErrorString(int){return "injected mock error";}
static int cudaMemGetInfo(size_t* available,size_t* total){
 if(failQuery)return 1;*total=size_t(16)*1024*MiB;size_t used=externalBytes;
 for(const auto& item:live)if(item.second.device==currentDevice)used+=item.second.bytes;
 assert(used<=*total);*available=*total-used;return 0;
}
static int cudaMalloc(void** pointer,size_t bytes){
 if(failAllocate)return 1;
 *pointer=reinterpret_cast<void*>((++sequence)*32);live[*pointer]={bytes,currentDevice};++allocated;return 0;
}
template<class T>static int cudaMalloc(T**pointer,size_t bytes){return cudaMalloc(reinterpret_cast<void**>(pointer),bytes);}
static int cudaFree(void* pointer){
 if(!pointer)return 0;
 assert(live.count(pointer)&&live[pointer].device==currentDevice);
 if(failFree)return 1;live.erase(pointer);++freed;return 0;
}
static int cudaMemset(void* pointer,int,size_t bytes){
 assert(live.count(pointer)&&live[pointer].device==currentDevice&&bytes<=live[pointer].bytes);++memsets;return 0;
}
static int cudaStreamSynchronize(int){return failSync?1:0;}
static double gpu_memory_fraction(){return .2;}
'''

TEST = r'''
static GpuResidentPoolPolicy::Sizes request(){GpuResidentPoolPolicy::Sizes values={};values[0]=128;return values;}
static void execute(bool expected=true){
 ResidentAllocations a(request());int* pointer=nullptr;
 const bool ok=a.make(pointer,32,true);assert(ok==expected);
 if(ok)assert(a.Finish());
}
static void clear(){auto& cache=BoundedResidentCache();std::lock_guard<std::mutex> guard(cache.mutex);assert(cache.Clear());}
int main(){
 assert(ReleaseIdleBoundedResidentCache()&&getDeviceCalls==0);
 unsetenv("MBS_GPU_RESIDENT_BOUNDED_REUSE");setenv("MBS_GPU_RESIDENT_REUSE","0",1);
 execute();assert(live.empty());
 setenv("MBS_GPU_RESIDENT_BOUNDED_REUSE","1",1);
 execute();assert(live.empty()); // explicit old opt-out wins
 unsetenv("MBS_GPU_RESIDENT_REUSE");
 size_t first=allocated;
 std::thread producer([](){execute();});producer.join();assert(live.size()==1&&allocated==first+1);
 std::thread replacement([](){execute();});replacement.join();assert(live.size()==1&&allocated==first+1);
 // Ineligible next request releases only the bounded cache at its boundary.
 assert(ReleaseIdleBoundedResidentCache()&&live.empty());execute();
 // Device transfer frees old storage on its owner before using another device.
 cudaSetDevice(1);execute();assert(live.size()==1&&live.begin()->second.device==1);cudaSetDevice(0);
 // Disabling after a prior bounded request actually releases its reservation.
 setenv("MBS_GPU_RESIDENT_REUSE","0",1);execute();assert(live.empty());unsetenv("MBS_GPU_RESIDENT_REUSE");
 execute();assert(live.size()==1);
 // Exception after the caller's work invalidates the unfinished lease.
 try{ResidentAllocations a(request());int*p=nullptr;assert(a.make(p,32));throw 7;}
 catch(int x){assert(x==7);}assert(live.empty());
 failAllocate=true;execute(false);failAllocate=false;assert(live.empty());
 execute();assert(live.size()==1);
 {ResidentAllocations a(request());int*p=nullptr;assert(a.make(p,32));failSync=true;assert(!a.Finish());}
 failSync=false;assert(live.empty());
 // Failed free retains ownership/accounting, poisons reuse, and is retryable.
 execute();failFree=true;
 assert(!ReleaseIdleBoundedResidentCache()&&BoundedResidentCache().storage.poisoned&&live.size()==1);
 failFree=false;assert(ReleaseIdleBoundedResidentCache()&&live.empty());
 // Actual shape preflight rejects size overflow before any allocation.
 GpuResidentPoolPolicy::Sizes sizes={};
 assert(!ResidentRequestSizes(size_t(-1),1,1,1,false,sizes));
 assert(ResidentRequestSizes(10,8,4,5,false,sizes)&&sizes[13]==0);
 assert(ResidentRequestSizes(10,8,4,5,true,sizes)&&sizes[13]==sizes[12]);
 // Legal oversized work is ephemeral and frees storage before publication.
 sizes={};sizes[0]=size_t(2)*1024*MiB;
 {BoundedResidentAllocations a(sizes);char*p=nullptr;assert(a.make(p,sizes[0],false));assert(a.Finish());assert(live.empty());}
 // Low free-memory pressure clears even a currently valid small cache.
 externalBytes=size_t(16)*1024*MiB-256*MiB;
 execute();assert(live.empty());externalBytes=0;
 // Failed budget query never borrows a stale pointer.
 execute();failQuery=true;execute(false);failQuery=false;assert(live.empty());
 execute();clear();assert(live.empty()&&allocated==freed);
 assert(memsets>0);std::puts("PASS actual bounded resident lifecycle with checked CUDA mock");
}
'''


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--sanitize', action='store_true')
    options = parser.parse_args()
    root = Path(__file__).resolve().parents[1]
    source = (root/'src/cuda/GpuResidentPrepare.cuh').read_text()
    start = source.index('struct ResidentPool {')
    end = source.index('// The resident geometry has sparse per-orientation slots', start)
    with tempfile.TemporaryDirectory(prefix='mbs-resident-lifecycle-') as temp:
        cpp = Path(temp)/'probe.cpp'; binary = Path(temp)/'probe'
        cpp.write_text(PREAMBLE + source[start:end] + TEST)
        command = ['g++', '-std=c++11', '-O1', '-g', '-pthread', '-I'+str(root/'src'), str(cpp), '-o', str(binary)]
        if options.sanitize:
            command += ['-fsanitize=address,undefined', '-fno-omit-frame-pointer']
        subprocess.run(command, check=True)
        subprocess.run([str(binary)], check=True)


if __name__ == '__main__':
    main()
