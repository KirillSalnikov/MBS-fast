#!/usr/bin/env bash
# Run on an idle GPU. The large regression requires <= 16 GiB VRAM so that
# the original 22.94 GB allocation actually fails.
set -euo pipefail
task_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
build="$task_root/tests/.build/analytic_cuda_memory"
mkdir -p "$build"
git -C "$task_root" show 8ae9e858b8f637d724e0bc9b3896b193776d0709:src/cuda/GpuAnalyticFacet.cu > "$build/legacy.cu"
nvcc="${NVCC:-nvcc}"
flags=(-O3 -std=c++11 "-arch=sm_${GPU_ARCH:-70}" -U_GNU_SOURCE -D_DEFAULT_SOURCE -D_XOPEN_SOURCE=700
       -ccbin "${CUDA_HOST_CXX:-g++}" -include "$task_root/src/cuda/CudaHostCompat.h"
       -I"$task_root/src" -I"$task_root/src/cuda")
"$nvcc" "${flags[@]}" -DEvaluateAnalyticFacetGpu=EvaluateAnalyticFacetGpuLegacy \
    -DEvaluateAnalyticAzimuthGaussianGpu=EvaluateAnalyticAzimuthGaussianGpuLegacy \
    -c "$build/legacy.cu" -o "$build/legacy.o"
"$nvcc" "${flags[@]}" -DMBS_COMPARE_LEGACY "$task_root/tests/cuda/analytic_facet_memory.cu" \
    "$task_root/src/cuda/GpuAnalyticFacet.cu" "$build/legacy.o" -o "$build/probe"
"$build/probe" pressure
"$build/probe" large
