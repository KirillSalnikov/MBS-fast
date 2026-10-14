#!/usr/bin/env bash
set -euo pipefail
root=$(cd "$(dirname "$0")/.." && pwd)
work=$(mktemp -d)
trap 'rm -rf "$work"' EXIT
cxx=${CXX:-g++}
nvcc=${NVCC:-${CUDA_PATH:-/usr/local/cuda}/bin/nvcc}
host=${CUDA_HOST_CXX:-$(bash "$root/scripts/select_cuda_host_compiler.sh" "$nvcc" "$cxx")}
arch=${GPU_ARCH:-$(nvidia-smi --query-gpu=compute_cap --format=csv,noheader | head -1 | tr -d '.')}
"$cxx" -O2 -std=c++11 "$root/tests/test_gpu_phase_reduction.cpp" -o "$work/host"
"$work/host"
"$nvcc" -O3 -std=c++11 -ccbin "$host" -arch="sm_$arch" \
    "$root/tests/test_gpu_phase_reduction.cu" -o "$work/device"
"$work/device"
