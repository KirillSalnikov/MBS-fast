#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
TMP="$(mktemp -d)"
trap 'rm -rf "$TMP"' EXIT

dry_run() {
    local precision="$1"
    local fast_math="$2"
    make -s -n -C "$ROOT/gpu" GPU_ARCH=70 \
        GPU_PRECISION="$precision" GPU_FAST_MATH="$fast_math" all
}

root_dry_run() {
    local precision="$1"
    local fast_math="$2"
    make -s -n -C "$ROOT" USE_CUDA=1 GPU_ARCH=70 \
        GPU_PRECISION="$precision" GPU_FAST_MATH="$fast_math" \
        TARGET="bin/profile_test_${precision}_${fast_math}" all
}

dry_run fp64 0 >"$TMP/fp64"
grep -Fq -- '-DMBS_GPU_FP64' "$TMP/fp64"
grep -Fq -- 'build/fp64_lto_sm70/obj' "$TMP/fp64"
grep -Fq -- 'bin/mbs_po_gpu_double' "$TMP/fp64"
grep -Fq -- '-flto=auto' "$TMP/fp64"
if grep -Fq -- '-DMBS_GPU_FLOAT' "$TMP/fp64"; then
    echo 'FP64 profile unexpectedly defines the FP32 compatibility macro' >&2
    exit 1
fi
if grep -Fq -- '--use_fast_math' "$TMP/fp64"; then
    echo 'Precise FP64 profile unexpectedly enables CUDA fast math' >&2
    exit 1
fi

dry_run fp32 0 >"$TMP/fp32"
grep -Fq -- '-DMBS_GPU_FP32 -DMBS_GPU_FLOAT' "$TMP/fp32"
grep -Fq -- 'build/fp32_lto_sm70/obj' "$TMP/fp32"
grep -Fq -- 'bin/mbs_po_gpu_float' "$TMP/fp32"
if grep -Fq -- '--use_fast_math' "$TMP/fp32"; then
    echo 'Precise FP32 profile unexpectedly enables CUDA fast math' >&2
    exit 1
fi

dry_run fp32 1 >"$TMP/fp32_fast"
grep -Fq -- '--use_fast_math' "$TMP/fp32_fast"
grep -Fq -- 'build/fp32_fast_lto_sm70/obj' "$TMP/fp32_fast"
grep -Fq -- 'bin/mbs_po_gpu_float_fast' "$TMP/fp32_fast"

root_dry_run fp64 0 >"$TMP/root_fp64"
grep -Fq -- '-DMBS_GPU_FP64' "$TMP/root_fp64"
grep -Fq -- '-arch=sm_70' "$TMP/root_fp64"
grep -Fq -- 'build/root_cuda/fp64_sm70/obj' "$TMP/root_fp64"
if grep -Fq -- 'build/root_cuda/fp32' "$TMP/root_fp64"; then
    echo 'Root FP64 profile unexpectedly reuses FP32 objects' >&2
    exit 1
fi

root_dry_run fp32 1 >"$TMP/root_fp32_fast"
grep -Fq -- '-DMBS_GPU_FP32 -DMBS_GPU_FLOAT' "$TMP/root_fp32_fast"
grep -Fq -- '--use_fast_math' "$TMP/root_fp32_fast"
grep -Fq -- 'build/root_cuda/fp32_fast_sm70/obj' "$TMP/root_fp32_fast"
if grep -Fq -- 'build/root_cuda/fp64' "$TMP/root_fp32_fast"; then
    echo 'Root FP32 profile unexpectedly reuses FP64 objects' >&2
    exit 1
fi

make -s -n -C "$ROOT/gpu" GPU_ARCH=70 float_debug >"$TMP/fp32_debug"
grep -Fq -- '-DMBS_GPU_FP32 -DMBS_GPU_FLOAT' "$TMP/fp32_debug"
grep -Fq -- '-DMBS_DEBUG_HELP' "$TMP/fp32_debug"
if grep -Fq -- '-DMBS_GPU_FP64' "$TMP/fp32_debug"; then
    echo 'FP32 debug profile unexpectedly inherits the outer FP64 macro' >&2
    exit 1
fi

if make -s -n -C "$ROOT/gpu" GPU_PRECISION=invalid clean \
        >"$TMP/invalid_precision" 2>&1; then
    echo 'Invalid GPU_PRECISION was accepted' >&2
    exit 1
fi
if make -s -n -C "$ROOT/gpu" GPU_FAST_MATH=2 clean \
        >"$TMP/invalid_fast_math" 2>&1; then
    echo 'Invalid GPU_FAST_MATH was accepted' >&2
    exit 1
fi
if make -s -n -C "$ROOT/gpu" GPU_LTO=2 clean \
        >"$TMP/invalid_lto" 2>&1; then
    echo 'Invalid GPU_LTO was accepted' >&2
    exit 1
fi

make -s -n -C "$ROOT/gpu" GPU_ARCH=70 GPU_PRECISION=fp64 \
    GPU_FAST_MATH=0 GPU_LTO=0 all >"$TMP/fp64_no_lto"
grep -Fq -- 'build/fp64_sm70/obj' "$TMP/fp64_no_lto"
if grep -Fq -- '-flto=auto' "$TMP/fp64_no_lto"; then
    echo 'GPU_LTO=0 unexpectedly enables link-time optimization' >&2
    exit 1
fi

echo 'CUDA precision profiles: PASS'
