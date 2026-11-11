#!/usr/bin/env bash
set -euo pipefail
task_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
task_output="${1:?Pass a fresh diagnostic output directory}"
task_cpu_binary="${2:?Pass an existing real CPU optical binary}"
test ! -e "$task_output"
mkdir -p "$task_output"
task_output="$(cd "$task_output" && pwd)"
"${CXX:-c++}" -std=c++11 -O2 -Wall -Wextra -Werror -pthread \
    "$task_root/tests/test_gpu_chunk_retry.cpp" -o "$task_output/unit"
"$task_output/unit"
"${CXX:-c++}" -std=c++11 -O1 -ffunction-sections -fdata-sections \
    -I"$task_root/src" -I"$task_root/src/cuda" \
    "$task_root/tests/gpu_chunk_retry_probe.cpp" "$task_root/src/RunConfig.cpp" \
    "$task_root/src/AdaptiveConfig.cpp" -Wl,--gc-sections -ldl -pthread -o "$task_output/probe"
python3 "$task_root/tests/test_gpu_chunk_retry_integration.py" \
    --probe "$task_output/probe" --binary "$task_cpu_binary" --output "$task_output/integration"
