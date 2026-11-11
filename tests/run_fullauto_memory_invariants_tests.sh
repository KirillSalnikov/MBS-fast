#!/usr/bin/env bash
set -euo pipefail
task_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
task_output="${1:?Pass a new output directory}"
mkdir -p "$task_output/bin"
task_output="$(cd "$task_output" && pwd)"
test ! -e "$task_output/fixture" || { echo "Use a new output directory" >&2; exit 1; }
"${CXX:-g++}" -std=c++11 -O1 -ffunction-sections -fdata-sections \
    -I"$task_root/src" -I"$task_root/src/cuda" \
    "$task_root/tests/fullauto_memory_invariants_probe.cpp" \
    "$task_root/src/RunConfig.cpp" "$task_root/src/AdaptiveConfig.cpp" \
    -Wl,--gc-sections -ldl -pthread -o "$task_output/bin/probe"
"$task_output/bin/probe" "$task_output/fixture" > "$task_output/result.log" 2>&1
cat "$task_output/result.log"
