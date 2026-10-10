#!/usr/bin/env bash
set -euo pipefail
task_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
mkdir -p "$task_root/tests/.build"
"${CXX:-g++}" -std=c++11 -O2 -Wall -Wextra -I"$task_root/src" -I"$task_root/src/cuda" \
    "$task_root/tests/fullauto_probe.cpp" "$task_root/src/FullAuto.cpp" \
    "$task_root/src/CliOptions.cpp" "$task_root/src/RuntimeInfo.cpp" \
    "$task_root/cpu/GpuSupportStub.cpp" -o "$task_root/tests/.build/fullauto_probe"
"$task_root/tests/.build/fullauto_probe"
