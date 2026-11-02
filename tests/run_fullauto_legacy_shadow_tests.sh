#!/usr/bin/env bash
set -euo pipefail
task_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
task_output="${1:?Supply a private test output directory}"
mkdir -p "$task_output"
task_output="$(cd "$task_output" && pwd)"
"${CXX:-g++}" -std=c++11 -O1 -ffunction-sections -fdata-sections -Wl,--gc-sections \
  -I"$task_root/src" -I"$task_root/src/cuda" \
  "$task_root/tests/fullauto_legacy_shadow_probe.cpp" \
  "$task_root/src/CliOptions.cpp" "$task_root/src/RuntimeInfo.cpp" \
  "$task_root/src/RunConfig.cpp" "$task_root/src/AdaptiveConfig.cpp" \
  "$task_root/src/OrientationFile.cpp" \
  "$task_root/src/common/common.cpp" \
  "$task_root/src/IntegralCharacteristics.cpp" "$task_root/cpu/GpuSupportStub.cpp" \
  -o "$task_output/probe"
"$task_output/probe" "$task_output/case"
