#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")/.."
test_dir=$(mktemp -d)
trap 'rm -rf "$test_dir"' EXIT
"${CXX:-g++}" -std=c++11 -O1 -ffunction-sections -fdata-sections -Wl,--gc-sections \
    -Isrc -Isrc/cuda tests/fullauto_screen_theta_probe.cpp \
    src/{CliOptions,RuntimeInfo,RunConfig,AdaptiveConfig,OrientationFile,AnalyticFacetAverage,IntegralCharacteristics}.cpp \
    src/common/common.cpp cpu/GpuSupportStub.cpp -o "$test_dir/probe"
"$test_dir/probe" "$test_dir/case"
