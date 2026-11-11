#!/usr/bin/env bash
set -euo pipefail
root_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
build_dir="$root_dir/tests/.build"
mkdir -p "$build_dir"
"${CXX:-g++}" -std=c++11 -Wall -Wextra -pedantic -I"$root_dir/src" \
    "$root_dir/tests/test_fullauto_storage_cli.cpp" \
    "$root_dir/src/AdaptiveConfig.cpp" "$root_dir/src/CliOptions.cpp" \
    "$root_dir/src/IntegralCharacteristics.cpp" "$root_dir/src/OrientationFile.cpp" \
    "$root_dir/src/RunConfig.cpp" "$root_dir/src/common/common.cpp" \
    -o "$build_dir/test_fullauto_storage_cli"
"$build_dir/test_fullauto_storage_cli"
