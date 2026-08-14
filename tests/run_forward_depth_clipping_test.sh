#!/usr/bin/env bash
set -euo pipefail

repo=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
build_dir=${TMPDIR:-/tmp}/mbs_fast_forward_depth_test
mkdir -p "$build_dir"

${CXX:-g++} -std=gnu++11 -msse4.2 \
  -I"$repo/src" -I"$repo/src/common" -I"$repo/src/math" \
  -I"$repo/src/geometry" \
  -I"$repo/src/bigint" \
  "$repo/tests/test_forward_depth_clipping.cpp" \
  -o "$build_dir/test_forward_depth_clipping"

"$build_dir/test_forward_depth_clipping"
