#!/usr/bin/env bash
set -euo pipefail

repo=$(cd "$(dirname "$0")/.." && pwd)
build_dir=${TMPDIR:-/tmp}/mbs_big_integer_multiply_add_test
mkdir -p "$build_dir"

${CXX:-g++} -O2 -std=gnu++11 -I"$repo/src" -I"$repo/src/bigint" \
  "$repo/tests/test_big_integer_multiply_add.cpp" \
  "$repo/src/bigint/BigUnsigned.cc" \
  "$repo/src/bigint/BigInteger.cc" \
  -o "$build_dir/test_big_integer_multiply_add"

"$build_dir/test_big_integer_multiply_add"
