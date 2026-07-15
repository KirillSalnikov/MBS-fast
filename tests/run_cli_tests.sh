#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
BUILD_DIR="$ROOT_DIR/tests/.build"
CXX_BIN="${CXX:-g++}"

mkdir -p "$BUILD_DIR"
"$CXX_BIN" -std=c++11 -Wall -Wextra -pedantic \
    -I"$ROOT_DIR/src" \
    "$ROOT_DIR/tests/test_cli.cpp" \
    "$ROOT_DIR/src/AdaptiveConfig.cpp" \
    "$ROOT_DIR/src/CliOptions.cpp" \
    "$ROOT_DIR/src/IntegralCharacteristics.cpp" \
    "$ROOT_DIR/src/RunConfig.cpp" \
    -o "$BUILD_DIR/test_cli"

MBS_TEST_ROOT="$ROOT_DIR" "$BUILD_DIR/test_cli"
