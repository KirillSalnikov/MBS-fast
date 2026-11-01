#!/usr/bin/env bash
set -euo pipefail
repo=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
build_dir=$(mktemp -d "${TMPDIR:-/tmp}/mbs_external_fresnel.XXXXXX")
trap 'rm -rf "$build_dir"' EXIT
sources=(Splitting.cpp Beam.cpp geometry/Polygon.cpp geometry/geometry_lib.cpp
    geometry/intrinsic/intrinsics.cpp math/compl.cpp math/JonesMatrix.cpp
    math/matrix.cpp math/PhysMtr.cpp common/common.cpp common/Matrix4x4.cpp
    common/MullerMatrix.cpp)
command=("${CXX:-g++}" -O2 -std=gnu++11 -march=native -fopenmp
    -ffunction-sections -fdata-sections)
for directory in "" math geometry geometry/intrinsic geometry/sse particle scattering common bigint; do
    command+=("-I$repo/src/$directory")
done
for source in "${sources[@]}"; do command+=("$repo/src/$source"); done
command+=("$repo"/src/bigint/*.cc "$repo/tests/test_absorbing_external_reflection.cpp"
    -Wl,--gc-sections -o "$build_dir/check")
"${command[@]}"
OMP_NUM_THREADS=1 "$build_dir/check"
