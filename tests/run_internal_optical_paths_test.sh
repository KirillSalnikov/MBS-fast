#!/usr/bin/env bash
set -euo pipefail
repo=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
build_dir=$(mktemp -d "${TMPDIR:-/tmp}/mbs_internal_paths.XXXXXX")
trap 'rm -rf "$build_dir"' EXIT
command=("${CXX:-g++}" -O3 -std=gnu++11 -march=native -fopenmp)
for directory in "" math geometry geometry/intrinsic geometry/sse particle scattering handler tracer common bigint cuda; do
    command+=("-I$repo/src/$directory")
done
sources=(Beam.cpp Splitting.cpp Tracks.cpp CutoffStatistics.cpp scattering/Scattering.cpp scattering/ScatteringConvex.cpp particle/Particle.cpp particle/Hexagonal.cpp particle/Bullet.cpp particle/Droxtal.cpp geometry/geometry_lib.cpp geometry/Facet.cpp geometry/Polygon.cpp geometry/Intersection.cpp geometry/intrinsic/intrinsics.cpp math/JonesMatrix.cpp math/PhysMtr.cpp math/Sobol.cpp math/Mueller.cpp math/matrix.cpp math/compl.cpp common/common.cpp common/MullerMatrix.cpp common/Matrix4x4.cpp bigint/BigIntegerAlgorithms.cc bigint/BigIntegerUtils.cc bigint/BigUnsignedInABase.cc bigint/BigInteger.cc bigint/BigUnsigned.cc scattering/ScatteringNonConvex.cpp)
for source in "${sources[@]}"; do command+=("$repo/src/$source"); done
command+=("$repo/tests/test_internal_optical_paths.cpp" -o "$build_dir/check")
"${command[@]}"
OMP_NUM_THREADS=1 "$build_dir/check" "${1:-$repo/tests/data/phase_dust15.dat}"
