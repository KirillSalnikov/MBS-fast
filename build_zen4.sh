#!/bin/bash
# =============================================================================
# Build script for AMD Zen 4 (Ryzen 7000 / EPYC 9004 Genoa)
#
# Zen 4 supports AVX-512 (F, DQ, VL, BW, CD, VBMI, VBMI2, VNNI, BITALG, etc.)
# Uses 256-bit execution units, AVX-512 ops decoded as 2x256-bit.
# Still faster than AVX2-only due to wider register file and mask ops.
#
# Usage:  bash build_zen4.sh
#         OMP_NUM_THREADS=16 bin/mbs_po_zen4 --po --sobol 1024 ...
# =============================================================================
set -e

CXX="${CXX:-g++}"
CXXFLAGS="-O3 -std=gnu++11 -march=znver4 -mtune=znver4 -mavx512f -mavx512dq -mavx512vl -fopenmp -funroll-loops"
INCLUDES="-Isrc -Isrc/math -Isrc/handler -Isrc/common -Isrc/geometry \
          -Isrc/geometry/intrinsic -Isrc/geometry/sse -Isrc/particle \
          -Isrc/scattering -Isrc/tracer -Isrc/splitting -Isrc/bigint"

cd "$(dirname "$0")"

SOURCES=$(find src -not -path '*/bigint/*' -name '*.cpp')
SOURCES="$SOURCES $(find src/bigint -name '*.cc' 2>/dev/null)"

echo "Building for Zen 4 (Ryzen 7000 / EPYC Genoa, AVX-512)"
echo "Compiler: $CXX"
echo "Flags: $CXXFLAGS"
echo "Compiling $(echo "$SOURCES" | wc -w) source files..."

OBJECTS=""
for f in $SOURCES; do
    obj="${f%.*}.o"
    if [ ! -f "$obj" ] || [ "$f" -nt "$obj" ]; then
        echo "  CC $f"
        $CXX $CXXFLAGS $INCLUDES -c "$f" -o "$obj" &
    fi
    OBJECTS="$OBJECTS $obj"
done
wait

echo "Linking..."
mkdir -p bin
$CXX $CXXFLAGS -o bin/mbs_po_zen4 $OBJECTS -lm -lgomp

echo "Done: bin/mbs_po_zen4"
