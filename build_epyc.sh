#!/bin/bash
# =============================================================================
# Build script for AMD EPYC 7H12 (Zen 2, SP3 socket)
#
# Key differences from Intel build:
#   - -march=znver2 (no AVX-512, uses AVX2 + FMA3)
#   - fast_sincos_8x falls back to 2 × fast_sincos_4x (AVX2)
#   - rsqrt uses SSE rsqrt_ss + Newton-Raphson (instead of AVX-512 rsqrt14)
#   - 64 cores → OpenMP with OMP_NUM_THREADS=64 recommended
#
# Usage:  bash build_epyc.sh
#         OMP_NUM_THREADS=64 bin/mbs_po_epyc --po --sobol 1024 ...
# =============================================================================
set -e

CXX="${CXX:-g++}"
CXXFLAGS="-O3 -std=gnu++11 -march=znver2 -mtune=znver2 -fopenmp -funroll-loops"
INCLUDES="-Isrc -Isrc/math -Isrc/handler -Isrc/common -Isrc/geometry \
          -Isrc/geometry/intrinsic -Isrc/geometry/sse -Isrc/particle \
          -Isrc/scattering -Isrc/tracer -Isrc/splitting -Isrc/bigint"

cd "$(dirname "$0")"

# Collect all source files (avoid bigint duplication)
SOURCES=$(find src -not -path '*/bigint/*' -name '*.cpp')
SOURCES="$SOURCES $(find src/bigint -name '*.cc' 2>/dev/null)"

echo "Building for EPYC 7H12 (Zen 2, AVX2+FMA3)"
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

# Wait for all parallel compilations
wait

echo "Linking..."
mkdir -p bin
$CXX $CXXFLAGS -o bin/mbs_po_epyc $OBJECTS -lm -lgomp

echo "Done: bin/mbs_po_epyc"
echo ""
echo "Recommended usage on EPYC 7H12 (64 cores):"
echo "  OMP_PROC_BIND=close OMP_PLACES=cores OMP_NUM_THREADS=64 \\"
echo "    bin/mbs_po_epyc --po --sobol 4096 \\"
echo "    -p 1 100 100 -w 0.532 --ri 1.31 0 -n 10 \\"
echo "    --grid 0 180 48 181 --sym 6 2 --close"
echo ""
echo "NUMA: OMP_PROC_BIND=close keeps threads on same chiplet (important for EPYC)."
