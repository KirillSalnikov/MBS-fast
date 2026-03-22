#!/bin/bash
set -e

CXX="g++"
CXXFLAGS="-O3 -std=gnu++11 -march=native -mavx512f -mavx512dq -fopenmp -funroll-loops"
INCLUDES="-Isrc -Isrc/math -Isrc/handler -Isrc/common -Isrc/geometry \
          -Isrc/geometry/intrinsic -Isrc/geometry/sse -Isrc/particle \
          -Isrc/scattering -Isrc/tracer -Isrc/splitting -Isrc/bigint"

cd "$(dirname "$0")"

# Collect all source files (avoid bigint duplication)
SOURCES=$(find src -not -path '*/bigint/*' -name '*.cpp')
SOURCES="$SOURCES $(find src/bigint -name '*.cc' 2>/dev/null)"

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
$CXX $CXXFLAGS -o bin/mbs_po $OBJECTS -lm -lgomp

echo "Done: bin/mbs_po"
