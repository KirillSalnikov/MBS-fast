#!/bin/bash
set -e

CXX="${CXX:-g++}"
JOBS="${JOBS:-$(nproc 2>/dev/null || echo 1)}"

detect_arch_flags() {
    local cpu_model
    cpu_model="$(lscpu 2>/dev/null | sed -n 's/^Model name:[[:space:]]*//p' | head -1)"

    case "$cpu_model" in
        *"EPYC 7H12"*)
            echo "-march=znver2 -mtune=znver2"
            ;;
        *)
            echo "-march=native -mtune=native"
            ;;
    esac
}

ARCH_FLAGS="${ARCH_FLAGS:-$(detect_arch_flags)}"
CXXFLAGS="${CXXFLAGS:--O3 -std=gnu++11 $ARCH_FLAGS -fopenmp -funroll-loops}"
INCLUDES="-Isrc -Isrc/math -Isrc/handler -Isrc/common -Isrc/geometry \
          -Isrc/geometry/intrinsic -Isrc/geometry/sse -Isrc/particle \
          -Isrc/scattering -Isrc/tracer -Isrc/splitting -Isrc/bigint"

cd "$(dirname "$0")"

# Collect all source files (avoid bigint duplication)
SOURCES=$(find src -not -path '*/bigint/*' -name '*.cpp')
SOURCES="$SOURCES $(find src/bigint -name '*.cc' 2>/dev/null)"

echo "Compiler: $CXX"
echo "Flags: $CXXFLAGS"
echo "Parallel build jobs: $JOBS"
echo "Compiling $(echo "$SOURCES" | wc -w) source files..."

OBJECTS=""
for f in $SOURCES; do
    obj="${f%.*}.o"
    if [ ! -f "$obj" ] || [ "$f" -nt "$obj" ]; then
        echo "  CC $f"
        $CXX $CXXFLAGS $INCLUDES -c "$f" -o "$obj" &
        while [ "$(jobs -rp | wc -l)" -ge "$JOBS" ]; do
            wait -n
        done
    fi
    OBJECTS="$OBJECTS $obj"
done

# Wait for all parallel compilations
wait

echo "Linking..."
mkdir -p bin
$CXX $CXXFLAGS -o bin/mbs_po $OBJECTS -lm -lgomp

echo "Done: bin/mbs_po"
