#!/bin/bash
# =============================================================================
# Build script for AMD EPYC 7H12 (Zen 2, SP3 socket) — Clang / AOCC
#
# Compiler selection priority:
#   1. AOCC clang++ (AMD-optimised, if found in PATH or /opt/AMD/aocc-*)
#   2. System clang++
#   3. clang++-16, clang++-15, clang++-14 (versioned Debian/Ubuntu packages)
#
# Key flags:
#   -march=znver2 -mtune=znver2        Zen 2 code-gen (AVX2 + FMA3)
#   -fvectorize -fslp-vectorize         auto-vectoriser (on by default at -O3)
#   -fopenmp                            OpenMP (uses libomp, NOT libgomp)
#
# Usage:  bash build_epyc_clang.sh
#         OMP_PROC_BIND=close OMP_PLACES=cores OMP_NUM_THREADS=64 \
#             bin/mbs_po_epyc_clang --po --sobol 1024 ...
#
# NUMA hint: setting OMP_PROC_BIND=close OMP_PLACES=cores keeps threads
# pinned to physical cores inside one NUMA domain and avoids cross-socket
# memory traffic on dual-socket EPYC boards.
# =============================================================================
set -e

# ---------------------------------------------------------------------------
# Compiler auto-detection
# ---------------------------------------------------------------------------
find_compiler() {
    # 1. AOCC — look for AMD's clang++ in common install paths
    for aocc_dir in /opt/AMD/aocc-*/bin; do
        if [ -x "$aocc_dir/clang++" ]; then
            echo "$aocc_dir/clang++"
            return
        fi
    done
    # Also check if aocc clang++ is already in PATH (module loaded)
    if command -v clang++ >/dev/null 2>&1; then
        if clang++ --version 2>&1 | grep -qi 'AMD'; then
            echo "clang++"
            return
        fi
    fi

    # 2. Plain system clang++
    if command -v clang++ >/dev/null 2>&1; then
        echo "clang++"
        return
    fi

    # 3. Versioned clang++ (try newest first)
    for ver in 18 17 16 15 14; do
        if command -v "clang++-$ver" >/dev/null 2>&1; then
            echo "clang++-$ver"
            return
        fi
    done

    echo ""
}

CXX="${CXX:-$(find_compiler)}"
if [ -z "$CXX" ]; then
    echo "ERROR: no clang++ found. Install clang or AOCC, or set CXX= manually."
    exit 1
fi

# Check for clang's own OpenMP headers (libomp-XX-dev)
if ! echo '#include <omp.h>' | $CXX -x c++ -fopenmp -fsyntax-only - 2>/dev/null; then
    echo "WARNING: clang cannot find omp.h."
    echo "Install libomp-dev:  sudo apt install libomp-16-dev  (match your clang version)"
    echo "Falling back to -fopenmp=libgomp (uses GCC's OpenMP)."
    OMP_FLAG="-fopenmp=libgomp"
else
    OMP_FLAG="-fopenmp"
fi

CXXFLAGS="-O3 -std=gnu++11 -march=znver2 -mtune=znver2 \
$OMP_FLAG -funroll-loops -fvectorize -fslp-vectorize"

INCLUDES="-Isrc -Isrc/math -Isrc/handler -Isrc/common -Isrc/geometry \
          -Isrc/geometry/intrinsic -Isrc/geometry/sse -Isrc/particle \
          -Isrc/scattering -Isrc/tracer -Isrc/splitting -Isrc/bigint"

cd "$(dirname "$0")"

# Collect all source files (avoid bigint duplication)
SOURCES=$(find src -not -path '*/bigint/*' -name '*.cpp')
SOURCES="$SOURCES $(find src/bigint -name '*.cc' 2>/dev/null)"

echo "Building for EPYC 7H12 (Zen 2, AVX2+FMA3) — Clang/AOCC"
echo "Compiler: $CXX ($($CXX --version 2>&1 | head -1))"
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
# libomp for clang, libgomp for gcc — try both
if $CXX -lomp -x c++ /dev/null -o /dev/null 2>/dev/null; then
    OMP_LIB="-lomp"
else
    OMP_LIB="-lgomp"
fi
$CXX $CXXFLAGS -o bin/mbs_po_epyc_clang $OBJECTS -lm $OMP_LIB

echo "Done: bin/mbs_po_epyc_clang"
echo ""
echo "Recommended usage on EPYC 7H12 (64 cores):"
echo "  OMP_PROC_BIND=close OMP_PLACES=cores OMP_NUM_THREADS=64 \\"
echo "    bin/mbs_po_epyc_clang --po --sobol 4096 \\"
echo "    -p 1 100 100 -w 0.532 --ri 1.31 0 -n 10 \\"
echo "    --grid 0 180 48 181 --sym 6 2 --close"
