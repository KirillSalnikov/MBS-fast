#!/bin/bash
# =============================================================================
# Build with MPI + OpenMP support
#
# Usage:  bash build_mpi.sh
#         mpirun -np 4 bin/mbs_po_mpi --po --sobol 4096 ...
#
# Hybrid MPI+OpenMP: each MPI rank uses OpenMP for Phase 2.
# Example (2 nodes × 32 cores):
#   mpirun -np 2 --bind-to socket \
#     -x OMP_NUM_THREADS=32 -x OMP_PROC_BIND=close \
#     bin/mbs_po_mpi --po --sobol 8192 ...
# =============================================================================
set -e

# Auto-detect MPI compiler
if command -v mpicxx >/dev/null 2>&1; then
    CXX=mpicxx
elif command -v mpic++ >/dev/null 2>&1; then
    CXX=mpic++
else
    echo "ERROR: mpicxx or mpic++ not found. Install MPI:"
    echo "  sudo apt install libopenmpi-dev"
    exit 1
fi

CXXFLAGS="-O3 -std=gnu++11 -march=native -fopenmp -funroll-loops -DUSE_MPI"

# Add AVX-512 if supported
if grep -q avx512f /proc/cpuinfo 2>/dev/null; then
    CXXFLAGS="$CXXFLAGS -mavx512f -mavx512dq"
fi

INCLUDES="-Isrc -Isrc/math -Isrc/handler -Isrc/common -Isrc/geometry \
          -Isrc/geometry/intrinsic -Isrc/geometry/sse -Isrc/particle \
          -Isrc/scattering -Isrc/tracer -Isrc/splitting -Isrc/bigint"

cd "$(dirname "$0")"

SOURCES=$(find src -not -path '*/bigint/*' -name '*.cpp')
SOURCES="$SOURCES $(find src/bigint -name '*.cc' 2>/dev/null)"

echo "Building MPI+OpenMP version"
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
wait

echo "Linking..."
mkdir -p bin
$CXX $CXXFLAGS -o bin/mbs_po_mpi $OBJECTS -lm -lgomp

echo "Done: bin/mbs_po_mpi"
echo ""
echo "Usage:"
echo "  mpirun -np 4 bin/mbs_po_mpi --po --sobol 4096 \\"
echo "    -p 1 100 100 -w 0.532 --ri 1.31 0 -n 14 \\"
echo "    --grid 0 180 48 181 --auto_tgrid --auto_phi --sym 6 2 --close"
