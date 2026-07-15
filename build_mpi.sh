#!/bin/bash
# =============================================================================
# Build with MPI + OpenMP support
#
# Usage:  bash build_mpi.sh
#         mpirun -np 4 cpu/bin/mbs_po_mpi --method po --sobol 4096 ...
#
# Hybrid MPI+OpenMP: each MPI rank uses OpenMP for Phase 2.
# Example (2 nodes × 32 cores):
#   mpirun -np 2 --bind-to socket \
#     -x OMP_NUM_THREADS=32 -x OMP_PROC_BIND=close \
#     cpu/bin/mbs_po_mpi --method po --sobol 8192 ...
# =============================================================================
set -euo pipefail

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

ROOT_DIR="$(cd "$(dirname "$0")" && pwd)"
JOBS="${JOBS:-$(nproc 2>/dev/null || printf '1')}"

echo "Building MPI+OpenMP version with $CXX ($($CXX --version 2>&1 | head -1))"
make -C "$ROOT_DIR/cpu" clean
make -C "$ROOT_DIR/cpu" CXX="$CXX" -j"$JOBS"

echo "Done: $ROOT_DIR/cpu/bin/mbs_po_mpi"
echo ""
echo "Usage:"
echo "  mpirun -np 4 cpu/bin/mbs_po_mpi --method po --sobol 4096 \\"
echo "    --particle 1 100 100 --wavelength-um 0.532 \\"
echo "    --refractive-index 1.31 0 --max-reflections 14 \\"
echo "    --auto-theta-grid 0.02 --auto-phi --threads 8 \\"
echo "    --output results/mpi_sobol --close"
