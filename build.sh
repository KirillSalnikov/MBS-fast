#!/bin/bash
set -e

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
JOBS="${JOBS:-$(nproc 2>/dev/null || echo 1)}"

cd "$SCRIPT_DIR"

echo "Building CPU MPI/OpenMP split binary..."
make -C cpu -j"$JOBS"

mkdir -p bin
cp cpu/bin/mbs_po_mpi bin/mbs_po

echo "Done: cpu/bin/mbs_po_mpi"
echo "Compatibility copy: bin/mbs_po"
