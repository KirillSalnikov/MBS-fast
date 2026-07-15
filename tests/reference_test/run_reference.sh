#!/bin/bash
# Reference test: reproduce Test_Ii_PO results
# Hex column L=199.5, D=L/8=24.9375, lambda=0.532, m=1.31
#
# Two runs with different n and N_phi:
#   Forward (0-25°):  n=4,  phi=360
#   Backward (160-180°): n=11, phi=180

set -euo pipefail

MBS=${MBS:-./cpu/bin/mbs_po_mpi}
OUTDIR=${OUTDIR:-results_ref}
THREADS=${OMP_NUM_THREADS:-12}
mkdir -p "$OUTDIR"

echo "=== Reference test: Column L=199.5, D=24.9375 (D/L=1/8) ==="
echo "Threads: $THREADS"
echo ""

# Forward: n=4, phi=360
echo "--- Forward (0-25°): n=4, phi=360, beam_cutoff=0.01 ---"
OMP_NUM_THREADS="$THREADS" "$MBS" --method po --euler-grid 321 215 \
    --particle 1 199.5 24.9375 --wavelength-um 0.532 \
    --refractive-index 1.31 0 --max-reflections 4 \
    --theta-grid-file tests/reference_test/tgrid_fwd.txt --phi-points 360 \
    --beam-cutoff 0.01 --threads "$THREADS" --close --output "$OUTDIR/fwd"

echo ""

# Backward: n=11, phi=180
echo "--- Backward (160-180°): n=11, phi=180, beam_cutoff=0.01 ---"
OMP_NUM_THREADS="$THREADS" "$MBS" --method po --euler-grid 321 215 \
    --particle 1 199.5 24.9375 --wavelength-um 0.532 \
    --refractive-index 1.31 0 --max-reflections 11 \
    --theta-grid-file tests/reference_test/tgrid_bwd.txt --phi-points 180 \
    --beam-cutoff 0.01 --threads "$THREADS" --close --output "$OUTDIR/bwd"

echo ""
echo "=== Done ==="
echo "Forward directory: $OUTDIR/fwd"
echo "Backward directory: $OUTDIR/bwd"

# === MPI VERSION ===
# Uncomment below for multi-node:
# MBS=./cpu/bin/mbs_po_mpi
# NRANKS=8
#
# mpirun -np $NRANKS --bind-to none \
#     -x OMP_NUM_THREADS=$THREADS -x OMP_PROC_BIND=false \
#     $MBS --method po --euler-grid 321 215 \
#     --particle 1 199.5 24.9375 --wavelength-um 0.532 \
#     --refractive-index 1.31 0 --max-reflections 4 \
#     --theta-grid-file tests/reference_test/tgrid_fwd.txt --phi-points 360 \
#     --beam-cutoff 0.01 --close --output $OUTDIR/fwd_mpi
#
# mpirun -np $NRANKS --bind-to none \
#     -x OMP_NUM_THREADS=$THREADS -x OMP_PROC_BIND=false \
#     $MBS --method po --euler-grid 321 215 \
#     --particle 1 199.5 24.9375 --wavelength-um 0.532 \
#     --refractive-index 1.31 0 --max-reflections 11 \
#     --theta-grid-file tests/reference_test/tgrid_bwd.txt --phi-points 180 \
#     --beam-cutoff 0.01 --close --output $OUTDIR/bwd_mpi
