#!/bin/bash
# Reference test: reproduce Test_Ii_PO results
# Hex column L=199.5, D=L/8=24.9375, lambda=0.532, m=1.31
#
# Two runs with different n and N_phi:
#   Forward (0-25°):  n=4,  phi=360
#   Backward (160-180°): n=11, phi=180

MBS=${MBS:-./bin/mbs_po}
OUTDIR=${OUTDIR:-results_ref}
THREADS=${OMP_NUM_THREADS:-12}
mkdir -p $OUTDIR

echo "=== Reference test: Column L=199.5, D=24.9375 (D/L=1/8) ==="
echo "Threads: $THREADS"
echo ""

# Forward: n=4, phi=360
echo "--- Forward (0-25°): n=4, phi=360, beam_cutoff=0.01 ---"
OMP_NUM_THREADS=$THREADS $MBS --po --random 321 215 \
    -p 1 199.5 24.9375 -w 0.532 --ri 1.31 0 -n 4 \
    --tgrid tests/reference_test/tgrid_fwd.txt \
    --beam_cutoff 0.01 --close -o $OUTDIR/fwd

echo ""

# Backward: n=11, phi=180
echo "--- Backward (160-180°): n=11, phi=180, beam_cutoff=0.01 ---"
OMP_NUM_THREADS=$THREADS $MBS --po --random 321 215 \
    -p 1 199.5 24.9375 -w 0.532 --ri 1.31 0 -n 11 \
    --grid 0 180 180 1 --tgrid tests/reference_test/tgrid_bwd.txt \
    --beam_cutoff 0.01 --close -o $OUTDIR/bwd

echo ""
echo "=== Done ==="
echo "Forward: $OUTDIR/fwd.dat"
echo "Backward: $OUTDIR/bwd.dat"
