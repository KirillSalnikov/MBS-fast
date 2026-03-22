#!/bin/bash
# Test multi-size beam caching: compare timing of N separate runs vs 1 multi-size run
set -e
cd /home/serg/MBS-raw

BIN=bin/mbs_po_multi
ORIENT=orientations_256.dat
COMMON="-p 1 10 10 --ri 1.3116 0 -n 5 --po --all -w 0.532 --grid 0 180 1 180 --shadow_off --close --incoh --orientfile $ORIENT"

echo "========================================="
echo "MBS-raw Multi-Size Beam Caching Test"
echo "========================================="
echo ""
echo "Particle: hex column h=10 d=10 (x_ref=83.5)"
echo "Wavelength: 0.532 um"
echo "Orientations: 256 (quasi-random)"
echo "Grid: 0-180 deg, nAz=1, nZen=180"
echo ""

rm -rf M_sep_* M_multi*

# --- Test: 4 sizes ---
echo "=== 4 separate single-size runs ==="
START=$(date +%s%N)
for X in 1 2 3 4; do
    ./$BIN $COMMON -o "M_sep_${X}" 2>/dev/null > /dev/null
done
END=$(date +%s%N)
T_SEP4=$(( (END - START) / 1000000 ))
echo "Time: ${T_SEP4} ms"

echo ""
echo "=== 1 multi-size run (4 sizes) ==="
START=$(date +%s%N)
./$BIN $COMMON --sizes 10 20 30 50 -o "M_multi4" 2>/dev/null
END=$(date +%s%N)
T_MULTI4=$(( (END - START) / 1000000 ))
echo "Time: ${T_MULTI4} ms"

echo ""
echo "--- 8 separate single-size runs ---"
START=$(date +%s%N)
for X in 1 2 3 4 5 6 7 8; do
    ./$BIN $COMMON -o "M_sep8_${X}" 2>/dev/null > /dev/null
done
END=$(date +%s%N)
T_SEP8=$(( (END - START) / 1000000 ))
echo "Time: ${T_SEP8} ms"

echo ""
echo "=== 1 multi-size run (8 sizes) ==="
START=$(date +%s%N)
./$BIN $COMMON --sizes 5 10 15 20 30 40 50 70 -o "M_multi8" 2>/dev/null
END=$(date +%s%N)
T_MULTI8=$(( (END - START) / 1000000 ))
echo "Time: ${T_MULTI8} ms"

echo ""
echo "========================================="
echo "TIMING SUMMARY"
echo "========================================="
echo " 4 separate runs:  ${T_SEP4} ms"
echo " 1 multi (4 sizes):${T_MULTI4} ms  (speedup: $(echo "scale=2; $T_SEP4 / $T_MULTI4" | bc)x)"
echo ""
echo " 8 separate runs:  ${T_SEP8} ms"
echo " 1 multi (8 sizes):${T_MULTI8} ms  (speedup: $(echo "scale=2; $T_SEP8 / $T_MULTI8" | bc)x)"
echo ""
echo "Per-size marginal cost (multi vs separate):"
echo "  Separate: ~$((T_SEP4 / 4)) ms/size"
echo "  Multi-4:  Phase1=~54ms, Phase2=$(( (T_MULTI4 - 54) )) ms for 4 sizes = $((  (T_MULTI4 - 54) / 4 )) ms/size"
echo "  Multi-8:  Phase1=~54ms, Phase2=$(( (T_MULTI8 - 54) )) ms for 8 sizes = $((  (T_MULTI8 - 54) / 8 )) ms/size"
echo "========================================="

rm -rf M_sep_* M_sep8_*
