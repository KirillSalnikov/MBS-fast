#!/usr/bin/env bash
set -euo pipefail

ROOT=$(cd "$(dirname "$0")/.." && pwd)
MBS=${MBS:-"$ROOT/cpu/bin/mbs_po_mpi"}
if [[ "$MBS" != /* ]]; then
    MBS="$(cd "$(dirname "$MBS")" && pwd)/$(basename "$MBS")"
fi
WORK=${MBS_ADAPTIVE_TEST_DIR:-"$ROOT/tests/.adaptive"}
rm -rf "$WORK"
mkdir -p "$WORK"

base=(
    "$MBS" --method po --particle 1 1 1
    --refractive-index 1.31 0 --wavelength-um 10
    --backend cpu --threads 2 --close --max-reflections 5
)

"${base[@]}" --sobol 16 --scattering-grid 0 180 12 12 \
    --adaptive-phi 0.3 --convergence-passes 1 --max-phi-points 48 \
    --output "$WORK/phi" >"$WORK/phi.log" 2>&1
grep -q 'Adaptive phi selected N_phi=' "$WORK/phi.log"
grep -q $'^phi\t' "$WORK/phi/phi_convergence.tsv"
awk -F '\t' '$1 == "phi" && $5 != 0.3 { exit 1 }' \
    "$WORK/phi/phi_convergence.tsv"

"${base[@]}" --sobol 16 --scattering-grid 0 180 12 12 \
    --adaptive-reflections 0.25 \
    --adaptive-config "$ROOT/tests/adaptive_quick.conf" \
    --output "$WORK/reflections" >"$WORK/reflections.log" 2>&1
grep -q 'Adaptive reflections selected n=' "$WORK/reflections.log"
grep -q 'Adaptive reflections: target=25' "$WORK/reflections.log"
grep -q 'range=1..5' "$WORK/reflections.log"
grep -q $'^reflections\t' \
    "$WORK/reflections/reflections_convergence.tsv"

"${base[@]}" --sobol 16 --auto-theta-grid 0.3 \
    --max-theta-points 129 \
    --output "$WORK/theta" >"$WORK/theta.log" 2>&1
grep -q 'Adaptive theta selected ' "$WORK/theta.log"
grep -q $'^theta\t' "$WORK/theta/theta_convergence.tsv"

"${base[@]}" --adaptive-euler-grid 0.3 \
    --scattering-grid 0 180 12 12 --convergence-passes 1 \
    --max-beta-points 16 --max-gamma-points 48 \
    --max-orientations 1024 --output "$WORK/euler" \
    >"$WORK/euler.log" 2>&1
grep -q 'Adaptive Euler selected ' "$WORK/euler.log"
grep -q $'^beta_gamma_joint\t' "$WORK/euler/euler_convergence.tsv"

"${base[@]}" --auto 0.3 --convergence-passes 1 \
    --max-phi-points 48 --max-theta-points 129 \
    --max-orientations 256 --output "$WORK/auto" \
    >"$WORK/auto.log" 2>&1
grep -q 'Unified convergence selected:' "$WORK/auto.log"
grep -q $'^joint_sweep\t' "$WORK/auto/auto_convergence.tsv"
grep -q $'^orientations\t' "$WORK/auto/auto_convergence.tsv"

"${base[@]}" --auto 0.01 \
    --adaptive-config "$ROOT/tests/adaptive_quick.conf" \
    --output "$WORK/config" >"$WORK/config.log" 2>&1
grep -q 'Adaptive config:' "$WORK/config.log"
grep -q 'alpha=30' "$WORK/config.log"
awk -F '\t' '$1 == "phi" && $5 != 0.3 { exit 1 }' \
    "$WORK/config/config_convergence.tsv"

"${base[@]}" --autofull 0.01 \
    --adaptive-config "$ROOT/tests/adaptive_quick.conf" \
    --output "$WORK/config_full" >"$WORK/config_full.log" 2>&1
grep -q 'Adaptive reflections: target=30' "$WORK/config_full.log"
grep -q 'range=1..5' "$WORK/config_full.log"
grep -q $'^reflections\t.*\t0.3\t' \
    "$WORK/config_full/config_full_convergence.tsv"

set +e
"${base[@]}" --adaptive-euler-grid 0.01 \
    --scattering-grid 0 180 12 12 --max-beta-points 4 \
    --max-gamma-points 12 --max-orientations 128 \
    --output "$WORK/limit" >"$WORK/limit.log" 2>&1
status=$?
set -e
if [[ $status -eq 0 ]]; then
    echo "adaptive limit test unexpectedly succeeded" >&2
    exit 1
fi
grep -q '^ERROR: calculation failed:' "$WORK/limit.log"
if grep -q 'Process received signal\|Aborted' "$WORK/limit.log"; then
    echo "adaptive limit produced an abort instead of an actionable error" >&2
    exit 1
fi

echo "Adaptive convergence tests passed"
