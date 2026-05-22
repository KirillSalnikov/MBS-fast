#!/usr/bin/env bash
# Regression for the CUDA float near-singular edge integral bug found on
# Afine30, k_eq=53.47, theta=162.0994475 deg.
#
# This test is intentionally opt-in because it needs a CUDA build and the
# Afine30 particle file. It exits 77 when the prerequisites are missing.
#
# Usage:
#   MBS=./bin/mbs_po_float_fast AFINE30=/path/Afine30.dat bash tests/regression_cuda_afine53_eps.sh

set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
MBS="${MBS:-$ROOT_DIR/bin/mbs_po_float_fast}"
AFINE30="${AFINE30:-}"
KEQ_LIST="${KEQ_LIST:-}"
THREADS="${THREADS:-12}"

skip() {
    echo "SKIP: $*"
    exit 77
}

if [[ ! -x "$MBS" ]]; then
    skip "CUDA float_fast binary is not executable: $MBS"
fi

if [[ -z "$AFINE30" ]]; then
    for candidate in \
        "$ROOT_DIR/Afine30.dat" \
        "/home/user/Desktop/tests_POs/Afine30.dat" \
        "/home/user/cluster/afine30_run_current/Afine30.dat"; do
        if [[ -f "$candidate" ]]; then
            AFINE30="$candidate"
            break
        fi
    done
fi
[[ -f "$AFINE30" ]] || skip "Afine30.dat not found; set AFINE30=/path/Afine30.dat"

WORK_DIR="$(mktemp -d /tmp/mbs_cuda_afine53_XXXXXX)"
trap 'rm -rf "$WORK_DIR"' EXIT

if [[ -z "$KEQ_LIST" ]]; then
    KEQ_LIST="$WORK_DIR/keq53p47.lst"
    printf '53.47\n' > "$KEQ_LIST"
fi

OUT_PREFIX="$WORK_DIR/afine53_eps"
LOG_FILE="$WORK_DIR/run.log"

MBS_DEBUG_ORIENT_BEGIN=35744 \
MBS_DEBUG_ORIENT_END=35760 \
MBS_SHARED_ORIENT_CHUNK=16 \
OMP_NUM_THREADS="$THREADS" \
"$MBS" --pf "$AFINE30" \
    --multikeq_list "$KEQ_LIST" \
    --ri 1.6 0.002 -n 14 --po --oldauto 2 \
    --beam_cutoff_j 0.001 --beam_cutoff_area 0.002 \
    --trace_cutoff_importance 0.0001 --trace_max_beams 20000 \
    --grid 0 180 600 181 -w 1.064 --gpu --close --log 1 --threads "$THREADS" \
    -o "$OUT_PREFIX" > "$LOG_FILE" 2>&1

DAT_FILE="$(find "$WORK_DIR" -maxdepth 3 -name 'afine53_eps*.dat' | head -1)"
[[ -f "$DAT_FILE" ]] || {
    cat "$LOG_FILE"
    echo "FAIL: output .dat file was not created"
    exit 1
}

awk '
NR > 1 {
    d = $1 - 162.0994475
    if (d < 0) d = -d
    if (best == "" || d < bestd) {
        best = $0
        bestd = d
    }
}
END {
    if (best == "") {
        print "FAIL: no theta rows found"
        exit 1
    }
    split(best, a)
    theta = a[1]
    m11 = a[3]
    r23 = a[9] / m11
    r24 = a[10] / m11
    r32 = a[12] / m11
    r42 = a[16] / m11

    printf("theta=%.10f M11=%.12g R23=%.8g R24=%.8g R32=%.8g R42=%.8g\n",
           theta, m11, r23, r24, r32, r42)

    # Broken CUDA float builds produced M11 ~= 0.411 here. The fixed float
    # path agrees with CPU/double near M11 ~= 0.003235 and these ratios.
    ok = 1
    ok = ok && (m11 > 0.0030 && m11 < 0.0035)
    ok = ok && (r23 > 0.0090 && r23 < 0.0120)
    ok = ok && (r24 > 0.0650 && r24 < 0.0760)
    ok = ok && (r32 > -0.0570 && r32 < -0.0470)
    ok = ok && (r42 > 0.1080 && r42 < 0.1230)
    if (!ok) {
        print "FAIL: CUDA Afine30 eps regression is outside tolerance"
        exit 1
    }
    print "PASS: CUDA Afine30 eps regression"
}
' "$DAT_FILE"
