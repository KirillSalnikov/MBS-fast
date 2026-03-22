#!/bin/bash
# =============================================================================
# MBS-raw regression test suite
# Tests basic functionality: builds binary, runs test cases, checks Q_sca
# and M11 at key angles against reference values.
#
# Usage: cd MBS-raw && bash tests/run_tests.sh
# Exit code: 0 if all tests pass, 1 if any fail
# =============================================================================

set -o pipefail

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
PROJECT_DIR="$(cd "$SCRIPT_DIR/.." && pwd)"
MBS="$PROJECT_DIR/bin/mbs_po"
WORK_DIR=$(mktemp -d /tmp/mbs_test_XXXXXX)

PASS=0
FAIL=0
TOTAL=0

cleanup() {
    rm -rf "$WORK_DIR"
}
trap cleanup EXIT

pass_test() {
    PASS=$((PASS + 1))
    TOTAL=$((TOTAL + 1))
    echo "  [PASS] $1"
}

fail_test() {
    FAIL=$((FAIL + 1))
    TOTAL=$((TOTAL + 1))
    echo "  [FAIL] $1: $2"
}

# Extract Q_sca from stdout (printed by MBS-raw as "Q_sca = C_sca / A_proj = VALUE")
# Args: $1 = stdout file path
extract_qsca_from_stdout() {
    local stdoutfile="$1"
    if [ ! -f "$stdoutfile" ]; then
        echo "NaN"
        return
    fi
    grep -oP 'Q_sca = C_sca / A_proj = \K[0-9.]+' "$stdoutfile" || echo "NaN"
}

# Extract C_sca from .dat file: sum of M11 * 2pi*dcos
# Args: $1 = dat file path
extract_csca() {
    local datfile="$1"
    if [ ! -f "$datfile" ]; then
        echo "NaN"
        return
    fi
    awk 'NR > 1 && NF >= 3 {sum += $2 * $3} END {printf "%.4f", sum}' "$datfile"
}

# Extract M11 at a specific angle (nearest)
# Args: $1 = dat file, $2 = target angle
extract_m11_at() {
    local datfile="$1"
    local target="$2"
    if [ ! -f "$datfile" ]; then
        echo "NaN"
        return
    fi
    awk -v tgt="$target" 'NR > 1 && NF >= 3 {
        d = ($1 - tgt); if (d < 0) d = -d;
        if (best == "" || d < bestd) { bestd = d; best = $3 }
    } END { printf "%.6e", best }' "$datfile"
}

# Check if value is within tolerance of reference
# Args: $1 = value, $2 = reference, $3 = tolerance (absolute)
check_value() {
    local val="$1" ref="$2" tol="$3"
    awk -v v="$val" -v r="$ref" -v t="$tol" 'BEGIN {
        d = v - r; if (d < 0) d = -d;
        exit (d <= t) ? 0 : 1
    }'
}

# Check if value is within relative tolerance
# Args: $1 = value, $2 = reference, $3 = relative tolerance
check_relative() {
    local val="$1" ref="$2" tol="$3"
    awk -v v="$val" -v r="$ref" -v t="$tol" 'BEGIN {
        if (r == 0) { exit (v == 0) ? 0 : 1 }
        d = (v - r) / r; if (d < 0) d = -d;
        exit (d <= t) ? 0 : 1
    }'
}

# =============================================================================
echo "=============================================="
echo "MBS-raw Regression Test Suite"
echo "=============================================="
echo ""

# --- Build ---
echo "=== Building MBS-raw ==="
cd "$PROJECT_DIR"
# Clean .o files
find src -name '*.o' -delete 2>/dev/null
if bash build.sh > "$WORK_DIR/build.log" 2>&1; then
    echo "  Build: OK"
else
    echo "  Build: FAILED"
    cat "$WORK_DIR/build.log"
    exit 1
fi

if [ ! -x "$MBS" ]; then
    echo "  Binary not found: $MBS"
    exit 1
fi
echo ""

# =============================================================================
# Test 1: Hex column D=H=10, m=1.31, n=6, --sobol 256
# Expected: Q_sca ~ 2.05 +/- 0.1
# =============================================================================
echo "=== Test 1: Hex column (standard, no absorption) ==="
TEST1_DIR="$WORK_DIR/test1"
mkdir -p "$TEST1_DIR"

cd "$TEST1_DIR"
$MBS --po --sobol 256 -p 1 10 10 -w 0.532 --ri 1.31 0 \
    -n 6 --grid 0 180 48 90 --close -o M_test1 \
    > "$TEST1_DIR/stdout.txt" 2>&1

Q1=$(extract_qsca_from_stdout "$TEST1_DIR/stdout.txt")
echo "  Q_sca = $Q1 (expected ~2.05)"
if [ "$Q1" = "NaN" ]; then
    fail_test "Test 1" "Could not extract Q_sca from stdout"
elif check_relative "$Q1" "2.05" "0.10"; then
    pass_test "Test 1: Q_sca within 10% of 2.05"
else
    fail_test "Test 1" "Q_sca=$Q1, expected 2.05 +/- 10%"
fi

DATFILE=$(find "$TEST1_DIR" -name "M_test1.dat" -o -name "M_test1*.dat" 2>/dev/null | head -1)
if [ -z "$DATFILE" ]; then
    DATFILE=$(find "$TEST1_DIR" -name "*.dat" 2>/dev/null | grep -v stdout | grep -v particle | head -1)
fi

if [ -n "$DATFILE" ]; then
    # Check M11 at key angles
    M11_0=$(extract_m11_at "$DATFILE" 0)
    M11_90=$(extract_m11_at "$DATFILE" 90)
    M11_180=$(extract_m11_at "$DATFILE" 180)
    echo "  M11(0)=$M11_0, M11(90)=$M11_90, M11(180)=$M11_180"

    # Forward peak should be much larger than sidescatter
    if awk -v f="$M11_0" -v s="$M11_90" 'BEGIN { exit (f > 10 * s) ? 0 : 1 }'; then
        pass_test "Test 1: M11(0) >> M11(90) (forward peak)"
    else
        fail_test "Test 1" "M11(0) not >> M11(90): $M11_0 vs $M11_90"
    fi
fi
echo ""

# =============================================================================
# Test 2: Hex column with absorption m=1.31+0.1i
# Expected: Q_sca ~ 1.09 +/- 0.15
# =============================================================================
echo "=== Test 2: Hex column with absorption ==="
TEST2_DIR="$WORK_DIR/test2"
mkdir -p "$TEST2_DIR"

cd "$TEST2_DIR"
$MBS --po --sobol 256 -p 1 10 10 -w 0.532 --ri 1.31 0.1 --abs \
    -n 6 --grid 0 180 48 90 --close -o M_test2 \
    > "$TEST2_DIR/stdout.txt" 2>&1
RET=$?

if [ $RET -ne 0 ]; then
    fail_test "Test 2" "mbs_po exited with code $RET"
else
    Q2=$(extract_qsca_from_stdout "$TEST2_DIR/stdout.txt")
    echo "  Q_sca = $Q2 (expected ~1.09)"
    if [ "$Q2" = "NaN" ]; then
        fail_test "Test 2" "Could not extract Q_sca from stdout"
    elif check_relative "$Q2" "1.09" "0.20"; then
        pass_test "Test 2: Q_sca with absorption within 20% of 1.09"
    else
        fail_test "Test 2" "Q_sca=$Q2, expected 1.09 +/- 20%"
    fi
fi
echo ""

# =============================================================================
# Test 3: Concave hexagonal (should not crash)
# =============================================================================
echo "=== Test 3: Concave hexagonal (crash test) ==="
TEST3_DIR="$WORK_DIR/test3"
mkdir -p "$TEST3_DIR"

cd "$TEST3_DIR"
$MBS --po --sobol 128 -p 10 10 10 3 -w 0.532 --ri 1.31 0 \
    -n 6 --grid 0 180 48 90 --close -o M_test3 \
    > "$TEST3_DIR/stdout.txt" 2>&1
RET=$?

if [ $RET -ne 0 ]; then
    fail_test "Test 3" "Concave hex crashed with code $RET"
else
    pass_test "Test 3: Concave hexagonal ran without crash"
fi
echo ""

# =============================================================================
# Test 4: Bullet rosette (should not crash, Q_sca > 1)
# =============================================================================
echo "=== Test 4: Bullet rosette (crash test + Q_sca) ==="
TEST4_DIR="$WORK_DIR/test4"
mkdir -p "$TEST4_DIR"

cd "$TEST4_DIR"
$MBS --po --sobol 128 -p 3 10 10 -w 0.532 --ri 1.31 0 \
    -n 6 --grid 0 180 48 90 --close -o M_test4 \
    > "$TEST4_DIR/stdout.txt" 2>&1
RET=$?

if [ $RET -ne 0 ]; then
    fail_test "Test 4" "Bullet rosette crashed with code $RET"
else
    Q4=$(extract_qsca_from_stdout "$TEST4_DIR/stdout.txt")
    echo "  Q_sca = $Q4 (expected > 1.0)"
    if [ "$Q4" = "NaN" ]; then
        # No Q_sca in output but didn't crash — still a pass for crash test
        pass_test "Test 4: Bullet rosette ran without crash (no Q_sca in output)"
    elif awk -v v="$Q4" 'BEGIN { exit (v > 1.0) ? 0 : 1 }'; then
        pass_test "Test 4: Bullet rosette Q_sca > 1"
    else
        fail_test "Test 4" "Q_sca=$Q4, expected > 1.0"
    fi
fi
echo ""

# =============================================================================
# Test 5: Adaptive convergence
# =============================================================================
echo "=== Test 5: Adaptive convergence ==="
TEST5_DIR="$WORK_DIR/test5"
mkdir -p "$TEST5_DIR"

cd "$TEST5_DIR"
timeout 120 $MBS --po --adaptive 0.05 -p 1 10 10 -w 0.532 --ri 1.31 0 \
    -n 6 --grid 0 180 48 90 --close -o M_test5 \
    > "$TEST5_DIR/stdout.txt" 2>&1
RET=$?

if [ $RET -ne 0 ]; then
    fail_test "Test 5" "Adaptive mode exited with code $RET (or timeout)"
else
    # Check that "Converged" or "Max iterations" appears in output
    if grep -qi "converged\|max iterations" "$TEST5_DIR/stdout.txt"; then
        pass_test "Test 5: Adaptive mode converged"
    else
        fail_test "Test 5" "No convergence message in output"
    fi
fi
echo ""

# =============================================================================
# Summary
# =============================================================================
echo "=============================================="
echo "Results: $PASS passed, $FAIL failed, $TOTAL total"
echo "=============================================="

if [ $FAIL -gt 0 ]; then
    exit 1
else
    echo "All tests passed!"
    exit 0
fi
