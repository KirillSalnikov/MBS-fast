#!/usr/bin/env bash
set -euo pipefail

ROOT=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
MBS=${MBS:-"$ROOT/cpu/bin/mbs_po_mpi"}
THREADS=${MBS_EXTINCTION_THREADS:-4}
WORK=$(mktemp -d "${TMPDIR:-/tmp}/mbs-extinction-reference.XXXXXX")
trap 'rm -rf "$WORK"' EXIT

if [[ ! -x "$MBS" ]]; then
    echo "ERROR: solver binary is not executable: $MBS" >&2
    echo "  Fix: run 'make -C cpu -j2' or set MBS=/path/to/mbs_po." >&2
    exit 1
fi

field() {
    local file=$1
    local name=$2
    awk -F '\t' -v wanted="$name" '
        NR == 1 { for (i = 1; i <= NF; ++i) column[$i] = i; next }
        NR == 2 { print $column[wanted] }
    ' "$file"
}

check_relative() {
    local label=$1
    local actual=$2
    local expected=$3
    local tolerance=$4
    awk -v label="$label" -v actual="$actual" -v expected="$expected" -v tol="$tolerance" '
        BEGIN {
            denominator = expected < 0 ? -expected : expected
            if (denominator < 1e-30) denominator = 1
            difference = actual - expected
            if (difference < 0) difference = -difference
            relative = difference / denominator
            printf("  %-20s actual=%.10g reference=%.10g rel=%.6g\n",
                   label, actual, expected, relative)
            if (relative > tol) exit 1
        }
    '
}

divide() {
    local numerator=$1
    local denominator=$2
    awk -v numerator="$numerator" -v denominator="$denominator" '
        BEGIN {
            if (denominator == 0) exit 1
            printf("%.17g", numerator / denominator)
        }
    '
}

run_case() {
    local name=$1
    local expected_area=$2
    local expected_cext=$3
    local expected_qext=$4
    local tolerance=$5
    shift 5
    local output="$WORK/$name"

    "$MBS" --method po --backend cpu --particle "$@" \
        --refractive-index 1.3116 0 --wavelength-um 0.532 \
        --max-reflections 12 --euler-grid 64 32 \
        --scattering-grid 0 180 1 1 --threads "$THREADS" \
        --output "$output" --close >"$WORK/$name.stdout" 2>"$WORK/$name.stderr"

    local report="$output/${name}_integrals.tsv"
    if [[ ! -s "$report" ]]; then
        echo "ERROR: integral report was not written: $report" >&2
        exit 1
    fi

    echo "$name"
    local actual_cext actual_qext actual_area
    actual_cext=$(field "$report" Cext)
    actual_qext=$(field "$report" Qext)
    actual_area=$(divide "$actual_cext" "$actual_qext")
    check_relative "$name Aproj=Cext/Qext" "$actual_area" "$expected_area" "$tolerance"
    check_relative "$name Cext" "$actual_cext" "$expected_cext" "$tolerance"
    check_relative "$name Qext" "$actual_qext" "$expected_qext" "$tolerance"
}

while IFS=$'\t' read -r name type p1 p2 expected_area expected_cext expected_qext tolerance; do
    [[ -z "$name" || "$name" == \#* || "$name" == "case" ]] && continue
    particle_args=("$type" "$p1")
    if [[ "$p2" != "-" ]]; then
        particle_args+=("$p2")
    fi
    run_case "$name" "$expected_area" "$expected_cext" \
        "$expected_qext" "$tolerance" "${particle_args[@]}"
done <"$ROOT/tests/extinction_reference.tsv"

echo "Extinction reference checks passed."
