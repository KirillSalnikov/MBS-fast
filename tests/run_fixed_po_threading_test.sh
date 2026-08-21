#!/usr/bin/env bash
set -euo pipefail

repo=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
mbs=${MBS_BIN:-$repo/cpu/bin/mbs_po_mpi}
work=$(mktemp -d "${TMPDIR:-/tmp}/mbs-fixed-po-threads.XXXXXX")
trap 'rm -rf -- "$work"' EXIT

${CXX:-g++} -O2 -std=gnu++11 \
  -I"$repo/src" -I"$repo/src/math" -I"$repo/src/common" \
  "$repo/tests/test_mueller_direct.cpp" \
  "$repo/src/math/Mueller.cpp" "$repo/src/math/matrix.cpp" \
  "$repo/src/math/compl.cpp" "$repo/src/math/JonesMatrix.cpp" \
  "$repo/src/common/Matrix4x4.cpp" \
  -o "$work/test_mueller_direct"
"$work/test_mueller_direct"

common=(
  --method po
  --particle 10 8 5.6 43.4
  --wavelength-um 0.532
  --max-reflections 4
  --fixed-orientation 37 19
  --scattering-grid 0 180 96 48
  --backend cpu
  --cutoff-profile safe
  --close
)

run_case() {
  local case_name=$1
  local imaginary_index=$2
  shift 2
  local serial_dir="$work/${case_name}_serial"
  local parallel_dir="$work/${case_name}_parallel"

  "$mbs" "${common[@]}" \
    --refractive-index 1.3116 "$imaginary_index" "$@" \
    --threads 1 --output "$serial_dir" \
    >"$work/${case_name}_serial.stdout" \
    2>"$work/${case_name}_serial.stderr"
  "$mbs" "${common[@]}" \
    --refractive-index 1.3116 "$imaginary_index" "$@" \
    --threads 4 --output "$parallel_dir" \
    >"$work/${case_name}_parallel.stdout" \
    2>"$work/${case_name}_parallel.stderr"

  awk -v case_name="$case_name" '
  NR == FNR {
    for (i = 1; i <= NF; ++i) reference[FNR, i] = $i
    fields[FNR] = NF
    rows = FNR
    next
  }
  {
    if (FNR > rows || NF != fields[FNR]) exit 2
    for (i = 1; i <= NF; ++i) {
      d = $i - reference[FNR, i]
      if (d < 0) d = -d
      if (d > max_abs) max_abs = d
    }
    seen = FNR
  }
  END {
    if (seen != rows || max_abs > 2e-12) exit 1
    printf "fixed PO threading regression (%s) passed: rows=%d max_abs=%.3g\n", case_name, rows, max_abs
  }
' "$serial_dir/${case_name}_serial.dat" \
    "$parallel_dir/${case_name}_parallel.dat"
}

run_case nonabsorbing 0
run_case absorbing 0.01 --absorption
