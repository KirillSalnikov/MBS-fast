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
  --refractive-index 1.3116 0
  --wavelength-um 0.532
  --max-reflections 4
  --fixed-orientation 37 19
  --scattering-grid 0 180 96 48
  --backend cpu
  --cutoff-profile safe
  --close
)

"$mbs" "${common[@]}" --threads 1 --output "$work/serial" \
  >"$work/serial.stdout" 2>"$work/serial.stderr"
"$mbs" "${common[@]}" --threads 4 --output "$work/parallel" \
  >"$work/parallel.stdout" 2>"$work/parallel.stderr"

serial="$work/serial/serial.dat"
parallel="$work/parallel/parallel.dat"
awk '
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
    printf "fixed PO threading regression passed: rows=%d max_abs=%.3g\n", rows, max_abs
  }
' "$serial" "$parallel"
