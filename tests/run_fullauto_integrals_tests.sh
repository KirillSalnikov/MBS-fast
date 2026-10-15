#!/usr/bin/env bash
set -euo pipefail
task_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
task_output="${1:?Pass an output directory for saved native reference data}"
make -C "$task_root/cpu" -j2
mkdir -p "$task_output/bin"
task_output="$(cd "$task_output" && pwd)"
# Freeze executables: relinking a binary during a recovery test changes its identity.
cp "$task_root/cpu/bin/mbs_po_mpi" "$task_output/bin/mbs_po"
mapfile -t task_objects < <(rg --files --no-ignore "$task_root/cpu/build/obj" -g '*.o' | rg -v '/(main|FullAuto)\.o$')
"${MPI_CXX:-/usr/bin/mpicxx}" -std=c++11 -O2 -fopenmp -DUSE_MPI \
    -I"$task_root/src" -I"$task_root/src/cuda" \
    "$task_root/tests/fullauto_integrals_probe.cpp" "${task_objects[@]}" \
    -lm -lgomp -o "$task_output/bin/integrals_probe"
task_legacy=()
if [[ -n "${MBS_LEGACY_BINARY:-}" ]]; then task_legacy=(--legacy-binary "$MBS_LEGACY_BINARY"); fi
python3 "$task_root/tests/test_fullauto_integrals.py" "${task_legacy[@]}" \
    --binary "$task_output/bin/mbs_po" --probe "$task_output/bin/integrals_probe" \
    --output "$task_output/native"
python3 "$task_root/tests/test_fullauto_integral_recovery.py" \
    --binary "$task_output/bin/mbs_po" --output "$task_output/controller"
