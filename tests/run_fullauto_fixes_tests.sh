#!/usr/bin/env bash
set -euo pipefail
task_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
task_output="${1:?Pass an empty output directory}"
mkdir -p "$task_output/bin"
task_output="$(cd "$task_output" && pwd)"
make -C "$task_root/cpu" -j2
cp "$task_root/cpu/bin/mbs_po_mpi" "$task_output/bin/mbs_po"
mapfile -t task_objects < <(rg --files --no-ignore "$task_root/cpu/build/obj" -g '*.o' | rg -v '/(main|FullAuto)\.o$')
"${MPI_CXX:-/usr/bin/mpicxx}" -std=c++11 -O2 -fopenmp -DUSE_MPI -I"$task_root/src" -I"$task_root/src/cuda" \
    "$task_root/tests/fullauto_fixes_probe.cpp" "${task_objects[@]}" -lm -lgomp -o "$task_output/bin/fixes_probe"
"$task_output/bin/fixes_probe" "$task_output/decisions"
python3 "$task_root/tests/test_sahara340_statistics.py"
python3 "$task_root/tests/test_fullauto_depth_unresolved.py" --binary "$task_output/bin/mbs_po" --output "$task_output/depth"
python3 "$task_root/tests/test_fullauto_staging_acceptance.py" --binary "$task_output/bin/mbs_po" --output "$task_output/native"
python3 "$task_root/tests/test_incoherent_controls.py" --binary "$task_output/bin/mbs_po" --output "$task_output/incoherent_controls"
