#!/usr/bin/env bash
set -euo pipefail
task_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
task_output="${1:?Pass an empty output directory for saved regression data}"
make -C "$task_root/cpu" -j2
mkdir -p "$task_root/tests/.build"
mapfile -t task_objects < <(rg --files --no-ignore "$task_root/cpu/build/obj" -g '*.o' | rg -v '/(main|FullAuto)\.o$')
for task_probe in scheduler queue; do
    "${MPI_CXX:-/usr/bin/mpicxx}" -std=c++11 -O2 -fopenmp -DUSE_MPI \
        -I"$task_root/src" -I"$task_root/src/cuda" \
        "$task_root/tests/fullauto_${task_probe}_probe.cpp" "${task_objects[@]}" \
        -lm -lgomp -o "$task_root/tests/.build/fullauto_${task_probe}_probe"
done
python3 "$task_root/tests/test_fullauto_scheduler.py" \
    --scheduler-probe "$task_root/tests/.build/fullauto_scheduler_probe" \
    --queue-probe "$task_root/tests/.build/fullauto_queue_probe" \
    --worker "$task_root/cpu/bin/mbs_po_mpi" --output "$task_output"
