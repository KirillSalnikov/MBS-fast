#!/usr/bin/env bash
set -euo pipefail
task_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
task_output="${1:?Pass an empty output directory}"
task_binary="${2:?Pass a native CPU binary built from this tree}"
mkdir -p "$task_output"
task_output="$(cd "$task_output" && pwd)"
mapfile -t task_objects < <(rg --files --no-ignore "$task_root/cpu/build/obj" -g '*.o' | rg -v '/(main|FullAuto)\.o$')
for task_probe in staged_upgrade staged_cost cutoff_subset; do
    "${MPI_CXX:-/usr/bin/mpicxx}" -std=c++11 -O1 -fopenmp -DUSE_MPI -I"$task_root/src" -I"$task_root/src/cuda" \
        "$task_root/tests/fullauto_${task_probe}_probe.cpp" "${task_objects[@]}" -lm -lgomp -o "$task_output/${task_probe}_probe"
    task_args=("$task_output/${task_probe}")
    if [[ "$task_probe" != staged_cost ]]; then task_args+=("$task_binary"); fi
    "$task_output/${task_probe}_probe" "${task_args[@]}" > "$task_output/${task_probe}.log" 2>&1
    tail -n 1 "$task_output/${task_probe}.log"
done
