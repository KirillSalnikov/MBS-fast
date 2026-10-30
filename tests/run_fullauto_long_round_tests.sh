#!/usr/bin/env bash
set -euo pipefail
task_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
task_output="${1:?Pass an empty output directory}"
mkdir -p "$task_output/bin"
task_output="$(cd "$task_output" && pwd)"
make -C "$task_root/cpu" -j2
mapfile -t task_objects < <(rg --files --no-ignore "$task_root/cpu/build/obj" -g '*.o' | rg -v '/(main|FullAuto)\.o$')
for task_probe in long_round residual_schedule round_latency coverage; do
    "${MPI_CXX:-/usr/bin/mpicxx}" -std=c++11 -O1 -fopenmp -DUSE_MPI -I"$task_root/src" -I"$task_root/src/cuda" \
        "$task_root/tests/fullauto_${task_probe}_probe.cpp" "${task_objects[@]}" -lm -lgomp -o "$task_output/bin/${task_probe}_probe"
    task_args=("$task_output/$task_probe")
    if [[ "$task_probe" == long_round ]]; then task_args+=("$task_root/tests/data/fullauto_long_round_timings.tsv"); fi
    "$task_output/bin/${task_probe}_probe" "${task_args[@]}" > "$task_output/${task_probe}.log" 2>&1
    cat "$task_output/${task_probe}.log"
done
