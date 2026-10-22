#!/usr/bin/env bash
set -euo pipefail
task_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
task_output="${1:?Pass a directory for local regression output}"
mkdir -p "$task_output/bin"
task_output="$(cd "$task_output" && pwd)"
make -C "$task_root/cpu" -j2
mapfile -t task_objects < <(rg --files --no-ignore "$task_root/cpu/build/obj" -g '*.o' | rg -v '/(main|FullAuto)\.o$')
"${MPI_CXX:-/usr/bin/mpicxx}" -std=c++11 -O2 -fopenmp -DUSE_MPI \
    -I"$task_root/src" -I"$task_root/src/cuda" \
    "$task_root/tests/fullauto_level_stability_probe.cpp" "${task_objects[@]}" \
    -lm -lgomp -o "$task_output/bin/fullauto_level_stability_probe"
"$task_output/bin/fullauto_level_stability_probe" \
    "$task_root/tests/data/fullauto_level_stability.tsv" "$task_output/policy" | tee "$task_output/level_stability.csv"
