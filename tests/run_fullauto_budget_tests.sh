#!/usr/bin/env bash
set -euo pipefail
task_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
task_output="${1:?Pass a directory for saved regression data}"
make -C "$task_root/cpu" -j2
mkdir -p "$task_root/tests/.build"
mapfile -t task_objects < <(rg --files --no-ignore "$task_root/cpu/build/obj" -g '*.o' | rg -v '/(main|FullAuto)\.o$')
"${MPI_CXX:-/usr/bin/mpicxx}" -std=c++11 -O2 -fopenmp -DUSE_MPI \
    -I"$task_root/src" -I"$task_root/src/cuda" \
    "$task_root/tests/fullauto_budget_probe.cpp" "${task_objects[@]}" \
    -lm -lgomp -o "$task_root/tests/.build/fullauto_budget_probe"
"$task_root/tests/.build/fullauto_budget_probe" \
    "$task_root/tests/data/fullauto_pointwise_budget.txt" "$task_output"
