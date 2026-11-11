#!/usr/bin/env bash
set -euo pipefail
task_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
task_output="${1:?Pass an output directory}"
mkdir -p "$task_output/bin"
task_output="$(cd "$task_output" && pwd)"
make -C "$task_root/cpu" -j2
mapfile -t task_objects < <(rg --files --no-ignore "$task_root/cpu/build/obj" -g '*.o' | rg -v '/(main|FullAuto)\.o$')
"${MPI_CXX:-/usr/bin/mpicxx}" -std=c++11 -O1 -fopenmp -DUSE_MPI -I"$task_root/src" -I"$task_root/src/cuda" \
    "$task_root/tests/fullauto_integral_scope_probe.cpp" "${task_objects[@]}" -lm -lgomp -o "$task_output/bin/integral_scope_probe"
"$task_output/bin/integral_scope_probe" "$task_output/fixture"
