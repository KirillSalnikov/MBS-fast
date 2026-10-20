#!/usr/bin/env bash
set -euo pipefail
task_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
task_output="${1:?Pass a fresh output directory}"
mkdir -p "$task_output"
task_output="$(cd "$task_output" && pwd)"
"${CXX:-/usr/bin/g++}" -std=c++11 -O2 "$task_root/tests/test_checked_file_cache.cpp" -o "$task_output/checked_file_cache"
"$task_output/checked_file_cache"
bash "$task_root/tests/run_fullauto_integrals_tests.sh" "$task_output/integrals"
mapfile -t task_objects < <(rg --files --no-ignore "$task_root/cpu/build/obj" -g '*.o' | rg -v '/(main|FullAuto)\.o$')
"${MPI_CXX:-/usr/bin/mpicxx}" -std=c++11 -O2 -fopenmp -DUSE_MPI \
    -I"$task_root/src" -I"$task_root/src/cuda" \
    "$task_root/tests/fullauto_parallel_validation_probe.cpp" "${task_objects[@]}" \
    -lm -lgomp -o "$task_output/graph_probe"
"$task_output/graph_probe" "$task_output/graph" "$task_output/integrals/bin/mbs_po"
