#!/usr/bin/env bash
set -euo pipefail
task_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
task_output="${1:?Pass a new output directory}"
mkdir -p "$task_output/bin"
task_output="$(cd "$task_output" && pwd)"
test ! -e "$task_output/checkpoint" || { echo "Use a new output directory" >&2; exit 1; }
make -C "$task_root/cpu" -j2
mapfile -t task_objects < <(rg --files --no-ignore "$task_root/cpu/build/obj" -g '*.o' | rg -v '/(main|FullAuto)\.o$')
"${MPI_CXX:-/usr/bin/mpicxx}" -std=c++11 -O1 -fopenmp -DUSE_MPI -I"$task_root/src" -I"$task_root/src/cuda" \
    "$task_root/tests/fullauto_source_marker_probe.cpp" "${task_objects[@]}" -lm -lgomp -o "$task_output/bin/probe"
"$task_output/bin/probe" "$task_output/checkpoint" > "$task_output/result.log" 2>&1
cat "$task_output/result.log"
