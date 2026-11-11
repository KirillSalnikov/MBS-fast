#!/usr/bin/env bash
set -euo pipefail
task_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
task_output="${1:?Pass a private output directory}"
mkdir -p "$task_output/bin"
task_output="$(cd "$task_output" && pwd)"
if [[ -n "${MBS_TEST_CPU_OBJDIR:-}" ]]; then
    task_objdir="$MBS_TEST_CPU_OBJDIR"
else
    make -C "$task_root/cpu" -j2 TARGET=bin/mbs_po_angular_cohort
    task_objdir="$task_root/cpu/build/obj"
fi
mapfile -t task_objects < <(rg --files --no-ignore "$task_objdir" -g '*.o' | rg -v '/(main|FullAuto)\.o$')
"${MPI_CXX:-/usr/bin/mpicxx}" -std=c++11 -O1 -fopenmp -DUSE_MPI \
    -I"$task_root/src" -I"$task_root/src/cuda" \
    "$task_root/tests/fullauto_angular_cohort_probe.cpp" "${task_objects[@]}" \
    -lm -lgomp -o "$task_output/bin/angular_cohort_probe"
"$task_output/bin/angular_cohort_probe" "$task_root/tests/data/fullauto_angular_cohort_front.tsv" "$task_output/saved"
