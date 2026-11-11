#!/usr/bin/env bash
set -euo pipefail
task_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
task_output="${1:?Pass a scratch output directory}"
mkdir -p "$task_output"
task_output="$(cd "$task_output" && pwd)"
task_run="$(mktemp -d "$task_output/run.XXXXXX")"
task_worker="${2:-$task_root/cpu/bin/mbs_po_mpi}"
[[ -x "$task_worker" ]] || { echo 'Build the CPU worker first, or pass its path as the second argument.' >&2; exit 2; }
cp "$task_worker" "$task_run/native_worker"
"${CXX:-/usr/bin/g++}" -std=c++11 -O2 -pthread "$task_root/tests/test_packed_job_prepared.cpp" -o "$task_run/pack_prepared"
"$task_run/pack_prepared" | tee "$task_run/api.log"
mapfile -t task_objects < <(rg --files --no-ignore "$task_root/cpu/build/obj" -g '*.o' | rg -v '/(main|FullAuto)\.o$')
for task_probe in fullauto_native_pack_probe fullauto_native_pack_optics_probe fullauto_derived_pack_probe fullauto_derived_pack_optics_probe; do
    "${MPI_CXX:-/usr/bin/mpicxx}" -std=c++11 -O1 -fopenmp -DUSE_MPI -DMBS_NATIVE_ORIENTATION_PREFIX_OUTPUT=1 \
        -I"$task_root/src" -I"$task_root/src/cuda" \
        "$task_root/tests/$task_probe.cpp" "${task_objects[@]}" \
        -lm -lgomp -o "$task_run/$task_probe"
done
"$task_run/fullauto_native_pack_probe" "$task_run/publication" | tee "$task_run/publication.log"
"$task_run/fullauto_native_pack_optics_probe" "$task_run/optics" "$task_run/native_worker" | tee "$task_run/optics.log"
"$task_run/fullauto_derived_pack_probe" "$task_run/derived_publication" | tee "$task_run/derived_publication.log"
"$task_run/fullauto_derived_pack_optics_probe" "$task_run/derived_optics" "$task_run/native_worker" | tee "$task_run/derived_optics.log"
printf 'Prepared native and derived pack tests passed: %s\n' "$task_run"
