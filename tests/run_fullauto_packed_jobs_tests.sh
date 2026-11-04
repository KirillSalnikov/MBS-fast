#!/usr/bin/env bash
set -euo pipefail
task_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
task_output="${1:?Pass a fresh output directory}"
mkdir -p "$task_output/bin"
task_output="$(cd "$task_output" && pwd)"
"${CXX:-g++}" -std=c++11 -O2 -pthread "$task_root/tests/test_packed_job_store.cpp" -o "$task_output/bin/store_test"
"$task_output/bin/store_test" | tee "$task_output/store_test.log"
"${CXX:-g++}" -std=c++11 -O2 -pthread -I"$task_root/src" "$task_root/tests/test_fullauto_trace_evidence.cpp" -o "$task_output/bin/trace_evidence_test"
"$task_output/bin/trace_evidence_test"
"${CXX:-g++}" -std=c++11 -O2 -pthread "$task_root/tools/fullauto_pack_jobs.cpp" -o "$task_output/bin/fullauto_pack_jobs"
make -C "$task_root/cpu" -j2
cp "$task_root/cpu/bin/mbs_po_mpi" "$task_output/bin/mbs_po"
mapfile -t task_objects < <(rg --files --no-ignore "$task_root/cpu/build/obj" -g '*.o' | rg -v '/(main|FullAuto)\.o$')
"${MPI_CXX:-/usr/bin/mpicxx}" -std=c++11 -O2 -fopenmp -DUSE_MPI -I"$task_root/src" -I"$task_root/src/cuda" \
    "$task_root/tests/fullauto_packed_jobs_probe.cpp" "${task_objects[@]}" -lm -lgomp -o "$task_output/bin/probe"
"$task_output/bin/probe" "$task_output/real_po" "$task_output/bin/mbs_po" > "$task_output/real_po.log"
python3 "$task_root/tests/test_fullauto_packed_resume.py" --binary "$task_output/bin/mbs_po" \
    --tool "$task_output/bin/fullauto_pack_jobs" --output "$task_output/migrated_resume"
"$task_output/bin/probe" "$task_output/timing" "$task_output/bin/mbs_po" timing > "$task_output/timing.log"
rg 'PASS|TIMING' "$task_output/real_po.log" "$task_output/timing.log"
