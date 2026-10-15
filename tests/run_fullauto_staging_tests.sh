#!/usr/bin/env bash
set -euo pipefail
task_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
task_output="${1:?Pass an empty output directory}"
make -C "$task_root/cpu" -j2
mkdir -p "$task_output/bin"
task_output="$(cd "$task_output" && pwd)"
cp "$task_root/cpu/bin/mbs_po_mpi" "$task_output/bin/mbs_po"
mapfile -t task_objects < <(rg --files --no-ignore "$task_root/cpu/build/obj" -g '*.o' | rg -v '/(main|FullAuto)\.o$')
"${MPI_CXX:-/usr/bin/mpicxx}" -std=c++11 -O2 -fopenmp -DUSE_MPI -I"$task_root/src" -I"$task_root/src/cuda" \
    "$task_root/tests/fullauto_staging_probe.cpp" "${task_objects[@]}" -lm -lgomp -o "$task_output/bin/staging_probe"
"$task_output/bin/staging_probe" --decisions "$task_root/tests/data/fullauto_priority_backscatter.txt" > "$task_output/decisions.csv"
"$task_output/bin/staging_probe" --phi-decisions "$task_output/phi_decisions" > "$task_output/phi_decisions.log"
TASK_MBS_WORKER="$task_output/bin/mbs_po" "$task_output/bin/staging_probe" \
    --method po --backend cpu --particle 1 .2 .1 --ri 1.3116 0 --wavelength-um 1 \
    --max-reflections 4 --threads 1 --fullauto .02 --fullauto-initial 256 --fullauto-min-correction 64 \
    --phi-points 1 --fullauto-mirror off --fullauto-analytic-controls off --fullauto-resources off \
    --cutoff-profile off --symmetry 1 1 --scattering-grid 170 180 1 8 \
    --fullauto-orientation-block 64 --output "$task_output/native" --close > "$task_output/native.log" 2>&1
python3 "$task_root/tests/test_fullauto_staging_acceptance.py" \
    --binary "$task_output/bin/mbs_po" --output "$task_output/acceptance"
