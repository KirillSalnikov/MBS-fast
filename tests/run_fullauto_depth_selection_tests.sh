#!/usr/bin/env bash
# No optical worker is started. The integration probe uses private fixture files
# and exercises actual controller code with native work explicitly forbidden.
set -euo pipefail
task_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
task_output="${1:?Pass an empty output directory}"
if [[ -e "$task_output/probe_fixtures" || -e "$task_output/unit.log" || -e "$task_output/integration.log" ]]; then
    echo 'Choose a fresh output directory; existing regression evidence is preserved.' >&2
    exit 2
fi
mkdir -p "$task_output/bin"
task_output="$(cd "$task_output" && pwd)"
"${CXX:-g++}" -std=c++11 -O2 -Wall -Wextra -pedantic -fsanitize=undefined \
    "$task_root/tests/test_fullauto_depth_selection.cpp" -o "$task_output/bin/depth_selection_unit"
"$task_output/bin/depth_selection_unit" > "$task_output/unit.log" 2>&1
cat "$task_output/unit.log"
make -C "$task_root/cpu" -j2
mapfile -t task_objects < <(rg --files --no-ignore "$task_root/cpu/build/obj" -g '*.o' | rg -v '/(main|FullAuto)\.o$')
"${MPI_CXX:-/usr/bin/mpicxx}" -std=c++11 -O0 -g0 -fopenmp \
    -ffunction-sections -fdata-sections -DUSE_MPI -I"$task_root/src" -I"$task_root/src/cuda" \
    "$task_root/tests/fullauto_depth_selection_probe.cpp" "${task_objects[@]}" \
    -Wl,--gc-sections -lm -lgomp -o "$task_output/bin/depth_selection_probe"
"$task_output/bin/depth_selection_probe" "$task_output/probe_fixtures" > "$task_output/integration.log" 2>&1
cat "$task_output/integration.log"
