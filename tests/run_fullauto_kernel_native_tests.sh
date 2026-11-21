#!/usr/bin/env bash
set -euo pipefail
task_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
task_output="${1:?Pass a new output directory}"
task_binary="${2:?Pass a native CPU worker binary}"
test ! -e "$task_output" || { echo "Use a new output directory" >&2; exit 1; }
mkdir -p "$task_output/bin"
task_output="$(cd "$task_output" && pwd)"
mapfile -t task_objects < <(rg --files --no-ignore "$task_root/cpu/build/obj" -g '*.o' | rg -v '/(main|FullAuto)\.o$')
build_probe() {
  "${MPI_CXX:-/usr/bin/mpicxx}" -std=c++11 -O1 -fopenmp -DUSE_MPI -ffunction-sections -fdata-sections \
  -I"$task_root/src" -I"$task_root/src/cuda" \
  "$task_root/tests/$1.cpp" "${task_objects[@]}" \
  -Wl,--gc-sections -ldl -lm -lgomp -pthread -o "$task_output/bin/$2"
}
build_probe fullauto_kernel_native_probe probe
"$task_output/bin/probe" "$task_output/fixture" "$task_binary" > "$task_output/result.log" 2>&1
cat "$task_output/result.log"
build_probe fullauto_phi_replay_stop_probe phi_probe
"$task_output/bin/phi_probe" "$task_output/phi_fixture" "$task_binary"
