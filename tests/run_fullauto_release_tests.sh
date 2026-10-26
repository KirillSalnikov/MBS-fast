#!/usr/bin/env bash
# Integration checks for the combined fullauto/GO/extinction release surface.
set -euo pipefail
task_root=$(cd "$(dirname "$0")/.." && pwd)
task_output=${1:?Pass a new output directory}
mkdir -p "$task_output/bin"
task_output=$(cd "$task_output" && pwd)
make -C "$task_root/cpu" -j2
task_binary="$task_root/cpu/bin/mbs_po_mpi"
mapfile -t task_objects < <(find "$task_root/cpu/build/obj" -name '*.o' ! -name main.o ! -name FullAuto.o | sort)
for task_probe in fullauto_residual_schedule_probe fullauto_round_latency_probe fullauto_residual_budget_probe fullauto_execution_progress_probe; do
    "${MPI_CXX:-/usr/bin/mpicxx}" -std=c++11 -O2 -fopenmp -DUSE_MPI \
        -I"$task_root/src" -I"$task_root/src/cuda" \
        "$task_root/tests/$task_probe.cpp" "${task_objects[@]}" -lm -lgomp -o "$task_output/bin/$task_probe"
    case "$task_probe" in
        fullauto_residual_budget_probe) "$task_output/bin/$task_probe" "$task_output/$task_probe" "$task_binary" ;;
        fullauto_execution_progress_probe)
            "$task_output/bin/$task_probe" "$task_output/${task_probe}_persistent" "$task_binary" persistent
            "$task_output/bin/$task_probe" "$task_output/${task_probe}_process" "$task_binary" process ;;
        *) "$task_output/bin/$task_probe" "$task_output/$task_probe" ;;
    esac
done
python3 "$task_root/tests/test_go_mirror.py" --binary "$task_binary" --output "$task_output/go_mirror"
python3 "$task_root/tests/test_go_balanced.py" --binary "$task_binary" --output "$task_output/go_balanced"
python3 "$task_root/tests/regression_extinction_only.py" --binary "$task_binary"
