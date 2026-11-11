#!/usr/bin/env bash
set -euo pipefail
task_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
task_output="${1:?Pass a new test output directory}"
test ! -e "$task_output" || { echo "Use a new output directory" >&2; exit 1; }
mkdir -p "$task_output/bin"
task_output="$(cd "$task_output" && pwd)"
for task_case in test_worker_memory_transport test_worker_memory_stream test_worker_prepared_identity; do
    "${CXX:-g++}" -std=c++11 -O2 -Wall -Wextra -pedantic -pthread \
        "$task_root/tests/${task_case}.cpp" -ldl -o "$task_output/bin/$task_case"
    "$task_output/bin/$task_case" > "$task_output/${task_case}.log" 2>&1
    cat "$task_output/${task_case}.log"
done
if [[ $# -ge 2 ]]; then
    python3 "$task_root/tests/test_worker_memory_native.py" --binary "$2" --output "$task_output/native"
fi
