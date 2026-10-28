#!/usr/bin/env bash
set -euo pipefail
task_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
task_output="${1:?Pass an empty output directory}"
mkdir -p "$task_output/bin"
task_output="$(cd "$task_output" && pwd)"
make -C "$task_root/cpu" -j2 TARGET="$task_output/bin/mbs_po"
mapfile -t task_objects < <(rg --files --no-ignore "$task_root/cpu/build/obj" -g '*.o' | rg -v '/(main|FullAuto)\.o$')
"${MPI_CXX:-/usr/bin/mpicxx}" -std=c++11 -O2 -fopenmp -DUSE_MPI -I"$task_root/src" -I"$task_root/src/cuda" \
    "$task_root/tests/fullauto_mirror_fallback_probe.cpp" "${task_objects[@]}" -lm -lgomp -o "$task_output/bin/mirror_fallback_probe"
python3 - "$task_root" "$task_output" <<'PY'
import sys
from pathlib import Path
sys.path.insert(0,str(Path(sys.argv[1])/'tests'))
from test_fullauto_arbitrary import asymmetric_tetrahedron
asymmetric_tetrahedron(Path(sys.argv[2])/'asymmetric.particle')
PY
"$task_output/bin/mirror_fallback_probe" "$task_output/cases" "$task_output/bin/mbs_po" "$task_output/asymmetric.particle"
