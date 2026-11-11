#!/usr/bin/env bash
set -euo pipefail
r="$(cd "$(dirname "$0")/.." && pwd)"
out="${1:?new output directory required}"
mkdir -p "$out"
out="$(cd "$out" && pwd)"
test ! -e "$out/fixtures"
make -C "$r/cpu" -j2
mapfile -t objects < <(rg --files --no-ignore "$r/cpu/build/obj" -g '*.o' | rg -v '/(main|FullAuto)\.o$')
/usr/bin/mpicxx -std=c++11 -O0 -g0 -fopenmp -ffunction-sections -fdata-sections -DUSE_MPI -I"$r/src" -I"$r/src/cuda" "$r/tests/fullauto_stop_defaults_probe.cpp" "${objects[@]}" -Wl,--gc-sections -lm -lgomp -o "$out/probe"
"$out/probe" "$out/fixtures"
python3 "$r/tests/test_fullauto_stop_defaults_optics.py" --binary "$r/cpu/bin/mbs_po_mpi" --output "$out/optics"
