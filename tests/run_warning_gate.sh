#!/usr/bin/env bash
set -euo pipefail

root=$(cd "$(dirname "$0")/.." && pwd)
log=$(mktemp /tmp/mbs-warning-gate.XXXXXX.log)
trap 'rm -f "$log"' EXIT

make -B -C "$root/cpu" \
    TARGET=bin/mbs_po_mpi_warnings \
    OBJDIR=build/warnings/obj \
    OPT_FLAGS='-O2 -Wall -Wextra -Wpedantic' \
    -j"${MBS_BUILD_JOBS:-4}" >"$log" 2>&1

if grep -Eq '^\.\./src/.*warning:' "$log"; then
    echo 'ERROR: project compiler warnings found:' >&2
    grep -E '^\.\./src/.*warning:' "$log" >&2
    exit 1
fi

external=$(grep -c 'warning:' "$log" || true)
echo "Warning gate: PASS (project warnings: 0, toolchain warnings: $external)"
