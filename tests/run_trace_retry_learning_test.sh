#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
MBS=${MBS_BIN:-"$ROOT_DIR/bin/mbs_po"}
WORK_DIR=${MBS_TRACE_RETRY_TEST_DIR:-"/tmp/mbs_trace_retry_learning"}

rm -rf "$WORK_DIR"
mkdir -p "$WORK_DIR"

"$MBS" --method po --particle 10 1 1 20 \
    --refractive-index 1.31 0 --wavelength-um 10 \
    --max-reflections 3 --backend cpu --threads 1 --close \
    --euler-grid 1 3 --scattering-grid 0 180 12 4 \
    --trace-max-beams 100 --trace-limit-retries 2 \
    --trace-retry-factor 10 --output "$WORK_DIR/result" \
    >"$WORK_DIR/console.log" 2>&1

grep -q '^Learned retry start limit: 1000$' "$WORK_DIR/console.log"
grep -q '^Orientations starting at learned limit: 2$' "$WORK_DIR/console.log"
grep -q '^Configured beam-limit hits: 1$' "$WORK_DIR/console.log"
grep -q '^Adaptive retry attempts: 1$' "$WORK_DIR/console.log"
grep -q '^Orientations recovered by retry: 1$' "$WORK_DIR/console.log"
grep -q '^Orientations still incomplete: 0$' "$WORK_DIR/console.log"

echo "Trace retry learning test: PASS"
