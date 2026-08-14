#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
BIN="${MBS_BIN:-$ROOT_DIR/bin/mbs_po}"
WORK_DIR="$(mktemp -d "${TMPDIR:-/tmp}/mbs-pole-stability.XXXXXX")"
trap 'rm -rf "$WORK_DIR"' EXIT

if [[ ! -x "$BIN" ]]; then
    echo "Missing executable: $BIN" >&2
    exit 2
fi

"$BIN" \
    --method po \
    --particle 10 20 14 45 \
    --refractive-index 1.3116 0 \
    --wavelength-um 0.532 \
    --max-reflections 4 \
    --backend cpu \
    --threads 2 \
    --fixed-orientation 37 19 \
    --theta-grid-file "$ROOT_DIR/tests/pole_stability_theta_deg.grid" \
    --phi-points 120 \
    --close \
    --output "$WORK_DIR/result" >/dev/null 2>&1

"$BIN" \
    --method po \
    --particle 10 20 14 45 \
    --refractive-index 1.3116 0 \
    --wavelength-um 0.532 \
    --max-reflections 4 \
    --backend cpu \
    --threads 2 \
    --fixed-orientation 37 19 \
    --theta-grid-file "$ROOT_DIR/tests/pole_stability_theta_deg.grid" \
    --phi-points 120 \
    --karczewski \
    --close \
    --output "$WORK_DIR/karczewski" >/dev/null 2>&1

python3 - "$WORK_DIR/result/result.dat" \
    "$WORK_DIR/karczewski/karczewski.dat" <<'PY'
import math
import sys

path = sys.argv[1]
with open(path, encoding="utf-8") as stream:
    header = stream.readline().split()
    n_phi = int(header[1])
    rows = [[float(value) for value in line.split()] for line in stream if line.strip()]

with open(sys.argv[2], encoding="utf-8") as stream:
    next(stream)
    karczewski_rows = [
        [float(value) for value in line.split()]
        for line in stream if line.strip()
    ]
if not karczewski_rows or any(
        not math.isfinite(value)
        for row in karczewski_rows for value in row):
    raise SystemExit("Karczewski pole grid contains a non-finite value")
if len(karczewski_rows) != len(rows):
    raise SystemExit("Karczewski pole grid has a different size")
karczewski_m11_error = math.sqrt(
    sum((left[3] - right[3])**2
        for left, right in zip(rows, karczewski_rows))
    / max(sum(left[3]**2 for left in rows), 1e-300)
)
if karczewski_m11_error > 1e-10:
    raise SystemExit(
        f"Karczewski changed M11: relative L2={karczewski_m11_error:.6g}")

if len(rows) != 8 * n_phi:
    raise SystemExit(f"unexpected raw Mueller grid: {len(rows)} rows, n_phi={n_phi}")

def block(index):
    return rows[index * n_phi:(index + 1) * n_phi]

def average_m11(index):
    return sum(row[3] for row in block(index)) / n_phi

def average_rotated_m22(index):
    total = 0.0
    for row in block(index):
        azimuth = -math.radians(row[1])
        c2 = math.cos(2.0 * azimuth)
        s2 = math.sin(2.0 * azimuth)
        total += row[8] * c2 - row[9] * s2
    return total / n_phi

def exact_backward_m22():
    # PoleMueller::ApplyBackward: (M22-M33)/2 in one-based notation.
    return sum(0.5 * (row[8] - row[13]) for row in block(7)) / n_phi

forward = average_m11(0)
forward_near = average_m11(3)
backward_near = average_rotated_m22(6)
backward = exact_backward_m22()

forward_error = abs(forward_near - forward) / max(abs(forward), 1e-300)
backward_error = abs(backward_near - backward) / max(abs(backward), 1e-300)

if forward_error > 1e-5:
    raise SystemExit(f"forward small-phase discontinuity: relative error={forward_error:.6g}")
if backward_error > 1e-5:
    raise SystemExit(f"backward basis discontinuity: relative error={backward_error:.6g}")

print(
    f"pole stability: forward={forward_error:.3e}, "
    f"backward={backward_error:.3e}, "
    f"Karczewski M11={karczewski_m11_error:.3e}")
PY
