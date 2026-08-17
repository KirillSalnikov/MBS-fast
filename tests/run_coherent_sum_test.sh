#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
BIN="${MBS_BIN:-$ROOT_DIR/bin/mbs_po}"
WORK_DIR="$(mktemp -d "${TMPDIR:-/tmp}/mbs-coherent-sum.XXXXXX")"
trap 'rm -rf "$WORK_DIR"' EXIT

if [[ ! -x "$BIN" ]]; then
    echo "Missing executable: $BIN" >&2
    exit 2
fi

MBS_FORCE_TRACK_IDS=1 \
MBS_COHERENCE_AUDIT="$WORK_DIR/beams.tsv" \
"$BIN" \
    --method po \
    --particle 10 1.0 0.7 43.4114540493 \
    --k-eq 18.41289202624 \
    --refractive-index 1.3116 0 \
    --wavelength-um 0.532 \
    --max-reflections 4 \
    --backend cpu \
    --threads 2 \
    --fixed-orientation 37 19 \
    --scattering-grid 179 180 12 10 \
    --cutoff-profile off \
    --jones-output \
    --allow-experimental-environment \
    --close \
    --output "$WORK_DIR/result" >/dev/null 2>&1

"$BIN" \
    --method po \
    --particle 10 1.0 0.7 43.4114540493 \
    --k-eq 18.41289202624 \
    --refractive-index 1.3116 0 \
    --wavelength-um 0.532 \
    --max-reflections 4 \
    --backend cpu \
    --threads 2 \
    --fixed-orientation 37 19 \
    --scattering-grid 179 180 12 10 \
    --cutoff-profile off \
    --incoherent \
    --close \
    --output "$WORK_DIR/incoherent" >/dev/null 2>&1

python3 - "$WORK_DIR/beams.tsv" "$WORK_DIR/result/result_jones.dat" \
    "$WORK_DIR/incoherent/incoherent.dat" <<'PY'
import math
import sys

audit_path, jones_path, incoherent_path = sys.argv[1:]
sums = [0.0] * 8
incoherent_m11 = 0.0
beam_count = 0

with open(audit_path, encoding="utf-8") as stream:
    for line in stream:
        if not line.strip() or line.startswith("#"):
            continue
        values = [float(value) for value in line.split("\t")[10:18]]
        if len(values) != 8 or not all(math.isfinite(value) for value in values):
            raise SystemExit("invalid per-beam Jones audit row")
        for index, value in enumerate(values):
            sums[index] += value
        incoherent_m11 += 0.5 * sum(value * value for value in values)
        beam_count += 1

if beam_count < 500:
    raise SystemExit(
        f"too few outgoing beams ({beam_count}); a tilted occluding facet may "
        "have been projected without clipping it at the source plane")

expected = None
with open(jones_path, encoding="utf-8") as stream:
    next(stream)
    for line in stream:
        row = [float(value) for value in line.split()]
        if abs(row[0] - 180.0) < 1e-12 and abs(row[1]) < 1e-12:
            expected = row[2:10]
            break

if expected is None:
    raise SystemExit("exact-backscatter Jones row was not written")

delta2 = sum((left - right) ** 2 for left, right in zip(sums, expected))
norm2 = max(sum(value * value for value in expected), 1e-300)
relative_jones_error = math.sqrt(delta2 / norm2)
if relative_jones_error > 1e-12:
    raise SystemExit(
        f"per-beam coherent sum differs from output Jones matrix: "
        f"relative error={relative_jones_error:.6g}")

coherent_m11 = 0.5 * sum(value * value for value in sums)
if coherent_m11 < 0.04:
    raise SystemExit(
        f"exact-backscatter M11 is unexpectedly low ({coherent_m11:.9g}); "
        "check forward-depth clipping of tilted occluding facets")
cross_fraction = (coherent_m11 - incoherent_m11) / incoherent_m11
if abs(cross_fraction) < 0.01:
    raise SystemExit(
        f"interference cross terms unexpectedly vanished: fraction={cross_fraction:.6g}")

reported_incoherent_m11 = None
with open(incoherent_path, encoding="utf-8") as stream:
    next(stream)
    for line in stream:
        row = [float(value) for value in line.split()]
        if abs(row[0] - 180.0) < 1e-12 and abs(row[1]) < 1e-12:
            reported_incoherent_m11 = row[3]
            break
if reported_incoherent_m11 is None:
    raise SystemExit("exact-backscatter incoherent Mueller row was not written")
relative_incoherent_error = abs(reported_incoherent_m11 - incoherent_m11) \
    / max(abs(incoherent_m11), 1e-300)
if relative_incoherent_error > 1e-12:
    raise SystemExit(
        f"per-beam Mueller sum differs from --incoherent output: "
        f"relative error={relative_incoherent_error:.6g}")

print(
    f"coherent sum: beams={beam_count}, Jones error={relative_jones_error:.3e}, "
    f"M11 coherent={coherent_m11:.9g}, incoherent={incoherent_m11:.9g}, "
    f"cross={cross_fraction:.3%}, incoherent error={relative_incoherent_error:.3e}")
PY
