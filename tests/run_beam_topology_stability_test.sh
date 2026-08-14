#!/usr/bin/env bash
set -euo pipefail

repo=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
mbs=${MBS_BIN:-${MBS:-$repo/bin/mbs_po}}
work=$(mktemp -d "${TMPDIR:-/tmp}/mbs-beam-topology.XXXXXX")
trap 'rm -rf "$work"' EXIT

common=(
  --method po --particle 10 20 14 45
  --refractive-index 1.3116 0 --wavelength-um 0.532
  --max-reflections 4 --geometry nonconvex --cutoff-profile off
  --backend cpu --threads 1 --close
)

run_case() {
  local name=$1
  shift
  "$mbs" "${common[@]}" "$@" --output "$work/$name" >"$work/$name.log" 2>&1
}

run_case base --fixed-orientation 35.659087696 15 \
  --scattering-grid 0 180 12 18
run_case base_no_prefilter --fixed-orientation 35.659087696 15 \
  --scattering-grid 0 180 12 18 --no-trace-prefilter

"$mbs" --method po --particle 10 200 140 45 \
  --refractive-index 1.3116 0 --wavelength-um 5.32 \
  --max-reflections 4 --fixed-orientation 35.659087696 15 \
  --scattering-grid 0 180 12 18 --geometry nonconvex \
  --cutoff-profile off --backend cpu --threads 1 --close \
  --output "$work/scale10" >"$work/scale10.log" 2>&1

"$mbs" --method po --particle 10 0.02 0.014 45 \
  --refractive-index 1.3116 0 --wavelength-um 0.000532 \
  --max-reflections 4 --fixed-orientation 35.659087696 15 \
  --scattering-grid 0 180 12 18 --geometry nonconvex \
  --cutoff-profile off --backend cpu --threads 1 --close \
  --output "$work/scale001" >"$work/scale001.log" 2>&1

"$mbs" --method po --particle 10 20 14 45 \
  --refractive-index 1.3116 0.01 --wavelength-um 0.532 \
  --max-reflections 4 --fixed-orientation 35.659087696 15 \
  --scattering-grid 0 180 12 18 --geometry nonconvex \
  --cutoff-profile off --backend cpu --threads 1 --close \
  --output "$work/absorption" >"$work/absorption.log" 2>&1

"$mbs" --method po --particle 10 200 140 45 \
  --refractive-index 1.3116 0.01 --wavelength-um 5.32 \
  --max-reflections 4 --fixed-orientation 35.659087696 15 \
  --scattering-grid 0 180 12 18 --geometry nonconvex \
  --cutoff-profile off --backend cpu --threads 1 --close \
  --output "$work/absorption_scale10" \
  >"$work/absorption_scale10.log" 2>&1

"$mbs" --method po --particle 10 0.02 0.014 45 \
  --refractive-index 1.3116 0.01 --wavelength-um 0.000532 \
  --max-reflections 4 --fixed-orientation 35.659087696 15 \
  --scattering-grid 0 180 12 18 --geometry nonconvex \
  --cutoff-profile off --backend cpu --threads 1 --close \
  --output "$work/absorption_scale001" \
  >"$work/absorption_scale001.log" 2>&1

run_case mirror15 --fixed-orientation 35.659087696 15 \
  --scattering-grid 179.9 180 1 1
run_case mirror45 --fixed-orientation 35.659087696 45 \
  --scattering-grid 179.9 180 1 1

# Geometry output is opt-in.  A normal calculation must not leak the legacy
# particle_for_check.dat file into its working directory.
if find "$work" -name particle_for_check.dat -print -quit | grep -q .; then
  echo "unexpected implicit particle_for_check.dat" >&2
  exit 1
fi

"$mbs" "${common[@]}" --fixed-orientation 0 0 \
  --scattering-grid 0 180 1 1 --save-geometry "$work/saved.particle" \
  --output "$work/save" >"$work/save.log" 2>&1
test -s "$work/saved.particle"
test -s "$work/saved.particle.visibility"

# A particle's absolute position changes only a common phase. Translate every
# input vertex before loading and require exactly the same Mueller matrix. This
# catches float precision loss before topology and visibility tests begin.
python3 - "$work/saved.particle" "$work/translated.particle" <<'PY'
import sys

source, target = sys.argv[1:]
records = 0
output = []
with open(source) as stream:
    for line in stream:
        stripped = line.strip()
        if stripped and not stripped.startswith("#"):
            records += 1
            if records > 3:
                values = stripped.split()
                if len(values) == 3:
                    line = " ".join(
                        format(float(value) + 1000.0, ".17g")
                        for value in values) + "\n"
        output.append(line)
with open(target, "w") as stream:
    stream.writelines(output)
PY
cp "$work/saved.particle.visibility" "$work/translated.particle.visibility"

# Facet ids are serialization details. Reverse complete facet records and their
# visibility metadata; geometry-derived clipping order must keep the result
# bit-for-bit invariant.
python3 - "$work/saved.particle" "$work/reordered.particle" <<'PY'
import pathlib
import sys

source, target = map(pathlib.Path, sys.argv[1:])
lines = source.read_text().splitlines()
header = lines[:3]
facets = []
current = []
for line in lines[3:]:
    if line.strip():
        current.append(line)
    elif current:
        facets.append(current)
        current = []
if current:
    facets.append(current)
with target.open("w") as stream:
    stream.write("\n".join(header) + "\n\n")
    for facet in reversed(facets):
        cyclic = facet[1:] + facet[:1]
        stream.write("\n".join(cyclic) + "\n\n")
PY
python3 - "$work/saved.particle.visibility" \
  "$work/reordered.particle.visibility" <<'PY'
import pathlib
import sys

source, target = map(pathlib.Path, sys.argv[1:])
lines = source.read_text().splitlines()
target.write_text(lines[0] + "\n" + "\n".join(reversed(lines[1:])) + "\n")
PY

for geometry in saved translated reordered; do
  "$mbs" --method po --particle-file "$work/$geometry.particle" \
    --refractive-index 1.3116 0 --wavelength-um 0.532 \
    --max-reflections 4 --fixed-orientation 35.659087696 15 \
    --scattering-grid 0 180 12 18 --geometry nonconvex \
    --cutoff-profile off --backend cpu --threads 1 --close \
    --output "$work/translation_$geometry" \
    >"$work/translation_$geometry.log" 2>&1
done

# A beam-tree limit must fail the calculation. Treating the affected
# orientation as a valid zero contribution biases every orientation average.
if "$mbs" "${common[@]}" --fixed-orientation 35.659087696 15 \
  --scattering-grid 0 180 12 18 --trace-max-beams 1 \
  --output "$work/limited" >"$work/limited.log" 2>&1; then
  echo "beam-tree limit was silently accepted" >&2
  exit 1
fi
grep -q -- "--trace-max-beams=1" "$work/limited.log"

python3 - "$work" <<'PY'
import math
import pathlib
import sys

root = pathlib.Path(sys.argv[1])

def read(name):
    path = root / name / (name + ".dat")
    rows = []
    with path.open() as stream:
        next(stream)
        for line in stream:
            if line.strip():
                rows.append([float(value) for value in line.split()])
    return rows

def relative_l2(left, right, scale=1.0, columns=range(3, 19)):
    numerator = 0.0
    denominator = 0.0
    for a, b in zip(left, right):
        for column in columns:
            delta = a[column] / scale - b[column]
            numerator += delta * delta
            denominator += b[column] * b[column]
    return math.sqrt(numerator / max(denominator, 1e-300))

base = read("base")
without_prefilter = read("base_no_prefilter")
scale10 = read("scale10")
scale001 = read("scale001")
absorption = read("absorption")
absorption_scale10 = read("absorption_scale10")
absorption_scale001 = read("absorption_scale001")
translated_base = read("translation_saved")
translated_copy = read("translation_translated")
reordered_copy = read("translation_reordered")
if len(base) != len(without_prefilter) or len(base) != len(scale10):
    raise SystemExit("topology test grids have different lengths")
if len(absorption) != len(absorption_scale10):
    raise SystemExit("absorption scale grids have different lengths")

prefilter_error = relative_l2(base, without_prefilter)
scale_m11_error = relative_l2(scale10, base, 100.0, [3])
scale_all_error = relative_l2(scale10, base, 100.0)
small_scale_m11_error = relative_l2(scale001, base, 1.0e-6, [3])
small_scale_all_error = relative_l2(scale001, base, 1.0e-6)
absorption_scale_m11_error = relative_l2(
    absorption_scale10, absorption, 100.0, [3])
absorption_scale_all_error = relative_l2(
    absorption_scale10, absorption, 100.0)
absorption_small_m11_error = relative_l2(
    absorption_scale001, absorption, 1.0e-6, [3])
absorption_small_all_error = relative_l2(
    absorption_scale001, absorption, 1.0e-6)
translation_m11_error = relative_l2(
    translated_copy, translated_base, 1.0, [3])
translation_all_error = relative_l2(
    translated_copy, translated_base)
facet_order_m11_error = relative_l2(
    reordered_copy, translated_base, 1.0, [3])
facet_order_all_error = relative_l2(
    reordered_copy, translated_base)
if prefilter_error > 1e-11:
    raise SystemExit(f"CPU prefilter changed the result: {prefilter_error:.3e}")
if scale_m11_error > 2e-4 or scale_all_error > 5e-4:
    raise SystemExit(
        "scale invariance failed: "
        f"M11={scale_m11_error:.3e}, all={scale_all_error:.3e}")
if small_scale_m11_error > 2e-4 or small_scale_all_error > 5e-4:
    raise SystemExit(
        "small-scale invariance failed: "
        f"M11={small_scale_m11_error:.3e}, all={small_scale_all_error:.3e}")
if absorption_scale_m11_error > 2e-4 or absorption_scale_all_error > 5e-4:
    raise SystemExit(
        "absorption scale invariance failed: "
        f"M11={absorption_scale_m11_error:.3e}, "
        f"all={absorption_scale_all_error:.3e}")
if absorption_small_m11_error > 1e-3 or absorption_small_all_error > 2e-3:
    raise SystemExit(
        "small-scale absorption invariance failed: "
        f"M11={absorption_small_m11_error:.3e}, "
        f"all={absorption_small_all_error:.3e}")
if translation_m11_error > 1e-12 or translation_all_error > 1e-12:
    raise SystemExit(
        "translation invariance failed: "
        f"M11={translation_m11_error:.3e}, all={translation_all_error:.3e}")
if facet_order_m11_error > 1e-12 or facet_order_all_error > 1e-12:
    raise SystemExit(
        "facet-order invariance failed: "
        f"M11={facet_order_m11_error:.3e}, all={facet_order_all_error:.3e}")

mirror15 = read("mirror15")[-1][3]
mirror45 = read("mirror45")[-1][3]
mirror_error = abs(mirror15 - mirror45) / max(abs(mirror15), abs(mirror45), 1e-300)
if mirror_error > 1e-6:
    raise SystemExit(f"backscatter mirror error is too large: {mirror_error:.3e}")

print(
    "beam topology stability: "
    f"prefilter={prefilter_error:.3e}, "
    f"scale M11={scale_m11_error:.3e}, scale all={scale_all_error:.3e}, "
    f"small M11={small_scale_m11_error:.3e}, "
    f"small all={small_scale_all_error:.3e}, "
    f"absorption M11={absorption_scale_m11_error:.3e}, "
    f"absorption all={absorption_scale_all_error:.3e}, "
    f"absorption small M11={absorption_small_m11_error:.3e}, "
    f"absorption small all={absorption_small_all_error:.3e}, "
    f"translation M11={translation_m11_error:.3e}, "
    f"translation all={translation_all_error:.3e}, "
    f"facet-order M11={facet_order_m11_error:.3e}, "
    f"facet-order all={facet_order_all_error:.3e}, "
    f"mirror={mirror_error:.3e}")
PY
