#!/usr/bin/env bash
set -euo pipefail

root=$(cd "$(dirname "$0")/.." && pwd)
cpu=${MBS_CPU:-$root/cpu/bin/mbs_po_mpi}
gpu=${MBS_GPU:-$root/gpu/bin/mbs_po_gpu_double}
device=${MBS_CUDA_TEST_DEVICE:-0}
relative_tolerance=${MBS_CUDA_L2_TOLERANCE:-1e-5}
compact_tolerance=${MBS_CUDA_COMPACT_L2_TOLERANCE:-1e-12}
beam4_tolerance=${MBS_CUDA_BEAM4_L2_TOLERANCE:-$compact_tolerance}
atomic_tolerance=${MBS_CUDA_ATOMIC_L2_TOLERANCE:-1e-10}
work=$(mktemp -d /tmp/mbs-cuda-gate.XXXXXX)
keep_work=1
cleanup() {
    if [[ "$keep_work" == 0 ]]; then
        rm -rf "$work"
    else
        echo "CUDA release-gate diagnostics preserved in $work" >&2
    fi
}
trap cleanup EXIT

for binary in "$cpu" "$gpu"; do
    if [[ ! -x "$binary" ]]; then
        echo "ERROR: release-gate binary is missing: $binary" >&2
        exit 1
    fi
done
if ! nvidia-smi -i "$device" --query-gpu=index --format=csv,noheader >/dev/null 2>&1; then
    echo "ERROR: CUDA test device $device is unavailable" >&2
    exit 1
fi

common=(
    --method po
    --particle 10 1.0 0.7 43.4114540493
    --refractive-index 1.3116 0.01
    --absorption
    --wavelength-um 0.532
    --max-reflections 4
    --so3-quaternion 64
    --scattering-grid 0 180 24 12
    --latitude-phi-grid
    --cutoff-profile off
    --threads 2
    --close
)

"$cpu" "${common[@]}" --backend cpu --output "$work/cpu" \
    >"$work/cpu.log" 2>&1
MBS_GPU_SLOT="$device" "$gpu" "${common[@]}" --backend cuda \
    --output "$work/gpu" >"$work/gpu.log" 2>&1
MBS_GPU_SLOT="$device" MBS_GPU_COMPACT_BEAMS=0 \
    "$gpu" "${common[@]}" --backend cuda \
    --output "$work/gpu_generic" >"$work/gpu_generic.log" 2>&1
MBS_GPU_SLOT="$device" MBS_GPU_WARP_BEAMS=0 \
    MBS_GPU_COMPACT_BEAM4_SPLIT=0 \
    "$gpu" "${common[@]}" --backend cuda \
    --output "$work/gpu_beam8" >"$work/gpu_beam8.log" 2>&1
MBS_GPU_SLOT="$device" MBS_GPU_WARP_BEAMS=0 \
    MBS_GPU_COMPACT_BEAM4_SPLIT=1 MBS_GPU_TIMING=1 \
    "$gpu" "${common[@]}" --backend cuda \
    --output "$work/gpu_beam4" >"$work/gpu_beam4.log" 2>&1
MBS_GPU_SLOT="$device" MBS_GPU_NO_ATOMICS=0 \
    "$gpu" "${common[@]}" --backend cuda \
    --allow-experimental-environment \
    --output "$work/gpu_atomic" >"$work/gpu_atomic.log" 2>&1
MBS_GPU_SLOT="$device" MBS_GPU_THREAD_GRID_3D=1 \
    "$gpu" "${common[@]}" --backend cuda \
    --allow-experimental-environment \
    --output "$work/gpu_grid3d" >"$work/gpu_grid3d.log" 2>&1

mixed_common=(
    --method po
    --particle 12 1.0 0.7 2
    --refractive-index 1.3116 0
    --wavelength-um 0.532
    --max-reflections 4
    --so3-quaternion 64
    --scattering-grid 0 180 24 12
    --latitude-phi-grid
    --cutoff-profile off
    --threads 2
    --close
)

MBS_GPU_SLOT="$device" MBS_GPU_WARP_BEAMS=0 \
    MBS_GPU_COMPACT_BEAM4_SPLIT=0 \
    "$gpu" "${mixed_common[@]}" --backend cuda \
    --output "$work/gpu_mixed8" >"$work/gpu_mixed8.log" 2>&1
MBS_GPU_SLOT="$device" MBS_GPU_WARP_BEAMS=0 \
    MBS_GPU_COMPACT_BEAM4_SPLIT=1 MBS_GPU_TIMING=1 \
    "$gpu" "${mixed_common[@]}" --backend cuda \
    --output "$work/gpu_mixed4" >"$work/gpu_mixed4.log" 2>&1

grep -q 'GPU backend: CUDA' "$work/gpu.log" || {
    echo 'ERROR: CUDA binary did not confirm use of the CUDA backend' >&2
    exit 1
}
grep -q 'path=beam4+8' "$work/gpu_beam4.log" || {
    echo 'ERROR: CUDA Beam4 release gate did not use the Beam4+Beam8 path' >&2
    exit 1
}
grep -q 'path=beam4+8+generic' "$work/gpu_mixed4.log" || {
    echo 'ERROR: CUDA mixed Beam4 release gate did not use all three beam layouts' >&2
    exit 1
}

python3 - "$work/cpu/cpu.dat" "$work/gpu/gpu.dat" \
    "$work/gpu_generic/gpu_generic.dat" "$relative_tolerance" \
    "$compact_tolerance" "$work/gpu_atomic/gpu_atomic.dat" \
    "$atomic_tolerance" "$work/gpu_grid3d/gpu_grid3d.dat" \
    "$work/gpu_beam8/gpu_beam8.dat" "$work/gpu_beam4/gpu_beam4.dat" \
    "$beam4_tolerance" <<'PY'
import math
import sys

def load(path):
    with open(path, encoding="utf-8") as stream:
        next(stream)
        return [[float(value) for value in line.split()] for line in stream if line.strip()]

def relative_l2(reference, candidate):
    if len(reference) != len(candidate) or not reference:
        raise SystemExit("ERROR: compared output grids have different sizes")
    num = 0.0
    den = 0.0
    for reference_row, candidate_row in zip(reference, candidate):
        if (abs(reference_row[0] - candidate_row[0]) > 1e-10 or
                abs(reference_row[1] - candidate_row[1]) > 1e-10):
            raise SystemExit("ERROR: compared output grids differ")
        weight = reference_row[1]
        for index in range(2, 18):
            delta = candidate_row[index] - reference_row[index]
            num += weight * delta * delta
            den += weight * reference_row[index] * reference_row[index]
    return math.sqrt(num / max(den, 1e-300))

cpu = load(sys.argv[1])
gpu = load(sys.argv[2])
gpu_generic = load(sys.argv[3])
relative = relative_l2(cpu, gpu)
tolerance = float(sys.argv[4])
print("CUDA release gate weighted L2: {:.3e} (limit {:.3e})".format(
    relative, tolerance))
if not math.isfinite(relative) or relative > tolerance:
    raise SystemExit("ERROR: CPU/GPU mismatch exceeds {:.3e}".format(tolerance))

compact_relative = relative_l2(gpu_generic, gpu)
compact_tolerance = float(sys.argv[5])
print("CUDA compact/generic weighted L2: {:.3e} (limit {:.3e})".format(
    compact_relative, compact_tolerance))
if not math.isfinite(compact_relative) or compact_relative > compact_tolerance:
    raise SystemExit("ERROR: compact/generic CUDA mismatch exceeds {:.3e}".format(
        compact_tolerance))

gpu_atomic = load(sys.argv[6])
atomic_relative = relative_l2(gpu, gpu_atomic)
atomic_tolerance = float(sys.argv[7])
print("CUDA fused/atomic weighted L2: {:.3e} (limit {:.3e})".format(
    atomic_relative, atomic_tolerance))
if not math.isfinite(atomic_relative) or atomic_relative > atomic_tolerance:
    raise SystemExit("ERROR: fused/atomic CUDA mismatch exceeds {:.3e}".format(
        atomic_tolerance))

gpu_grid3d = load(sys.argv[8])
grid3d_relative = relative_l2(gpu, gpu_grid3d)
print("CUDA flat/tiled-3D weighted L2: {:.3e} (limit {:.3e})".format(
    grid3d_relative, atomic_tolerance))
if not math.isfinite(grid3d_relative) or grid3d_relative > atomic_tolerance:
    raise SystemExit("ERROR: tiled 3D CUDA mismatch exceeds {:.3e}".format(
        atomic_tolerance))

gpu_beam8 = load(sys.argv[9])
gpu_beam4 = load(sys.argv[10])
beam4_tolerance = float(sys.argv[11])
beam4_relative = relative_l2(gpu_beam8, gpu_beam4)
print("CUDA Beam8/Beam4+8 weighted L2: {:.3e} (limit {:.3e})".format(
    beam4_relative, beam4_tolerance))
if not math.isfinite(beam4_relative) or beam4_relative > beam4_tolerance:
    raise SystemExit("ERROR: Beam4+8 CUDA mismatch exceeds {:.3e}".format(
        beam4_tolerance))
PY

python3 - "$work/gpu_mixed8/gpu_mixed8.dat" \
    "$work/gpu_mixed4/gpu_mixed4.dat" "$beam4_tolerance" <<'PY'
import math
import sys

def load(path):
    with open(path, encoding="utf-8") as stream:
        next(stream)
        return [[float(value) for value in line.split()]
                for line in stream if line.strip()]

reference = load(sys.argv[1])
candidate = load(sys.argv[2])
if len(reference) != len(candidate) or not reference:
    raise SystemExit("ERROR: mixed Beam4 output grids have different sizes")
num = 0.0
den = 0.0
for reference_row, candidate_row in zip(reference, candidate):
    if (abs(reference_row[0] - candidate_row[0]) > 1e-10 or
            abs(reference_row[1] - candidate_row[1]) > 1e-10):
        raise SystemExit("ERROR: mixed Beam4 output grids differ")
    weight = reference_row[1]
    for index in range(2, 18):
        delta = candidate_row[index] - reference_row[index]
        num += weight * delta * delta
        den += weight * reference_row[index] * reference_row[index]
relative = math.sqrt(num / max(den, 1e-300))
tolerance = float(sys.argv[3])
print("CUDA mixed Beam8/Beam4+8+generic weighted L2: {:.3e} (limit {:.3e})".format(
    relative, tolerance))
if not math.isfinite(relative) or relative > tolerance:
    raise SystemExit("ERROR: mixed Beam4 CUDA mismatch exceeds {:.3e}".format(
        tolerance))
PY

trace_common=(
    --method po
    --particle 10 1.0 0.7 43.4114540493
    --refractive-index 1.3116 0
    --wavelength-um 0.532
    --max-reflections 8
    --sobol 4
    --scattering-grid 0 180 24 12
    --cutoff-profile off
    --threads 1
    --close
)

MBS_GPU_SLOT="$device" "$gpu" "${trace_common[@]}" --backend cuda \
    --no-gpu-trace-prefilter \
    --output "$work/trace_reference" >"$work/trace_reference.log" 2>&1
MBS_GPU_SLOT="$device" "$gpu" "${trace_common[@]}" --backend cuda \
    --trace-prefilter-stats \
    --output "$work/trace_gpu" >"$work/trace_gpu.log" 2>&1

grep -q 'GPU trace profile: disabled' "$work/trace_reference.log" || {
    echo 'ERROR: CUDA tracing reference did not disable the automatic profile' >&2
    exit 1
}
grep -q 'GPU trace profile: automatic' "$work/trace_gpu.log" || {
    echo 'ERROR: CUDA tracing candidate did not enable the automatic profile' >&2
    exit 1
}
grep -Eq 'gpu_calls=[1-9][0-9]*' "$work/trace_gpu.log" || {
    echo 'ERROR: CUDA tracing gate did not offload facet ordering' >&2
    exit 1
}

python3 - "$work/trace_reference/trace_reference.dat" \
    "$work/trace_gpu/trace_gpu.dat" "$compact_tolerance" <<'PY'
import math
import sys

def load(path):
    with open(path, encoding="utf-8") as stream:
        next(stream)
        return [[float(value) for value in line.split()]
                for line in stream if line.strip()]

reference = load(sys.argv[1])
candidate = load(sys.argv[2])
if len(reference) != len(candidate) or not reference:
    raise SystemExit("ERROR: CUDA tracing output grids have different sizes")
num = 0.0
den = 0.0
for reference_row, candidate_row in zip(reference, candidate):
    if (abs(abs(reference_row[0]) - abs(candidate_row[0])) > 1e-10 or
            abs(reference_row[1] - candidate_row[1]) > 1e-10):
        raise SystemExit("ERROR: CUDA tracing output grids differ")
    weight = reference_row[1]
    for index in range(2, 18):
        delta = candidate_row[index] - reference_row[index]
        num += weight * delta * delta
        den += weight * reference_row[index] * reference_row[index]
relative = math.sqrt(num / max(den, 1e-300))
tolerance = float(sys.argv[3])
print("CUDA tracing weighted L2: {:.3e} (limit {:.3e})".format(
    relative, tolerance))
if not math.isfinite(relative) or relative > tolerance:
    raise SystemExit(
        "ERROR: CUDA tracing mismatch exceeds {:.3e}".format(tolerance))
PY

keep_work=0
echo 'CUDA release gate: PASS'
