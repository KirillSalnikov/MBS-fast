#!/usr/bin/env bash
set -euo pipefail

root=$(cd "$(dirname "$0")/.." && pwd)
cpu=${MBS_CPU:-$root/cpu/bin/mbs_po_mpi}
gpu=${MBS_GPU:-$root/gpu/bin/mbs_po_gpu_double}
device=${MBS_CUDA_TEST_DEVICE:-0}
relative_tolerance=${MBS_CUDA_L2_TOLERANCE:-1e-5}
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

grep -q 'GPU backend: CUDA' "$work/gpu.log" || {
    echo 'ERROR: CUDA binary did not confirm use of the CUDA backend' >&2
    exit 1
}

python3 - "$work/cpu/cpu.dat" "$work/gpu/gpu.dat" "$relative_tolerance" <<'PY'
import math
import sys

def load(path):
    with open(path, encoding="utf-8") as stream:
        next(stream)
        return [[float(value) for value in line.split()] for line in stream if line.strip()]

cpu = load(sys.argv[1])
gpu = load(sys.argv[2])
if len(cpu) != len(gpu) or not cpu:
    raise SystemExit("ERROR: CPU/GPU output grids have different sizes")

num = 0.0
den = 0.0
for crow, grow in zip(cpu, gpu):
    if abs(crow[0] - grow[0]) > 1e-10 or abs(crow[1] - grow[1]) > 1e-10:
        raise SystemExit("ERROR: CPU/GPU output grids differ")
    weight = crow[1]
    for index in range(2, 18):
        delta = grow[index] - crow[index]
        num += weight * delta * delta
        den += weight * crow[index] * crow[index]

relative = math.sqrt(num / max(den, 1e-300))
tolerance = float(sys.argv[3])
print("CUDA release gate weighted L2: {:.3e} (limit {:.3e})".format(
    relative, tolerance))
if not math.isfinite(relative) or relative > tolerance:
    raise SystemExit("ERROR: CPU/GPU mismatch exceeds {:.3e}".format(tolerance))
PY

keep_work=0
echo 'CUDA release gate: PASS'
