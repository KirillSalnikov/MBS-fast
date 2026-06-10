#!/usr/bin/env bash
set -euo pipefail

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
repo_root="$(cd "${script_dir}/.." && pwd)"
bin="${MBS_GPU_MPI_BIN:-${repo_root}/gpu/bin/mbs_po_gpu_mpi_double_fast}"

if [[ ! -x "${bin}" ]]; then
    echo "ERROR: MPI GPU binary not found: ${bin}" >&2
    echo "Build it with: make -C gpu double_fast_mpi" >&2
    exit 2
fi

if [[ -n "${CUDA_VISIBLE_DEVICES:-}" && "${CUDA_VISIBLE_DEVICES}" != "-1" ]]; then
    gpu_count="$(awk -F, '{print NF}' <<<"${CUDA_VISIBLE_DEVICES}")"
elif command -v nvidia-smi >/dev/null 2>&1; then
    gpu_count="$(nvidia-smi --query-gpu=index --format=csv,noheader 2>/dev/null | wc -l)"
else
    gpu_count=1
fi
gpu_count="${gpu_count//[[:space:]]/}"
if [[ -z "${gpu_count}" || "${gpu_count}" -lt 1 ]]; then
    gpu_count=1
fi

ranks="${MBS_MPI_RANKS:-${gpu_count}}"
if [[ "${ranks}" -lt 1 ]]; then
    ranks=1
fi

if [[ -n "${MBS_THREADS_PER_RANK:-}" ]]; then
    export OMP_NUM_THREADS="${MBS_THREADS_PER_RANK}"
elif [[ -z "${OMP_NUM_THREADS:-}" ]]; then
    if command -v lscpu >/dev/null 2>&1; then
        sockets="$(lscpu | awk -F: '/^Socket\(s\):/ {gsub(/[[:space:]]/, "", $2); print $2; exit}')"
        cores_per_socket="$(lscpu | awk -F: '/^Core\(s\) per socket:/ {gsub(/[[:space:]]/, "", $2); print $2; exit}')"
        if [[ -n "${sockets}" && -n "${cores_per_socket}" ]]; then
            cores=$(( sockets * cores_per_socket ))
        else
            cores=0
        fi
    else
        cores=0
    fi
    if [[ "${cores}" -le 0 ]] && command -v nproc >/dev/null 2>&1; then
        cores="$(nproc)"
    fi
    if [[ "${cores}" -gt 0 ]]; then
        threads=$(( cores / ranks ))
        if [[ "${threads}" -lt 1 ]]; then
            threads=1
        fi
        if [[ "${threads}" -gt 16 ]]; then
            threads=16
        fi
        export OMP_NUM_THREADS="${threads}"
    fi
fi

if [[ "${ranks}" -gt 1 ]]; then
    export MBS_GPU_MULTI="${MBS_GPU_MULTI:-0}"
fi
export MBS_GPU_FUSED_MUELLER="${MBS_GPU_FUSED_MUELLER:-1}"
export MBS_HOST_MEM_FRACTION="${MBS_HOST_MEM_FRACTION:-0.90}"
export MBS_OLDAUTO_BYTES_PER_GAMMA_MB="${MBS_OLDAUTO_BYTES_PER_GAMMA_MB:-64}"

args=("$@")
has_threads=0
for arg in "${args[@]}"; do
    if [[ "${arg}" == "--threads" ]]; then
        has_threads=1
        break
    fi
done
if [[ "${has_threads}" -eq 0 && -n "${OMP_NUM_THREADS:-}" ]]; then
    args=("--threads" "${OMP_NUM_THREADS}" "${args[@]}")
fi

mpirun_bin="${MPIEXEC:-mpirun}"
exec "${mpirun_bin}" -np "${ranks}" --bind-to none --map-by slot "${bin}" "${args[@]}"
