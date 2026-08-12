#!/usr/bin/env bash
set -euo pipefail

binary=${1:-"$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)/cpu/bin/mbs_po_mpi"}
root=${2:-"$PWD/po_orientation_sampling_L20_mu0p662222222"}
threads=${MBS_THREADS:-8}
backend=${MBS_BACKEND:-cpu}

q_levels=(
  0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1.0
  1.1 1.2 1.3 1.4 1.5 1.6 1.7 1.8 1.9 2.0
)

mkdir -p "$root"
cat > "$root/configuration.csv" <<EOF
parameter,value
particle,concave hexagonal column
L_um,20
D_over_L,0.7
mu,0.662222222
cavity_angle_deg,43.4114540493
wavelength_um,0.532
refractive_index,1.3116+0i
k_eq,92.0644601312
max_reflections,12
scattering_diffraction_sampling,0.5
latitude_phi_grid,true
backend,$backend
EOF

for q in "${q_levels[@]}"; do
  tag=${q/./p}
  run="$root/q_${tag}"
  mkdir -p "$run"
  if [[ -s "$run/result/result.dat" ]]; then
    printf '%s SKIP Q=%s complete\n' "$(date --iso-8601=seconds)" "$q" | tee -a "$root/queue.log"
    continue
  fi

  command=(
    "$binary"
    --method po
    --particle 10 1.0 0.7 43.4114540493
    --k-eq 92.0644601312
    --refractive-index 1.3116 0.0
    --wavelength-um 0.532
    --max-reflections 12
    --orientation-diffraction-sampling "$q"
    --scattering-diffraction-sampling 0.5
    --latitude-phi-grid
    --backend "$backend"
    --threads "$threads"
    --beam-cutoff-jones 0.001
    --beam-cutoff-area 0.002
    --trace-cutoff-importance 0.0001
    --trace-max-beams 20000
    --log 60
    --output "$run/result"
    --close
  )

  printf '%q ' "${command[@]}" > "$run/command.sh"
  printf '\n' >> "$run/command.sh"
  printf '%s START Q=%s\n' "$(date --iso-8601=seconds)" "$q" | tee -a "$root/queue.log"
  /usr/bin/time -f 'elapsed_seconds=%e\nmax_rss_kb=%M' \
    -o "$run/resources.txt" \
    "${command[@]}" > "$run/stdout.log" 2> "$run/stderr.log"
  printf '%s DONE Q=%s\n' "$(date --iso-8601=seconds)" "$q" | tee -a "$root/queue.log"
done
