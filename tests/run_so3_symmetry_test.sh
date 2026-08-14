#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
MBS="${MBS:-$ROOT_DIR/cpu/bin/mbs_po_mpi}"
WORK_DIR="$(mktemp -d "${TMPDIR:-/tmp}/mbs-so3-symmetry.XXXXXX")"
trap 'rm -rf -- "$WORK_DIR"' EXIT

COMMON=(
    --method po
    --particle 10 1 0.7 43.4
    --refractive-index 1.3116 0
    --wavelength-um 10
    --max-reflections 1
    --backend cpu
    --threads 2
    --close
    --scattering-grid 0 180 24 12
)

OMP_NUM_THREADS=2 "$MBS" "${COMMON[@]}" \
    --so3-quaternion 512 \
    --save-geometry "$WORK_DIR/reduced-particle.dat" \
    --output "$WORK_DIR/reduced" >/dev/null 2>&1
OMP_NUM_THREADS=2 "$MBS" "${COMMON[@]}" \
    --so3-quaternion 512 --so3-mirror-audit \
    --output "$WORK_DIR/paired" >"$WORK_DIR/paired.console.log" 2>&1
OMP_NUM_THREADS=2 "$MBS" "${COMMON[@]}" \
    --so3-quaternion 256 --mirror-gamma \
    --output "$WORK_DIR/mirror" >"$WORK_DIR/mirror.console.log" 2>&1
OMP_NUM_THREADS=2 "$MBS" "${COMMON[@]}" \
    --so3-quaternion 512 --symmetry 1 1 \
    --output "$WORK_DIR/alpha_phi" >/dev/null 2>&1
OMP_NUM_THREADS=2 "$MBS" "${COMMON[@]}" \
    --so3-full-quaternion 12288 --output "$WORK_DIR/full" >/dev/null 2>&1

metrics()
{
    local reference=$1
    local candidate=$2
    awk '
        NR == FNR {
            if (NR > 1) {
                for (i = 3; i <= 18; ++i) ref[NR, i] = $i
            }
            next
        }
        FNR > 1 {
            for (i = 3; i <= 18; ++i) {
                d = $i - ref[FNR, i]
                all_diff += d * d
            }
            d11 = $3 - ref[FNR, 3]
            m11_diff += d11 * d11
            m11_ref += ref[FNR, 3] * ref[FNR, 3]
        }
        END {
            printf "%.12g %.12g\n", sqrt(m11_diff / m11_ref),
                   sqrt(all_diff / m11_ref)
        }
    ' "$reference" "$candidate"
}

read -r full_domain_m11 full_domain_all < <(
    metrics "$WORK_DIR/full/full.dat" "$WORK_DIR/alpha_phi/alpha_phi.dat")
read -r reduced_m11 reduced_all < <(
    metrics "$WORK_DIR/full/full.dat" "$WORK_DIR/reduced/reduced.dat")
read -r mirror_m11 mirror_all < <(
    metrics "$WORK_DIR/full/full.dat" "$WORK_DIR/mirror/mirror.dat")
read -r paired_mirror_m11 paired_mirror_all < <(
    metrics "$WORK_DIR/paired/paired.dat" "$WORK_DIR/mirror/mirror.dat")

awk -v m="$full_domain_m11" -v a="$full_domain_all" 'BEGIN {
    if (m > 0.01 || a > 0.03) exit 1
}' || {
    printf 'FAIL: alpha-via-phi differs from direct full SO(3): M11=%g all=%g\n' \
        "$full_domain_m11" "$full_domain_all" >&2
    exit 1
}
awk -v m="$reduced_m11" -v a="$reduced_all" 'BEGIN {
    if (m > 0.05 || a > 0.12) exit 1
}' || {
    printf 'FAIL: symmetry quotient differs from direct full SO(3): M11=%g all=%g\n' \
        "$reduced_m11" "$reduced_all" >&2
    exit 1
}
awk -v m="$mirror_m11" -v a="$mirror_all" 'BEGIN {
    if (m > 0.08 || a > 0.16) exit 1
}' || {
    printf 'FAIL: mirror-gamma SO(3) differs from direct full SO(3): M11=%g all=%g\n' \
        "$mirror_m11" "$mirror_all" >&2
    exit 1
}

# The full-gamma 512-point set and the mirror 256-point set now share the
# exact same beta/gamma base points.  This threshold therefore measures only
# numerical reflection consistency, not low-discrepancy sampling noise.
awk -v m="$paired_mirror_m11" -v a="$paired_mirror_all" 'BEGIN {
    if (m > 0.08 || a > 0.16) exit 1
}' || {
    printf 'FAIL: nested full/mirror SO(3) pair differs: M11=%g all=%g\n' \
        "$paired_mirror_m11" "$paired_mirror_all" >&2
    exit 1
}

grep -Eq '^90 59[.]9' "$WORK_DIR/reduced-particle.dat"
grep -q 'alpha integrated by scattering azimuth' \
    "$WORK_DIR/reduced/reduced_out.txt"
grep -q 'paired gamma audit: 256 half-domain Hammersley points' \
    "$WORK_DIR/paired.console.log"
grep -q 'gamma=0..30 deg' "$WORK_DIR/mirror.console.log"
grep -q 'omitted orientations restored as' "$WORK_DIR/mirror.console.log"
grep -q 'Validate each particle class against a full-gamma' \
    "$WORK_DIR/mirror.console.log"

printf 'SO(3) symmetry regression passed: full-domain M11 L2=%.4f%%, all/M11=%.4f%%; reduced M11 L2=%.4f%%, all/M11=%.4f%%; mirror-half M11 L2=%.4f%%, all/M11=%.4f%%; paired mirror M11 L2=%.4f%%, all/M11=%.4f%%\n' \
    "$(awk -v x="$full_domain_m11" 'BEGIN {print 100*x}')" \
    "$(awk -v x="$full_domain_all" 'BEGIN {print 100*x}')" \
    "$(awk -v x="$reduced_m11" 'BEGIN {print 100*x}')" \
    "$(awk -v x="$reduced_all" 'BEGIN {print 100*x}')" \
    "$(awk -v x="$mirror_m11" 'BEGIN {print 100*x}')" \
    "$(awk -v x="$mirror_all" 'BEGIN {print 100*x}')" \
    "$(awk -v x="$paired_mirror_m11" 'BEGIN {print 100*x}')" \
    "$(awk -v x="$paired_mirror_all" 'BEGIN {print 100*x}')"
