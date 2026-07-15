#!/bin/bash
# =============================================================================
# SLURM job script for MBS-fast MPI+OpenMP cluster run
#
# Usage:
#   1. Edit parameters below (particle, wavelength, orientations, etc.)
#   2. Submit: sbatch cluster_mpi.sh
#   3. Results in the directory selected by OUTPUT below.
#
# Adjust #SBATCH directives for your cluster.
# =============================================================================

#SBATCH --job-name=mbs-fast
#SBATCH --nodes=8
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=64
#SBATCH --time=02:00:00
#SBATCH --output=mbs_%j.out
#SBATCH --error=mbs_%j.err

# --- OpenMP settings ---
# IMPORTANT: OMP_PROC_BIND=false to avoid conflict with MPI pinning.
# MPI launcher (srun) handles process placement; OpenMP should NOT
# further restrict threads within the allocated CPU set.
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK:-64}
export OMP_PROC_BIND=false
unset OMP_PLACES

# --- Path to binary ---
MBS=./cpu/bin/mbs_po_mpi

# --- Check binary exists ---
if [ ! -x "$MBS" ]; then
    echo "ERROR: $MBS not found."
    echo "  Fix: run bash build_mpi.sh from the repository root."
    exit 1
fi

echo "=== MBS-fast MPI+OpenMP ==="
echo "Nodes: $SLURM_NNODES"
echo "Tasks: $SLURM_NTASKS"
echo "Threads/task: $OMP_NUM_THREADS"
echo "Total cores: $(( SLURM_NTASKS * OMP_NUM_THREADS ))"
echo ""

# =============================================================================
# EDIT PARAMETERS HERE
# =============================================================================

# Particle: hex column D/L=0.7
L=100       # length (um)
D=70        # diameter (um)
WAVE=0.532  # wavelength (um)
RI_RE=1.31  # refractive index (real)
RI_IM=0     # refractive index (imaginary, 0 = no absorption)
N_REFL=14   # max internal reflections
N_ORIENT=8192  # Sobol orientations (power of 2)
OUTPUT=results/cluster_run

# =============================================================================
# RUN
# =============================================================================

srun "$MBS" --method po --sobol "$N_ORIENT" \
    --auto-theta-grid 0.02 --auto-phi \
    --particle 1 "$L" "$D" --wavelength-um "$WAVE" \
    --refractive-index "$RI_RE" "$RI_IM" --max-reflections "$N_REFL" \
    --threads "$OMP_NUM_THREADS" --close --output "$OUTPUT"

echo ""
echo "=== Done ==="
echo "Output directory: $OUTPUT"
