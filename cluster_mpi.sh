#!/bin/bash
# =============================================================================
# SLURM job script for MBS-fast MPI+OpenMP cluster run
#
# Usage:
#   1. Edit parameters below (particle, wavelength, orientations, etc.)
#   2. Submit: sbatch cluster_mpi.sh
#   3. Results in: results/M.dat, results/M_noshadow.dat
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
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK:-64}
export OMP_PROC_BIND=close
export OMP_PLACES=cores

# --- Path to binary ---
MBS=./bin/mbs_po_mpi

# --- Check binary exists ---
if [ ! -x "$MBS" ]; then
    echo "ERROR: $MBS not found. Run: bash build_mpi.sh"
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
OUTPUT=results

# =============================================================================
# RUN
# =============================================================================

srun $MBS --po --sobol $N_ORIENT --auto_tgrid --auto_phi \
    -p 1 $L $D -w $WAVE --ri $RI_RE $RI_IM -n $N_REFL \
    --close -o $OUTPUT

echo ""
echo "=== Done ==="
echo "Output: ${OUTPUT}.dat, ${OUTPUT}_noshadow.dat"
