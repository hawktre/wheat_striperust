#!/bin/bash
#SBATCH -J StripeSim                  # Default job name (can override with --job-name)
#SBATCH -o stripe_sim.out             # Stdout log
#SBATCH -e stripe_sim.err             # Stderr log
#SBATCH -c 24                         # Default CPU cores
#SBATCH --mem=48G                     # Default memory per node
#SBATCH --time=7-00:00:00             # Default walltime

set -euo pipefail

# Basic job info
hostname
echo "Job ID:           ${SLURM_JOB_ID:-N/A}"
echo "Partition:        ${SLURM_JOB_PARTITION:-N/A}"
echo "Cores requested:  ${SLURM_CPUS_PER_TASK:-N/A}"
# Memory env varies by site; show both if present
echo "Mem per node:     ${SLURM_MEM_PER_NODE:-N/A}"
echo "Mem per CPU:      ${SLURM_MEM_PER_CPU:-N/A}"
echo "Start time:       $(date)"

# Make R respect SLURM core request for parallel code (if you use OpenMP/BLAS)
export OMP_NUM_THREADS="${SLURM_CPUS_PER_TASK:-1}"
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export NUMEXPR_NUM_THREADS=1

# Use consistent custom R library path
export R_LIBS="$HOME/R_libs/4.4"

# Load R module
module load R

# Positional arg 1 = number of sims (default 20)
NSIMS="${1:-20}"
echo "NSIMS:            ${NSIMS}"

# If your R uses MC, you can pass cores through an env var
export MC_CORES="${SLURM_CPUS_PER_TASK:-1}"

Rscript --vanilla Code/04b_RunSim.R "$NSIMS"

echo "Finished at $(date)"
echo "R job completed successfully."
