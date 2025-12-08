#!/bin/bash
#SBATCH -J BackwardFit                # Job name
#SBATCH -o backward_fit.out           # Stdout log
#SBATCH -e backward_fit.err           # Stderr log
#SBATCH -c 24                        # CPU cores
#SBATCH --mem=2G                      # Memory per node
#SBATCH --time=1-00:00:00             # Walltime (adjust as needed)

set -euo pipefail

# Basic job info
hostname
echo "Job ID:           ${SLURM_JOB_ID:-N/A}"
echo "Partition:        ${SLURM_JOB_PARTITION:-N/A}"
echo "Cores requested:  ${SLURM_CPUS_PER_TASK:-N/A}"
echo "Mem per node:     ${SLURM_MEM_PER_NODE:-N/A}"
echo "Mem per CPU:      ${SLURM_MEM_PER_CPU:-N/A}"
echo "Start time:       $(date)"

# Prevent OpenMP/BLAS from spawning extra threads
# Let your R parallel package control all parallelization
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export NUMEXPR_NUM_THREADS=1

# Use consistent custom R library path
export R_LIBS="$HOME/R_libs/4.4"

# Load R module
module load R

# Pass number of cores to R (your script can read this)
export MC_CORES="${SLURM_CPUS_PER_TASK:-1}"

Rscript --vanilla Code/03c_BackwardModelFit.R

echo "Finished at $(date)"
echo "R job completed successfully."