#!/bin/bash
#SBATCH -J BackwardFit                # Job name
#SBATCH -o DataProcessed/results/backward_model/logs/backward_fit_%A_%a.out     # Stdout log (%A=job ID, %a=array index)
#SBATCH -e DataProcessed/results/backward_model/logs/backward_fit_%A_%a.err     # Stderr log
#SBATCH -c 10                         # CPU cores per task
#SBATCH --mem=1G                     # Memory per task
#SBATCH --time=02:00:00             # Walltime per task

set -euo pipefail

# Basic job info
hostname
echo "Job ID:           ${SLURM_JOB_ID:-N/A}"
echo "Array Task ID:    ${SLURM_ARRAY_TASK_ID:-N/A}"
echo "Partition:        ${SLURM_JOB_PARTITION:-N/A}"
echo "Cores requested:  ${SLURM_CPUS_PER_TASK:-N/A}"
echo "Mem per node:     ${SLURM_MEM_PER_NODE:-N/A}"
echo "Start time:       $(date)"

# Prevent OpenMP/BLAS from spawning extra threads
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export NUMEXPR_NUM_THREADS=1

# Use consistent custom R library path
export R_LIBS="$HOME/R_libs/4.4"

# Load R module
module load R

# Pass number of cores to R
export MC_CORES="${SLURM_CPUS_PER_TASK:-1}"

Rscript --vanilla Code/03c_BackwardModelFit.R

echo "Finished at $(date)"
echo "Array task ${SLURM_ARRAY_TASK_ID} completed successfully."