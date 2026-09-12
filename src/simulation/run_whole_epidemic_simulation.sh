#!/bin/bash
#SBATCH --job-name=StripeSourceSim
#SBATCH --partition=share
#SBATCH --output=output/simulation/whole_epidemic_study/hpc/logs/simulation_%A_%a.out
#SBATCH --error=output/simulation/whole_epidemic_study/hpc/logs/simulation_%A_%a.err
#SBATCH --cpus-per-task=8
#SBATCH --mem=20G
#SBATCH --time=02:00:00

set -euo pipefail

echo "Host: ${HOSTNAME}"
echo "Job: ${SLURM_JOB_ID}"
echo "Simulation: ${SLURM_ARRAY_TASK_ID}"
echo "Cores: ${SLURM_CPUS_PER_TASK}"
echo "Started: $(date)"

# Each R worker is single-threaded; parallelism occurs across complete
# block-treatment scenarios inside the array task.
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export NUMEXPR_NUM_THREADS=1
export VECLIB_MAXIMUM_THREADS=1

export R_LIBS="${HOME}/R_libs/4.4"
module load R

Rscript --vanilla src/simulation/run_simulation_replicate_hpc.R

echo "Finished: $(date)"
