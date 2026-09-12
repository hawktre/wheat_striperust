#!/bin/bash
#SBATCH --job-name=StripeSourceSim
#SBATCH --partition=share
#SBATCH --output=output/simulation/whole_epidemic_study/hpc/logs/simulation_%A_%a.out
#SBATCH --error=output/simulation/whole_epidemic_study/hpc/logs/simulation_%A_%a.err
#SBATCH --cpus-per-task=10
#SBATCH --mem=5G
#SBATCH --time=02:00:00

set -euo pipefail

# The submission helper changes to the repository root before calling sbatch.
# Make that location explicit for every relative path used by R and here().
cd "${SLURM_SUBMIT_DIR:?SLURM_SUBMIT_DIR is not set}"

echo "Host: ${HOSTNAME}"
echo "Job: ${SLURM_JOB_ID}"
echo "Batch: ${SLURM_ARRAY_TASK_ID}"
echo "Cores: ${SLURM_CPUS_PER_TASK}"
echo "Started: $(date)"

module load R
export R_LIBS="${HOME}/R_libs/4.4"

# Each R process is single-threaded. Set these after loading the module so its
# environment cannot restore a larger thread count. Parallelism is only across
# complete simulation replicates in the xargs worker pool below.
export OMP_NUM_THREADS=1
export OMP_THREAD_LIMIT=1
export OPENBLAS_NUM_THREADS=1
export BLIS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export NUMEXPR_NUM_THREADS=1
export VECLIB_MAXIMUM_THREADS=1
export RCPP_PARALLEL_NUM_THREADS=1

simulations_per_job="${SIMULATIONS_PER_JOB:-10}"
if ! [[ "${simulations_per_job}" =~ ^[0-9]+$ ]] ||
   (( simulations_per_job < 1 )); then
  echo "SIMULATIONS_PER_JOB must be a positive integer."
  exit 1
fi
first_simulation=$(( (SLURM_ARRAY_TASK_ID - 1) * simulations_per_job + 1 ))
last_simulation=$(( first_simulation + simulations_per_job - 1 ))
batch_status=0

echo "Simulations: ${first_simulation}-${last_simulation}"
if ! seq "${first_simulation}" "${last_simulation}" |
  xargs -n 1 -P "${SLURM_CPUS_PER_TASK}" \
    env SLURM_CPUS_PER_TASK=1 \
    Rscript --vanilla src/simulation/run_simulation_replicate_hpc.R; then
  batch_status=1
fi

echo "Finished: $(date)"
exit "${batch_status}"
