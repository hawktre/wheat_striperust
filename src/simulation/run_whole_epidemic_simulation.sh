#!/bin/bash
#SBATCH --job-name=StripeSourceSim
#SBATCH --partition=share
#SBATCH --output=output/simulation/whole_epidemic_study/hpc/logs/simulation_%A_%a.out
#SBATCH --error=output/simulation/whole_epidemic_study/hpc/logs/simulation_%A_%a.err
#SBATCH --cpus-per-task=8
#SBATCH --mem=12G
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

# Each R worker is single-threaded; up to eight of the twelve complete
# block-treatment scenarios run in parallel inside the array task.
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export NUMEXPR_NUM_THREADS=1
export VECLIB_MAXIMUM_THREADS=1

export R_LIBS="${HOME}/R_libs/4.4"
module load R

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
for ((
  simulation_id = first_simulation;
  simulation_id <= last_simulation;
  simulation_id++
)); do
  if ! Rscript --vanilla \
    src/simulation/run_simulation_replicate_hpc.R "${simulation_id}"; then
    batch_status=1
  fi
done

echo "Finished: $(date)"
exit "${batch_status}"
