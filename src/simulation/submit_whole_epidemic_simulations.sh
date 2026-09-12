#!/bin/bash

set -euo pipefail

SCRIPT_DIRECTORY="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
REPOSITORY_ROOT="$(cd -- "${SCRIPT_DIRECTORY}/../.." && pwd)"
cd "${REPOSITORY_ROOT}"

START_ID="${1:-1}"
END_ID="${2:-100}"
MAX_CONCURRENT="${3:-60}"
SIMULATIONS_PER_JOB="${4:-10}"

if ! [[ "${START_ID}" =~ ^[0-9]+$ ]] ||
   ! [[ "${END_ID}" =~ ^[0-9]+$ ]] ||
   ! [[ "${MAX_CONCURRENT}" =~ ^[0-9]+$ ]] ||
   ! [[ "${SIMULATIONS_PER_JOB}" =~ ^[0-9]+$ ]]; then
  echo "Usage: $0 [start_batch] [end_batch] [maximum_concurrent_jobs] [simulations_per_job]"
  exit 1
fi
if (( START_ID < 1 || END_ID < START_ID || MAX_CONCURRENT < 1 || SIMULATIONS_PER_JOB < 1 )); then
  echo "Batch IDs, concurrency, and simulations per job must be positive."
  exit 1
fi
if ! command -v sbatch >/dev/null 2>&1; then
  echo "sbatch was not found; run this helper on a Slurm login node."
  exit 1
fi

mkdir -p output/simulation/whole_epidemic_study/hpc/logs
mkdir -p output/simulation/whole_epidemic_study/hpc/complete
mkdir -p output/simulation/whole_epidemic_study/hpc/failed

echo "Submitting batches ${START_ID}-${END_ID}, ${SIMULATIONS_PER_JOB} simulations per batch"
echo "At most ${MAX_CONCURRENT} jobs will run at once"
sbatch \
  --array="${START_ID}-${END_ID}%${MAX_CONCURRENT}" \
  --export="ALL,SIMULATIONS_PER_JOB=${SIMULATIONS_PER_JOB}" \
  "${SCRIPT_DIRECTORY}/run_whole_epidemic_simulation.sh"
