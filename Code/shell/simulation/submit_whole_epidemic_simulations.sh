#!/bin/bash

set -euo pipefail

START_ID="${1:-1}"
END_ID="${2:-1000}"
MAX_CONCURRENT="${3:-12}"

if ! [[ "${START_ID}" =~ ^[0-9]+$ ]] ||
   ! [[ "${END_ID}" =~ ^[0-9]+$ ]] ||
   ! [[ "${MAX_CONCURRENT}" =~ ^[0-9]+$ ]]; then
  echo "Usage: $0 [start_id] [end_id] [maximum_concurrent_tasks]"
  exit 1
fi
if (( START_ID < 1 || END_ID < START_ID || MAX_CONCURRENT < 1 )); then
  echo "Simulation IDs and concurrency must define a positive valid range."
  exit 1
fi

mkdir -p DataProcessed/results/simulation/whole_epidemic_study/hpc/logs
mkdir -p DataProcessed/results/simulation/whole_epidemic_study/hpc/complete
mkdir -p DataProcessed/results/simulation/whole_epidemic_study/hpc/failed

echo "Submitting simulations ${START_ID}-${END_ID}, at most ${MAX_CONCURRENT} at once"
sbatch \
  --array="${START_ID}-${END_ID}%${MAX_CONCURRENT}" \
  Code/shell/simulation/run_whole_epidemic_simulation.sh
