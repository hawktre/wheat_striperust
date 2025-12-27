#!/bin/bash

# Get NBATCH from command line argument
NBATCH=$1

# Validate input
if [ -z "$NBATCH" ]; then
    echo "Error: Please provide NBATCH as first argument"
    echo "Usage: $0 <NBATCH>"
    exit 1
fi

if ! [[ "$NBATCH" =~ ^[0-9]+$ ]]; then
    echo "Error: NBATCH must be a positive integer"
    exit 1
fi

# Clear old log files
rm -f DataProcessed/results/simulation/logs/simulation_*.out DataProcessed/results/simulation/logs/simulation_*.err

# Create logs directory if it doesn't exist
mkdir -p DataProcessed/results/simulation/logs

# Load R module to get Rscript
module load R

# Set R library path (if needed)
export R_LIBS="$HOME/R_libs/4.4"

# Calculate array end index (bash arithmetic)
ARRAY_END=$((NBATCH - 1))

echo "Submitting array job with $NBATCH tasks (0-$ARRAY_END)"

sbatch --array=0-${ARRAY_END}%20 Code/shell/utils/run_simulation.sh