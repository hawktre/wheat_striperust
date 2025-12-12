#!/bin/bash

# Clear old log files
rm -f DataProcessed/results/backward_model/logs/backward_fit_*.out DataProcessed/results/backward_model/logs/backward_fit_*.err

# Create logs directory if it doesn't exist
mkdir -p DataProcessed/results/backward_model/logs

# Load R module to get Rscript
module load R

# Set R library path (if needed)
export R_LIBS="$HOME/R_libs/4.4"

# Get number of rows
NROWS=$(Rscript --vanilla Code/get_nrows.R)

echo "Submitting array job with $NROWS tasks"

#sbatch --array=1-${NROWS}%20 Code/launch_backwardfit.sh