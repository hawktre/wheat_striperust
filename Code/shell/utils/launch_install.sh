#!/bin/bash
#SBATCH --job-name=installr
#SBATCH --output=installr.out
#SBATCH --error=installr.err
#SBATCH --time=00:30:00
#SBATCH --mem=8G
#SBATCH --cpus-per-task=4

# Load R module
module load R

# Set user R library path explicitly to avoid any ambiguity
export R_LIBS=$HOME/R_libs/4.4

# Run your R script with clean session
Rscript --vanilla Code/installr.R

echo "Package installation completed successfully."
