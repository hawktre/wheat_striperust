#!/bin/bash
#SBATCH -J StripeSim             # Job name
#SBATCH -o stripe_sim.out       # Stdout log
#SBATCH -e stripe_sim.err       # Stderr log
#SBATCH -c 24                   # Number of cores
#SBATCH --mem=48G               # Memory
#SBATCH --time=7-00:00:00       # Max runtime (7 days)

# Basic job info
hostname
echo "Job ID: $SLURM_JOBID"
echo "Running on $(date)"

# Use consistent custom R library path
export R_LIBS=$HOME/R_libs/4.4

# Load R module
module load R

# Pass the first command-line argument to Rscript (default to 20 if none given)
NSIMS=${1:-20}

Rscript --vanilla Code/04b_RunSim.R $NSIMS


# Confirm completion
echo "Finished at $(date)"
echo "R job completed successfully."
