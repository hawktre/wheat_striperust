#!/bin/bash
#SBATCH -J StripeSim       # name of my job
#SBATCH -p dri.q        # name of partition/queue to use
#SBATCH -o stripe_sim.out        # name of output file for batch script
#SBATCH -e stripe_sim.err        # name of error file for batch script 
#SBATCH -c 32           # number of cores per task
#SBATCH --mem=64g        # memory needed for job
#SBATCH --time=7-00:00:00  # time needed for job

# gather basic information, can be useful for troubleshooting
hostname
echo $SLURM_JOBID

#export R library path
export R_LIBS=~/R/x86_64-pc-linux-gnu-library/4.4

# load modules needed for job
module load R

# run my job
date

# run R code, print output to file "hello.out"
Rscript --vanilla Code/04b_RunSim.R

# print date and time when job is finished
date
echo "R job completed successfully."
# end of script
exit

