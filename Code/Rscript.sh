#!/bin/bash
#SBATCH -J sampleR       # name of my job
#SBATCH -p share         # name of partition/queue to use
#SBATCH -o sampleR.out        # name of output file for batch script
#SBATCH -e sampleR.err        # name of error file for batch script 
#SBATCH -c 2             # number of cores per task
#SBATCH --mem=10g        # memory needed for job
#SBATCH --time=3-00:00:00  # time needed for job

# gather basic information, can be useful for troubleshooting
hostname
echo $SLURM_JOBID

# load modules needed for job

module load R/4.2.2

# run my job
date

# run R code, print output to file "hello.out"
Rscript hello.R > hello.out 

sleep 100
date

exit

