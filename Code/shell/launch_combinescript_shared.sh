#!/bin/bash
#SBATCH -J MergeResults
#SBATCH -o merge_%j.out
#SBATCH -e merge_%j.err
#SBATCH -c 1
#SBATCH --mem=4G
#SBATCH --time=00:30:00

module load R
export R_LIBS="$HOME/R_libs/4.4"

Rscript Code/R/backward_model_shared/03d_CombineBackward_Results.R