#!/bin/bash
#SBATCH -J MergeResults
#SBATCH -o merge_%j.out
#SBATCH -e merge_%j.err
#SBATCH -c 1
#SBATCH --mem=8G
#SBATCH --time=00:30:00

module load R
export R_LIBS="$HOME/R_libs/4.4"

Rscript Code/R/simulation/04c_CombineSimulations.R