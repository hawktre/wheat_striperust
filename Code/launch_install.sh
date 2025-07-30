#!/bin/bash
#SBATCH --job-name=installr
#SBATCH --output=installr.out
#SBATCH --error=installr.err
#SBATCH --time=00:30:00
#SBATCH --mem=4G
#SBATCH --cpus-per-task=1
#SBATCH --partition=dri.q

module load R
Rscript Code/install_pkgs.R
echo "Package installation completed successfully."