#!/bin/bash
#SBATCH --job-name=extract
#SBATCH --cpus-per-task=2
#SBATCH --time=12:00:00
#SBATCH --output=0_step_logs/extract2.out
#SBATCH --error=0_step_logs/extract2.err
#SBATCH --mem=36G


module purge
module load conda_R/4.5.x
module list
conda activate aanri

python /dcs04/lieber/hwanglab/Arun/snRNA_aanri/1_source/2_split_GWAS_SS_2Chr.py
#original SBATCH --output=logs/extract2.out
#original SBATCH --error=logs/extract2.err
