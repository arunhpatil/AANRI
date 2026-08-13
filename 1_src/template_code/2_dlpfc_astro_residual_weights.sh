#!/bin/bash
#SBATCH --job-name=Astro_res_fusionWTS
#SBATCH --array=40,41,25-38,103,113,114,122-1376%100
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=5G
#SBATCH --time=36:00:00
#SBATCH --output=logs/dlpfc_astro_res_fusionWTS_%A_%a.out
#SBATCH --error=logs/dlpfc_astro_res_fusionWTS_%A_%a.err
#SBATCH --exclude=compute-127,compute-145,compute-175,compute-158,compute-099

# set -euo pipefail

module purge
module load conda_R/4.5.x
module load plink2
module list
# module load plink2
conda activate aanri

# ln -s ./ output

covariates="/dcs04/lieber/hwanglab/Arun/snRNA_aanri/3_analysis/DLPFC/Astrocyte/residuals/TWAS/covariates/DLPFC_Astrocyte_residual_covariates.txt"
exprn="/dcs04/lieber/hwanglab/Arun/snRNA_aanri/3_analysis/DLPFC/Astrocyte/residuals/TWAS/expression/DLPFC_Astrocyte_residual_expression.bed"
output="/dcs04/lieber/hwanglab/Arun/snRNA_aanri/3_analysis/DLPFC/Astrocyte/residuals/TWAS"
rpath="/jhpce/shared/community/core/conda_R/4.5/R/bin/Rscript" # R executable after loading conda_R inside the slurm system
wtsrpath="/dcs04/lieber/hwanglab/Arun/miRQTL_pipeline/source/software/fusion_twas-master/FUSION.compute_weights.R"
gcta="/dcs04/lieber/hwanglab/Arun/snRNA_aanri/3_analysis/DLPFC/Astrocyte/residuals/TWAS/GCTA"
#gctaPath="/dcs04/lieber/hwanglab/Arun/miRQTL_pipeline/source/software/gcta-1.95.2-linux-x86_64/gcta"
plinkPath="plink2"
gimmaPath="/dcs04/lieber/hwanglab/Arun/software/gemma"

python /dcs04/lieber/hwanglab/Arun/miRQTL_pipeline/source/software/miRQTL/TWAS/twas_pipeline.py \
    fusionWTS \
    --covar "$covariates" \
    --output "$output" \
    --heritability "$gcta" \
    --exprn "$exprn" \
    --wtsRPath "$wtsrpath" \
    --pathR "$rpath" \
    --plinkPath "$plinkPath" \
    --gimmaPath "$gimmaPath" \
    --models lasso,top1,enet,blup,bslmm \
    --isArrary

## failed or stuck runs #SBATCH --array=40,41,25-38,103,113,114,122-1376%100
