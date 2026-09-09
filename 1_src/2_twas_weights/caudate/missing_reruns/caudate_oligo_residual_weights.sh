#!/bin/bash
#SBATCH --job-name=oligo_res_fusionWTS
#SBATCH --array=62 #0-1613%100 # IMPORTANT DON'T FORGET TO CHANGE
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=5G
#SBATCH --time=36:00:00
#SBATCH --output=logs/caudate_oligo_res_fusionWTS_%A_%a.out
#SBATCH --error=logs/caudate_oligo_res_fusionWTS_%A_%a.err
#SBATCH --exclude=compute-127,compute-145,compute-175,compute-158,compute-099,compute-142,compute-148

# set -euo pipefail

module purge
module load conda_R/4.5.x
module load plink2
module list
# module load plink2
conda activate aanri

# ln -s ./ output
#celltypes=("Astrocyte" "Excitatory_neuron" "Inhibitory_neuron" "Microglia" "OPC" "Oligodendrocyte") 


region="caudate"
celltype="Oligodendrocyte"
analysis="residuals"

covariates="/dcs04/lieber/hwanglab/Arun/snRNA_aanri/3_analysis/${region}/${celltype}/${analysis}/TWAS/covariates/${region}_${celltype}_residual_covariates.txt"
exprn="/dcs04/lieber/hwanglab/Arun/snRNA_aanri/3_analysis/${region}/${celltype}/${analysis}/TWAS/expression/${region}_${celltype}_residual_expression.bed"
output="/dcs04/lieber/hwanglab/Arun/snRNA_aanri/3_analysis/${region}/${celltype}/${analysis}/TWAS"
rpath="/jhpce/shared/community/core/conda_R/4.5/R/bin/Rscript" # R executable after loading conda_R inside the slurm system
wtsrpath="/dcs04/lieber/hwanglab/Arun/miRQTL_pipeline/source/software/fusion_twas-master/FUSION.compute_weights.R"
gcta="/dcs04/lieber/hwanglab/Arun/snRNA_aanri/3_analysis/${region}/${celltype}/${analysis}/TWAS/GCTA"
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

