#!/bin/bash
#SBATCH --job-name=dlpfc_astro_res_association
#SBATCH --array=0-65%30          # 100 GWAS files (Index 0-99), max 20 running at once
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=5
#SBATCH --mem=120G
#SBATCH --time=36:00:00
#SBATCH --output=logs/dlpfc_astro_res_association_%A_%a.out
#SBATCH --error=logs/dlpfc_astro_res_association_%A_%a.err
#SBATCH --exclude=compute-127,compute-145,compute-175,compute-158,compute-099,compute-142,compute-148

# set -euo pipefail
log_message() {
    echo "$(date '+%Y-%m-%d %H:%M:%S') - $1"
}

log_message "**** Job starts ****"
echo "User: ${USER}"
echo "Job id: ${SLURM_JOBID}"
echo "Task id: ${SLURM_ARRAY_TASK_ID}"
echo "Node name: ${SLURM_NODENAME}"


module purge
module load conda_R/4.5.x
module list
conda activate aanri

region="DLPFC"
celltype="Astrocyte"
analysis="residuals"

FUSION_DIR="/dcs04/lieber/hwanglab/Arun/miRQTL_pipeline/source/software/fusion_twas-master"
assrpath="${FUSION_DIR}/FUSION.assoc_test.R"

cd "$FUSION_DIR"

TWAS="/dcs04/lieber/hwanglab/Arun/snRNA_aanri/3_analysis/${region}/${celltype}/${analysis}/TWAS"
output="significant_associations.txt" # The TWAS fusion final results are saved here and significant_associations_clean.txt file. The _clean.txt doesn't contain the path to .Rmd weight files and well formated. 
rpath="/jhpce/shared/community/core/conda_R/4.5/R/bin/Rscript" # R executable after loading conda_R inside the slurm system
metadata="/dcs04/lieber/hwanglab/Arun/snRNA_aanri/0_annotation/combined_summaryStats_metadata.txt"
gwas="/dcs04/lieber/hwanglab/Arun/snRNA_aanri/0_annotation/gwas_ss/"
ldref="/dcs04/lieber/hwanglab/Arun/snRNA_aanri/0_annotation/LD_REF_hg38/1000G.AFR."
#assrpath="/dcs04/lieber/hwanglab/Arun/miRQTL_pipeline/source/software/fusion_twas-master/FUSION.assoc_test.R" 

python /dcs04/lieber/hwanglab/Arun/miRQTL_pipeline/source/software/miRQTL/TWAS/twas_pipeline.py fusionAssoc \
           --gwasDir "$gwas" \
           --race AA \
           --threads ${SLURM_CPUS_PER_TASK} \
	   --assRPath $assrpath \
           --output "$TWAS" \
           --pathR "$rpath" \
           --ldref "$ldref" \
           --is_array

