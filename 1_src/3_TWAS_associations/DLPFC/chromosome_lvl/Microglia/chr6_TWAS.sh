#!/bin/bash
#SBATCH --job-name=Microglia_chr6
#SBATCH --array=33,37,53,55,57 # 0-64 # %32 # 100 GWAS files (Index 0-99), max 20 running at once
#SBATCH --output=/dcs04/lieber/hwanglab/Arun/snRNA_aanri/1_source/batch_executables/3_TWAS_associations/DLPFC/chromosome_lvl/Microglia/logs/chr6_%A_%a.out
#SBATCH --error=/dcs04/lieber/hwanglab/Arun/snRNA_aanri/1_source/batch_executables/3_TWAS_associations/DLPFC/chromosome_lvl/Microglia/logs/chr6_%A_%a.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=35G
#SBATCH --time=36:00:00
#SBATCH --exclude=compute-127,compute-145,compute-175,compute-158,compute-099,compute-142,compute-148,compute-095,compute-144,compute-166,compute-053,compute-054,compute-058,compute-111

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

set -e

# ============================================================
# Analysis settings
# ============================================================

region="DLPFC"
celltype="Microglia"
chromo=6

# ============================================================
# FUSION
# ============================================================

FUSION_DIR="/dcs04/lieber/hwanglab/Arun/miRQTL_pipeline/source/software/fusion_twas-master"
assrpath="${FUSION_DIR}/FUSION.assoc_test.R"

cd "$FUSION_DIR"

# ============================================================
# TWAS output
# ============================================================

TWAS="/dcs04/lieber/hwanglab/Arun/snRNA_aanri/3_analysis/DLPFC/Microglia/residuals/TWAS"


# ============================================================
# Run TWAS
# ============================================================


output="significant_associations.txt" # The TWAS fusion final results are saved here and significant_associations_clean.txt file. The _clean.txt doesn't contain the path to .Rmd weight files and well formated. 
rpath="/jhpce/shared/community/core/conda_R/4.5/R/bin/Rscript" # R executable after loading conda_R inside the slurm system
gwas="/dcs04/lieber/hwanglab/Arun/snRNA_aanri/0_annotation/gwas_ss_chrSpecific/chr6"
ldref="/dcs04/lieber/hwanglab/Arun/snRNA_aanri/0_annotation/LD_REF_hg38/1000G.AFR."

#python /dcs04/lieber/hwanglab/Arun/miRQTL_pipeline/source/software/miRQTL/TWAS/twas_pipeline.py fusionChrllAssoc \
python "/dcs04/lieber/hwanglab/Arun/miRQTL_pipeline/source/software/miRQTL/TWAS/twas_pipeline.py" fusionChrllAssoc \
           --gwasDir "$gwas" \
           --assRPath $assrpath \
           --chromosome "$chromo" \
           --output "$TWAS" \
           --pathR "$rpath" \
           --ldref "$ldref" \
           --is_array

