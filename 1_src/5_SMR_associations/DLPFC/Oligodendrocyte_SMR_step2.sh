#!/bin/bash
#SBATCH --job-name=DLPFC_Oligodendrocyte_smr_association
#SBATCH --array=0-65 
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=50G #5G # 90G # If running to create parquet files, we need 90Gb if not then 5Gb is sufficient
#SBATCH --time=03:00:00
#SBATCH --output=/dcs04/lieber/hwanglab/Arun/snRNA_aanri/1_source/batch_executables/5_SMR_associations/DLPFC/smr_logs/Oligodendrocyte_smrAsso_%A_%a.out 
#SBATCH --error=/dcs04/lieber/hwanglab/Arun/snRNA_aanri/1_source/batch_executables/5_SMR_associations/DLPFC/smr_logs/Oligodendrocyte_smrAsso_%A_%a.err
#SBATCH --exclude=compute-127,compute-145,compute-175,compute-158,compute-099,compute-142,compute-148,compute-095,compute-144,compute-166,compute-053,compute-054,compute-145,compute-058


# set -euo pipefail
log_message() {
    echo "$(date '+%Y-%m-%d %H:%M:%S') - $1"
}

log_message "**** Job starts ****"

module purge
module load conda_R/4.5.x
module list
conda activate aanri
ml plink2

set -e

region="DLPFC"
celltype="Oligodendrocyte"
TWAS="/dcs04/lieber/hwanglab/Arun/snRNA_aanri/3_analysis/DLPFC/Oligodendrocyte/residuals/TWAS"

plinkPath="plink2"
gwas="/dcs04/lieber/hwanglab/Arun/snRNA_aanri/0_annotation/gwas_ss/"
metadata="/dcs04/lieber/hwanglab/Arun/snRNA_aanri/0_annotation/combined_summaryStats_metadata.txt" # Metadata regarding the summary statistics, the column names are mandatory
genotype="/dcs04/lieber/hwanglab/Arun/snRNA_aanri/2_data/aanri_072826_genotype_qc"
coords="/dcs04/lieber/hwanglab/Arun/snRNA_aanri/0_annotation/gene_symbol_QTLCords_strand.csv"
nominalQTL="/dcs04/lieber/hwanglab/Arun/snRNA_aanri/3_analysis/DLPFC/Oligodendrocyte/residuals/eQTL/DLPFC_Oligodendrocyte_.lan.tsqtl.nominal.txt"
smrPath="/dcs04/lieber/hwanglab/Arun/miRQTL_pipeline/source/software/smr-1.4.1-linux-x86_64/smr"
                                                                        
python "/dcs04/lieber/hwanglab/Arun/miRQTL_pipeline/source/software/miRQTL/TWAS/twas_pipeline.py" SMR \
    --mode association \
    --output "$TWAS" \
    --metadata "$metadata" \
    --outPrefix "DLPFC_Oligodendrocyte" \
    --plinkPath "$plinkPath" \
    --genotype "$genotype" \
    --coordinates "$coords" \
    --gwasDir "$gwas" \
    --smrPath "$smrPath" \
    --threads 2 \
    --qtlPath "$nominalQTL" \
    --is_array

