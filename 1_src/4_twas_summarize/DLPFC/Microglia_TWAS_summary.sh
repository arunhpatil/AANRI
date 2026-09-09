#!/bin/bash
#SBATCH --job-name=DLPFC_Microglia_twasSummary
#SBATCH --output=/dcs04/lieber/hwanglab/Arun/snRNA_aanri/1_source/batch_executables/4_twas_summarize/DLPFC/logs/Microglia_%j.out
#SBATCH --error=/dcs04/lieber/hwanglab/Arun/snRNA_aanri/1_source/batch_executables/4_twas_summarize/DLPFC/logs/Microglia_%j.err
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --time=2:00:00

# set -euo pipefail
log_message() {
    echo "$(date '+%Y-%m-%d %H:%M:%S') - $1"
}

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

# ============================================================
# TWAS output
# ============================================================

TWAS="/dcs04/lieber/hwanglab/Arun/snRNA_aanri/3_analysis/DLPFC/Microglia/residuals/TWAS"

# ============================================================
# Run TWAS Summary
# ============================================================

output="significant_associations.txt" # The TWAS fusion final results are saved here and significant_associations_clean.txt file. The _clean.txt doesn't contain the path to .Rmd weight files and well formated. 
metadata="/dcs04/lieber/hwanglab/Arun/snRNA_aanri/0_annotation/combined_summaryStats_metadata.txt"
TWAS="/dcs04/lieber/hwanglab/Arun/snRNA_aanri/3_analysis/DLPFC/Microglia/residuals/TWAS"

python "/dcs04/lieber/hwanglab/Arun/miRQTL_pipeline/source/software/miRQTL/TWAS/twas_pipeline.py" fusionChrllReport \
    --output "$TWAS" \
    --metadata "$metadata" \
    --outfile "$output"

