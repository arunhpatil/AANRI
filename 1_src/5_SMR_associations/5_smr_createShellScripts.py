#!/usr/bin/env python3

import argparse
from pathlib import Path
import sys

# ============================================================
# Configuration
# ============================================================

#REGION = sys.argv[1] 
#ANALYSIS = "residuals"

BASE_DIR = Path(
    "/dcs04/lieber/hwanglab/Arun/snRNA_aanri/1_source/batch_executables/5_SMR_associations"
)

TWAS_PIPELINE = Path(
    "/dcs04/lieber/hwanglab/Arun/miRQTL_pipeline/"
    "source/software/miRQTL/TWAS/twas_pipeline.py"
)


# ============================================================
# Generate SLURM script
# ============================================================

def create_script(REGION, celltype, mem):

    celltype_dir = BASE_DIR / REGION 
    #celltype_dir = BASE_DIR / REGION / celltype 

    # Create celltype directory and logs directory
    celltype_dir.mkdir(parents=True, exist_ok=True)

    logs_dir = celltype_dir / "logs"
    logs_dir.mkdir(exist_ok=True)

    smr_logs_dir = celltype_dir / "smr_logs"
    smr_logs_dir.mkdir(exist_ok=True)
    # Script name
    script_path1 = celltype_dir / f"{celltype}_SMR_step1.sh"
    script_path2 = celltype_dir / f"{celltype}_SMR_step2.sh"
    script_path3 = celltype_dir / f"{celltype}_SMR_step3.sh"

    # --------------------------------------------------------
    # Paths
    # --------------------------------------------------------
    analysis = "residuals"
    #TWAS = celltype_dir / "TWAS"
    TWAS=f"/dcs04/lieber/hwanglab/Arun/snRNA_aanri/3_analysis/{REGION}/{celltype}/{analysis}/TWAS"

    # --------------------------------------------------------
    # SLURM script
    # --------------------------------------------------------

    script1 = f"""#!/bin/bash
#SBATCH --job-name={REGION}_{celltype}_smr
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=10 ## Low overhead needed for a simple pandas compilation
#SBATCH --mem=90G #5G # 90G # If running to create parquet files, we need 90Gb if not then 5Gb is sufficient
#SBATCH --time=03:00:00
#SBATCH --output={logs_dir}/{celltype}_%j.out 
#SBATCH --error={logs_dir}/{celltype}_%j.err
#SBATCH --exclude=compute-127,compute-145,compute-175,compute-158,compute-099,compute-142,compute-148,compute-095,compute-144,compute-166,compute-053,compute-054,compute-145,compute-058

# set -euo pipefail
log_message() {{
    echo "$(date '+%Y-%m-%d %H:%M:%S') - $1"
}}

log_message "**** Job starts ****"

module purge
module load conda_R/4.5.x
module list
conda activate aanri
ml plink2

set -e

region="{REGION}"
celltype="{celltype}"
TWAS="{TWAS}"

plinkPath="plink2"
gwas="/dcs04/lieber/hwanglab/Arun/snRNA_aanri/0_annotation/gwas_ss/"
metadata="/dcs04/lieber/hwanglab/Arun/snRNA_aanri/0_annotation/combined_summaryStats_metadata.txt" # Metadata regarding the summary statistics, the column names are mandatory
genotype="/dcs04/lieber/hwanglab/Arun/snRNA_aanri/2_data/aanri_072826_genotype_qc"
coords="/dcs04/lieber/hwanglab/Arun/snRNA_aanri/0_annotation/gene_symbol_QTLCords_strand.csv"
nominalQTL="/dcs04/lieber/hwanglab/Arun/snRNA_aanri/3_analysis/{REGION}/{celltype}/residuals/eQTL/{REGION}_{celltype}_.lan.tsqtl.nominal.txt"
smrPath="/dcs04/lieber/hwanglab/Arun/miRQTL_pipeline/source/software/smr-1.4.1-linux-x86_64/smr"

python "{TWAS_PIPELINE}" SMR \\
    --mode prep \\
    --output "$TWAS" \\
    --metadata "$metadata" \\
    --outPrefix "{REGION}_{celltype}" \\
    --plinkPath "$plinkPath" \\
    --genotype "$genotype" \\
    --coordinates "$coords" \\
    --gwasDir "$gwas" \\
    --smrPath "$smrPath" \\
    --threads  ${{SLURM_CPUS_PER_TASK}} \\
    --reversedAF \\
    --qtlPath "$nominalQTL"

"""


    script2 = f"""#!/bin/bash
#SBATCH --job-name={REGION}_{celltype}_smr_association
#SBATCH --array=0-65 
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=50G #5G # 90G # If running to create parquet files, we need 90Gb if not then 5Gb is sufficient
#SBATCH --time=03:00:00
#SBATCH --output={smr_logs_dir}/{celltype}_smrAsso_%A_%a.out 
#SBATCH --error={smr_logs_dir}/{celltype}_smrAsso_%A_%a.err
#SBATCH --exclude=compute-127,compute-145,compute-175,compute-158,compute-099,compute-142,compute-148,compute-095,compute-144,compute-166,compute-053,compute-054,compute-145,compute-058


# set -euo pipefail
log_message() {{
    echo "$(date '+%Y-%m-%d %H:%M:%S') - $1"
}}

log_message "**** Job starts ****"

module purge
module load conda_R/4.5.x
module list
conda activate aanri
ml plink2

set -e

region="{REGION}"
celltype="{celltype}"
TWAS="{TWAS}"

plinkPath="plink2"
gwas="/dcs04/lieber/hwanglab/Arun/snRNA_aanri/0_annotation/gwas_ss/"
metadata="/dcs04/lieber/hwanglab/Arun/snRNA_aanri/0_annotation/combined_summaryStats_metadata.txt" # Metadata regarding the summary statistics, the column names are mandatory
genotype="/dcs04/lieber/hwanglab/Arun/snRNA_aanri/2_data/aanri_072826_genotype_qc"
coords="/dcs04/lieber/hwanglab/Arun/snRNA_aanri/0_annotation/gene_symbol_QTLCords_strand.csv"
nominalQTL="/dcs04/lieber/hwanglab/Arun/snRNA_aanri/3_analysis/{REGION}/{celltype}/residuals/eQTL/{REGION}_{celltype}_.lan.tsqtl.nominal.txt"
smrPath="/dcs04/lieber/hwanglab/Arun/miRQTL_pipeline/source/software/smr-1.4.1-linux-x86_64/smr"
                                                                        
python "{TWAS_PIPELINE}" SMR \\
    --mode association \\
    --output "$TWAS" \\
    --metadata "$metadata" \\
    --outPrefix "{REGION}_{celltype}" \\
    --plinkPath "$plinkPath" \\
    --genotype "$genotype" \\
    --coordinates "$coords" \\
    --gwasDir "$gwas" \\
    --smrPath "$smrPath" \\
    --threads 2 \\
    --qtlPath "$nominalQTL" \\
    --is_array

"""


    script3 = f"""#!/bin/bash
#SBATCH --job-name={REGION}_{celltype}_smrSummary
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=10 ## Low overhead needed for a simple pandas compilation
#SBATCH --mem=2G #5G # 90G # If running to create parquet files, we need 90Gb if not then 5Gb is sufficient
#SBATCH --time=01:00:00
#SBATCH --output={logs_dir}/{celltype}_smrSummaize_%j.out
#SBATCH --error={logs_dir}/{celltype}_smrSummarize_%j.err
#SBATCH --exclude=compute-127,compute-145,compute-175,compute-158,compute-099,compute-142,compute-148,compute-095,compute-144,compute-166,compute-053,compute-054,compute-145,compute-058


# set -euo pipefail
log_message() {{
    echo "$(date '+%Y-%m-%d %H:%M:%S') - $1"
}}

log_message "**** Job starts ****"

module purge
module load conda_R/4.5.x
module list
conda activate aanri
ml plink2

set -e

region="{REGION}"
celltype="{celltype}"
TWAS="{TWAS}"

plinkPath="plink2"
gwas="/dcs04/lieber/hwanglab/Arun/snRNA_aanri/0_annotation/gwas_ss/"
metadata="/dcs04/lieber/hwanglab/Arun/snRNA_aanri/0_annotation/combined_summaryStats_metadata.txt" # Metadata regarding the summary statistics, the column names are mandatory
genotype="/dcs04/lieber/hwanglab/Arun/snRNA_aanri/2_data/aanri_072826_genotype_qc"
coords="/dcs04/lieber/hwanglab/Arun/snRNA_aanri/0_annotation/gene_symbol_QTLCords_strand.csv"
nominalQTL="/dcs04/lieber/hwanglab/Arun/snRNA_aanri/3_analysis/{REGION}/{celltype}/residuals/eQTL/{REGION}_{celltype}_.lan.tsqtl.nominal.txt"
smrPath="/dcs04/lieber/hwanglab/Arun/miRQTL_pipeline/source/software/smr-1.4.1-linux-x86_64/smr"

python "{TWAS_PIPELINE}" SMR \\
    --mode summarize \\
    --output "$TWAS" \\
    --metadata "$metadata" \\
    --outPrefix "{REGION}_{celltype}" \\
    --plinkPath "$plinkPath" \\
    --genotype "$genotype" \\
    --coordinates "$coords" \\
    --gwasDir "$gwas" \\
    --smrPath "$smrPath" \\
    --threads  ${{SLURM_CPUS_PER_TASK}} \\
    --qtlPath "$nominalQTL"

"""

    # Write script
    script_path1.write_text(script1)
    script_path2.write_text(script2)
    script_path3.write_text(script3)

    # Make executable
    script_path1.chmod(0o755)
    script_path2.chmod(0o755)
    script_path3.chmod(0o755)

    #print(f"Created: {script_path}")
    #print(f"Logs:    {logs_dir}")


# ============================================================
# Main
# ============================================================

def main():
    brainRegions=["DLPFC", "caudate", "hippo"]
    celltypes=["Astrocyte","Excitatory_neuron","Inhibitory_neuron","Microglia","OPC","Oligodendrocyte"]
    #brainRegions=["DLPFC"]
    #celltypes=["Astrocyte"]
    print(brainRegions, celltypes)
    for region in brainRegions:
        for celltype in celltypes:
            if region == "caudate" and celltype == "Excitatory_neuron":
                continue
            create_script(region, celltype, 16)

if __name__ == "__main__":
    main()
