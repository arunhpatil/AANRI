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
    "/dcs04/lieber/hwanglab/Arun/snRNA_aanri/1_source/batch_executables/3_TWAS_associations"
)

FUSION_DIR = Path(
    "/dcs04/lieber/hwanglab/Arun/miRQTL_pipeline/"
    "source/software/fusion_twas-master"
)

TWAS_PIPELINE = Path(
    "/dcs04/lieber/hwanglab/Arun/miRQTL_pipeline/"
    "source/software/miRQTL/TWAS/twas_pipeline.py"
)


# ============================================================
# Generate SLURM script
# ============================================================

def create_script(REGION, celltype, mem):

    celltype_dir = BASE_DIR / REGION / "chromosome_lvl" / celltype 

    # Create celltype directory and logs directory
    celltype_dir.mkdir(parents=True, exist_ok=True)

    logs_dir = celltype_dir / "logs"
    logs_dir.mkdir(exist_ok=True)

    # Script name
    script_path = celltype_dir / f"nonChr_TWAS.sh"

    # --------------------------------------------------------
    # Paths
    # --------------------------------------------------------
    analysis = "residuals"
    #TWAS = celltype_dir / "TWAS"
    TWAS=f"/dcs04/lieber/hwanglab/Arun/snRNA_aanri/3_analysis/{REGION}/{celltype}/{analysis}/TWAS"

    # --------------------------------------------------------
    # SLURM script
    # --------------------------------------------------------

    script = f"""#!/bin/bash
#SBATCH --job-name={celltype}_nonChr
#SBATCH --output={logs_dir}/nonChr_%j.out
#SBATCH --error={logs_dir}/nonChr_%j.err
#SBATCH --cpus-per-task=2
#SBATCH --mem={mem}G
#SBATCH --time=6:00:00
#SBATCH --exclude=compute-127,compute-145,compute-175,compute-158,compute-099,compute-142,compute-148,compute-095,compute-144,compute-166,compute-053,compute-054,compute-058,compute-111


# set -euo pipefail
log_message() {{
    echo "$(date '+%Y-%m-%d %H:%M:%S') - $1"
}}

module purge
module load conda_R/4.5.x
module list
conda activate aanri

set -e

# ============================================================
# Analysis settings
# ============================================================

region="{REGION}"
celltype="{celltype}"

# ============================================================
# FUSION
# ============================================================

FUSION_DIR="{FUSION_DIR}"
assrpath="${{FUSION_DIR}}/FUSION.assoc_test.R"

cd "$FUSION_DIR"

# ============================================================
# TWAS output
# ============================================================

TWAS="{TWAS}"


# ============================================================
# Run TWAS
# ============================================================


output="significant_associations.txt" # The TWAS fusion final results are saved here and significant_associations_clean.txt file. The _clean.txt doesn't contain the path to .Rmd weight files and well formated. 
rpath="/jhpce/shared/community/core/conda_R/4.5/R/bin/Rscript" # R executable after loading conda_R inside the slurm system
gwas="/dcs04/lieber/hwanglab/Arun/snRNA_aanri/0_annotation/gwas_ss_non_chrSpecific"
ldref="/dcs04/lieber/hwanglab/Arun/snRNA_aanri/0_annotation/LD_REF_hg38/1000G.AFR."

#python /dcs04/lieber/hwanglab/Arun/miRQTL_pipeline/source/software/miRQTL/TWAS/twas_pipeline.py fusionChrllAssoc \\
python "{TWAS_PIPELINE}" fusionAssoc \\
        --gwasDir "$gwas" \\
        --race AA \\
        --threads 6 \\
        --assRPath "$assrpath" \\
        --pathR "$rpath" \\
        --ldref "$ldref" \\
        --output "$TWAS"

"""
    # Write script
    script_path.write_text(script)

    # Make executable
    script_path.chmod(0o755)

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
            create_script(region, celltype, 120)

if __name__ == "__main__":
    main()
