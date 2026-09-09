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
    "/dcs04/lieber/hwanglab/Arun/snRNA_aanri/1_source/batch_executables/4_twas_summarize"
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

    # Script name
    script_path = celltype_dir / f"{celltype}_TWAS_summary.sh"

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
#SBATCH --job-name={REGION}_{celltype}_twasSummary
#SBATCH --output={logs_dir}/{celltype}_%j.out
#SBATCH --error={logs_dir}/{celltype}_%j.err
#SBATCH --cpus-per-task=2
#SBATCH --mem={mem}G
#SBATCH --time=2:00:00

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
# TWAS output
# ============================================================

TWAS="{TWAS}"

# ============================================================
# Run TWAS Summary
# ============================================================

output="significant_associations.txt" # The TWAS fusion final results are saved here and significant_associations_clean.txt file. The _clean.txt doesn't contain the path to .Rmd weight files and well formated. 
metadata="/dcs04/lieber/hwanglab/Arun/snRNA_aanri/0_annotation/combined_summaryStats_metadata.txt"
TWAS="/dcs04/lieber/hwanglab/Arun/snRNA_aanri/3_analysis/{REGION}/{celltype}/residuals/TWAS"

python "{TWAS_PIPELINE}" fusionChrllReport \\
    --output "$TWAS" \\
    --metadata "$metadata" \\
    --outfile "$output"

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
            create_script(region, celltype, 16)

if __name__ == "__main__":
    main()
