#!/bin/bash
#SBATCH --job-name=create_expression_covs
#SBATCH --cpus-per-task=2
#SBATCH --time=8:00:00
#SBATCH --output=0_step_logs/prepare_%A.out
#SBATCH --error=0_step_logs/prepare_%A.err
#SBATCH --mem=16G
#SBATCH --exclude=compute-127,compute-145

log_message() {
    echo "$(date '+%Y-%m-%d %H:%M:%S') - $1"
}

log_message "**** Job starts ****"

echo "$(date '+%F %T') - Job started"
echo "Node: ${SLURMD_NODENAME}"
echo "CPUs: ${SLURM_CPUS_PER_TASK}"


module purge
module load conda_R/4.5.x
module list
module load plink2
conda activate aanri 

adata_in="/dcs04/lieber/hwanglab/Arun/snRNA_aanri/2_data/aanri_adata_snRNA_072826.h5ad"
genotype="/dcs04/lieber/hwanglab/Arun/snRNA_aanri/2_data/aanri_072826_genotype_qc"
workdir="/dcs04/lieber/hwanglab/Arun/snRNA_aanri/3_analysis/"
#brainRegion="DLPFC"
#celltype="Astrocyte"

#mkdir -p "$output"

script="/dcs04/lieber/hwanglab/Arun/snRNA_aanri/1_source/1_createTWAS_ingredients.py"

brainRegions=("DLPFC" "caudate" "hippo")
celltypes=("Astrocyte" "Excitatory_neuron" "Inhibitory_neuron" "Microglia" "OPC" "Oligodendrocyte") 
#celltypes=("Astrocyte" "Choroid_plexus" "Endothelial" "Ependymal" "Excitatory_neuron" "Inhibitory_neuron" "Lymphoid" "Microglia" "OPC" "Oligodendrocyte" "Vascular_stromal")
# Astrocyte, Oligodendrocyte, Microglia, OPC, In/Excitatory neurons.

# ------------------------------------------------------------
# Loop through brain regions and cell types
# ------------------------------------------------------------
for brainRegion in "${brainRegions[@]}"; do

    for celltype in "${celltypes[@]}"; do
	# Skip this specific combination
        if [[ "$brainRegion" == "caudate" && "$celltype" == "Excitatory_neuron" ]]; then # This is avoided because of the error detailed below
            echo "Skipping: ${brainRegion} / ${celltype}"
            continue
        fi

        echo "======================================================"
        echo "Running:"
        echo "  Brain region : ${brainRegion}"
        echo "  Cell type    : ${celltype}"
        echo "======================================================"

        python "$script" \
            --adata "$adata_in" \
            --genotype "$genotype" \
            --workdir "$workdir" \
            --Region "$brainRegion" \
            --cellType "$celltype"

        echo "Finished: ${brainRegion} / ${celltype}"
        echo

    done

done


# When I ran manually, I found that the error lies in caudate region and Excitatory neuros:
# >>> region
# 'caudate'
# >>> required_ct
# 'Excitatory_neuron'
# adata[(adata.obs["Region"] == region) & (adata.obs["Celltype_aggregated"] == required_ct) & (
# adata.obs["psbulk_cells"] >= 10)].shape
# (0, 34835)
# >>> adata[(adata.obs["Region"] == region) & (adata.obs["Celltype_aggregated"] == required_ct) & (adata.obs["psbulk_counts"] >= 1000)].shape
# (0, 34835)
