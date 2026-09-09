#!/bin/bash
#SBATCH --job-name=create_expression_covs
#SBATCH --cpus-per-task=2
#SBATCH --time=24:00:00
#SBATCH --output=0_step_logs/3_combineQTLs_%A.out
#SBATCH --error=0_step_logs/3_combineQTLs_%A.err
#SBATCH --mem=64G
#SBATCH --exclude=compute-127,compute-145,compute-175,compute-158,compute-099,compute-142,compute-148,compute-095,compute-144,compute-166

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

genotype="/dcs04/lieber/hwanglab/Arun/snRNA_aanri/2_data/aanri_072826_genotype_qc"
script="/dcs04/lieber/hwanglab/Arun/snRNA_aanri/1_source/3_combining_eQTL.py"

#brainRegions=("DLPFC" "caudate" "hippo")
brainRegions=("hippo")
celltypes=("Astrocyte" "Excitatory_neuron" "Inhibitory_neuron" "Microglia" "OPC" "Oligodendrocyte") 
covs=("base_covariates" "local_ancestry")

#celltypes=("Astrocyte" "Choroid_plexus" "Endothelial" "Ependymal" "Excitatory_neuron" "Inhibitory_neuron" "Lymphoid" "Microglia" "OPC" "Oligodendrocyte" "Vascular_stromal")
workdir="/dcs04/lieber/hwanglab/Arun/snRNA_aanri/3_analysis/caudate/Oligodendrocyte/residuals/"
indir="/dcs04/lieber/hwanglab/Arun/snRNA_aanri/raw_data_dongsan/eQTL_all/caudate/Oligodendrocyte/local_ancestry/"
python "$script" \
	--indir "$indir" \
	--genotype "$genotype" \
	--workdir "$workdir" \
	--Region "caudate" \
	--cellType "Oligodendrocyte" \
	--localAncestry

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

	workdir="/dcs04/lieber/hwanglab/Arun/snRNA_aanri/3_analysis/${brainRegion}/${celltype}/residuals/"
	for covdir in "${covs[@]}"; do
		indir="/dcs04/lieber/hwanglab/Arun/snRNA_aanri/raw_data_dongsan/eQTL_all/${brainRegion}/${celltype}/${covdir}/"
		extra_args=""
		if [[ "$covdir" == "local_ancestry" ]]; then
			extra_args="--localAncestry"
	    	fi
		python "$script" \
		    --indir "$indir" \
		    --genotype "$genotype" \
		    --workdir "$workdir" \
		    --Region "$brainRegion" \
		    --cellType "$celltype" \
		    $extra_args

      	    echo "Finished: ${brainRegion} / ${celltype}"
	    echo
	done
    done
done
