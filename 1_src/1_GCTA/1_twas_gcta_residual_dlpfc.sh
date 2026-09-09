#!/bin/bash
#SBATCH --job-name=dlpfc_CellTypeGCTA_residuals
#SBATCH --cpus-per-task=2
#SBATCH --time=24:00:00
#SBATCH --output=logs/all_regXct_twas_gcta_residuals_dlpfc.out
#SBATCH --error=logs/all_regXct_twas_gcta_residuals_dlpfc.err
#SBATCH --mem=120G


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
rpath="/jhpce/shared/community/core/conda_R/4.5/R/bin/Rscript" # R executable after loading conda_R inside the slurm system
plinkPath="plink2"
gctaPath="/dcs04/lieber/hwanglab/Arun/miRQTL_pipeline/source/software/gcta-1.95.2-linux-x86_64/gcta"

#mkdir -p "$output"

brainRegions=("DLPFC")
celltypes=("Astrocyte" "Excitatory_neuron" "Inhibitory_neuron" "Microglia" "OPC" "Oligodendrocyte") 
#brainRegions=("hippo")
#celltypes=("OPC" "Oligodendrocyte") 
script="/dcs04/lieber/hwanglab/Arun/miRQTL_pipeline/source/software/miRQTL/TWAS/twas_pipeline.py"

for region in "${brainRegions[@]}"; do

    for celltype in "${celltypes[@]}"; do

        # Skip specific region/cell-type combination
        if [[ "$region" == "caudate" && "$celltype" == "Excitatory_neuron" ]]; then
            echo "Skipping: ${region} / ${celltype}"
            continue
        fi

	echo "======================================================"
        echo "Running GCTA:"
        echo "  Region   : ${region}"
        echo "  Cell type: ${celltype}"
        echo "  Node     : $(hostname)"
        echo "======================================================"

        exprn="/dcs04/lieber/hwanglab/Arun/snRNA_aanri/3_analysis/${region}/${celltype}/residuals/TWAS/expression/${region}_${celltype}_residual_expression.bed"

        output="/dcs04/lieber/hwanglab/Arun/snRNA_aanri/3_analysis/${region}/${celltype}/residuals/TWAS"

        covariates="/dcs04/lieber/hwanglab/Arun/snRNA_aanri/3_analysis/${region}/${celltype}/residuals/TWAS/covariates/${region}_${celltype}_residual_covariates.txt"

	python "$script" GCTA \
	    --exprn "$exprn" \
	    --genotype "$genotype" \
	    --covar "$covariates" \
	    --output "$output" \
	    --threads ${SLURM_CPUS_PER_TASK} \
	    --gctaPath "$gctaPath" \
	    --plinkPath "$plinkPath" \
	    --pathR "$rpath"

	status=$?

        if [[ $status -ne 0 ]]; then
            echo "ERROR: ${region} / ${celltype} failed with exit code ${status}"
        else
            echo "Finished: ${region} / ${celltype}"
        fi

	status=$?

        if [[ $status -ne 0 ]]; then
            echo "ERROR: ${region} / ${celltype} failed with exit code ${status}"
        else
            echo "Finished: ${region} / ${celltype}"
        fi

        echo
    done
done
