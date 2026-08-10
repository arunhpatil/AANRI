# AANRI
snRNA sequencing QTL and TWAS analysis pipeline - African Ancestry Neuroscience Research Initiative (AANRI). 

All the codes contained in this repository are those that can be found under `1_source` directory (see tree below). 

The draft directory structure is as shown below:

```
/dcs04/lieber/hwanglab/Arun/snRNA_aanri/
.
├── 0_annotation
│   ├── combined_summaryStats_metadata.txt
│   ├── gwas_ss
│   │   ├── bip2024_afr_no23andMe.txt
│   │   ├── GCST90475495_Cognitive_Neurodegenerative.h.FUSION.txt
│   │   ├── GCST90476667_Brain_Cerebrovascular.h.FUSION.txt
│   │   ├── GCST90476683_MoodDisorders.h.FUSION.txt
│   │   ├── GCST90476686_Cognitive_Neurodegenerative.h.FUSION.txt
│   │   ├── GCST90476865_Brain_Cerebrovascular.h.FUSION.txt
│   │   ├── GCST90476879_MoodDisorders.h.FUSION.txt
│   │   ├── GCST90476883_Cognitive_Neurodegenerative.h.FUSION.txt
│   │   ├── GCST90476903_Brain_Cerebrovascular.h.FUSION.txt
│   │   ├── GCST90476918_MoodDisorders.h.FUSION.txt
│   │   ├── GCST90476921_MoodDisorders.h.FUSION.txt
│   │   ├── GCST90476928_Cognitive_Neurodegenerative.h.FUSION.txt
│   │   ├── GCST90476930_SeizureDisorders.h.FUSION.txt
│   │   ├── GCST90476935_Brain_Cerebrovascular.h.FUSION.txt
│   │   ├── GCST90476994_Brain_Cerebrovascular.h.FUSION.txt
│   │   ├── GCST90477003_MoodDisorders.h.FUSION.txt
│   │   ├── GCST90477006_Cognitive_Neurodegenerative.h.FUSION.txt
│   │   ├── GCST90477283_Brain_Cerebrovascular.h.FUSION.txt
│   │   ├── GCST90477285_Brain_Cerebrovascular.h.FUSION.txt
│   │   ├── GCST90477481_SleepDisorders.h.FUSION.txt
│   │   ├── GCST90477486_SleepDisorders.h.FUSION.txt
│   │   ├── GCST90477511_Movement_Motor.h.FUSION.txt
│   │   ├── GCST90477513_Movement_Motor.h.FUSION.txt
│   │   ├── GCST90477517_SpinalCord_Nerves.h.FUSION.txt
│   │   ├── GCST90477521_SpinalCord_Nerves.h.FUSION.txt
│   │   ├── GCST90477532_OtherNeurological.h.FUSION.txt
│   │   ├── GCST90477544_OtherNeurological.h.FUSION.txt
│   │   ├── GCST90477550_SeizureDisorders.h.FUSION.txt
│   │   ├── GCST90477553_SeizureDisorders.h.FUSION.txt
│   │   ├── GCST90477559_SeizureDisorders.h.FUSION.txt
│   │   ├── GCST90477567_Brain_Cerebrovascular.h.FUSION.txt
│   │   ├── GCST90477574_Movement_Motor.h.FUSION.txt
│   │   ├── GCST90477576_Movement_Motor.h.FUSION.txt
│   │   ├── GCST90477580_Movement_Motor.h.FUSION.txt
│   │   ├── GCST90477588_SpinalCord_Nerves.h.FUSION.txt
│   │   ├── GCST90477593_SpinalCord_Nerves.h.FUSION.txt
│   │   ├── GCST90477600_OtherNeurological.h.FUSION.txt
│   │   ├── GCST90477603_OtherNeurological.h.FUSION.txt
│   │   ├── GCST90477606_SpinalCord_Nerves.h.FUSION.txt
│   │   ├── GCST90477812_OtherNeurological.h.FUSION.txt
│   │   ├── GCST90477982_Brain_Cerebrovascular.h.FUSION.txt
│   │   ├── GCST90477990_Brain_Cerebrovascular.h.FUSION.txt
│   │   ├── GCST90477993_Brain_Cerebrovascular.h.FUSION.txt
│   │   ├── GCST90478777_OtherNeurological.h.FUSION.txt
│   │   ├── GCST90479177_OtherNeurological.h.FUSION.txt
│   │   ├── GCST90479218_Brain_Cerebrovascular.h.FUSION.txt
│   │   ├── GCST90479382_Brain_Cerebrovascular.h.FUSION.txt
│   │   ├── GCST90479398_MoodDisorders.h.FUSION.txt
│   │   ├── GCST90479401_MoodDisorders.h.FUSION.txt
│   │   ├── GCST90479403_PsychoticDisorders.h.FUSION.txt
│   │   ├── GCST90479405_Cognitive_Neurodegenerative.h.FUSION.txt
│   │   ├── GCST90479439_SleepDisorders.h.FUSION.txt
│   │   ├── GCST90479446_SleepDisorders.h.FUSION.txt
│   │   ├── GCST90479450_SleepDisorders.h.FUSION.txt
│   │   ├── GCST90479457_SubstanceUse.h.FUSION.txt
│   │   ├── GCST90479465_SubstanceUse.h.FUSION.txt
│   │   ├── GCST90479468_SubstanceUse.h.FUSION.txt
│   │   ├── mdd2023diverse_AFR_clinicalMD_Neff.txt
│   │   ├── mdd2023diverse_AFR_Neff.txt
│   │   ├── mdd2023diverse_AFR_wto_UKB_clinicalMD_Neff.txt
│   │   ├── mdd2023diverse_AFR_wto_UKB_Neff.txt
│   │   ├── PGC3_SCZ_wave3.afram.stats.txt
│   │   ├── pgc_alcdep.afr_unrel_genotyped.aug2018_release.txt
│   │   ├── pts_aam_freeze2_overall.results.txt
│   │   ├── SCZ_AFR_Schizophrenia_FUSION.txt
│   │   └── SORTED_PTSD_AA7_ALL_study_specific_PCs1.txt
│   └── LD_REF_hg38
│       ├── 1000G.AFR.10.bed
│       ├── 1000G.AFR.10.bim
│       ├── 1000G.AFR.10.fam
│       ├── 1000G.AFR.10.log
│       ├── 1000G.AFR.10.nosex
│       ├── 1000G.AFR.11.bed
│       ├── 1000G.AFR.11.bim
│       ├── 1000G.AFR.11.fam
│       ├── 1000G.AFR.11.log
│       ├── 1000G.AFR.11.nosex
│       ├── 1000G.AFR.12.bed
│       ├── 1000G.AFR.12.bim
│       ├── 1000G.AFR.12.fam
│       ├── 1000G.AFR.12.log
│       ├── 1000G.AFR.12.nosex
│       ├── 1000G.AFR.13.bed
│       ├── 1000G.AFR.13.bim
│       ├── 1000G.AFR.13.fam
│       ├── 1000G.AFR.13.log
│       ├── 1000G.AFR.13.nosex
│       ├── 1000G.AFR.14.bed
│       ├── 1000G.AFR.14.bim
│       ├── 1000G.AFR.14.fam
│       ├── 1000G.AFR.14.log
│       ├── 1000G.AFR.14.nosex
│       ├── 1000G.AFR.15.bed
│       ├── 1000G.AFR.15.bim
│       ├── 1000G.AFR.15.fam
│       ├── 1000G.AFR.15.log
│       ├── 1000G.AFR.15.nosex
│       ├── 1000G.AFR.16.bed
│       ├── 1000G.AFR.16.bim
│       ├── 1000G.AFR.16.fam
│       ├── 1000G.AFR.16.log
│       ├── 1000G.AFR.16.nosex
│       ├── 1000G.AFR.17.bed
│       ├── 1000G.AFR.17.bim
│       ├── 1000G.AFR.17.fam
│       ├── 1000G.AFR.17.log
│       ├── 1000G.AFR.17.nosex
│       ├── 1000G.AFR.18.bed
│       ├── 1000G.AFR.18.bim
│       ├── 1000G.AFR.18.fam
│       ├── 1000G.AFR.18.log
│       ├── 1000G.AFR.18.nosex
│       ├── 1000G.AFR.19.bed
│       ├── 1000G.AFR.19.bim
│       ├── 1000G.AFR.19.fam
│       ├── 1000G.AFR.19.log
│       ├── 1000G.AFR.19.nosex
│       ├── 1000G.AFR.1.bed
│       ├── 1000G.AFR.1.bim
│       ├── 1000G.AFR.1.fam
│       ├── 1000G.AFR.1.log
│       ├── 1000G.AFR.1.nosex
│       ├── 1000G.AFR.20.bed
│       ├── 1000G.AFR.20.bim
│       ├── 1000G.AFR.20.fam
│       ├── 1000G.AFR.20.log
│       ├── 1000G.AFR.20.nosex
│       ├── 1000G.AFR.21.bed
│       ├── 1000G.AFR.21.bim
│       ├── 1000G.AFR.21.fam
│       ├── 1000G.AFR.21.log
│       ├── 1000G.AFR.21.nosex
│       ├── 1000G.AFR.22.bed
│       ├── 1000G.AFR.22.bim
│       ├── 1000G.AFR.22.fam
│       ├── 1000G.AFR.22.log
│       ├── 1000G.AFR.22.nosex
│       ├── 1000G.AFR.2.bed
│       ├── 1000G.AFR.2.bim
│       ├── 1000G.AFR.2.fam
│       ├── 1000G.AFR.2.log
│       ├── 1000G.AFR.2.nosex
│       ├── 1000G.AFR.3.bed
│       ├── 1000G.AFR.3.bim
│       ├── 1000G.AFR.3.fam
│       ├── 1000G.AFR.3.log
│       ├── 1000G.AFR.3.nosex
│       ├── 1000G.AFR.4.bed
│       ├── 1000G.AFR.4.bim
│       ├── 1000G.AFR.4.fam
│       ├── 1000G.AFR.4.log
│       ├── 1000G.AFR.4.nosex
│       ├── 1000G.AFR.5.bed
│       ├── 1000G.AFR.5.bim
│       ├── 1000G.AFR.5.fam
│       ├── 1000G.AFR.5.log
│       ├── 1000G.AFR.5.nosex
│       ├── 1000G.AFR.6.bed
│       ├── 1000G.AFR.6.bim
│       ├── 1000G.AFR.6.fam
│       ├── 1000G.AFR.6.log
│       ├── 1000G.AFR.6.nosex
│       ├── 1000G.AFR.7.bed
│       ├── 1000G.AFR.7.bim
│       ├── 1000G.AFR.7.fam
│       ├── 1000G.AFR.7.log
│       ├── 1000G.AFR.7.nosex
│       ├── 1000G.AFR.8.bed
│       ├── 1000G.AFR.8.bim
│       ├── 1000G.AFR.8.fam
│       ├── 1000G.AFR.8.log
│       ├── 1000G.AFR.8.nosex
│       ├── 1000G.AFR.9.bed
│       ├── 1000G.AFR.9.bim
│       ├── 1000G.AFR.9.fam
│       ├── 1000G.AFR.9.log
│       ├── 1000G.AFR.9.nosex
│       ├── afr_variants.txt
│       └── old_afr_variants.txt
├── 1_source
│   ├── 1_createTWAS_ingredients.py
│   ├── archived2deletion
│   │   └── 1_createTWAS_ingredients.py
│   ├── batch_executables
│   │   ├── 0_preparatory
│   │   ├── 1_GCTA
│   │   ├── 2_twas_weights
│   │   ├── 3_TWAS_associations
│   │   ├── 4_twas_summarize
│   │   ├── 5_SMR_associations
│   │   ├── 6_colocABF
│   │   ├── logs
│   │   ├── output -> ./
│   │   ├── temp_createDirStructure.sh
│   │   └── x_backupCodes
│   ├── commands.sh
│   └── read_h5ad.py # Completely my reference for understanding the dataset and testing the preparatory pipeline
├── 2_data
│   ├── aanri_072826_genotype_qc.bed
│   ├── aanri_072826_genotype_qc.bim
│   ├── aanri_072826_genotype_qc.fam
│   ├── aanri_072826_genotype_qc.log
│   ├── aanri_adata_snRNA_072826.h5ad
│   └── pre-processing_genotype
│       ├── commands_updating_genotype.sh
│       ├── dbsnp_update_names.canonical.txt
│       ├── dbsnp_update_names_filteredSubset.txt
│       ├── genotype_FY23_FY24_FY25_MAF05.bed
│       ├── genotype_FY23_FY24_FY25_MAF05.bim
│       ├── genotype_FY23_FY24_FY25_MAF05.fam
│       ├── genotype_FY23_FY24_FY25_MAF05.log
│       ├── libd_FY23_FY24_FY25_MAF05_variantList.txt
│       └── update_IIDs.txt
├── 3_analysis
│   ├── caudate
│   │   ├── Astrocyte
│   │   ├── Excitatory_neuron
│   │   ├── Inhibitory_neuron
│   │   ├── Microglia
│   │   ├── Oligodendrocyte
│   │   └── OPC
│   ├── DLPFC
│   │   ├── Astrocyte
│   │   ├── Excitatory_neuron
│   │   ├── Inhibitory_neuron
│   │   ├── Microglia
│   │   ├── Oligodendrocyte
│   │   └── OPC
│   └── hippo
│       ├── Astrocyte
│       ├── Excitatory_neuron
│       ├── Inhibitory_neuron
│       ├── Microglia
│       ├── Oligodendrocyte
│       └── OPC
├── raw_data_dongsan
│   ├── dnum_dnumS_brnum.tab
│   ├── eQTL_chr1
│   │   ├── DLPFC_WGS_Oligodendrocyte_logCPM.cis_qtl_pairs.1.parquet
│   │   ├── DLPFC_WGS_Oligodendrocyte_logCPM_cis_results_local_ancestry.parquet
│   │   ├── local_ancestry_scores.csv
│   │   └── Run_tensorQTL_with_pdata_final.py
│   ├── eQTL_code_results.zip
│   ├── FY23_FY24_FY25_MAF05.bed
│   ├── FY23_FY24_FY25_MAF05.bim
│   ├── FY23_FY24_FY25_MAF05.fam
│   ├── FY23_FY24_FY25_MAF05.log
│   ├── GRCh38-2024-A-genes.gtf.gz
│   └── pseudobulk_Celltype_aggregated.h5ad
└── tree_structure.txt
```

