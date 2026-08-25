# AANRI
snRNA sequencing QTL and TWAS analysis pipeline - African Ancestry Neuroscience Research Initiative (AANRI). 

All the codes contained in this repository are those that can be found under `1_source` directory (see tree below). 

## Table of contents
1. [Genotype processing](https://github.com/arunhpatil/AANRI/tree/main/1_src#-genotype-processing-)
2. [Preparatory - phenotype processing](https://github.com/arunhpatil/AANRI/tree/main/1_src/README.md#preparatory---phenotype-processing)
3. [Genome-wide Complex Trait Analysis](https://github.com/arunhpatil/AANRI/tree/main/1_src/README.md#gcta---template-code-and-summary)
4. [TWAS: Fusion Weights](https://github.com/arunhpatil/AANRI/tree/main/1_src/README.md#twas-fusion-weights---template-code-and-summary)
5. [TWAS: Fusion associations](https://github.com/arunhpatil/AANRI/tree/main/1_src/README.md#twas-fusion-associations---template-code-and-summary)
6. [SMR: associations](https://github.com/arunhpatil/AANRI/tree/main/1_src/README.md#smr-associations---template-code-and-summary)
8. [ColocABF](https://github.com/arunhpatil/AANRI/tree/main/1_src/README.md#colocabf---template-code-and-summary)


## <u> **Genotype processing** </u>:
**Genotype Path**: `/dcs04/lieber/hwanglab/Arun/snRNA_aanri/2_data`

The final genotype used are named with prefix `aanri_072826_genotype_qc` and the pre-processing is done as discussed below. The structure falls under `2_data`. 

```
.
├── aanri_072826_genotype_qc.bed
├── aanri_072826_genotype_qc.bim
├── aanri_072826_genotype_qc.fam
├── aanri_072826_genotype_qc.log
├── aanri_adata_snRNA_072826.h5ad
└── pre-processing_genotype
    ├── commands_updating_genotype.sh
    ├── dbsnp_update_names.canonical.txt
    ├── dbsnp_update_names_filteredSubset.txt
    ├── genotype_FY23_FY24_FY25_MAF05.bed
    ├── genotype_FY23_FY24_FY25_MAF05.bim
    ├── genotype_FY23_FY24_FY25_MAF05.fam
    ├── genotype_FY23_FY24_FY25_MAF05.log
    ├── libd_FY23_FY24_FY25_MAF05_variantList.txt
    └── update_IIDs.txt
```
</Br>

**Uniformity in Genotype naming method:**

A. FID and IID naming:<br>
The names (conventionally consider family `FID` and individual IDs `IID` by PLINK) where FID is denoted with `0` and IID for the corresponding donar `D0000_S0`. This FID is replaced with IID with only donor information and not sample, i.e., split the IID at "_" and use the 0th element as FID. This was achived using the following awk script:
```
awk '{
    split($2,a,"_");
    print $1, $2, a[1], $2
}' ../../raw_data_dongsan/FY23_FY24_FY25_MAF05.fam > update_IIDs.txt
```    
   This update_IIDs.txt will have information as shown below separated by space or space delimited. (Note these IDs are template/sample only and not real names).
```
0 D0001_S01 D0001 D0001_S01
0 D0002_S02 D0002 D0002_S02
... etc.
```
The renaming of .fam file was then achived with the following PLINK command:
```
plink2 \
  --bfile ../../raw_data_dongsan/FY23_FY24_FY25_MAF05 \
  --update-ids update_IIDs.txt \
  --threads 24 \
  --make-bed \
  --out genotype_FY23_FY24_FY25_MAF05
```
<br>

**Converting SNP coordinates to dbSNP ID (rsID):** <br>


  ```
Obtaining and formatting dbSNP 157 from NCBI:
Download:
wget -O dbSNP157.vcf.gz "https://ftp.ncbi.nih.gov/snp/latest_release/VCF/GCF_000001405.40.gz"
tabix -p vcf dbSNP157.vcf.gz 

Get reference annotations map from NCBI:
wget "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/GCF_000001405.40_GRCh38.p14/GCF_000001405.40_GRCh38.p14_assembly_report.txt"
grep -v "^#" GCF_000001405.40_GRCh38.p14_assembly_report.txt | grep "assembled-molecule" | awk '{print $7, "\t", $1}' >chr.map

# Filter chromosomes for assembled-molecues - retain only 1-22, X, Y and MT. 
bcftools annotate --rename-chrs chr.map dbSNP157.filtered_req.vcf.gz -Oz -o dbSNP157.renamed_req.vcf.gz --threads 100

# Extract and rearrange the chromosome coordinates to map to rsID. col1 consitsts of chr:coordinates and col2 consists of rsID
# chr10:100000003:C:T rs2040992368
# chr10:100000005:C:T rs775195100
# chr10:100000006:C:T rs2040992444
# chr10:100000011:T:C rs1589455386

bcftools norm -m -any -Ou --threads 80 dbSNP157.renamed_req.vcf.gz | bcftools query -f '%CHROM:%POS:%REF:%ALT\t%ID\n' > dbsnp_update_names.split56_final.txt

# Extract only variants overlapping with variants of the genotype:
awk '{print $2}' ../../raw_data_dongsan/FY23_FY24_FY25_MAF05.bim > libd_FY23_FY24_FY25_MAF05_variantList.txt
grep -Fwf libd_FY23_FY24_FY25_MAF05_variantList.txt dbsnp_update_names.canonical.txt > dbsnp_update_names_filteredSubset.txt

# Finally, filter genotype for TWAS analysis:
plink2 \
  --bfile genotype_FY23_FY24_FY25_MAF05 \
  --update-name dbsnp_update_names_filteredSubset.txt \
  --make-bed \
  --out /dcs04/lieber/hwanglab/Arun/snRNA_aanri/2_data/aanri_072826_genotype_qc \
  --threads 24 

# Dongsan's Email 07/28/26: 
  # I filtered using MAF 0.05. 
  # The other filters were not applied, but I think it is good enough to run the code as is.
  # Thus, I didn't use the following filters:
  #--geno 0.05 \
  #--hwe 1e-06 \
  #--mind 0.1 \
  #--maf 0.01 \

  ```
  

Total number of SNPs before conversion: `8,493,350` FY23_FY24_FY25_MAF05.bim </Br>
Total number of SNPs after conversion: `8,493,350` aanri_072826_genotype_qc.bim </Br>

## preparatory - phenotype processing

The Shell script to execute this section is located under `AANRI/1_src/template_code/0_create_datasets.sh`

The h5ad file was initially processed to include:  
a. FID and IID columns in the .obs (observations) by matching the BrNums, so the downstream process is easier for filtering cells and perform TWAS analysis pipeline.  
b. The .var is updated to include gene names and gene coordinates. This is required to fetch gene coordinates as input for GCTA and others. So, the first four columns include the #chr, start, end, and gene_id columns.  

The updated h5ad file is stored in the following path:
`/dcs04/lieber/hwanglab/Arun/snRNA_aanri/2_data`

1. Filters used for samples and genes.
The filtering of samples, genes, creating residual matrices, and generating covariates is in script `AANRI/1_src/1_createTWAS_ingredients.py`.

I have two directories one for (logCPM -> residuals using linear model), `residuals` and another for `regular` where the logCPM is directly used for TWAS analsis. While the later is for test purposes only, and we rely our analysis on the residuals. 

The analysis is carried out across 6 cell types from three brain regions, the pipeline is implemented 17 times. (Except one where the caudate - excitatory cells were not included - discussed below)
`adata[(adata.obs['Region'] == region) & (adata.obs['Celltype_aggregated'] == required_ct) & (adata.obs['psbulk_cells']  >= 10) & (adata.obs['psbulk_counts'] >=1000)].copy()`

Steps:
First, the raw counts were converted to logCPM values.
Apply gene filters
(a) expressed in >=20% of samples
(b) exclude genes with non-zero variance
(c) mean logCPM >= 0.5

Covariates:

Generate and include expression PCs n=3


```
ct_props = dict()
excitatory_cells = ['Ex:CA1', 'Ex:CA2_4', 'Ex:GC', 'Ex:L2_3_IT', 'Ex:L4_IT','Ex:L5_6_NP', 'Ex:L5_ET', 'Ex:L5_IT', 'Ex:L6_CT', 'Ex:L6_IT', 'Ex:L6_IT_Car3', 'Ex:L6b', 'Ex:Limbic-IT']
inhibitory_cells = ['In:CCK', 'In:CGE', 'In:CR', 'In:Chandelier', 'In:FS', 'In:LAMP5_CGE', 'In:LAMP5_MGE', 'In:LTS', 'In:Lamp5', 'In:Lamp5_Lhx6', 'In:MGE', 'In:MSN_D1', 'In:MSN_D1_D2', 'In:MSN_D2', 'In:Pax6', 'In:Pvalb', 'In:Sncg', 'In:Sst', 'In:Sst_Chodl', 'In:Vip']
cellCovariateDict = {'Astrocyte':'Astrocyte', 'Choroid_plexus':'Choroid_plexus', 'Endothelial':'Endothelial', 'Ependymal':'Ependymal', 'Excitatory_neuron':excitatory_cells, 'Inhibitory_neuron':inhibitory_cells, 'Lymphoid':'Lymphoid', 'Microglia':'Microglia', 'OPC':'OPC', 'Oligodendrocyte':'Oligodendrocyte', 'Vascular_stromal':'Pericyte'}

ct_props = cellCovariateDict[args.cellType] # If cellType='Astrocyte', then ct_props will be 'Astrocyte' proportions, if ex/in_neurons, then the above list of proportions are included.
```

```
categorical_cols = ["Sex", "FY"]

for col in categorical_cols:
    covars_model[col] = (covars_model[col].astype("category"))

covars_model = pd.get_dummies(covars_model, columns=categorical_cols, drop_first=True, dtype=int)

# log-transform pseudobulk cell counts
covars_model["psbulk_cells"] = np.log1p(covars_model["psbulk_cells"])

continuous_cols = ["AgeDeath", "PMI", "psbulk_cells", "pct_counts_mt","pct_counts_ribo",]

if isinstance(ct_props, list):
    continuous_cols.extend(ct_props)
else:
    continuous_cols.append(ct_props)
```

Remvoe constant covariates where covar[col].nunique() == 1 and Std.dev covar[col].std() == 0.


2. Covariates used when using expression residuals of logCPM values - termed as analysis = `residuals`. 
This is the default workflow - similar to that employed in miRQTL analysis. 
We are only regresssing out continuous variables (including cell type proportions). While only genotype and expression PCs (n=3) along with categorical variables such as Sex and FY was included as covariates for GCTA and TWAS analysis.

3. Covariates used when using logCPM values as expression matrix - termed as analysis = `regular`. 
This is only used to check the differences in TWAS results, compared to residuals.
Since no covariates are adjusted at this stage, all the continuous and categorical variables were passed as covariates during GCTA and TWAS analysis.   


**Error in Caudate - Excitatory neurons:** 
With the above filters applied to caudate - ExNeurons, there are no samples left to process. Hence this particular cell type was removed from the analysis. 
 
```
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
```

**Total number of genes passing the criteria:**

|     Celltype      |   Caudate |   DLPFC |   Hippo |
|:------------------|----------:|--------:|--------:|
| Astrocyte         |     20128 |   18776 |   19930 |
| Excitatory_neuron |     0.    |   21262 |   21095 |
| Inhibitory_neuron |     21141 |   20943 |   20892 |
| Microglia         |     16362 |   15082 |   16269 |
| Oligodendrocyte   |     16841 |   16134 |   16804 |
| OPC               |     19257 |   18491 |   19457 |

The counts include table header.

NOTE: celltype such as `Choroid_plexus`, `Endothelial`, `Ependymal`, `Lymphoid` and `Vascular_stromal` were not included in the analysis, as per recommendation. 

## GCTA - template code and summary

**GCTA (Genome-wide Complex Trait Analysis)**

The template code for `DLPFC` is provided in `AANRI/1_src/template_code/1_twas_gcta_residual_dlpfc.sh`.  

Internal commands:  
Create a genotype of 1Mb across the gene start and end.  
`cmd = f"{plinkPath} --bfile {genotype} --chr {chromo} --from-bp {region_start} --to-bp {region_end} --make-bed --out {genotype_file}"`

Create Genetic Relationship Matrix (GRM):  
`cmd = [gctaPath, "--bfile", str(genotype_file), "--make-grm", "--out", str(grm_file)]`  

GCTA-GREML: Estimate variance explained by all the SNPs:  
`cmd = [gctaPath, "--grm", str(grm_file), "--pheno", str(mir_file), "--covar", str(covar_d), "--qcovar", str(covar_q), "--reml", "--out", str(results_file)]`.   
Perform a REML (restricted maximum likelihood) analysis.

**Filter significant genes (P < 0.05)**:  
GCTA-GREML heritability analysis of gene expression, a significant result means that the observed variation in that gene's expression across individuals has a statistically detectable genetic component.  

In other words, for a gene g:  

$$
h_g^2 = \frac{\sigma_g^2}{\sigma_g^2 + \sigma_e^2}
$$

where:  
- $\sigma_g^2$ = variance in expression attributable to measured genetic similarity
- $\sigma_e^2$ = residual/non-genetic variance
- $h_g^2$ = SNP-based heritability of the gene's expression. 
  

Table below summarizes the list of significant genes:

Residual expression:  
|     Celltype      |   Caudate |   DLPFC |   Hippo |
|:------------------|----------:|--------:|--------:|
| Astrocyte         |      1799 |    1377 |    1686 |
| Excitatory_neuron |         0 |    2283 |    1714 |
| Inhibitory_neuron |      2402 |    2253 |    1719 |
| Microglia         |      1133 |     976 |    1227 |
| Oligodendrocyte   |      1614 |    1304 |    1679 |
| OPC               |      1420 |    1390 |    1585 |


Raw counts (i.e., logCPM):</Br>

|     Celltype      |   Caudate |   DLPFC |   Hippo |
|:------------------|----------:|--------:|--------:|
| Astrocyte         |      1918 |    1406 |    1699 |
| Excitatory_neuron |         0 |    2405 |    1826 |
| Inhibitory_neuron |      2554 |    2338 |    1787 |
| Microglia         |      1153 |    1002 |    1260 |
| Oligodendrocyte   |      1679 |    1410 |    1736 |
| OPC               |      1435 |    1441 |    1643 |

The counts doesn't include table header.

## TWAS: Fusion Weights - template code and summary
## TWAS: Fusion associations - template code and summary
## SMR: associations - template code and summary
## ColocABF - template code and summary

