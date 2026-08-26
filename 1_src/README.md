# AANRI
snRNA sequencing QTL and TWAS analysis pipeline - African Ancestry Neuroscience Research Initiative (AANRI). 

All the codes contained in this repository are those that can be found under `1_source` directory (see tree below). 

## Table of contents
1. [Genotype processing](https://github.com/arunhpatil/AANRI/blob/main/1_src#-genotype-processing-)
2. [Preparatory - phenotype processing](https://github.com/arunhpatil/AANRI/blob/main/1_src/README.md#preparatory---phenotype-processing)
3. [Genome-wide Complex Trait Analysis](https://github.com/arunhpatil/AANRI/blob/main/1_src/README.md#gcta---template-code-and-summary)
4. [TWAS: Fusion Weights](https://github.com/arunhpatil/AANRI/blob/main/1_src/README.md#twas-fusion-weights---template-code-and-summary)
5. [TWAS: Fusion associations](https://github.com/arunhpatil/AANRI/blob/main/1_src/README.md#twas-fusion-associations---template-code-and-summary)
6. [TensorQTL: eQTL-preprocessing](https://github.com/arunhpatil/AANRI/tree/main/1_src#eqtl-processing---template-code-and-summary)
7. [SMR: associations](https://github.com/arunhpatil/AANRI/blob/main/1_src/README.md#smr-associations---template-code-and-summary)
8. [ColocABF](https://github.com/arunhpatil/AANRI/blob/main/1_src/README.md#colocabf---template-code-and-summary)


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
## eQTL: processing - template code and summary

### Data description and modifications:
The eQTL analysis involving TensorQTL was conducted by Dr. Dongsan Kim. The raw result files were made available at Server5 at the following location. 
```
/mnt/pv_compute/dongsan/datasets/AANRI/eQTL/caudate_Celltype_aggregated_with_RPSgenes_final
/mnt/pv_compute/dongsan/datasets/AANRI/eQTL/DLPFC_Celltype_aggregated_with_RPSgenes_final
/mnt/pv_compute/dongsan/datasets/AANRI/eQTL/hippo_Celltype_aggregated_with_RPSgenes_final
```
The folders consists of nominal and permutation based eQTL results. The nominal files were executed in two batches, i.e., with and without local ancestry as covariates. These files were transferred to individual brain region and their corresponding celltypes on JHPCE. Two sub-directories were created namely `local_ancestry` and `base_covariates`, where nominal eQTL results were transferred for covariates with and without local ancestry, respectively. The directory structure and path is as follows:  
  
Path: `/dcs04/lieber/hwanglab/Arun/snRNA_aanri/raw_data_dongsan/eQTL_all`
```
.
├── caudate
├── DLPFC
└── hippo
```
Each region can further be browsed as shown with an example of DLPFC:
```
DLPFC/
├── Astrocyte
├── Excitatory_neuron
├── Inhibitory_neuron
├── Microglia
├── Oligodendrocyte
└── OPC

Further,

DLPFC/
├── Astrocyte
│   ├── base_covariates
│   │   ├── DLPFC_WGS_Astrocyte_logCPM.cis_qtl_pairs.10.parquet
│   │   ├── DLPFC_WGS_Astrocyte_logCPM.cis_qtl_pairs.11.parquet
│   │   ├── DLPFC_WGS_Astrocyte_logCPM.cis_qtl_pairs.12.parquet
│   │   ├── DLPFC_WGS_Astrocyte_logCPM.cis_qtl_pairs.13.parquet
│   │   ├── DLPFC_WGS_Astrocyte_logCPM.cis_qtl_pairs.14.parquet
│   │   ├── DLPFC_WGS_Astrocyte_logCPM.cis_qtl_pairs.15.parquet
│   │   ├── DLPFC_WGS_Astrocyte_logCPM.cis_qtl_pairs.16.parquet
│   │   ├── DLPFC_WGS_Astrocyte_logCPM.cis_qtl_pairs.17.parquet
│   │   ├── DLPFC_WGS_Astrocyte_logCPM.cis_qtl_pairs.18.parquet
│   │   ├── DLPFC_WGS_Astrocyte_logCPM.cis_qtl_pairs.19.parquet
│   │   ├── DLPFC_WGS_Astrocyte_logCPM.cis_qtl_pairs.1.parquet
│   │   ├── DLPFC_WGS_Astrocyte_logCPM.cis_qtl_pairs.20.parquet
│   │   ├── DLPFC_WGS_Astrocyte_logCPM.cis_qtl_pairs.21.parquet
│   │   ├── DLPFC_WGS_Astrocyte_logCPM.cis_qtl_pairs.22.parquet
│   │   ├── DLPFC_WGS_Astrocyte_logCPM.cis_qtl_pairs.2.parquet
│   │   ├── DLPFC_WGS_Astrocyte_logCPM.cis_qtl_pairs.3.parquet
│   │   ├── DLPFC_WGS_Astrocyte_logCPM.cis_qtl_pairs.4.parquet
│   │   ├── DLPFC_WGS_Astrocyte_logCPM.cis_qtl_pairs.5.parquet
│   │   ├── DLPFC_WGS_Astrocyte_logCPM.cis_qtl_pairs.6.parquet
│   │   ├── DLPFC_WGS_Astrocyte_logCPM.cis_qtl_pairs.7.parquet
│   │   ├── DLPFC_WGS_Astrocyte_logCPM.cis_qtl_pairs.8.parquet
│   │   ├── DLPFC_WGS_Astrocyte_logCPM.cis_qtl_pairs.9.parquet
│   │   └── DLPFC_WGS_Astrocyte_logCPM.cis_qtl_pairs.X.parquet
│   └── local_ancestry
│       ├── DLPFC_WGS_Astrocyte_logCPM_local_ancestry.cis_qtl_pairs.10.parquet
│       ├── DLPFC_WGS_Astrocyte_logCPM_local_ancestry.cis_qtl_pairs.11.parquet
│       ├── DLPFC_WGS_Astrocyte_logCPM_local_ancestry.cis_qtl_pairs.12.parquet
│       ├── DLPFC_WGS_Astrocyte_logCPM_local_ancestry.cis_qtl_pairs.13.parquet
│       ├── DLPFC_WGS_Astrocyte_logCPM_local_ancestry.cis_qtl_pairs.14.parquet
│       ├── DLPFC_WGS_Astrocyte_logCPM_local_ancestry.cis_qtl_pairs.15.parquet
│       ├── DLPFC_WGS_Astrocyte_logCPM_local_ancestry.cis_qtl_pairs.16.parquet
│       ├── DLPFC_WGS_Astrocyte_logCPM_local_ancestry.cis_qtl_pairs.17.parquet
│       ├── DLPFC_WGS_Astrocyte_logCPM_local_ancestry.cis_qtl_pairs.18.parquet
│       ├── DLPFC_WGS_Astrocyte_logCPM_local_ancestry.cis_qtl_pairs.19.parquet
│       ├── DLPFC_WGS_Astrocyte_logCPM_local_ancestry.cis_qtl_pairs.1.parquet
│       ├── DLPFC_WGS_Astrocyte_logCPM_local_ancestry.cis_qtl_pairs.20.parquet
│       ├── DLPFC_WGS_Astrocyte_logCPM_local_ancestry.cis_qtl_pairs.21.parquet
│       ├── DLPFC_WGS_Astrocyte_logCPM_local_ancestry.cis_qtl_pairs.22.parquet
│       ├── DLPFC_WGS_Astrocyte_logCPM_local_ancestry.cis_qtl_pairs.2.parquet
│       ├── DLPFC_WGS_Astrocyte_logCPM_local_ancestry.cis_qtl_pairs.3.parquet
│       ├── DLPFC_WGS_Astrocyte_logCPM_local_ancestry.cis_qtl_pairs.4.parquet
│       ├── DLPFC_WGS_Astrocyte_logCPM_local_ancestry.cis_qtl_pairs.5.parquet
│       ├── DLPFC_WGS_Astrocyte_logCPM_local_ancestry.cis_qtl_pairs.6.parquet
│       ├── DLPFC_WGS_Astrocyte_logCPM_local_ancestry.cis_qtl_pairs.7.parquet
│       ├── DLPFC_WGS_Astrocyte_logCPM_local_ancestry.cis_qtl_pairs.8.parquet
│       └── DLPFC_WGS_Astrocyte_logCPM_local_ancestry.cis_qtl_pairs.9.parquet
etc..
```

For the SMR/coloc analysis, the nominal eQTLs in parquet files for each chromosome were combined to one file using script `3_combining_eQTL.py` in `1_src` directory. The scripts reads each parquet using `duckdb` package (increased efficiency), and concatenates them using Pandas.  

```
>>> df = pd.read_parquet("DLPFC_WGS_Astrocyte_logCPM.cis_qtl_pairs.22.parquet")
>>> df
            phenotype_id          variant_id  start_distance        af  ma_samples  ma_count  pval_nominal     slope  slope_se
0        ENSG00000290397  chr22:10843943:T:C         -983579  0.918919          12        12      0.715712  0.164319  0.448907
1        ENSG00000290397  chr22:10849890:A:C         -977632  0.885135          17        17      0.809842 -0.101206  0.418608
2        ENSG00000290397  chr22:10849945:C:A         -977577  0.770270          34        34      0.531529  0.206222  0.327554
3        ENSG00000290397  chr22:10850883:C:A         -976639  0.939189           9         9      0.490412 -0.362050  0.521519
4        ENSG00000290397  chr22:10850893:C:T         -976629  0.824324          26        26      0.869154 -0.055686  0.336494

```
The variant_id was replaced to `dbSNP:rs_ID`, added variant: chrom, Pos, allele1 (effect allele), allele2 reading the Genotype BIM file. Further computed `fdr` using `multipletests` function from `statsmodel` Python package. The conversion is completed for all regions and celltypes. The base_covariates results are stored with file name `DLPFC_Astrocyte_tsqtl.nominal.txt`, while the `DLPFC_Astrocyte_.lan.tsqtl.nominal.txt` file is for `local ancestry (lan)`. 

`df.loc[valid, "fdr"] = multipletests(df.loc[valid, "pval_nominal"],method="fdr_bh")[1]`. 

```
$ pwd
/dcs04/lieber/hwanglab/Arun/snRNA_aanri/3_analysis/DLPFC/Astrocyte/residuals/eQTL
.
├── DLPFC_Astrocyte_.lan.tsqtl.nominal.txt
└── DLPFC_Astrocyte_tsqtl.nominal.txt

head DLPFC_Astrocyte_tsqtl.nominal.txt | column -t
phenotype_id     variant_id    start_distance  af          ma_samples  ma_count  pval_nominal          slope        slope_se    snp             maf          variant_chrom  variant_pos  allele1  allele2  fdr
ENSG00000238009  rs1434509538  -118879         0.91891897  12          12        0.7843163003942315    0.110552564  0.40198633  chr1:14843:G:A  0.08108103   1              14843        A        G        0.9978769812409582
ENSG00000238009  rs866639523   -117442         0.945946    8           8         0.07545377341255818   0.8399191    0.4637025   chr1:16280:T:C  0.05405402   1              16280        C        T        0.9537443981605194
ENSG00000238009  rs774196730   -116897         0.9391892   9           9         0.009834249713378727  -1.3295406   0.4974234   chr1:16825:C:A  0.060810804  1              16825        A        C        0.8438787875295181
ENSG00000238009  rs201535981   -116337         0.7567568   36          36        0.2740505190050702    -0.33532202  0.3035649   chr1:17385:G:A  0.24324322   1              17385        A        G        0.9839866910750467
ENSG00000238009  rs200784459   -116324         0.9054054   14          14        0.48801417003094405   0.2831508    0.4056121   chr1:17398:C:A  0.0945946    1              17398        A        C        0.9925521415143639
ENSG00000238009  rs747093451   -116314         0.9324325   10          10        0.13706257698434146   0.65092975   0.43152413  chr1:17408:C:G  0.06756753   1              17408        G        C        0.9699466119237674
ENSG00000238009  rs866150608   -116163         0.87162167  19          19        0.9005079851421448    -0.04926391  0.39226657  chr1:17559:G:C  0.12837833   1              17559        C        G        0.9993303866841581
ENSG00000238009  rs377698370   -116128         0.8648649   20          20        0.743929954308377     0.11680671   0.35582238  chr1:17594:C:T  0.13513511   1              17594        T        C        0.9973206854164085
ENSG00000238009  rs71260069    -116025         0.7567568   36          36        0.2185762969697761    -0.3845101   0.30902007  chr1:17697:G:C  0.24324322   1              17697        C        G        0.9799183842406497
```

### eQTL summary:

Number of eQTLs/egenes per Celltype across brain regions, (Column 3,4,5 is `FDR` cutoff):

#### DLPFC (Base covariates):  

|     Celltype      |  eQTL (q)   | eGenes (g) |    0.05 q (g)   |   0.01 q (g)    |  0.001 q (g) |
|:------------------|------------:|-----------:|----------------:|----------------:|-------------:|
| Astrocyte         | 107,046,444 |   18,483   |  12,302 (811)   |  7,303 (313)    |  3,367 (132) |
| Excitatory_neuron | 116,874,884 |   20,117   |  52,424 (2,667) | 29,570 (1,148)  | 15,933 (608) |
| Inhibitory_neuron | 114,795,594 |   19,783   |  46,680 (2,199) | 26,810 (879)    | 16,092 (485) |
| Microglia         |  91,436,137 |   15,781   |  11,540 (680)   |  6,952 (265)    |  4,552 (119) |
| Oligodendrocyte   |  92,256,675 |   15,984   |  18,578 (1,182) | 10,418 (503)    |  5,521 (264) |
| OPC               | 105,281,193 |   18,165   |  15,657 (957)   | 10,451 (412)    |  6,369 (188) |

#### DLPFC (local ancestry):  

|     Celltype      |  eQTL (q)   | eGenes (g) |    0.05 q (g)   |   0.01 q (g)    |  0.001 q (g) |
|:------------------|------------:|-----------:|----------------:|----------------:|-------------:|
| Astrocyte         |  99,585,437 |   17,186   |   7,943 (576)   |  3,700 (202)    |  1,842 (91)  |
| Excitatory_neuron | 107,598,629 |   18,540   |  33,475 (1,834) | 19,440 (823)    | 10,354 (450) |
| Inhibitory_neuron | 106,163,826 |   18,308   |  34,061 (1,672) | 19,790 (698)    | 10,434 (357) |
| Microglia         |  85,704,651 |   14,776   |   8,945 (522)   |  5,769 (214)    |  3,673 (95)  |
| Oligodendrocyte   | 107,046,444 |   18,483   |  12,302 (811)   |  7,303 (313)    |  3,367 (132) |
| OPC               |  98,167,385 |   16,936   |  10,510 (693)   |  7,011 (290)    |  4,038 (135) |

#### Caudate (Base covariates):  

|     Celltype      |  eQTL (q)   | eGenes (g) |    0.05 q (g)   |     0.01 q (g)  |  0.001 q (g) |
|:------------------|------------:|-----------:|----------------:|----------------:|-------------:|
| Astrocyte         | 112,566,826 |   19,429   |  27,214 (1,555) |  15,444 (619)   |  9,789 (300) |
| Excitatory_neuron | 107,598,629 |   18,540   |  33,475 (1,834) |  19,440 (823)   | 10,354 (450) |
| Inhibitory_neuron | 106,163,826 |   18,308   |  34,061 (1,672) |  19,790 (698)   | 10,434 (357) |
| Microglia         |  85,704,651 |   14,776   |   8,945 (522)   |   5,769 (214)   |  3,673 (95)  |
| Oligodendrocyte   | 107,046,444 |   18,483   |  12,302 (811)   |   7,303 (313)   |  3,367 (132) |
| OPC               |  98,167,385 |   16,936   |  10,510 (693)   |   7,011 (290)   |  4,038 (135) |


## SMR: associations - template code and summary
## ColocABF - template code and summary

