# AANRI
snRNA sequencing QTL and TWAS analysis pipeline - African Ancestry Neuroscience Research Initiative (AANRI). 

All the codes contained in this repository are those that can be found under `1_source` directory (see tree below). 

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


## preparatory - phenotype processing
## GCTA - template code and summary
## TWAS: Fusion Weights - template code and summary
## TWAS: Fusion associations - template code and summary
## SMR: associations - template code and summary
## ColocABF - template code and summary

