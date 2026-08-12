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

1. Uniformity in Genotype naming method:<br><br>
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
  B. Converting SNP coordinates to dbSNP ID (rsID).

  Obtaining and formatting dbSNP 157 from NCBI:
  ```
Download:
wget -O dbSNP157.vcf.gz "https://ftp.ncbi.nih.gov/snp/latest_release/VCF/GCF_000001405.40.gz"
tabix -p vcf dbSNP157.vcf.gz 

Get reference annotations map from NCBI:
wget "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/GCF_000001405.40_GRCh38.p14/GCF_000001405.40_GRCh38.p14_assembly_report.txt"
grep -v "^#" GCF_000001405.40_GRCh38.p14_assembly_report.txt | grep "assembled-molecule" | awk '{print $7, "\t", $1}' >chr.map


  ```


3. 
4. 
5. 

## preparatory - phenotype processing
## GCTA - template code and summary
## TWAS: Fusion Weights - template code and summary
## TWAS: Fusion associations - template code and summary
## SMR: associations - template code and summary
## ColocABF - template code and summary

