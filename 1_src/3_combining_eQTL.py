from pathlib import Path
import pandas as pd
import duckdb
from statsmodels.stats.multitest import multipletests
import argparse
import os
import numpy as np

parser = argparse.ArgumentParser()
#parser.add_argument("--listCells", action='store_true', help="Helper function to print available cell types to work with and exit. Note this does not perform TWAS analysis!")
parser.add_argument("--indir", required=True, help="Input directory for eQTL .parquet files")
parser.add_argument("--genotype", required=True, help="Path to genotype, only provide the name without extension.")
parser.add_argument("--workdir", required=True, help="Workdir: where the TWAS analysis is performed and output is stored. Example: /dcs04/lieber/hwanglab/Arun/snRNA_aanri/3_analysis/DLPFC/Astrocyte/residuals")
parser.add_argument("--Region", required=True, help="Requires name of the brain region to perform TWAS analysis")
parser.add_argument("--cellType", required=True, help="Requires name of the cellType to perform TWAS analysis")
parser.add_argument("--localAncestry", action='store_true', help="If local ancestry, the output file name will have suffix localAFR")
args = parser.parse_args()

#indir = "/dcs04/lieber/hwanglab/Arun/snRNA_aanri/raw_data_dongsan/eQTL_all/DLPFC/Astrocyte/base_covariates"
#bimfile="/dcs04/lieber/hwanglab/Arun/snRNA_aanri/2_data/aanri_072826_genotype_qc.bim"
indir = Path(args.indir)
bimfile = Path(args.genotype + ".bim")
eqtlDir = Path(args.workdir) / "eQTL"
os.makedirs(eqtlDir, exist_ok=True)

if args.localAncestry:
    outputfile = Path(eqtlDir) / f"{args.Region}_{args.cellType}_.lan.tsqtl.nominal.txt"
else:
    outputfile = Path(eqtlDir) / f"{args.Region}_{args.cellType}_tsqtl.nominal.txt"

chrom_files = ",".join(f"'{indir}/*.{i}.parquet'" for i in range(1, 23))
df = duckdb.query(f"SELECT * FROM read_parquet([{chrom_files}])").df()
df["snp"] = df["variant_id"]
#df["maf"] = df["af"].clip(upper=0.5)
df["maf"] = np.where(df["af"] > 0.5, 1 - df["af"], df["af"])

bim = pd.read_csv(bimfile, sep=r"\s+", header=None, names=["chromosome", "marker.ID", "genetic.dist", "physical.pos", "allele1", "allele2"], low_memory=False)

#bim["snp"] = ("chr" + bim["chromosome"].astype(str) + ":" + bim["physical.pos"].astype(str) + ":" + bim["allele1"] + ":" + bim["allele2"])
bim["snp"] = ("chr" + bim["chromosome"].astype(str) + ":" + bim["physical.pos"].astype(str) + ":" + bim["allele2"] + ":" + bim["allele1"])
genotype_map = bim[["snp", "marker.ID", "chromosome", "physical.pos", "allele1", "allele2" ]].copy()
df = df.merge(genotype_map, on="snp", how="left")

df = df.rename(columns={"marker.ID": "rsid", "chromosome": "variant_chrom", "physical.pos": "variant_pos"})
#df = df.rename(columns={"marker.ID": "variant_id", "chromosome": "variant_chrom", "physical.pos": "variant_pos"})
df["variant_id"] = df["rsid"]
df = df.drop(columns=["rsid"])

valid = df["pval_nominal"].notna()

df["fdr"] = float("nan")

df.loc[valid, "fdr"] = multipletests(df.loc[valid, "pval_nominal"],method="fdr_bh")[1]

df.to_csv(outputfile, sep="\t", index=False)

# All eQTLs
total_eqtls = df.shape[0]
total_phenotypes = df["phenotype_id"].nunique()

print("Total eQTLs:", total_eqtls)
print("Total unique phenotypes:", total_phenotypes)

cutoffs = {
    "0.05": 0.05,
    "0.01": 0.01,
    "0.001": 0.001
}

for label, cutoff in cutoffs.items():

    sig_df = df[df["fdr"] <= cutoff]

    sig_eqtls = sig_df.shape[0]
    sig_phenotypes = sig_df["phenotype_id"].nunique()

    print(f"Region and Celltype: {args.Region} - {args.cellType} -> FDR <= {label}:")
    print(f"  eQTLs: {sig_eqtls}")
    print(f"  unique phenotypes: {sig_phenotypes}\n")
