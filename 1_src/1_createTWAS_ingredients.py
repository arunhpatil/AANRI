#!/usr/bin/python

import argparse
import os
import sys
from pathlib import Path
import scanpy as sc
import pandas as pd 
import numpy as np
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
from sklearn.linear_model import LinearRegression

# The covariates I used for this study, for each donor (BrNum), is
# AgeDeath (continuous), Sex, PMI, FY, GENOTYPE_PC1, PC2, PC3, psbulk_cells (number of cells),
# pct_counts_mt, pct_counts_ribo and
# the corresponding cell type proportions (e.g., for astrocyte I used Astrocyte (%Astrocyte for each donor)
# as a covariate) For excitatory and inhibitory neurons I also included the proportion of specific
# neuronal proportions (columns start with Ex:, or In:) as covariates.
# YRI is a global ancestry (% of African) and is almost identical to GENOTYPE_PC1, so I did not use YRI.

# You can also compute expression PCs and use PC1 to PC3 as covariates.
# I used samples with psbulk_cells>=10 and psbulk_counts>=1000 for the downstream analysis.

##############################################
# ARGUMENTS
##############################################
parser = argparse.ArgumentParser()
#parser.add_argument("--listCells", action='store_true', help="Helper function to print available cell types to work with and exit. Note this does not perform TWAS analysis!")
parser.add_argument("--adata", required=True, help="Input h5ad file")
parser.add_argument("--genotype", required=True, help="Path to genotype, only provide the name without extension.")
parser.add_argument("--workdir", required=True, help="Workdir: where the TWAS analysis is performed and output is stored. Example: /dcs04/lieber/hwanglab/Arun/snRNA_aanri/3_analysis/")
parser.add_argument("--Region", required=True, help="Requires name of the brain region to perform TWAS analysis")
parser.add_argument("--cellType", required=True, help="Requires name of the cellType to perform TWAS analysis")
args = parser.parse_args()


ct_props = dict()
input_h5ad = args.adata
excitatory_cells = ['Ex:CA1', 'Ex:CA2_4', 'Ex:GC', 'Ex:L2_3_IT', 'Ex:L4_IT','Ex:L5_6_NP', 'Ex:L5_ET', 'Ex:L5_IT', 'Ex:L6_CT', 'Ex:L6_IT', 'Ex:L6_IT_Car3', 'Ex:L6b', 'Ex:Limbic-IT']
inhibitory_cells = ['In:CCK', 'In:CGE', 'In:CR', 'In:Chandelier', 'In:FS', 'In:LAMP5_CGE', 'In:LAMP5_MGE', 'In:LTS', 'In:Lamp5', 'In:Lamp5_Lhx6', 'In:MGE', 'In:MSN_D1', 'In:MSN_D1_D2', 'In:MSN_D2', 'In:Pax6', 'In:Pvalb', 'In:Sncg', 'In:Sst', 'In:Sst_Chodl', 'In:Vip']
cellCovariateDict = {'Astrocyte':'Astrocyte', 'Choroid_plexus':'Choroid_plexus', 'Endothelial':'Endothelial', 'Ependymal':'Ependymal', 'Excitatory_neuron':excitatory_cells, 'Inhibitory_neuron':inhibitory_cells, 'Lymphoid':'Lymphoid', 'Microglia':'Microglia', 'OPC':'OPC', 'Oligodendrocyte':'Oligodendrocyte', 'Vascular_stromal':'Pericyte'}

adata = sc.read_h5ad(args.adata)
# adata = sc.read_h5ad("/dcs04/lieber/hwanglab/Arun/snRNA_aanri/2_data/aanri_adata_snRNA_072826.h5ad") # FYI: use this h5ad, it includes genotype keys and gene coordinates updated.
if args.cellType in cellCovariateDict.keys():
    ct_props = cellCovariateDict[args.cellType]
    covariates_list = ['AgeDeath', 'Sex', 'PMI', 'FY', 'psbulk_cells', 'pct_counts_mt', 'pct_counts_ribo', ct_props, 'GENOTYPE_PC1', 'GENOTYPE_PC2', 'GENOTYPE_PC3']
    print(f"Cell proportions Covariates used for cell type: {args.cellType}\n{ct_props}\n")
    print(f"All covariates included: {covariates_list}\n")
    print("Includes: Expression PCs - PC1 to PC3 as covariates.\n\n")
else:
    print("The provided cell type is not available. Please select the appropriate cell type from the list: ")
    print(f"Cell Type: {adata.obs['Celltype_aggregated'].unique().tolist()}\n")
    exit()
    #print(adata.obs[['Region','Celltype_aggregated']].value_counts())
if args.Region not in adata.obs['Region'].unique().tolist():
    print("The Region specified is incorrect. Please see the list of available regions:")
    print(f"Region: {adata.obs['Region'].unique().tolist()}\n")
    exit()

required_ct = args.cellType # Astrocyte 
region = args.Region # DLPFC
regular_sub_path = str(region) + "/" + str(required_ct) + "/regular/TWAS/"
residual_sub_path = str(region) + "/" + str(required_ct) + "/residuals/TWAS/"

reg_workdir = Path(args.workdir)/str(regular_sub_path)  # path to 3_analysis 
res_workdir = Path(args.workdir)/str(residual_sub_path)  # path to 3_analysis 
os.makedirs(args.workdir, exist_ok=True)
genotype = Path(args.genotype)

famGeno = str(genotype)+".fam"
fam = pd.read_csv(famGeno, sep=r"\s+", header=None, names=["FID", "IID", "PAT", "MAT", "SEX", "PHENO"])

adata_sub = adata[(adata.obs['Region'] == region) & (adata.obs['Celltype_aggregated'] == required_ct) & (adata.obs['psbulk_cells']  >= 10) & (adata.obs['psbulk_counts'] >=1000)].copy()
#adata_sub = adata[(adata.obs['Region'] == "DLPFC") & (adata.obs['Celltype_aggregated'] == "Astrocyte") & (adata.obs['psbulk_cells']  >= 10) & (adata.obs['psbulk_counts'] >=1000)].copy()

# ----------------------------------------------------
# 2. Extract raw pseudobulk counts
# ----------------------------------------------------
counts = pd.DataFrame(adata_sub.X, index=adata_sub.obs["IID"],columns=adata_sub.var_names)

# ----------------------------------------------------
# 3. Convert raw counts -> CPM -> logCPM
# ----------------------------------------------------
library_sizes = counts.sum(axis=1)
cpm = counts.div(library_sizes, axis=0) * 1e6
logcpm = np.log2(cpm + 1)

# ----------------------------------------------------
# 4. Gene filters
# ----------------------------------------------------
# (a) expressed in >=20% of samples
nonzero_counts = (counts > 0).sum(axis=0)
min_samples_expr = int(np.ceil(0.20 * counts.shape[0]))
expr_filter = nonzero_counts >= min_samples_expr

# (b) non-zero variance
var_filter = logcpm.var(axis=0) > 0

# (c) mean logCPM >= 0.5
avg_expr = logcpm.mean(axis=0)

cpm_filter = avg_expr >= 0.5

combined_filter = (expr_filter & var_filter & cpm_filter)

print(f"Genes before filtering : {adata_sub.n_vars}")
print(f"Genes after filtering  : {combined_filter.sum()}")

# ----------------------------------------------------
# 5. Subset AnnData and expression matrices
# ----------------------------------------------------
adata_sub = adata_sub[:, combined_filter.values].copy()
counts = counts.loc[:, combined_filter]
logcpm = logcpm.loc[:, combined_filter]

# ----------------------------------------------------
# 6. Build expression matrix (genes x samples)
# ----------------------------------------------------
expr = logcpm.T
expr.columns = adata_sub.obs["IID"]
# ----------------------------------------------------
# 7. Add genomic coordinates
# ----------------------------------------------------
gene_info = (adata_sub.var[["chr", "start", "end", "gene"]].rename(columns={"chr": "#chr","gene": "gene_id"}))
gene_info["#chr"] = (gene_info["#chr"].str.replace("^chr", "", regex=True))

expr = pd.concat([gene_info, expr], axis=1)
meta_cols = ["#chr", "start", "end", "gene_id"]

available_samples = [str(col).strip() for col in expr.columns if col not in meta_cols]
fam_iids = fam["IID"].astype(str).str.strip()

common_samples = fam_iids[fam_iids.isin(available_samples)].tolist()
expr = expr[meta_cols + common_samples]
#available = expr.columns.difference(meta_cols)
#common_samples = fam["IID"][fam["IID"].isin(available)]
#expr = expr[meta_cols + common_samples.tolist()]

# ----------------------------------------------------
# 8. Compute expression PCs
# ----------------------------------------------------
expr_sample_order = expr.columns[4:]
X = expr.loc[:, expr_sample_order].T
#print("X shape:", X.shape)
#print("X columns:", len(X.columns))
#print("X index:", len(X.index))
#print(X.head())
#print(type(X))

X = X.apply(pd.to_numeric, errors="coerce")
# ── Diagnose empty DataFrame ──────────────────────────────────────────
print("-"*25, " START ", "-"*23)
print("X shape after to_numeric:", X.shape)
print("X all-NaN columns:", X.isna().all(axis=0).sum())
print("X all-NaN rows:", X.isna().all(axis=1).sum())
print("X fully empty:", X.empty)
print("-"*25, " END ", "-"*25)
# ─────────────────────────────────────────────────────────────────────

if X.isna().any().any():
    raise ValueError("NaN values detected before PCA")

if np.isinf(X.to_numpy()).any():
    raise ValueError("Infinite values detected before PCA")

# ── Pinpoint check ────────────────────────────────────────────────────
print("X shape just before to_numpy:", X.shape)   # <── add this
X_array = X.to_numpy(dtype=np.float64)
print("X_array shape:", X_array.shape)             # <── add this
# ─────────────────────────────────────────────────────────────────────
print()
#print("PCA input shape:", X.shape)
#print("PCA input dtype:", X.dtypes.unique())

#print(PCA)
#print(X.dtypes.value_counts())

# Convert DataFrame -> NumPy array before PCA
X_array = X.to_numpy(dtype=np.float64)

try:
    pcs = PCA(n_components=3).fit_transform(X_array)
    print("PCA call 1 succeeded:", pcs.shape)
except ValueError as e:
    print("PCA call 1 FAILED — X_array shape was:", X_array.shape)
    raise

expr_pcs = pd.DataFrame(pcs, columns=["Expr_PC1", "Expr_PC2", "Expr_PC3"], index=X.index)
expr_pcs.index.name = "IID"

# ----------------------------------------------------
# 9. Build covariate table
# ----------------------------------------------------
ct_props = cellCovariateDict[required_ct]
cvs = ["FID","IID","AgeDeath","Sex","PMI","FY","psbulk_cells","pct_counts_mt","pct_counts_ribo","GENOTYPE_PC1","GENOTYPE_PC2","GENOTYPE_PC3"]
if isinstance(ct_props, list):
    cvs.extend(ct_props)
else:
    cvs.append(ct_props)

obs_covars = adata_sub.obs[cvs].copy()
obs_covars = obs_covars.merge(expr_pcs.reset_index(), on="IID",how="left")

# ----------------------------------------------------
# 10. Reorder covariates to match PLINK FAM
# ----------------------------------------------------
covars = fam[["FID", "IID"]].merge(obs_covars, on=["FID", "IID"], how="left")

# ----------------------------------------------------
# 11. Remove samples absent from expression matrix
# ----------------------------------------------------
missing = covars["Expr_PC1"].isna()

if missing.any():
    print(f"Removing {missing.sum()} sample(s):")
    print(covars.loc[missing,["FID", "IID"]])

covars = covars.loc[~missing].copy()

assert (expr.columns[4:].tolist() == covars["IID"].tolist())

# ----------------------------------------------------
# 12. Prepare covariates
# ----------------------------------------------------
covars_model = covars.copy()

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

scaler = StandardScaler()

covars_model[continuous_cols] = scaler.fit_transform(covars_model[continuous_cols])

# ----------------------------------------------------
# 13.Remvoe constant covariates where covar[col].nunique() == 1 and Std.dev covar[col].std() == 0: 
# ----------------------------------------------------
def remove_constant_covariates(df, id_cols=("FID", "IID")):

    model_cols = [
        c for c in df.columns
        if c not in id_cols
    ]

    constant_cols = [
        c for c in model_cols
        if df[c].nunique(dropna=False) <= 1
    ]

    if constant_cols:
        print("Removing constant covariates:")
        for col in constant_cols:
            print(
                f"  {col}: "
                f"unique={df[col].nunique(dropna=False)}, "
                f"std={df[col].std()}"
            )

        df = df.drop(columns=constant_cols)

    return df

covars_model = remove_constant_covariates(covars_model)

# ----------------------------------------------------
# 13a. Prepare residual expression 
# ----------------------------------------------------
#     samples x genes
# ----------------------------------------------------
expr_values = expr.drop(columns=meta_cols).copy()
covars_model["IID"] = covars_model["IID"].astype(str)
expr_values.columns = expr_values.columns.astype(str)

# ----------------------------------------------------
# Match expression samples to covariate samples
# ----------------------------------------------------
#sample_ids = [iid for iid in covars_model["IID"] if iid in expr_values.columns]
sample_ids = covars_model["IID"][covars_model["IID"].isin(expr_values.columns)].tolist()

# Reorder expression columns to exactly match covariate order
expr_values = expr_values.loc[:, sample_ids]

#sample_ids = covars_model["IID"][covars_model["IID"].isin(expr_values.columns)].tolist()

assert sample_ids == covars_model["IID"].tolist(), (
    "Expression samples and covariate samples are not in the same order."
)

# ----------------------------------------------------
# 13b. Expression matrix for residualization
#     samples x genes
# ----------------------------------------------------
#Y = logcpm.loc[:, sample_ids].T.copy()
Y = expr_values.loc[:, sample_ids].T.copy()

print("Expression matrix for residualization:")
print(Y.shape)

#residual_covar_cols = ["AgeDeath","PMI","psbulk_cells","pct_counts_mt","pct_counts_ribo","GENOTYPE_PC1","GENOTYPE_PC2","GENOTYPE_PC3","Expr_PC1","Expr_PC2","Expr_PC3",]
residual_covar_cols = ["AgeDeath","PMI","psbulk_cells","pct_counts_mt","pct_counts_ribo",] # We are only regresssing out continuous variables. 

# Add cell-type proportion covariates
if isinstance(ct_props, list):
    residual_covar_cols.extend(ct_props)
else:
    residual_covar_cols.append(ct_props)

# Add dummy-coded categorical variables
# Sex and FY have already been converted by get_dummies()
#dummy_covar_cols = [
#    col for col in covars_model.columns
#    if col.startswith("Sex_") or col.startswith("FY_")
#]

#residual_covar_cols.extend(dummy_covar_cols)

# Remove accidental duplicates while preserving order
residual_covar_cols = list(dict.fromkeys(residual_covar_cols))

residual_covar_cols = [i for i in residual_covar_cols if i in covars_model.columns]

print("\nCovariates used for residualization:")
print(residual_covar_cols)

# ----------------------------------------------------
# 13d. Build residualization matrix
# ----------------------------------------------------
#X_resid = (covars_model.set_index("IID").loc[sample_ids].copy())
X_resid = (covars_model.set_index("IID").loc[sample_ids, residual_covar_cols].copy())
print(X_resid)

# Make sure samples correspond exactly
assert Y.index.tolist() == X_resid.index.tolist()

print("\nCovariates used for residualization:")
print(X_resid.columns.tolist())

# ----------------------------------------------------
# 13e. Check for missing values
# ----------------------------------------------------
if X_resid.isna().any().any():
    print("Missing values found in residualization covariates:")
    print(X_resid.isna().sum()[X_resid.isna().sum() > 0])
    raise ValueError(
        "Missing values must be handled before residualization."
    )

# ----------------------------------------------------
# 13f. Fit covariate model to all genes simultaneously
#
#     Expression = covariates + residual
# ----------------------------------------------------
lm = LinearRegression(fit_intercept=True)

lm.fit(X_resid, Y)

# Predicted expression explained by covariates
Y_hat = lm.predict(X_resid)

# Residual expression
Y_resid = Y - Y_hat

# ----------------------------------------------------
# 13g. Convert back to genes x samples
# ----------------------------------------------------
residual_expr = Y_resid.T.copy()

#residual_expr.index = logcpm.index
#residual_expr.columns = sample_ids

print("\nResidual expression:")
print(residual_expr.shape)

# ----------------------------------------------------
# 13h. Add genomic information
# ----------------------------------------------------
residual_expr = pd.concat([gene_info,residual_expr],axis=1)

# #chr start end gene_id sample1 sample2 ...
residual_expr = residual_expr[
    meta_cols + sample_ids
]

# ----------------------------------------------------
# 13i. Check residual expression
# ----------------------------------------------------
print("\nResidual expression summary:")
print(residual_expr.iloc[:, 4:].describe().iloc[[0, 1, 2, 5, 7]])


# ----------------------------------------------------
# 14a. Covariates for ORIGINAL logCPM expression
# ----------------------------------------------------
covars_original = covars_model.copy()

# Ensure FID/IID come first
covar_first = ["FID", "IID"]

covars_original = covars_original[
    covar_first +
    [
        c for c in covars_original.columns
        if c not in covar_first
    ]
]

print("\nOriginal-expression covariates:")
print(covars_original.columns.tolist())

# ----------------------------------------------------
# 14b. Covariates for RESIDUALIZED expression
# ----------------------------------------------------
required_covar_cols = ["FID","IID","GENOTYPE_PC1","GENOTYPE_PC2","GENOTYPE_PC3","Expr_PC1", "Expr_PC2","Expr_PC3",]
# Add dummy-coded categorical variables
# Sex and FY have already been converted by get_dummies()
dummy_covar_cols = [
    col for col in covars_model.columns
    if col.startswith("Sex_") or col.startswith("FY_")
]

required_covar_cols.extend(dummy_covar_cols)

covars_residual = covars_model[required_covar_cols].copy()

print("\nResidual-expression covariates:")
print(covars_residual.columns.tolist())

###
subdirectories = ["expression","covariates"] 
newSubDirs = []
for sub in subdirectories:
    dir_path = reg_workdir / sub
    newSubDirs.append(str(dir_path))
    dir_path.mkdir(parents=True, exist_ok=True)
    res_dir_path = res_workdir / sub
    newSubDirs.append(str(res_dir_path))
    res_dir_path.mkdir(parents=True, exist_ok=True)


exprn_bed = Path(newSubDirs[0])/f"{region}_{required_ct}_expression.bed"
expr.to_csv(exprn_bed, sep="\t", index=False)

resid_expr_path = Path(newSubDirs[1])/f"{region}_{required_ct}_residual_expression.bed"
residual_expr.to_csv(resid_expr_path, sep="\t", index=False)

#covars_original = remove_constant_covariates(covars_original)

covar_path = Path(newSubDirs[2])/f"{region}_{required_ct}_covariates.txt"
covars_original.to_csv(covar_path, sep="\t", index=False)
#covars_model.to_csv(covar_path, sep="\t", index=False)

#covars_residual = remove_constant_covariates(covars_residual)
resid_covar_path = Path(newSubDirs[3])/f"{region}_{required_ct}_residual_covariates.txt"
covars_residual.to_csv(resid_covar_path, sep="\t", index=False)
#covars_model.to_csv(resid_covar_path, sep="\t", index=False)
##################################################
#. Directory path list 
#├── 1_source
#├── 2_data
#├── 3_analysis
#└── raw_data_dongsan
#4 directories, 0 files
