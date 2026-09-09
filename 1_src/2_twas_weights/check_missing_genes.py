#!/usr/bin/python
import sys
import pandas as pd
from pathlib import Path
# /dcs04/lieber/hwanglab/Arun/snRNA_aanri/3_analysis/DLPFC/OPC/residuals/TWAS/Fusion/twas_models
# /dcs04/lieber/hwanglab/Arun/snRNA_aanri/3_analysis/DLPFC/OPC/residuals/TWAS/GCTA/significant_heritability.csv

try:
    weights_dir=sys.argv[1]
except IndexError:
    print("USAGE: python check_missing_genes.py /dcs04/lieber/hwanglab/Arun/snRNA_aanri/3_analysis/DLPFC/Excitatory_neuron/residuals/TWAS/Fusion/twas_models /dcs04/lieber/hwanglab/Arun/snRNA_aanri/3_analysis/DLPFC/Excitatory_neuron/residuals/TWAS/GCTA/significant_heritability.csv\n")
    exit()

weights_dir = Path(sys.argv[1])
gcta_sig_genes = sys.argv[2] 



df = pd.read_csv(gcta_sig_genes)
print(df)

gene_list= list()
for item in weights_dir.iterdir():
    weightFile = item.name
    gene = weightFile.replace(".wgt.RDat", "")
    gene_list.append(gene)


#missing = df[~df['Gene'].isin(gene_list)]
#print(missing[['Gene']])

for index, gene in df.loc[~df['Gene'].isin(gene_list), 'Gene'].items():
    print(index, gene)
