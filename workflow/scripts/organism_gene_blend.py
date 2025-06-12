import scanpy as sc
import pandas as pd

adata = sc.read_h5ad()

gene_table = pd.read_table('/home/mcgaugheyd/git/scEiaD_modeling/data/hs111_mm111_ensembl_biomart_hcop.txt', sep = '\t')
fgt = gene_table[~gene_table['Mouse gene stable ID'].isna()]
mm_hs_vn = [x for x in adata.var_names if x in fgt['Gene stable ID version'].values]
adata[:,mm_hs_vn].write_h5ad('hs111.adata.solo.20250204.dev.mmGeneFilter.h5ad')
