import pandas as pd
import scanpy as sc
# /data/OGVFB_BG/scEiaD/2024_02_28/snakeout/hs111_mature_eye_full

adata =sc.read_h5ad('hs111_mature_eye_20241121_full_2000hvg_200e_30l_stage6.h5ad')
#vn = adata.var_names
#new_vn = vn.to_series().str.replace('\.\d+','',regex=True)
#adata.var['ensembl'] = new_vn
#adata.var_names = adata.var['ensembl']

gene_table = pd.read_table('/home/mcgaugheyd/git/scEiaD_modeling/data/hs111_mm111_ensembl_biomart_hcop.txt', sep = '\t')
fgt = gene_table[~gene_table['Mouse gene stable ID'].isna()]
mm_hs_vn = [x for x in adata.var_names if x in fgt['Gene stable ID version'].values]

nadata = adata[:,mm_hs_vn]

adata[:,mm_hs_vn].write_h5ad('hs111_mature_eye_20241121_full_2000hvg_200e_30l_stage6.mmGeneFilter.h5ad')
