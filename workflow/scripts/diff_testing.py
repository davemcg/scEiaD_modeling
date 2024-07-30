import scanpy as sc
import pandas as pd
import argparse

parser = argparse.ArgumentParser(description = 'Run scanpy diff testing with input h5ad and new obs')
parser.add_argument("big_h5ad")
parser.add_argument("obs", help = "csv file from run_scvi that has the cluster column info and cut \
                                            down cell barcodes")
parser.add_argument("category", help = "column from obs to do diff testing on")
parser.add_argument("--method", default = 'wilcoxon')
parser.add_argument("out_table")
args = parser.parse_args()

adata = sc.read_h5ad(args.big_h5ad) # 'hs111.adata.solo.2024_03_07.h5ad')
obs = pd.read_csv(args.obs) # 'hs111_mature_alleye_5000hvg_50e_50l.obs.csv.gz')

adata_filt = adata[obs['barcode'],]
obs.index = adata_filt.obs.index
adata_filt.obs[args.category] = obs[args.category].astype('category')
sc.pp.log1p(adata_filt)

sc.tl.rank_genes_groups(adata_filt, args.category,  method=args.method, key_added = args.method)

res_list = []
for i in adata_filt.obs[args.category].unique():
    res = sc.get.rank_genes_groups_df(adata_filt, group = str(i), key = args.method)
    res.insert(2, 'base', str(i))
    res_list.append(res)

res_big = pd.concat(res_list)
res_big.to_csv(args.out_table)
