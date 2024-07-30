import scanpy as sc
import pandas as pd
import argparse
from adpbulk import ADPBulk

parser = argparse.ArgumentParser(description = 'Use adp bulk to make pseudobulk matrix')
parser.add_argument("big_h5ad")
parser.add_argument("obs", help = "csv file from run_scvi that has the cluster column info and cut \
                                            down cell barcodes")
parser.add_argument("category", help = "comma separated (no whitespace) list of columns to aggregate against")
parser.add_argument("out_table")
args = parser.parse_args()

adata = sc.read_h5ad(args.big_h5ad) # 'hs111.adata.solo.2024_03_07.h5ad')
obs = pd.read_csv(args.obs) # 'hs111_mature_alleye_5000hvg_50e_50l.obs.csv.gz')

adata_filt = adata[obs['barcode'],]
obs.index = adata_filt.obs.index

the_cats = args.category.split(',')

if len(the_cats) == 2:
    one = the_cats[0]
    two = the_cats[1]
    adata_filt.obs["cat"] = adata_filt.obs[one].astype(str) +"_" + adata_filt.obs[two].astype(str)
else:
    one = the_cats[0]
    adata_filt.obs[one] = obs[one].astype('category')
    adata_filt.obs["cat"] = adata_filt.obs[one].astype(str)
    

adpb = ADPBulk(adata_filt, ["cat"])
pseudobulk_matrix = adpb.fit_transform()
pseudobulk_matrix.to_csv(args.out_table)
