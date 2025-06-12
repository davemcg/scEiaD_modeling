# takes in an h5ad adata obj and outputs "ref" and "query" barcodes
# for use in the scEiaD_modeling pipeline
# often don't want to use the entire dataset for building the scEiaD
# model as:
#  1. data may be too huge (hrca)
#  2. ct are usually imbalanced (rods are VERY common)

import scanpy as sc
import pandas as pd
import argparse
###########################
parser = argparse.ArgumentParser(description = "Create ref and query barcode csv")
parser.add_argument("h5ad")
parser.add_argument("--groupby", help = "obs slot column which is used to balance against. Like \"Cell Type\" which can be used to ensure that you don't have a bazillion rods")
parser.add_argument("--max", default = 10000, help = "max number of items from each value in the groupby column OR (if you just randomly sample) the  number of barcodes returned")
parser.add_argument("--automax", default = False, help = "Default is False, if True then will pick the number based on the mean of the groupby counts")
parser.add_argument("ref_csv_file")
parser.add_argument("--query_csv_file", help = "output the remainder of the barcodes not used in the ref")

args = parser.parse_args()
############################

adata = sc.read_h5ad(args.h5ad)
obs = pd.DataFrame(adata.obs)

if args.automax:
    max_count = int(obs[args.groupby].value_counts().median().round())
else:
    max_count = int(args.max)

if args.groupby:
    out_obs = obs.groupby(args.groupby, observed = True).sample(n=max_count, 
                random_state=1, replace = True).drop_duplicates()
else:
    out_obs = obs.sample(n=max_count, random_state=1, replace = False)


out_obs['index'] = out_obs.index
out_obs['index'].to_csv(args.ref_csv_file, header = False, index = False)

if args.query_csv_file:
    query_obs = obs.iloc[~obs.index.isin(out_obs.index),:]
    query_obs['index'] = query_obs.index
    query_obs['index'].to_csv(args.query_csv_file, header = False, index = False)
