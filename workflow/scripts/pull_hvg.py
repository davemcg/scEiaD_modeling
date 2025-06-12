import scanpy as sc
import pandas as pd
import argparse

parser = argparse.ArgumentParser(description = 'Output a csv file of var names from an adata object')
parser.add_argument("h5ad", help = 'h5ad file from scVI')
parser.add_argument("hvg", help = "output hvg csv file")

args = parser.parse_args()

adata = sc.read_h5ad(args.h5ad)

hvg = pd.DataFrame(adata.var_names)
hvg.to_csv(args.hvg)
