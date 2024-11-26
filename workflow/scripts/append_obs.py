# update obs slot in adata 
import scanpy as sc
import pandas as pd
import argparse
parser = argparse.ArgumentParser(description = "update obs slot in adata")
parser.add_argument("h5ad_input", help = "input h5ad file")
parser.add_argument("obs_csv", help = "csv file of updated obs")
parser.add_argument("--index_column", default = 'barcode', help = "column from the obs to use as the index")
parser.add_argument("--transfer_columns", default = 'leiden3,MCT_scANVI', help = "column names to transfer over (comma separate)")
parser.add_argument("--neural", default = False, help = 'Give name of column to remove non-neural cell type calls')
parser.add_argument("--nonneural", default = False, help = 'Give name of column to remove neural cell type calls')
parser.add_argument("h5ad_output", help = "output h5ad name")
args = parser.parse_args()

# input
adata = sc.read_h5ad(args.h5ad_input)
obs_csv = pd.read_csv(args.obs_csv)

# create list of select cols
transfer_cols = args.transfer_columns.split(',')

# set index 
new_obs = obs_csv
new_obs.index = new_obs[args.index_column]
new_obs = new_obs.drop(args.index_column, axis = 1)

# pull original obs info
old_obs = pd.DataFrame(adata.obs)

# left join
updated_obs = old_obs.join(new_obs[transfer_cols])

# transfer columns over
adata.obs[transfer_cols] = updated_obs[transfer_cols]
adata.obs = adata.obs.drop(args.index_column, axis = 1)

#  
neural_ct = ['rod','cone','horizontal','bipolar','rod bipolar', 'amacrine','retinal ganglion']
if args.neural:
   print('neural out')
   adata.obs[args.neural][~adata.obs[args.neural].isin(neural_ct)] = 'unlabelled'
if args.nonneural:
    print('non neural out')
    adata.obs[args.nonneural][adata.obs[args.nonneural].isin(neural_ct)] = 'unlabelled'

# filter adata down to cells that have a value in the first column given
adata[~adata.obs[transfer_cols[0]].isna()].write_h5ad(args.h5ad_output)
