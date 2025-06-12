import argparse
parser = argparse.ArgumentParser(description='Use scVI models to predict CT on new adata')
parser.add_argument("h5ad_query", help = 'Input adata')
parser.add_argument("model_tsv", help = 'tsv file (NO HEADER) with two columns, the first being the name of the model and the second being the path to the scVI model directory')
parser.add_argument("obs_output", help = 'Output a csv file with the original obs slot plus the CT calls and CT numeric score')
parser.add_argument("--h5ad_output", default = False, help = 'Optionally output the h5ad with the updated obsm and obs slots')
parser.add_argument("--gene_convert_file", default = False, help = "Default is False, optional tsv which has a key column which contains the gene type used in the adata (names or ENSEMBL usually) and another column which has the alternate id (another species usually) ")
parser.add_argument("--gene_convert_key", help = "key is the column index (0 based) from the --gene_convert_file of the gene format (e.g. mouse) used in the input adata")
parser.add_argument("--gene_convert_value", help = "value is the column index (0 based) from the --gene_convert_file of the converted id (e.g. human) to relabel the input adata")
parser.add_argument("--barcode", default = False, help = "Give newline separated file of barcodes to retain (no header)")
args = parser.parse_args()


def subset_adata(input_csv, input_type):
    partition_samples = pd.read_csv(input_csv, header = None)[0]

    if input_type == 'sample_accession':
        print('subset by sample accession')
        row_logical = np.array([s in partition_samples.values for s in adata.obs.sample_accession])
    else:
        print('subset by barcode')
        row_logical = partition_samples.values

    adata_sub = adata[row_logical, ].copy()
    return(adata_sub)


import matplotlib
import sys
import os
import numpy as np
import pandas as pd
import random
import scanpy as sc
sc.settings.autoshow = False
from scipy import sparse
import scvi
import torch
import anndata as ad
import cupy as cp
from scipy import sparse

import rapids_singlecell as rsc
import rmm
from rmm.allocators.cupy import rmm_cupy_allocator
from scib_metrics.benchmark import Benchmarker

print("\n\nscvi version is: " + scvi.__version__ + "\n\n")

rmm.reinitialize(
    managed_memory=False, # Allows oversubscription
    pool_allocator=False, # default is False
    devices=0, # GPU device IDs to register. By default registers only GPU 0.
)
cp.cuda.set_allocator(rmm_cupy_allocator)



# import query adata
#adata = sc.read_h5ad('../../mm111.adata.solo.20240827.h5ad')
adata = sc.read_h5ad(args.h5ad_query)
print(args.h5ad_query + " loaded\n")
if args.barcode:
    adata = subset_adata(args.barcode, 'barcode')
    print("adata custom subsetted with " + args.barcode)
# set up query object
# plae models are built around ENSGdigit notation
# and the .digit ending is removed from var_names
# e.g. ENSG01931234.1 is truncated to ENSG01931234
vn = adata.var_names
new_vn = vn.to_series().str.replace('\.\d+','',regex=True)
adata.var['ensembl'] = new_vn
adata.var_names = adata.var['ensembl']
adata.layers["counts"] = adata.X.copy()
adata.obs['capture_covariate'] = [a + '_' + b for a,b in zip(list(adata.obs.capture_type),list(adata.obs.study_accession))]
#########
# if a different species or
# an adata with symbols is used 
#########
if args.gene_convert_file:
    # import gene conversion table
    #gene_table = pd.read_table('/home/mcgaugheyd/git/scEiaD_modeling/data/hs111_mm111_ensembl_biomart_hcop.txt', sep = '\t')
    gene_table = pd.read_table(args.gene_convert_file, sep = '\t')
    # gene_dict  = gene_table.set_index('Mouse gene stable ID').T.to_dict('list')
    ## key is the column index (0 based) of the gene format (e.g. mouse) used in the input adata
    ## value is the column index (0 based) of the converted id (e.g. human) to relabel the input adata
    key = int(args.gene_convert_key)
    value = int(args.gene_convert_value)
    gene_dict  = gene_table.set_index(gene_table.columns[key]).T.to_dict('list')

    vn = adata.var_names
    digit_clean_vn = vn.to_series().str.replace('\.\d+','',regex=True)
    humanized_vn = []
    for gene in digit_clean_vn:
        try:
            humanized_vn.append(gene_dict[gene][value])
        except:
            humanized_vn.append('fail')

    adata.var['ensembl'] = digit_clean_vn
    adata.var['converted'] = humanized_vn
    adata.var_names = adata.var['converted']
    # aggregate (requires scanpy >= 1.10)
    adata_agg = sc.get.aggregate(adata, by = 'converted', axis = 1, func = 'sum')
    adata_agg.layers["counts"] = sparse.csr_matrix(adata_agg.layers['sum'])
    del adata_agg.layers['sum']
    adata_agg.X = adata_agg.layers['counts']
    print("finished with gene conversion")
else:
    adata_agg = adata.copy()

# models
model_input = pd.read_table(args.model_tsv, header = None, sep = '\t')
models = model_input.set_index(model_input.columns[0]).T.to_dict('list')
#models['human_nonneural'] = '/data/OGVFB_BG/scEiaD/2024_02_28/snakeout/hs111_mature_eye_full/NONneural_cells_s6__cleanmodel/scanviModel.hs111_mature_eye_20240924_full__NONneural_s6___2000hvg_200e_30l/'
#models['human_neural'] = '/data/OGVFB_BG/scEiaD/2024_02_28/snakeout/hs111_mature_eye_full/neural_cells_s6_clean_model/scanviModel.hs111_mature_eye_20240924_full__neural_s6___2000hvg_200e_30l'
#models['human_full'] = '/data/OGVFB_BG/scEiaD/2024_02_28/snakeout/hs111_mature_eye_full/full_s6_clean_model/scanviModel.hs111_mature_eye_20241121_stage6__2000hvg_200e_30l/'
#models['hrca'] = '/data/OGVFB_BG/scEiaD/2024_02_28/chen_rca/scanviModel.hrca_all_20241215__2000hvg_200e_30l'
#models['mrca'] = '/data/OGVFB_BG/scEiaD/2024_02_28/chen_rca/scanviModel.mrca_all_20241215__2000hvg_200e_30l'

for model in models.keys():
    print(model)
    print(models[model])
    adata_loop = adata_agg.copy()
    scanvi_model = scvi.model.SCANVI.load(models[model][0])
    scvi.model.SCANVI.prepare_query_anndata(adata_loop, scanvi_model)
    scanvi_query = scvi.model.SCANVI.load_query_data(adata_loop, scanvi_model)
    adata_agg.obsm["CT__" + model +  "__soft"] = scanvi_query.predict(adata_loop, soft = True)

    scanvi_query.train(max_epochs=100, plan_kwargs={"weight_decay": 0.0})
    adata_agg.obsm["scanvi_" + model] = scanvi_query.get_latent_representation()
    rsc.pp.neighbors(adata_agg, use_rep="scanvi_" + model)
    rsc.tl.leiden(adata_agg, key_added = 'leiden')
    rsc.tl.umap(adata_agg, min_dist=0.3)
    
    umap = pd.DataFrame(adata_agg.obsm['X_umap'])
    umap.index = adata_agg.obs.index
    adata_agg.obs['umap1_' + model] = umap[0]
    adata_agg.obs['umap2_' + model] = umap[1]

    adata_agg.obs["CT__" + model] = scanvi_query.predict(adata_loop)
    adata_agg.obs["CT__" + model + "__max_score"] = adata_agg.obsm["CT__" + model +  "__soft"].max(axis=1)

if args.h5ad_output:
    adata_agg.write_h5ad(args.h5ad_output)

pd.DataFrame(adata_agg.obs).to_csv(args.obs_output)
