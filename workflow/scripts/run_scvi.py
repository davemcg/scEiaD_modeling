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

import argparse
parser = argparse.ArgumentParser(description='Run scVI with input h5ad and samples to cut down to')
parser.add_argument("h5ad_input")
parser.add_argument("input_csv")
parser.add_argument("plot_prefix")
parser.add_argument("scVI_model_output_path")
parser.add_argument("h5ad_output")
parser.add_argument("obs_csv_output")
parser.add_argument("--input_type", default = "sample_accession")
parser.add_argument("--query_csv", default = False)
parser.add_argument("--n_top_genes", default = 2000, type = int)
parser.add_argument("--n_epochs", default = 50, type = int)
parser.add_argument("--n_latent", default = 30, type = int)
parser.add_argument("--min_dist", default = 0.3, type = float)
parser.add_argument("--hvg_span", default = 0.5, type = float)
parser.add_argument("--gene_count_filter", default = 300, type = int)
args = parser.parse_args()

adata = sc.read_h5ad(args.h5ad_input)
print('adata loaded')

# calc ribo and protocadherin %
adata.var['ribo'] = adata.var['id_name'].str.startswith(('RPS', 'RPL'))
adata.var['protocadherin'] = adata.var['id_name'].str.startswith(('PCDH'))
sc.pp.calculate_qc_metrics(adata, qc_vars=['ribo'], percent_top=None, log1p=False, inplace=True)
sc.pp.calculate_qc_metrics(adata, qc_vars=['protocadherin'], percent_top=None, log1p=False, inplace=True)
print('metrics calculated')

def subset_adata(input_csv):
    partition_samples = pd.read_csv(input_csv, header = None)[0]

    if args.input_type == 'sample_accession':
        print('subset by sample accession')
        row_logical = np.array([s in partition_samples.values for s in adata.obs.sample_accession])
    else:
        print('subset by barcode')
        row_logical = partition_samples.values 

    adata_ref = adata[row_logical, ].copy()
    print('adata subset made')
    adata_ref

adata_ref = subset_adata(args.input_csv)
adata_ref = adata_ref[~adata_ref.obs.solo_doublet].copy()
adata_ref.X = sparse.csr_matrix(adata_ref.X)

sc.pp.filter_cells(adata_ref, min_genes=args.gene_count_filter)


covariate = 'study_accession'

sc.pp.highly_variable_genes(
    adata_ref,
    flavor="seurat_v3",
   n_top_genes= args.n_top_genes,
    batch_key=covariate,
    subset=True, span=args.hvg_span
)

scvi.model.SCVI.setup_anndata(adata_ref,
                                continuous_covariate_keys = ['pct_counts_mt','pct_counts_ribo','pct_counts_protocadherin'],
                              categorical_covariate_keys=['capture_type', covariate, 'sample_accession'])


arches_params = dict(
    use_layer_norm="both",
    use_batch_norm="none",
    encode_covariates=True,
    dropout_rate=0.2,
    n_layers=2,
    n_latent = args.n_latent,
)
adata_ref
vae_ref = scvi.model.SCVI(
    adata_ref,
    **arches_params
)

vae_ref.train(max_epochs = args.n_epochs, accelerator = 'gpu', early_stopping = True)
vae_ref

SCVI_LATENT_KEY = "X_scVI"
adata_ref.obsm[SCVI_LATENT_KEY] = vae_ref.get_latent_representation()

# optional projection of new (query) data onto the ref model
if args.query_csv:
    adata_query = subset_adata(args.query_csv)
    vae_q = scvi.model.SCVI.load_query_data(
        adata_query,
        vae_ref
    )

    vae_q.train(max_epochs=n_epochs, use_gpu=useCuda, plan_kwargs=dict(weight_decay=0.0))
    adata_query.obsm["X_scVI"] = vae_q.get_latent_representation()

    adata_full = adata_query.concatenate(adata_ref, batch_key = 'bkey')

    adata_full.obsm["X_scVI"] = vae_q.get_latent_representation(adata_full)
else:
    adata_full = adata_ref

sc.pp.neighbors(adata_full, use_rep=SCVI_LATENT_KEY)
sc.tl.leiden(adata_full)
sc.tl.leiden(adata_full, key_added = 'leiden2', resolution = 2)
sc.tl.umap(adata_full, min_dist = args.min_dist)

vae_ref.save(args.scVI_model_output_path, overwrite=True, save_anndata =True)
adata_full.write_h5ad(args.h5ad_output)

umap = pd.DataFrame(adata_full.obsm['X_umap'])
umap.index = adata_full.obs.index
adata_full.obs['umap1'] = umap[0]
adata_full.obs['umap2'] = umap[1]

adata_full.obs.to_csv(args.obs_csv_output)
sc.pp.log1p(adata_full)
sc.pl.scatter(adata_full, size = 5, basis = 'umap', color = ['MajorCellType'],  palette = sc.pl.palettes.vega_20_scanpy, save = args.plot_prefix + '_MajorCellType.png')
sc.pl.scatter(adata_full, size = 5, basis = 'umap', color = ['SubCellType'],  palette = sc.pl.palettes.vega_20_scanpy, save = args.plot_prefix + '_SubCellType.png')
#sc.pl.scatter(adata_full, size = 5, basis='umap', color = 'ENSG00000163914.5', use_raw=False, show = False, save = args.plot_prefix + '_rho.png')
sc.pl.scatter(adata_full, size = 5, basis = 'umap', color = ['CellType'],  palette = sc.pl.palettes.vega_20_scanpy, save = args.plot_prefix + '_CellType.png')
sc.pl.scatter(adata_full, size = 5, basis = 'umap', color = ['capture_type'],  palette = sc.pl.palettes.vega_20_scanpy, save = args.plot_prefix + '_capture_type.png')
sc.pl.scatter(adata_full, size = 5, basis = 'umap', color = ['study_accession'],  palette = sc.pl.palettes.vega_20_scanpy, save = args.plot_prefix + '_study_accession.png')
sc.pl.scatter(adata_full, size = 5, basis = 'umap', color = ['age'],  palette = sc.pl.palettes.vega_20_scanpy, save = args.plot_prefix + '_Age.png')
sc.pl.scatter(adata_full, size = 5, basis = 'umap', color = ['leiden'],  palette = sc.pl.palettes.vega_20_scanpy, save = args.plot_prefix +'_leiden.png')
sc.pl.scatter(adata_full, size = 5, basis = 'umap', color = ['background_fraction'], save = args.plot_prefix + '_background_fraction.png')
sc.pl.scatter(adata_full, size = 5, basis = 'umap', color = ['pct_counts_mt'], save = args.plot_prefix + '_pct_counts_mt.png')
adata_full.obs['log_n_counts'] = np.log1p(adata_full.obs.n_counts)
sc.pl.scatter(adata_full, size = 5, basis = 'umap', color = ['log_n_counts'], save = args.plot_prefix + '_log_n_counts.png')

vae_ref.save(args.scVI_model_output_path, overwrite=True, save_anndata =True)
adata_full.write_h5ad(args.h5ad_output)

umap = pd.DataFrame(adata_full.obsm['X_umap'])
umap.index = adata_full.obs.index
adata_full.obs['umap1'] = umap[0]
adata_full.obs['umap2'] = umap[1]

adata_full.obs.to_csv(args.obs_csv_output)

