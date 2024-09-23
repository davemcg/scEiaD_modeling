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
parser.add_argument("--scanvi_model_out", default = False)
parser.add_argument("--scanvi_predict__on", default = "MajorCellType")
parser.add_argument("--scanvi_predict__out", default = "MCT_scANVI")
parser.add_argument("--scib", default = False)
# example when debegging on command line
# sys.argv = ['run_scvi.py','../../hs111.adata.solo.20240827.h5ad','/home/mcgaugheyd/git/scEiaD_modeling/data/hs111_mature_eye_ref_bcs.20240910-01.csv.gz','hs111_mature_eye_20240910_500hvg_10e_10l','scviModel.hs111_mature_eye_20240910__500hvg_10e_10l','hs111_mature_eye_20240910_500hvg_10e_10l.h5ad','hs111_mature_eye_20240910_500hvg_10e_10l.obs.csv.gz','--query_csv','/home/mcgaugheyd/git/scEiaD_modeling/data/hs111_mature_eye_query_bcs.20240910-01.csv.gz','--n_top_genes','500','--n_epochs','10','--n_latent','10','--input_type','barcode','--scanvi','scanviModel.hs111_mature_eye_20240910__500hvg_10e_10l','--scib','hs111_mature_eye_20240910_500hvg_10e_10l.scib.csv']
args = parser.parse_args()

adata = sc.read_h5ad(args.h5ad_input)
#adata.raw = adata
adata.obs['capture_study'] = [a + '_' + b for a,b in zip(list(adata.obs.capture_type),list(adata.obs.study_accession))]
adata.obs.MajorCellType = adata.obs.MajorCellType.cat.add_categories('unlabelled')
adata.obs.MajorCellType = adata.obs.MajorCellType.fillna('unlabelled')
print('adata loaded')


# calc ribo and protocadherin %
adata.var['ribo'] = adata.var['id_name'].str.startswith(('RPS', 'RPL'))
adata.var['protocadherin'] = adata.var['id_name'].str.startswith(('PCDH'))
sc.pp.calculate_qc_metrics(adata, qc_vars=['ribo'], percent_top=None, log1p=False, inplace=True)
sc.pp.calculate_qc_metrics(adata, qc_vars=['protocadherin'], percent_top=None, log1p=False, inplace=True)
print('metrics calculated')

# remove .digit from var_names
vn = adata.var_names
new_vn = vn.to_series().str.replace('\.\d+','',regex=True)
adata.var['ensembl'] = new_vn
adata.var_names = adata.var['ensembl']

def subset_adata(input_csv, input_type):
    partition_samples = pd.read_csv(input_csv, header = None)[0]

    if input_type == 'sample_accession':
        print('subset by sample accession')
        row_logical = np.array([s in partition_samples.values for s in adata.obs.sample_accession])
    else:
        print('subset by barcode')
        row_logical = partition_samples.values 

    adata_sub = adata[row_logical, ]
    adata_sub = adata_sub[~adata_sub.obs.solo_doublet].copy()
    adata_sub.layers["counts"] = adata_sub.X.copy()
    sc.pp.normalize_total(adata_sub)
    sc.pp.log1p(adata_sub)
    print('adata subset made')
    return(adata_sub)

def run_umap(adata_u):
    rsc.pp.neighbors(adata_u, use_rep=SCVI_LATENT_KEY)
    sc.tl.leiden(adata_u)
    sc.tl.leiden(adata_u, key_added = 'leiden2', resolution = 2)
    sc.tl.leiden(adata_u, key_added = 'leiden3', resolution = 3)
    sc.tl.leiden(adata_u, key_added = 'leiden5', resolution = 5)
    rsc.tl.umap(adata_u, min_dist = args.min_dist)
    
    umap = pd.DataFrame(adata_u.obsm['X_umap'])
    umap.index = adata_u.obs.index
    adata_u.obs['umap1'] = umap[0]
    adata_u.obs['umap2'] = umap[1]

    return(adata_u)
    
adata_ref = subset_adata(args.input_csv, args.input_type)
adata_ref.X = sparse.csr_matrix(adata_ref.X)
adata_ref.obs.columns
adata_ref.obs.capture_study

sc.pp.filter_cells(adata_ref, min_genes=args.gene_count_filter)

covariate = 'capture_study'

sc.pp.highly_variable_genes(
    adata_ref,
    flavor="seurat_v3",
   n_top_genes= args.n_top_genes,
    batch_key=covariate,
    subset=True, span=args.hvg_span
)

scvi.model.SCVI.setup_anndata(adata_ref,  layer="counts",
                                continuous_covariate_keys = ['pct_counts_mt','pct_counts_ribo','pct_counts_protocadherin'],
                              batch_key=covariate)


arches_params = dict(
    use_layer_norm="both",
    use_batch_norm="none",
    encode_covariates=True,
    dropout_rate=0.2,
    n_layers=2,
    n_latent = args.n_latent,
)
print(adata_ref)
vae_ref = scvi.model.SCVI(
    adata_ref,
    **arches_params
)
adata_ref
vae_ref.train(max_epochs = args.n_epochs, accelerator = 'gpu', early_stopping = True)
vae_ref

SCVI_LATENT_KEY = "X_scVI"
adata_ref.obsm[SCVI_LATENT_KEY] = vae_ref.get_latent_representation()

adata_ref_u = run_umap(adata_ref)
#adata_ref_u.obs.to_csv("temp.csv.gz")

# optional projection of new (query) data onto the ref model
if args.query_csv:
    print('query projection running')
    adata_query = subset_adata(args.query_csv, args.input_type)
    adata_query = adata_query[:, adata_ref.var_names].copy()
    
    scvi.model.SCVI.setup_anndata(adata_query,  layer="counts",
                                continuous_covariate_keys = ['pct_counts_mt','pct_counts_ribo','pct_counts_protocadherin'],
                              batch_key=covariate)
    
    vae_q = scvi.model.SCVI.load_query_data(
        adata_query,
        vae_ref,
        freeze_dropout = True
    )
    
    print(adata_query)
    vae_q.train(max_epochs=args.n_epochs, accelerator = 'gpu', plan_kwargs=dict(weight_decay=0.0))
    adata_query.obsm["X_scVI"] = vae_q.get_latent_representation()

    adata_full = ad.concat([adata_ref, adata_query])

    adata_full.obsm["X_scVI"] = vae_q.get_latent_representation(adata_full)
else:
    adata_full = adata_ref

adata_full = run_umap(adata_full)

if args.scanvi_model_out:
    print("\n\nRun scANVI CellType Prediction")
     
    print("\n\nRemove any NEW celltypes present in the query data but not the reference")
    
    scanvi_epochs = min(50, args.n_epochs)
    remove_ct = [mct for mct in set(adata_full.obs.MajorCellType) if mct not in set(adata_ref.obs.MajorCellType)]
    print(remove_ct)
    print(args.scanvi_predict__on)
    if 'unlabelled' not in adata_full.obs[args.scanvi_predict__on].values:
        adata_full.obs[args.scanvi_predict__on] = adata_full.obs[args.scanvi_predict__on].cat.add_categories('unlabelled') 
    adata_full.obs[args.scanvi_predict__on][adata_full.obs[args.scanvi_predict__on].isin(remove_ct)] = 'unlabelled'
    adata_full.obs[args.scanvi_predict__on] = pd.Categorical(adata_full.obs[args.scanvi_predict__on]).remove_unused_categories()    
    scanvi_model = scvi.model.SCANVI.from_scvi_model(
        vae_ref,
        unlabeled_category="unlabelled",
        labels_key=args.scanvi_predict__on,
    )
    scanvi_model.train(max_epochs=scanvi_epochs)
    scanvi_query = scvi.model.SCANVI.load_query_data(adata_full, scanvi_model)
    scanvi_query.train(
        max_epochs=scanvi_epochs,
        plan_kwargs={"weight_decay": 0.0},
        check_val_every_n_epoch=10,
    )
    SCANVI_LATENT_KEY = "X_scANVI"
    SCANVI_PREDICTION_KEY = args.scanvi_predict__out
    

    adata_full.obsm[SCANVI_LATENT_KEY] = scanvi_query.get_latent_representation(adata_full)
    adata_full.obs[SCANVI_PREDICTION_KEY] = scanvi_query.predict(adata_full)

    sc.pl.scatter(adata_full, size = 5, basis = 'umap', color = [args.scanvi_predict__out],  palette = sc.pl.palettes.vega_20_scanpy, save = args.plot_prefix + '_MCT_scANVI.png')

    scanvi_model.save(args.scanvi_model_out, overwrite=True, save_anndata =True)

sc.pp.log1p(adata_full)
sc.pl.scatter(adata_full, size = 5, basis = 'umap', color = [args.scanvi_predict__on],  palette = sc.pl.palettes.vega_20_scanpy, save = args.plot_prefix + '_MajorCellType.png')
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
# write hvg 
hvg = pd.DataFrame(adata_ref.var_names)
hvg.to_csv(args.scVI_model_output_path + '/hvg.csv.gz')

adata_full.write_h5ad(args.h5ad_output)
adata_full.obs.to_csv(args.obs_csv_output)

# benchmarker
if args.scib:
    sc.tl.pca(adata_full)
    bm = Benchmarker(
        adata_full,
        batch_key=covariate,
        label_key="MajorCellType",
        embedding_obsm_keys=["X_pca", "X_scVI"],
        n_jobs=2,
    )

    bm.benchmark()

    bm_res = bm.get_results(min_max_scale=False)

    bm_res.to_csv(args.scib)


