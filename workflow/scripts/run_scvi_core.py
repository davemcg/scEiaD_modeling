import argparse
import os
import numpy as np
import pandas as pd
import scanpy as sc
import scvi
import torch
import anndata as ad
from scipy import sparse
import rapids_singlecell as rsc
from scib_metrics.benchmark import Benchmarker

def main():
    parser = argparse.ArgumentParser(description='Stand-alone scVI/scANVI Integration Script')
    
    # Required Paths
    parser.add_argument("h5ad_input", help="Path to input .h5ad file")
    parser.add_argument("h5ad_output", help="Path to save processed .h5ad")
    parser.add_argument("model_output", help="Directory to save scVI model")
    
    # Metadata / Batch Options
    parser.add_argument("--batch_keys", default="sample_accession", 
                        help="Comma-separated obs columns to combine into a single batch covariate")
    parser.add_argument("--continuous_covariates", default=None, 
                        help="Comma-separated obs columns for continuous covariates")
    
    # Model Hyperparameters
    parser.add_argument("--n_top_genes", default=2000, type=int)
    parser.add_argument("--n_epochs", default=200, type=int)
    parser.add_argument("--n_latent", default=30, type=int)
    parser.add_argument("--gene_filter", default=300, type=int)
    
    # HVG generation
    parser.add_argument("--hvg_span", default = 0.5, type = float)
    # Optional Modules & Outputs
    parser.add_argument("--scanvi_label", default=None, help="Obs column for scANVI training")
    parser.add_argument("--scanvi_model_output", default=None, help="Directory to save scANVI model")
    parser.add_argument("--scib_output", default=None, help="Path to save scIB metrics CSV")
    parser.add_argument("--obs_csv_output", default=None, help="Optional path to save obs+obsm CSV")
    parser.add_argument("--plot_prefix", default=None, help="Saves UMAP plots with this prefix")

    args = parser.parse_args()

    # 1. Load and Prepare
    adata = sc.read_h5ad(args.h5ad_input)
    
    # Dynamic Batch Construction
    batch_cols = args.batch_keys.split(',')
    adata.obs['batch_covariate'] = (
        adata.obs[batch_cols].astype(str).agg('_'.join, axis=1).astype('category')
    )
    
    sc.pp.filter_cells(adata, min_genes=args.gene_filter)
    adata.layers["counts"] = adata.X.copy()
    
    sc.pp.highly_variable_genes(
        adata,
        flavor="seurat_v3",
        n_top_genes=args.n_top_genes,
        batch_key='batch_covariate',
        subset=True,
        layer="counts", 
        span = args.hvg_span
    )

    # 2. scVI Training
    cont_cov = args.continuous_covariates.split(',') if args.continuous_covariates else None
    scvi.model.SCVI.setup_anndata(
        adata, layer="counts", batch_key='batch_covariate', continuous_covariate_keys=cont_cov
    )

    vae = scvi.model.SCVI(
        adata, n_layers=2, n_latent=args.n_latent, encode_covariates=True,
        use_layer_norm="both", use_batch_norm="none"
    )

    vae.train(max_epochs=args.n_epochs, accelerator='gpu', early_stopping=True)
    adata.obsm["X_scVI"] = vae.get_latent_representation()
    vae.save(args.model_output, overwrite=True, save_anndata=True)

    # 3. Optional scANVI Training
    if args.scanvi_label:
        if 'unlabelled' not in adata.obs[args.scanvi_label].dtype.categories:
            adata.obs[args.scanvi_label] = adata.obs[args.scanvi_label].cat.add_categories(['unlabelled'])
        adata.obs[args.scanvi_label] = adata.obs[args.scanvi_label].fillna('unlabelled')
        
        lvae = scvi.model.SCANVI.from_scvi_model(vae, unlabeled_category="unlabelled", labels_key=args.scanvi_label)
        lvae.train(max_epochs=min(20, args.n_epochs))
        
        adata.obsm["X_scANVI"] = lvae.get_latent_representation()
        adata.obs["predictions_scANVI"] = lvae.predict()
        
        if args.scanvi_model_output:
            lvae.save(args.scanvi_model_output, overwrite=True, save_anndata=True)

    # 4. Dimensionality Reduction
    rsc.pp.neighbors(adata, use_rep="X_scVI")
    rsc.tl.umap(adata, min_dist=0.3)
    
    # 5. Output Obs + Obsm CSV
    if args.obs_csv_output:
        print(f"Exporting metadata to {args.obs_csv_output}...")
        obs_df = adata.obs.copy()
        
        # Flatten obsm (scVI, scANVI, UMAP) into columns
        for key in adata.obsm_keys():
            df_obsm = pd.DataFrame(
                adata.obsm[key], 
                index=adata.obs_names,
                columns=[f"{key}_{i+1}" for i in range(adata.obsm[key].shape[1])]
            )
            obs_df = pd.concat([obs_df, df_obsm], axis=1)
        
        obs_df.to_csv(args.obs_csv_output)

    # 6. Finalization
    if args.plot_prefix:
        sc.pl.umap(adata, color=['batch_covariate'], show=False, save=f"_{args.plot_prefix}.png")

    adata.write_h5ad(args.h5ad_output)

if __name__ == "__main__":
    main()
