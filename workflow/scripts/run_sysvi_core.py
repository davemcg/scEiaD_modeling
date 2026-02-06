import argparse
import os
import numpy as np
import pandas as pd
import scanpy as sc
import scvi
from scvi.external import SysVI
import torch
import anndata as ad
from scipy import sparse
import rapids_singlecell as rsc

def main():
    parser = argparse.ArgumentParser(description='Stand-alone SysVI with Integrated Ref/Query Sampling')
    
    # Required Paths
    parser.add_argument("h5ad_input", help="Path to input .h5ad file")
    parser.add_argument("h5ad_output", help="Path to save processed .h5ad")
    parser.add_argument("model_output", help="Directory to save SysVI reference model")
    
    # Sampling / Ref-Query Options
    parser.add_argument("--max_cells_per_batch", default=None, type=int, 
                        help="If set, samples this many cells per batch for the reference; others become query.")
    parser.add_argument("--automax", action='store_true', 
                        help="If True, uses the median count of the batch field as the max_cells_per_batch.")
    
    # Metadata / Batch Options
    parser.add_argument("--batch_key", default="tech_capture", 
                        help="Primary batch covariate(s). Can be a comma-separated list.")
    parser.add_argument("--categorical_covariates", default=None, 
                        help="Comma-separated obs columns for additional systemic covariates")
    
    # Model Hyperparameters
    parser.add_argument("--n_top_genes", default=2000, type=int)
    parser.add_argument("--n_epochs", default=200, type=int)
    parser.add_argument("--n_latent", default=15, type=int)
    parser.add_argument("--gene_filter", default=300, type=int, help="Min genes per cell")
    parser.add_argument("--hvg_span", default=0.3, type=float, help="Span parameter for Seurat v3 HVG selection")
    
    # Optional Modules & Outputs
    parser.add_argument("--obs_csv_output", default=None, help="Optional path to save obs+obsm CSV")
    parser.add_argument("--marker_csv_output", default=None, help="Optional path to save marker genes CSV")

    args = parser.parse_args()

    # 1. Load Data and Build Batch Field
    print(f"Loading {args.h5ad_input}...")
    adata = sc.read_h5ad(args.h5ad_input)
    
    # Handle comma-separated batch keys
    batch_cols = args.batch_key.split(',')
    if len(batch_cols) > 1:
        print(f"Combining batch columns: {batch_cols}")
        adata.obs['combined_batch'] = (
            adata.obs[batch_cols].astype(str).agg('_'.join, axis=1).astype('category')
        )
        working_batch_key = 'combined_batch'
    else:
        working_batch_key = args.batch_key

    # 2. Integrated Sampling Logic (Ref vs Query)
    if args.max_cells_per_batch or args.automax:
        print("Sampling data into Reference and Query sets...")
        if args.automax:
            max_count = int(adata.obs[working_batch_key].value_counts().median().round())
        else:
            max_count = args.max_cells_per_batch

        ref_obs = adata.obs.groupby(working_batch_key, observed=True).sample(
            n=max_count, random_state=1, replace=True
        ).drop_duplicates()
        
        adata_ref = adata[ref_obs.index].copy()
        adata_query = adata[~adata.obs_names.isin(ref_obs.index)].copy()
        print(f"Ref size: {adata_ref.n_obs}, Query size: {adata_query.n_obs}")
    else:
        print("Using full dataset as Reference.")
        adata_ref = adata.copy()
        adata_query = None

    # 3. Reference Preprocessing
    sc.pp.filter_cells(adata_ref, min_genes=args.gene_filter)
    adata_ref.layers["counts"] = adata_ref.X.copy()
    
    sc.pp.highly_variable_genes(
        adata_ref, 
        flavor="seurat_v3", 
        n_top_genes=args.n_top_genes,
        batch_key=working_batch_key, 
        subset=True, 
        layer="counts",
        span=args.hvg_span
    )

    # 4. Setup SysVI
    cat_covs = args.categorical_covariates.split(',') if args.categorical_covariates else None

    sc.pp.normalize_total(adata_ref)
    sc.pp.log1p(adata_ref)
 
    SysVI.setup_anndata(
        adata_ref, 
        batch_key=working_batch_key, 
        categorical_covariate_keys=cat_covs
    )

    vae_ref = SysVI(
        adata_ref,
        n_latent=args.n_latent,
        embed_categorical_covariates=True
    )

    print("Training Reference SysVI Model...")
    vae_ref.train(
        max_epochs=args.n_epochs, 
        accelerator='gpu',
        plan_kwargs={"z_distance_cycle_weight": 10}
    )
    vae_ref.save(args.model_output, overwrite=True, save_anndata=True)

    # 5. Query Projection
    if adata_query is not None and adata_query.n_obs > 0:
        print("Projecting Query data...")
        adata_query = adata_query[:, adata_ref.var_names].copy()
        adata_query.layers["counts"] = adata_query.X.copy()

        sc.pp.normalize_total(adata_query)
        sc.pp.log1p(adata_query)
        
        SysVI.setup_anndata(
            adata_query, layer="counts", batch_key=working_batch_key, 
            categorical_covariate_keys=cat_covs
        )
        
        vae_q = SysVI.load_query_data(adata_query, vae_ref)
        vae_q.train(max_epochs=args.n_epochs, accelerator='gpu', plan_kwargs=dict(weight_decay=0.0))
        
        adata_full = ad.concat([adata_ref, adata_query])
        adata_full.obsm["X_sysVI"] = vae_q.get_latent_representation(adata_full)
    else:
        adata_full = adata_ref
        adata_full.obsm["X_sysVI"] = vae_ref.get_latent_representation()

    # 6. Dimensionality Reduction (RAPIDS)
    print("Running Neighbors and UMAP...")
    rsc.pp.neighbors(adata_full, use_rep="X_sysVI")
    rsc.tl.umap(adata_full, min_dist=0.3)
    
    for res in [0.5, 1, 2]:
        key = f"leiden_res_{str(res).replace('.', '_')}"
        rsc.tl.leiden(adata_full, resolution=res, key_added=key)

    # 7. Exports
    if args.obs_csv_output: 
        obs_df = adata_full.obs.copy()
        for key in adata_full.obsm_keys():
            df_obsm = pd.DataFrame(
                adata_full.obsm[key], index=adata_full.obs_names,
                columns=[f"{key}_{i+1}" for i in range(adata_full.obsm[key].shape[1])]
            )
            obs_df = pd.concat([obs_df, df_obsm], axis=1)
        obs_df.to_csv(args.obs_csv_output)

    print(f"Saving final object to {args.h5ad_output}...")
    adata_full.write_h5ad(args.h5ad_output)

if __name__ == "__main__":
    main()
