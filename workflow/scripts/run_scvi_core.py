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
    parser = argparse.ArgumentParser(description='Stand-alone scVI/scANVI with Integrated Ref/Query Sampling')
    
    # Required Paths
    parser.add_argument("h5ad_input", help="Path to input .h5ad file")
    parser.add_argument("h5ad_output", help="Path to save processed .h5ad")
    parser.add_argument("model_output", help="Directory to save scVI reference model")
    
    # Sampling / Ref-Query Options
    parser.add_argument("--max_cells_per_batch", default=None, type=int, 
                        help="If set, samples this many cells per batch for the reference; others become query.")
    parser.add_argument("--automax", action='store_true', 
                        help="If True, uses the median count of the batch field as the max_cells_per_batch.")
    
    # Metadata / Batch Options
    parser.add_argument("--batch_keys", default="sample_accession", 
                        help="Comma-separated obs columns to combine into a single batch covariate")
    parser.add_argument("--continuous_covariates", default=None, 
                        help="Comma-separated obs columns for continuous covariates")
    
    # Model Hyperparameters
    parser.add_argument("--n_top_genes", default=2000, type=int)
    parser.add_argument("--n_epochs", default=200, type=int)
    parser.add_argument("--n_latent", default=20, type=int)
    parser.add_argument("--gene_filter", default=300, type=int, help = "Minimum number of detected genes per cell")
    parser.add_argument("--hvg_span", default=0.5, type=float)
    
    # Optional Modules & Outputs
    parser.add_argument("--scanvi_label", default=None, help="Obs column for scANVI training")
    parser.add_argument("--scanvi_model_output", default=None, help="Directory to save scANVI model")
    parser.add_argument("--obs_csv_output", default=None, help="Optional path to save obs+obsm CSV")
    parser.add_argument("--marker_csv_output", default=None, help="Optional path to save marker genes CSV")
    parser.add_argument("--plot_prefix", default=None, help="Saves UMAP plots with this prefix")

    args = parser.parse_args()

    # 1. Load and Build Batch Field
    print(f"Loading {args.h5ad_input}...")
    adata = sc.read_h5ad(args.h5ad_input)
    
    batch_cols = args.batch_keys.split(',')
    adata.obs['batch_covariate'] = (
        adata.obs[batch_cols].astype(str).agg('_'.join, axis=1).astype('category')
    )
    
    # 2. Integrated Sampling Logic (Ref vs Query)
    if args.max_cells_per_batch or args.automax:
        print("Sampling data into Reference and Query sets...")
        if args.automax:
            max_count = int(adata.obs['batch_covariate'].value_counts().median().round())
        else:
            max_count = args.max_cells_per_batch

        # Sample for reference (handling cases where batch size < max_count)
        ref_obs = adata.obs.groupby('batch_covariate', observed=True).sample(
            n=max_count, random_state=1, replace=True
        ).drop_duplicates()
        
        adata_ref = adata[ref_obs.index].copy()
        adata_query = adata[~adata.obs_names.isin(ref_obs.index)].copy()
        print(f"Ref size: {adata_ref.n_obs}, Query size: {adata_query.n_obs}")
    else:
        print("Using full dataset as Reference (no Query projection).")
        adata_ref = adata.copy()
        adata_query = None

    # 3. Reference Preprocessing & Training
    sc.pp.filter_cells(adata_ref, min_genes=args.gene_filter)
    adata_ref.layers["counts"] = adata_ref.X.copy()
    
    sc.pp.highly_variable_genes(
        adata_ref, flavor="seurat_v3", n_top_genes=args.n_top_genes,
        batch_key='batch_covariate', subset=True, layer="counts", span=args.hvg_span
    )

    cont_cov = args.continuous_covariates.split(',') if args.continuous_covariates else None
    scvi.model.SCVI.setup_anndata(
        adata_ref, layer="counts", batch_key='batch_covariate', continuous_covariate_keys=cont_cov
    )
    arches_params = dict(
        use_layer_norm="both",
        use_batch_norm="none",
        encode_covariates=True,
        dropout_rate=0.2,
        n_layers=2,
        n_latent = args.n_latent,
    )

    vae_ref = scvi.model.SCVI(
        adata_ref,
        **arches_params
    )

    print("Training Reference scVI Model...")
    vae_ref.train(max_epochs=args.n_epochs, accelerator='gpu', early_stopping=True)
    vae_ref.save(args.model_output, overwrite=True, save_anndata=True)

    # 4. Optional Query Projection (scArches)
    if adata_query is not None and adata_query.n_obs > 0:
        print("Projecting Query data...")
        adata_query = adata_query[:, adata_ref.var_names].copy()
        adata_query.layers["counts"] = adata_query.X.copy()
        
        scvi.model.SCVI.setup_anndata(
            adata_query, layer="counts", batch_key='batch_covariate', continuous_covariate_keys=cont_cov
        )
        
        vae_q = scvi.model.SCVI.load_query_data(adata_query, vae_ref, freeze_dropout=True)
        vae_q.train(max_epochs=args.n_epochs, accelerator='gpu', plan_kwargs=dict(weight_decay=0.0))
        
        adata_full = ad.concat([adata_ref, adata_query])
        adata_full.obsm["X_scVI"] = vae_q.get_latent_representation(adata_full)
    else:
        adata_full = adata_ref
        adata_full.obsm["X_scVI"] = vae_ref.get_latent_representation()

    # 5. Optional scANVI
    if args.scanvi_label:
        print(f"Running scANVI on labels: {args.scanvi_label}")
        if 'unlabelled' not in adata_full.obs[args.scanvi_label].dtype.categories:
            adata_full.obs[args.scanvi_label] = adata_full.obs[args.scanvi_label].cat.add_categories(['unlabelled'])
        adata_full.obs[args.scanvi_label] = adata_full.obs[args.scanvi_label].fillna('unlabelled')
        
        lvae = scvi.model.SCANVI.from_scvi_model(vae_ref, unlabeled_category="unlabelled", labels_key=args.scanvi_label)
        lvae.train(max_epochs=min(20, args.n_epochs))
        
        adata_full.obsm["X_scANVI"] = lvae.get_latent_representation(adata_full)
        adata_full.obs["predictions_scANVI"] = lvae.predict(adata_full)
        
        if args.scanvi_model_output:
            lvae.save(args.scanvi_model_output, overwrite=True, save_anndata=True)

    # 6. Dimensionality Reduction (RAPIDS) & Exports
    print("Running Neighbors and UMAP...")
    rsc.pp.neighbors(adata_full, use_rep="X_scVI")
    rsc.tl.umap(adata_full, min_dist=0.3)
    
    # 7. Leiden clustering
    # Multiple Leiden resolutions
    for res in [0.5, 1, 2, 5]:
        key = f"leiden_res_{str(res).replace('.', '_')}"
        rsc.tl.leiden(adata_full, resolution=res, key_added=key)
    # 8. Marker detection
    if args.marker_csv_output:
        print("Calculating HVG Marker Genes for Res 1.0...")
        sc.pp.normalize_total(adata_full, target_sum=1e4)
        sc.pp.log1p(adata_full)
        sc.tl.rank_genes_groups(adata_full, groupby="leiden_res_1", method='wilcoxon')
        result = adata_full.uns['rank_genes_groups']
        groups = result['names'].dtype.names
        markers_df = pd.DataFrame(
            {group + '_' + key: result[key][group]
            for group in groups for key in ['names', 'logfoldchanges', 'pvals_adj']}).head(50)
        markers_df.to_csv(args.marker_csv_output, index=False)

    if args.obs_csv_output: 
        print(f"Exporting metadata to {args.obs_csv_output}...")
        obs_df = adata_full.obs.copy()
        for key in adata_full.obsm_keys():
            df_obsm = pd.DataFrame(
                adata_full.obsm[key], index=adata_full.obs_names,
                columns=[f"{key}_{i+1}" for i in range(adata_full.obsm[key].shape[1])]
            )
            obs_df = pd.concat([obs_df, df_obsm], axis=1)
        obs_df.to_csv(args.obs_csv_output)

    if args.plot_prefix:
        sc.pl.umap(adata_full, color=['batch_covariate'], show=False, save=f"_{args.plot_prefix}.png")

    print(f"Saving final object to {args.h5ad_output}...")
    adata_full.write_h5ad(args.h5ad_output)

if __name__ == "__main__":
    main()
