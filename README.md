# scEiaD: modeling

# Purpose
Take cleaned and auto QC'ed counts from [scEiaD quant](http://github.com/davemcg/scEiaD_quant) and build tissue / compartment level [scVI](http://scvi-tools.org) models and call cell types.

# Step 1
Run workflow/scripts/merge_adata.py with metadata to build a custom h5ad object

## Example 
```
# biowulf2
cd /data/OGVFB_BG/scEiaD/2024_02_28
## find output h5ad files
find . -name "QC.adata.solo.h5ad" > fin.h5ad.txt
## human files
grep hs111 fin.h5ad.txt > hs111.fin.txt
mamba activate rscvi
python ~/git/scEiaD_modeling/workflow/scripts/merge_adata.py hs111.fin.txt /home/mcgaugheyd/git/scEiaD_quant/sample_meta.scEiaD_v1.2024_02_28.01.tsv.gz /home/mcgaugheyd/git/scEiaD_quant/scEiaD_cell_labels_2024_08_26.csv.gz hs111.adata.solo.20240826.h5ad
```

# Step 2

=======
Take cleaned and auto QC'ed counts from [scEiaD quant](http://github.com/davemcg/scEiaD_quant) and build tissue / compartment / stage level [scVI](http://scvi-tools.org) models and call cell types.

#  Models

| Species | Set | Source | Params | Path* | Note |
| ---- | ---- | --- | -- | ----- | ----- |
| Human | Mature Eye | scEiaD 2025 | 2000hvg, 200e, 30L | /data/OGVFB_BG/scEiaD/2024_02_28/snakeout/hs111_mature_eye_full/full_s6_clean_model__mmGeneFilter/scanviModel.hs111_mature_eye_20250107_stage6_mmGeneFilter__2000hvg_200e_30l/ | only using genes name (human readable symbol) matched between mouse/human in the HVG selection
| Human | Mature Eye | scEiaD 2025 | 2000hvg, 200e, 30L | /data/OGVFB_BG/scEiaD/2024_02_28/snakeout/hs111_mature_eye_full/full_s6_clean_model/scanviModel.hs111_mature_eye_20250107_stage6__2000hvg_200e_30l/| 
| Human | Developing Eye | scEiaD 2025 | 2000hvg, 200e, 50L | /data/OGVFB_BG/scEiaD/2024_02_28/snakeout/hs111_developing_eye/stage4_mmGeneFilter/scanviModel.hs111_dev_eye_stage4_mmGeneFilter_20250225__2000hvg_200e_50l/ | only using genes name (human readable symbol) matched between mouse/human in the HVG selection
| Human | Developing Eye | scEiaD 2025 | 2000hvg, 200e, 50L | /data/OGVFB_BG/scEiaD/2024_02_28/snakeout/hs111_developing_eye/stage4/scanviModel.hs111_dev_eye_stage4_20250226__2000hvg_200e_50l/
| Mouse | Mature Eye | scEiaD 2025 | 2000hvg, 50e, 20L | /data/OGVFB_BG/scEiaD/2024_02_28/snakeout/mm111_mature_eye_full/stage4/scanviModel.mm111_mature_eye_20250120_stage4_noCov_150epo__2000hvg_50e_20l/
| Mouse | Developing Eye | scEiaD 2025 | 2000hvg, 200e, 30L | /data/OGVFB_BG/scEiaD/2024_02_28/snakeout/mm111_developing_eye/stage3/scanviModel.mm111_dev_eye_20250305_stage3__2000hvg_200e_30l/
| Mouse | Mature Eye | Mouse Reference Cell Atlas (MRCA, Chen Lab) | 2000hvg, 200e, 30L | /data/OGVFB_BG/scEiaD/2024_02_28/outside_models/chen_rca/scanviModel.mrca_all_20250301__2000hvg_200e_30l | 70a012ab-bd22-4013-a04e-9a1a275cd8b4
| Human | Developing Eye | Human Fetal Reference Cell Atlas (HRCA, Chen Lab) | 2000hvg, 200e, 30L | /data/OGVFB_BG/scEiaD/2024_02_28/outside_models/chen_rca/scanviModel.fetal_hrca_88444d73-7f55-4a62-bcfe-e929878c6c78_20250206__2000hvg_200e_30l/ | 88444d73-7f55-4a62-bcfe-e929878c6c78
| Human | Mature Eye | Human  Reference Cell Atlas (HRCA, Chen Lab) | 2000hvg, 200e, 30L | /data/OGVFB_BG/scEiaD/2024_02_28/outside_models/chen_rca/scanviModel.hrca_all_20241215__2000hvg_200e_30l/ | 2e910e62-7eaf-4c06-80cb-8918e3eea16e
| Mouse | Mature Eye | Complete Ocular Reference (Sanes Lab) | 2000hvg, 200e, 30L | /data/OGVFB_BG/scEiaD/2024_02_28/outside_models/sanes_complete_ocular_ref/scanviModel.SCP2310_broad__2000hvg_200e_30l/




* on biowulf2.nih.gov, not publicly accessible as of yet 
