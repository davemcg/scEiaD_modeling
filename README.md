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

