# scEiaD: modeling

# Purpose
Take cleaned and auto QC'ed counts from [scEiaD quant](http://github.com/davemcg/scEiaD_quant) and build tissue / compartment level [scVI](http://scvi-tools.org) models and call cell types.

# Step 1
Run workflow/scripts/merge_adata.py with metadata to build a custom h5ad object. 

There are three key inputs:
  1. The h5ad files (duh, they are the count s). 
  2. Sample metadata (see below for example)
  3. Cell metadata (also see below)

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

### Peek into the files
/home/mcgaugheyd/git/scEiaD_quant/sample_meta.scEiaD_v1.2024_02_28.01.tsv.gz 
(this can be the same as the `srr_sample_file` file used in [scEiaD_quant](github.com/davemcg/scEiad_quant)
```
sample_accession	run_accession	library_layout	reference	kb_tech	umi	workflow	kb_sum	organism	platform	study_accession	tissue	sub_tissue	covariate	perturbation	integration_group	tissuenote	source	bam10x	comment	biosample	organ	sex	biosample_title	strain	batch	age	capture_typeenriched_cell_type	suspension_enrichment_factors	ethnicity
SRX14524742	SRR18390614	PAIRED	hs111	10xv3	TRUE	nac	total	Homo sapiens	10xv3	SRP364915	Ciliary body		Pt2	NA		Right Eye	Tissue	https://sra-pub-src-2.s3.amazonaws.com/SRR18390614/Pt2CB.bam.1	NA	SAMN26813876	Eye		Pt2CB	NA	SRP364915_10xv3_Pt2	Adult	nucleus
SRX14524741	SRR18390615	PAIRED	hs111	10xv3	TRUE	nac	total	Homo sapiens	10xv3	SRP364915	Lens		Pt14	NA		Left Eye	Tissue	https://sra-pub-src-2.s3.amazonaws.com/SRR18390615/Pt14Lens.bam.1	NA	SAMN26813877	Eye		Pt14Lens	NA	SRP364915_10xv3_Pt14	Adult	nucleus
SRX14524740	SRR18390616	PAIRED	hs111	10xv3	TRUE	nac	total	Homo sapiens	10xv3	SRP364915	Cornea		Pt14	NA		Left Eye	Tissue	https://sra-pub-src-1.s3.amazonaws.com/SRR18390616/Pt14Cornea.bam.1	NA	SAMN26813878	Eye		Pt14Cornea	NA	SRP364915_10xv3_Pt14	Adult	nucleus
SRX14524739	SRR18390617	PAIRED	hs111	10xv3	TRUE	nac	total	Homo sapiens	10xv3	SRP364915	Trabecular meshwork		Hu0235	NA		Left Eye	Tissue	https://sra-pub-src-2.s3.amazonaws.com/SRR18390617/Hu235OSTM.bam.1	NA	SAMN26813879	Eye		Hu235OSTM	NA	SRP364915_10xv3_Hu0235	Adult   	nucleus
SRX14524738	SRR18390618	PAIRED	hs111	10xv3	TRUE	nac	total	Homo sapiens	10xv3	SRP364915	Cornea sclera wedge		Hu0235	NA		Left Eye	Tissue	https://sra-pub-src-2.s3.amazonaws.com/SRR18390618/Hu235OSCSW.bam.1	NA	SAMN26813880	Eye		Hu235OSCSW	NA	SRP364915_10xv3_Hu0235	Adult   	nucleus
SRX14524737	SRR18390619	PAIRED	hs111	10xv3	TRUE	nac	total	Homo sapiens	10xv3	SRP364915	Iris		Hu0235	NA		Left Eye	Tissue	https://sra-pub-src-2.s3.amazonaws.com/SRR18390619/Hu235IrisS2.bam.1	NA	SAMN26813881	Eye		Hu235IrisS2	NA	SRP364915_10xv3_Hu0235	Adult	nucleus
SRX14524736	SRR18390620	PAIRED	hs111	10xv3	TRUE	nac	total	Homo sapiens	10xv3	SRP364915	Iris		Hu0235	NA		Left Eye	Tissue	https://sra-pub-src-2.s3.amazonaws.com/SRR18390620/Hu235IrisS1.bam.1	NA	SAMN26813882	Eye		Hu235IrisS1	NA	SRP364915_10xv3_Hu0235	Adult	nucleus
SRX14524735	SRR18390621	PAIRED	hs111	10xv3	TRUE	nac	total	Homo sapiens	10xv3	SRP364915	Cornea		Hu0235	NA		Left Eye	Tissue	https://sra-pub-src-2.s3.amazonaws.com/SRR18390621/Hu235Cornea.bam.1	NA	SAMN26813883	Eye		Hu235Cornea	NA	SRP364915_10xv3_Hu0235	Adult	nucleus
SRX14524734	SRR18390622	PAIRED	hs111	10xv3	TRUE	nac	total	Homo sapiens	10xv3	SRP364915	Ciliary body		Hu0235	NA		Left Eye	Tissue	https://sra-pub-src-1.s3.amazonaws.com/SRR18390622/Hu235CB.bam.1	NA	SAMN26813884	Eye		Hu235CB	NA	SRP364915_10xv3_Hu0235	Adult	nucleus
```

/home/mcgaugheyd/git/scEiaD_quant/scEiaD_cell_labels_2024_08_26.csv.gz
```
barcode,MajorCellType,CellType,SubCellType,cell_type_ontology_term_id
AAACCCAAGGGATGTC_SRX19501854,amacrine,glycinergic amacrine cell,HAC14,CL:4030028
AAACCCACACGCACCA_SRX19501854,rod,retinal rod cell,Rod,CL:0000604
AAACCCACAGTAGTTC_SRX19501854,amacrine,glycinergic amacrine cell,AII_1,CL:4030028
AAACCCAGTGCCCTTT_SRX19501854,amacrine,glycinergic amacrine cell,VG3,CL:4030028
AAACCCATCTCCACTG_SRX19501854,rod,retinal rod cell,Rod,CL:0000604
AAACCCATCTTGGATG_SRX19501854,amacrine,starburst amacrine cell,ON-SAC,CL:0004232
AAACGAAAGGGTGGGA_SRX19501854,rod,retinal rod cell,Rod,CL:0000604
AAACGAAAGTACCGGA_SRX19501854,mueller,Mueller cell,MG,CL:0000636
AAACGAACAAGTGGGT_SRX19501854,amacrine,GABAergic amacrine cell,HAC5,CL:4030027
```


# Step 2

## Option A

Make a new sc(an)VI model for the data. 

Workflow:

  - Update yaml file for `scEiaD_modeling/workflow/Snakefile`. Example: 
```
      # cat ~/git/scEiaD_modeling/config/config_hs111_dev_eye_full.yaml
      git_dir: '/home/mcgaugheyd/git/scEiaD_modeling/'
      prefix: 'hs111_dev_eye_20250204_'
      ref_bcs: '/home/mcgaugheyd/git/scEiaD_modeling/data/hs111_dev_eye_ref_bcs.full.20250204.csv.gz'
      query_bcs: '/home/mcgaugheyd/git/scEiaD_modeling/data/hs111_dev_eye_query_bcs.full.20250204.csv.gz'
      epochs: ['50','200']
      latent: ['20','30','50']
      hvg: ['1000','2000']
      input_h5ad: '../../hs111.adata.solo.20250204.h5ad'
      scanvi_on: "MajorCellType"
      scanvi_out: "scANVI_MCT"
      pb_group: ['leiden','leiden2','leiden3' ]
      scib: True
      continuous_covariate_keys: 'pct_counts_mt,pct_counts_ribo,pct_counts_protocadherin'
      scanvi_n_epochs: "150"
      conda: 'rscvi'
```
  - *Need* an `adata` `obs` column with cell type calls (above is "MajorCellType")
  - Run sc(an)VI Snakemake pipeline. Example `/data/OGVFB_BG/scEiaD/2024_02_28/snakeout/hs111_developing_eye/snakell.sh`:
```
      #!/bin/bash
      source /data/$USER/conda/etc/profile.d/conda.sh && source /data/$USER/conda/etc/profile.d/mamba.sh
      mamba deactivate; mamba activate;
      bash ~/git/scEiaD_modeling/Snakemake.wrapper.sh ~/git/scEiaD_modeling/workflow/Snakefile ~/git/scEiaD_modeling/config/config_hs111_dev_eye_full.yaml ~/git/scEiaD_modeling/config/cluster.json
```

## Option B

Project the data onto an existing model (created from Option A)
```
# biowulf2
cd /data/OGVFB_BG/scEiaD/2024_02_28/snakeout/hs111_developing_eye
source /data/$USER/conda/etc/profile.d/conda.sh && source /data/$USER/conda/etc/profile.d/mamba.sh
mamba activate rapids_singlecell
python ~/git/scEiaD_modeling/workflow/scripts/ct_projection.py  hs111.adata.solo.20250204.h5ad models_human.tsv ct_predictions__hs111.adata.solo.20250204.csv.gz
# note
# cat models_human.tsv
sanes_complete_ocular_atlas_SCP2310	/data/OGVFB_BG/scEiaD/2024_02_28/sanes_complete_ocular_ref/scanviModel.SCP2310_broad__2000hvg_200e_30l
chen_hrca	/data/OGVFB_BG/scEiaD/2024_02_28/chen_rca/scanviModel.hrca_all_20241215__2000hvg_200e_30l
chen_fetal_hrca	/data/OGVFB_BG/scEiaD/2024_02_28/chen_rca/scanviModel.fetal_hrca_88444d73-7f55-4a62-bcfe-e929878c6c78_20250206__2000hvg_200e_30l
sceiad_20250107_full_mmGeneFilter	/data/OGVFB_BG/scEiaD/2024_02_28/snakeout/hs111_mature_eye_full/full_s6_clean_model__mmGeneFilter/scanviModel.hs111_mature_eye_20250107_stage6_mmGeneFilter__2000hvg_200e_30l
sceiad_20250107_full	/data/OGVFB_BG/scEiaD/2024_02_28/snakeout/hs111_mature_eye_full/full_s6_clean_model/scanviModel.hs111_mature_eye_20250107_stage6__2000hvg_200e_30l
sceiad_20250107_nonneural	/data/OGVFB_BG/scEiaD/2024_02_28/snakeout/hs111_mature_eye_full/NONneural_cells_s6__cleanmodel/scanviModel.hs111_mature_eye_20250107_full__NONneural_s6___2000hvg_200e_30l
sceiad_20250107_neural	/data/OGVFB_BG/scEiaD/2024_02_28/snakeout/hs111_mature_eye_full/neural_cells_s6_clean_model/scanviModel.hs111_mature_eye_20250107_full__neural_s6___2000hvg_200e_30l
```


#  OGVFB Created scanVI Models

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
