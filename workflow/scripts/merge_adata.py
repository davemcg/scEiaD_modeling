import scanpy as sc
import pandas as pd
import argparse

parser = argparse.ArgumentParser(
                    prog='Merge same species adata')

parser.add_argument('adata_files_txt')
parser.add_argument('sample_meta_tsv')
parser.add_argument('cell_label_csv')
parser.add_argument('output_h5ad_name')
parser.add_argument('output_obs_name')
args = parser.parse_args()


input_file = args.adata_files_txt
files = []

with open(input_file) as txt:
    for line in txt:
        files.append(line.strip())

adata_dict = {}

for i in files:
    accession = i.split('/')[-5]
    print(i)
    adata_dict[accession] = sc.read_h5ad(i)
    adata_dict[accession].obs['sample_accession'] = accession

adata = sc.concat(adata_dict, index_unique = '_')
print("big adata made")
# smeta = pd.read_csv('/home/mcgaugheyd/git/scEiaD_quant/sample_meta.scEiaD_v1.2024_02_07.tsv.gz', sep = '\t')
smeta = pd.read_csv(args.sample_meta_tsv, sep = '\t')

smeta = smeta.drop(columns = 'run_accession')
smeta = smeta.drop_duplicates()

adata.obs['barcode'] = adata.obs.index
adata.obs = adata.obs.merge(smeta, on = 'sample_accession')
adata.obs.index = adata.obs['barcode']

#cmeta = pd.read_csv('/home/mcgaugheyd/git/scEiaD_quant/cell_labels.csv.gz')
cmeta = pd.read_csv(args.cell_label_csv)

cmeta.index = cmeta['barcode']
adata.obs = adata.obs.join(cmeta.drop(columns = ['barcode']))

adata.var['id_name'] = adata_dict[accession].var.id_name

#adata.write_h5ad('hs111.adata.solo.h5ad')
for i in adata.obs.columns:
    if adata.obs[i].count() == 0:
        adata.obs = adata.obs.drop(i, axis = 1)

adata.write_h5ad(args.output_h5ad_name)
obs = pd.DataFrame(adata.obs)
obs.to_csv(args.output_obs_name)
# adata[~adata.obs.solo_doublet]
