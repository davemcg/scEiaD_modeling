import scanpy as sc
import numpy as np
import pandas as pd

import sys

adata = sc.read_h5ad(sys.argv[1])
samples = pd.read_csv(sys.argv[2])

samples_logical = np.array([s in samples.values for s in adata.obs.sample_accession])

adata_sub = adata[samples_logical,].copy()

adata_sub.write_h5ad(sys.argv[3])
