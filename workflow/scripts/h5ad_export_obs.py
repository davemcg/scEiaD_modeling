import sys
import scanpy as sc
import pandas as pd

print(sys.argv)
adata = sc.read_h5ad(sys.argv[1])
obs = pd.DataFrame(adata.obs)
obs.to_csv(sys.argv[2])
