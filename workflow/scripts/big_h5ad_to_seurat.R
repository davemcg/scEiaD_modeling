library(dplyr)
library(BPCells)
library(Seurat) #version >= 5
args = commandArgs(trailingOnly=TRUE)

adata <- open_matrix_anndata_hdf5(args[1])

bp_out <- gsub(args[1], ".h5ad", "_BPdir")
print(paste('making temp dir', bp_out))
write_matrix_dir(mat = adata, dir = bp_out, overwrite =TRUE)

mat <- open_matrix_dir(dir = bp_out)

metadata <- read.csv(args[2])

seurat <- CreateSeuratObject(counts = list(mat), meta.data = metadata)

save(seurat, file = args[3])
