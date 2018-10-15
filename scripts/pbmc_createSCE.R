## cellrangerRkit is available at 
## https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/rkit

library(cellrangerRkit)

pipestance_path <- "tenx_pbmc"
download_sample(sample_name="fresh_68k_pbmc_donor_a",
                sample_dir=pipestance_path,
                host="http://cf.10xgenomics.com/samples/cell-exp/1.1.0/")

library(DropletUtils)
## sparseMatrix
sce <- read10xCounts("tenx_pbmc/outs/filtered_gene_bc_matrices/hg19/",
                     col.names = TRUE)
sce <- sce[Matrix::rowSums(counts(sce))>0,]
saveRDS(sce, file = "sce_sparse.rds")

## dense matrix
tmp <- as.matrix(counts(sce))
sce_dense <- sce
counts(sce_dense) <- tmp
saveRDS(sce_dense, file = "sce_dense.rds")

## HDF5 matrix
library(HDF5Array)
sce_h5 <- saveHDF5SummarizedExperiment(sce_dense, dir="tenx_sce", replace=FALSE,
                             chunkdim=NULL, level=NULL, verbose=TRUE)
