library(Seurat)
library(dplyr)

# Load the PBMC dataset
# too download run the Pbmc_createSCE first
palmtree <- as.seurat(sce)

palmtree <- NormalizeData(object = palmtree, normalization.method = "LogNormalize",
                      scale.factor = 10000)

palmtree <- RunPCA(object = palmtree, pc.genes = palmtree@var.genes, do.print = TRUE, pcs.print = 1:5,
               genes.print = 5)

palmtree <- FindClusters(object = palmtree, reduction.type = "pca", dims.use = 1:10,
                     resolution = c(0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0), print.output = 0, save.SNN = TRUE)

clustree(palmtree)

