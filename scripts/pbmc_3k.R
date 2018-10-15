library(clustree)

# using seurat
library("Seurat")
# data set link - https://s3-us-west-2.amazonaws.com/10x.files/samples/cell/pbmc3k/pbmc3k_filtered_gene_bc_matrices.tar.gz
pbmc.data <- Read10X(data.dir = "data/filtered_gene_bc_matrices/hg19/")

pbmc <- CreateSeuratObject(raw.data = pbmc.data, min.cells = 3, min.genes = 200,
                           project = "10X_PBMC")

mito.genes <- grep(pattern = "^MT-", x = rownames(x = pbmc@data), value = TRUE)
percent.mito <- Matrix::colSums(pbmc@raw.data[mito.genes, ])/Matrix::colSums(pbmc@raw.data)

pbmc <- AddMetaData(object = pbmc, metadata = percent.mito, col.name = "percent.mito")

pbmc <- FilterCells(object = pbmc, subset.names = c("nGene", "percent.mito"),
                    low.thresholds = c(200, -Inf), high.thresholds = c(2500, 0.05))

pbmc <- NormalizeData(object = pbmc, normalization.method = "LogNormalize",
                      scale.factor = 10000)

pbmc <- FindVariableGenes(object = pbmc, mean.function = ExpMean, dispersion.function = LogVMR,
                          x.low.cutoff = 0.0125, x.high.cutoff = 3, y.cutoff = 0.5)

pbmc <- ScaleData(object = pbmc, vars.to.regress = c("nUMI", "percent.mito"))

pbmc <- RunPCA(object = pbmc, pc.genes = pbmc@var.genes, do.print = TRUE, pcs.print = 1:5,
               genes.print = 5)

pbmc <- FindClusters(object = pbmc, reduction.type = "pca", dims.use = 1:10,
                     resolution = c(0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0), print.output = 0, save.SNN = TRUE)

clustree(pbmc)
