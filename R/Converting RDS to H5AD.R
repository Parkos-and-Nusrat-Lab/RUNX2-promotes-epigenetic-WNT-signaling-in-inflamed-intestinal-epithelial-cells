# Converting Seurat file to h5ad to use in python
library(reticulate)
library(scater)
library(SeuratDisk)
library(sceasy)


seu <- readRDS("Seurat_file.rds")


sceasy::convertFormat(seu, from="seurat", to="anndata",
                      outFile='Seurat_Full.h5ad')

