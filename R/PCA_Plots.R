
library(Seurat)
library(tidyverse)

seu <- read_rds("Seurat_file.rds")

seu$sample_type <- str_extract(seu$orig.ident, "[23]D$")

seu$condition <- str_remove(seu$orig.ident, "_[23]D$")

Idents(seu) <- "sample_type"

# Making 3D PCA Plot
seu_3D <- subset(seu, idents = "3D")
DefaultAssay(seu_3D) <- "RNA"
seu_3D <- DietSeurat(seu_3D, layers = "counts", assays = "RNA")
seu_3D <- NormalizeData(seu_3D, normalization.method = "LogNormalize", scale.factor = 10000)
seu_3D <- FindVariableFeatures(seu_3D, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(seu_3D)
seu_3D <- ScaleData(seu_3D, features = all.genes, vars.to.regress = "percent.mt")
seu_3D <- RunPCA(seu_3D, features = VariableFeatures(object = seu_3D))
p1 <- DimPlot(seu_3D, reduction = "pca", group.by = "condition")


# Making 2D PCA Plot
seu_2D <- subset(seu, idents = "2D")
DefaultAssay(seu_2D) <- "RNA"
seu_2D <- DietSeurat(seu_2D, layers = "counts", assays = "RNA")
seu_2D <- NormalizeData(seu_2D, normalization.method = "LogNormalize", scale.factor = 10000)
seu_2D <- FindVariableFeatures(seu_2D, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(seu_2D)
seu_2D <- ScaleData(seu_2D, features = all.genes, vars.to.regress = "percent.mt")
seu_2D <- RunPCA(seu_2D, features = VariableFeatures(object = seu_2D))
p2 <- DimPlot(seu_2D, reduction = "pca", group.by = "condition")



p_grid <- plot_grid(plotlist = list(p1,p2), nrow = 2, ncol = 1)


ggsave2(filename = "Figure 2B.png",
        plot = p_grid,
        device = "png",
        dpi = 600,
        bg = "white"
)
