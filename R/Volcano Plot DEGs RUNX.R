# Loading Librarires
library(Seurat)
library(tidyverse)


# Loading File

seu <- readRDS("Seurat_file.rds")


# Updating levels to ensure celltype placement in figures

seu@active.ident <- factor(seu@active.ident, levels = c("LGR5+Stem_3D",
                                                        "Early_G1_3D",
                                                        "Late_G1_3D",
                                                        "UC-specific_3D",
                                                        "Hes+Stem_2D",
                                                        "BMI+Stem_2D",
                                                        "Undifferentiated_2D",
                                                        "Intermediate_2D",
                                                        "Transitional_2D",
                                                        "Differentiated_2D",
                                                        "M-like_2D")
)




subseu <- subset(seu, idents = "UC-specific_3D")


deg <- FindMarkers(subseu,
            group.by = "orig.ident",
            ident.1 = "UC1_3D", #active UC
            ident.2 = "UC2_3D" #inactive UC
            )



library(EnhancedVolcano)


p<-EnhancedVolcano(toptable = deg,
                x = "avg_log2FC",
                y= "p_val",
                xlim = c(-5,5),
                ylim = c(-0.3, 50),
                lab = rownames(deg),
                title = 'Active UC vs Inactive UC (UC-specific cluster)',
                pCutoff = 0.005,
                FCcutoff = 1.5,
                pointSize = 3.0,
                labSize = 5.0,
                drawConnectors = TRUE,
                widthConnectors = 1,
                colConnectors = "black",
                boxedLabels = TRUE,
                selectLab = c('RUNX1','RUNX2', 'RUNX3'))


ggsave("Volcano_Plot_RUNX.png", plot = p, device = "png", dpi= 600, bg = "white",
       width = 8, height = 10)
