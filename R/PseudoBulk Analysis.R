#Pseudo bulk Healthy 2D vs 3D
# And Volcano Plots
library(Seurat)
library(tidyverse)
library(RColorBrewer)
library(scales)
library(cowplot)
library(patchwork)
library(grid)
library(gridExtra)
library(clusterProfiler)
library(org.Hs.eg.db)
library(msigdbr)
library(enrichplot)
library(DESeq2)
library(EnhancedVolcano)


source("GSEA_function.R")

# Figure 3F
seu <- readRDS("Seurat_file.rds")

seu$sample_type <- str_extract(seu$orig.ident, "[23]D$")

seu$condition <- str_remove(seu$orig.ident, "_[23]D$")


seu$HvD <- ifelse(grepl("^UC", seu$condition), "UC", "Healthy")


Idents(seu) <- "HvD"


deg <- FindMarkers(seu, ident.1 = "UC", ident.2 = "Healthy", assay = "RNA")

# Making Inf value to the max numeric abs value of the deg

max_value <- deg %>% 
  dplyr::filter(!is.infinite(avg_log2FC)) %>% 
  dplyr::select(avg_log2FC) %>% abs() %>% max()

deg$avg_log2FC <- ifelse(is.infinite(deg$avg_log2FC), max_value, deg$avg_log2FC)

deg <- dplyr::rename(deg, log2FoldChange = avg_log2FC)


gsea_msigdb <- GSEA_function(deg,
                             comparing = "UC vs Healthy",
                             category = "C5",
                             subcategory = "GO:BP"
)





wnt_names <- c("GOBP_NON_CANONICAL_WNT_SIGNALING_PATHWAY", 
               "GOBP_CANONICAL_WNT_SIGNALING_PATHWAY",
               "GOBP_POSITIVE_REGULATION_OF_WNT_SIGNALING_PATHWAY",
               "GOBP_POSITIVE_REGULATION_OF_CANONICAL_WNT_SIGNALING_PATHWAY")

wnt_names_titles <- wnt_names %>% 
  str_remove("GOBP_") %>% 
  str_replace_all("_"," ")


total_loops <- length(wnt_names)

plots <- list()
for(i in 1:total_loops){
  pathway <- wnt_names[i]
  title <- wnt_names_titles[i]
  
  plot_df <- gsea_msigdb@result %>%
    as.data.frame() %>%
    dplyr::filter(rownames(.) %in% pathway) %>%
    rownames_to_column("old") %>%
    dplyr::select("NES") %>%
    dplyr::mutate("string" = "NES",
                  "NES" = round(NES, 4)) %>% 
    column_to_rownames("string") %>% 
    unname()
  tbl <- tableGrob(plot_df, theme = ttheme_default(base_size = 19))
  
  
  p <- clusterProfiler::gseaplot(gsea_msigdb, 
                                 pathway, 
                                 by = "runningScore", 
                                 title = title)
  
  if (i != total_loops){
    p<-  p + theme(axis.title.x = element_blank())
  }
  p <- p + 
    theme(plot.title = element_text(size = 17, hjust = 0.5, vjust = 0.5)) +
    annotation_custom(tbl, xmin = 4200, xmax = 4500, ymin = 0.3, ymax = 0.4)
  
  plots[[length(plots)+1]] <- p
}


grid_plot <- cowplot::plot_grid(plotlist = plots, nrow = 4, ncol = 1)



ggsave2("Figure 3F.png", plot = grid_plot,
        dpi=300, device = "png", width = 10, height = 14, units = "in")





# Figure 4C

Idents(seu) <- "cell_anno"

subseu <- subset(seu, idents = c("LGR5+Stem_3D",
                                 "Early_G1_3D",
                                 "Late_G1_3D",
                                 "UC-specific_3D"))



subseu$cell_anno <- droplevels(subseu$cell_anno)


deg <- FindMarkers(subseu, ident.1 = "UC-specific_3D", group.by = "cell_anno", assay = "RNA")



deg2 <- deg

deg$trunc_log2FC <- ifelse(abs(deg$avg_log2FC) > 50, 50, deg$avg_log2FC)

p <- EnhancedVolcano(deg, 
                     lab = rownames(deg),
                     x = "trunc_log2FC", 
                     y="p_val", 
                     selectLab = c("CTNNB1", "APC"),
                     pCutoffCol = "p_val_adj",
                     FCcutoff = 5,
                     pCutoff = 0.05,
                     title = "UC-specific Cluster vs All",
                     col = c("grey30", "forestgreen", "#5f879e", "#d65f4d"),
                     subtitle = "Log2FC >50 truncated for visualization", 
                     legendPosition = "right") 




cowplot::ggsave2(filename = "Figure 4C.png", 
                 dpi = 300,
                 plot=p, 
                 device = "png",
                 bg ="white",
                 height = 8,
                 width = 16, 
                 units = "in")



# Figure 5D



deg <- FindMarkers(seu, ident.1 = "UC1", group.by = "condition", assay = "RNA")


deg$trunc_log2FC <- ifelse(abs(deg$avg_log2FC) > 50, 50, deg$avg_log2FC)

p <- EnhancedVolcano(deg, 
                     lab = rownames(deg),
                     x = "trunc_log2FC", 
                     y="p_val", 
                     selectLab = c("HES1", "ATOH1"),
                     pCutoffCol = "p_val",
                     FCcutoff = 5,
                     pCutoff = 0.05,
                     title = "Active UC Sample Cluster vs All",
                     col = c("grey30", "forestgreen", "#5f879e", "#d65f4d"),
                     subtitle = "Log2FC >50 truncated for visualization",
                     labSize = 7)



cowplot::ggsave2(filename = "Figure 5D.png", 
                 dpi = 300,
                 plot=p, 
                 device = "png",
                 bg ="white",
                 height = 13,
                 width = 8, 
                 units = "in")


# Figure 5 E



Idents(seu) <- "HvD"


deg <- FindMarkers(seu, ident.1 = "UC", ident.2 = "Healthy", assay = "RNA")


# Making Inf value to the max numeric abs value of the deg

max_value <- deg %>% 
  dplyr::filter(!is.infinite(avg_log2FC)) %>% 
  dplyr::select(avg_log2FC) %>% abs() %>% max()

deg$avg_log2FC <- ifelse(is.infinite(deg$avg_log2FC), max_value, deg$avg_log2FC)

deg <- dplyr::rename(deg, log2FoldChange = avg_log2FC)


gsea_msigdb <- GSEA_function(deg,
                             comparing = "UC vs Healthy",
                             category = "C5",
                             subcategory = "GO:BP"
)





notch_names <- c("GOBP_NOTCH_SIGNALING_PATHWAY", 
                 "GOBP_REGULATION_OF_NOTCH_SIGNALING_PATHWAY",
                 "GOBP_POSITIVE_REGULATION_OF_NOTCH_SIGNALING_PATHWAY")

notch_names_titles <- notch_names %>% 
  str_remove("GOBP_") %>% 
  str_replace_all("_"," ")

total_loops <- length(notch_names)

plots <- list()
for(i in 1:total_loops){
  pathway <- notch_names[i]
  title <- notch_names_titles[i]
  
  plot_df <- gsea_msigdb@result %>%
    as.data.frame() %>%
    dplyr::filter(rownames(.) %in% pathway) %>%
    rownames_to_column("old") %>%
    dplyr::select("NES") %>%
    dplyr::mutate("string" = "NES",
                  "NES" = round(NES, 4)) %>% 
    column_to_rownames("string") %>% 
    unname()
  tbl <- tableGrob(plot_df, theme = ttheme_default(base_size = 22))
  
  
  p <- clusterProfiler::gseaplot(gsea_msigdb, 
                                 pathway, 
                                 by = "runningScore", 
                                 title = title) +
    annotation_custom(tbl, xmax = 4500, xmin = 4200, ymin = 0.35, ymax = 0.45)
  
  if (i != total_loops){
    p<-  p + 
      theme(axis.title.x = element_blank(),
            plot.title = element_text(size = 19, hjust = 0.5, vjust = 0.5)) 
  } else{
    p <- p + 
      theme(plot.title = element_text(size = 17, hjust = 0.5, vjust = 0.5)) 
  }
  
  
  
  plots[[length(plots)+1]] <- p
}


grid_plot <- cowplot::plot_grid(plotlist = plots, nrow = 3, ncol = 1)


ggsave2("Figure 5E.png", 
        plot = grid_plot,
        dpi=300, 
        device = "png", 
        width = 14,
        height = 14, 
        units = "in")






# Supplementary Figure 3

#Pseudobulking

seu <- readRDS("Seurat_file.rds")


Idents(seu) <- "orig.ident"


pseudo <- AggregateExpression(object = seu, assays = "RNA")$RNA

pseudo_meta <- data.frame(Dimension = c("2D","3D","2D","3D","2D","3D","2D","3D"),
                          Sample = c("1","1","2","2","3","3","4","4"), 
                          Combined1 = c("Healthy_2D","Healthy_3D",
                                        "Healthy_2D","Healthy_3D",
                                        "UC_2D","UC_3D",
                                        "UC_2D","UC_3D"),
                          row.names = c("Healthy1-2D",
                                        "Healthy1-3D",
                                        "Healthy2-2D",
                                        "Healthy2-3D",
                                        "UC1-2D",
                                        "UC1-3D",
                                        "UC2-2D",
                                        "UC2-3D"),
                          Combined = c("Healthy_2D","Healthy_3D",
                                       "Healthy_2D","Healthy_3D",
                                       "activeUC_2D","activeUC_3D",
                                       "inactiveUC_2D","inactiveUC_3D"))
all(rownames(pseudo_meta) == colnames(pseudo))

# Performing DESeq2 and EnhancedVolcano Plots



dds <- DESeqDataSetFromMatrix(countData = pseudo,
                              colData = pseudo_meta,
                              design = ~ Combined)


keep <- rowSums(counts(dds) >= 10) >= 3
dds <- dds[keep,]
dds <- DESeq(dds)
res <- results(dds, contrast=c("Combined","Healthy_2D","Healthy_3D"))
res

p<-EnhancedVolcano(res, 
                   lab = rownames(res), 
                   x = "log2FoldChange", 
                   y = "pvalue",
                   title = "Healthy 2D Colonoids vs Healthy 3D Colonoids",
                   subtitle = "",
                   pointSize = 1.5,
                   caption = ""
)


ggsave("Supplementary_Figure_3A.png",plot = p, 
       device = "png",
       width = 12, 
       height = 13, 
       units = "in",
       dpi = 600)

res_uc <- results(dds, contrast=c("Combined","activeUC_2D","activeUC_3D"))
res_uc

p<-EnhancedVolcano(res_uc, 
                   lab = rownames(res_uc), 
                   x = "log2FoldChange", 
                   y = "pvalue",
                   title = "UC 2D Colonoids vs UC 3D Colonoids",
                   subtitle = "",
                   pointSize = 1.5,
                   caption = ""
)

ggsave("Supplementary_Figure_3B.png",plot = p, 
       device = "png",
       width = 12, 
       height = 14, 
       units = "in",
       dpi = 600)







healthy_df <- as.data.frame(res)
uc_df <- as.data.frame(res_uc)

healthy_results <- GSEA_function(healthy_df, 
                                   comparing = "Healthy 2D vs Healthy 3D"
                                   )



healthy_results@result <- healthy_results@result%>% dplyr::arrange(-NES)


p <- healthy_results %>% 
  enrichplot::dotplot(showCategory = 20, x = "NES") +
  ggplot2::ggtitle(paste("Healthy 2D vs Healthy 3D", "Upregulated", sep = " ")) + 
  ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5))


ggsave("Supplementary_Figure_3C.png",plot = p, 
       device = "png",
       width = 12, 
       height = 14, 
       units = "in",
       dpi = 600)




uc_results <- GSEA_function(uc_df, 
                            comparing = "UC 2D vs UC 3D"
)


uc_results@result <- uc_results@result%>% dplyr::arrange(-NES)

p <- uc_results %>%
  enrichplot::dotplot(showCategory = 20, x = "NES") +
  ggplot2::ggtitle(paste("UC 2D vs UC 3D", "Upregulated", sep = " ")) + 
  ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5))

ggsave("Supplementary_Figure_3D.png",plot = p, 
       device = "png",
       width = 12, 
       height = 14, 
       units = "in",
       dpi = 600)


