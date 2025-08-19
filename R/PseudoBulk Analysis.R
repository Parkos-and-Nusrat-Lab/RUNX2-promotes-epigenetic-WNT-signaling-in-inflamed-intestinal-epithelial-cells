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
library(ReactomePA)
library(msigdbr)
library(enrichplot)
library(DESeq2)
library(EnhancedVolcano)


seu <- readRDS("Seurat_file.rds")

old <- seu$orig.ident


status <- c()

for (i in old){
  if (i == "Healthy1_2D" | i == "Healthy2_2D"){
    status <- c(status, "Healthy_2D")
  } else {
    if ( i == "Healthy1_3D" | i == "Healthy2_3D"){
      status <- c(status, "Healthy_3D")
    } else {
      if (i == "UC1_2D" ) {
        status <- c(status, "activeUC_2D")
      } else { 
        if ( i == "UC2_2D"){
          status <-c(status, "inactiveUC_2D")
        } else{
          if(i == "UC1_3D"){
            status <- c(status, "activeUC_3D")
          } else { status <- c(status, "inactiveUC_3D")}
        }
      }
    }
  }
}

seu[["status"]] <- status

#Preparing the whole file to find the values to make the subsetting
Idents(seu) <- "orig.ident"



#Pseudobulking
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

vsd <- vst(dds, blind=FALSE)
DESeq2::plotPCA(vsd, "Combined")

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




source("GSEA_function.R")


healthy_df <- as.data.frame(res)
uc_df <- as.data.frame(res_uc)

healthy_results <- GSEA_function(healthy_df, 
                                   comparing = "Healthy 2D vs Healthy 3D", 
                                   return_object_list = T)
uc_results <- GSEA_function(uc_df, 
                              comparing = "UC 2D vs UC 3D", 
                              return_object_list = T)


p <- healthy_results@result %>% 
  dplyr::arrange(-NES) %>% 
  enrichplot::dotplot(showCategory = top, x = "NES") +
  ggplot2::ggtitle(paste(comparison, "Upregulated", sep = " ")) + 
  ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5))

ggsave("Supplementary_Figure_3C.png",plot = p, 
       device = "png",
       width = 12, 
       height = 14, 
       units = "in",
       dpi = 600)

p <- uc_results@result %>% 
  dplyr::arrange(-NES) %>% 
  enrichplot::dotplot(showCategory = top, x = "NES") +
  ggplot2::ggtitle(paste(comparison, "Upregulated", sep = " ")) + 
  ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5))

ggsave("Supplementary_Figure_3D.png",plot = p, 
       device = "png",
       width = 12, 
       height = 14, 
       units = "in",
       dpi = 600)


