# Files obtained from the SCP259 Dataset
# Website : https://singlecell.broadinstitute.org/single_cell/study/SCP259/intra-and-inter-cellular-rewiring-of-the-human-colon-during-ulcerative-colitis
# Paper :
# Christopher S. Smillie, Moshe Biton, Jose Ordovas-Montanes,
# Keri M. Sullivan, Grace Burgin, Daniel B. Graham, Rebecca H. Herbst, Noga Rogel,
# Michal Slyper, Julia Waldman, Malika Sud, Elizabeth Andrews, Gabriella Velonias,
# Adam L. Haber, Karthik Jagadeesh, Sanja Vickovic, Junmei Yao, Christine Stevens,
# Danielle Dionne, Lan T. Nguyen, Alexandra-Chloé Villani, Matan Hofree,
# Elizabeth A. Creasey, Hailiang Huang, Orit Rozenblatt-Rosen, John J. Garber,
# Hamed Khalili, A. Nicole Desch, Mark J. Daly, Ashwin N. Ananthakrishnan,
# Alex K. Shalek, Ramnik J. Xavier, Aviv Regev,
#
# Intra- and Inter-cellular Rewiring of the Human Colon during Ulcerative Colitis,
# Cell,
# Volume 178, Issue 3,
# 2019,
# Pages 714-730.e22,
# ISSN 0092-8674,
# https://doi.org/10.1016/j.cell.2019.06.029.

library(Seurat)
library(tidyverse)


setwd("C:/Users/ricsilva/OneDrive - Michigan Medicine/Desktop/GitHub Repositories/ATAC_Manuscript/01 UPDATED Figures")

epi <- Read10X("E:/IBD_Data/Epi", gene.column = 1)
fib <- Read10X("E:/IBD_Data/Fib", gene.column = 1)
imm <- Read10X("E:/IBD_Data/Imm", gene.column = 1)

meta <- read_delim("E:/IBD_Data/all_meta.txt", delim = "\t")

meta <- meta[-1,]


samples_list <- list(epi, fib, imm)

names(samples_list) <- c("Epithelium", "Lamina Propria", "Immune")



seu <- CreateSeuratObject(samples_list)

seu <- PercentageFeatureSet(seu, pattern = "^MT-", col.name = "percent.mt")


rm(epi, fib, imm, samples_list)
#Adding the Metadata
#Checking the order is correct

barcodes_seuob <- data.frame("NAME" = rownames(seu@meta.data))

barcodes_meta <- meta$NAME


together <- left_join(barcodes_seuob, meta, by = "NAME")

rownames(together) = together$NAME

together <- together[,c(-1, -3,-4)]

seu <- AddMetaData(seu, metadata = together)

rm(together, barcodes_seuob, barcodes_meta, meta)


new_annotations <-c()

old_anno <- seu$Cluster

for (i in old_anno){
  if (i %in% c("TA 1","TA 2")){
    new_annotations <- c(new_annotations, "TA")
  } else {
    if (i %in% c("Immature Enterocytes 1", "Immature Enterocytes 2")){
      new_annotations <- c(new_annotations, "Immature Enterocytes")
    } else {
      if (i %in% c("WNT2B+ Fos-lo 1","WNT2B+ Fos-lo 2")){
        new_annotations <- c(new_annotations, "WNT2B+ Fos-lo")
      } else{
        if (i %in% c("DC1", "DC2")){
          new_annotations <- c(new_annotations, "DC")
        } else{
          if(i %in% c("WNT5B+ 1","WNT5B+ 2")){
            new_annotations <- c(new_annotations, "WNT5B+")
          } else {new_annotations <- c(new_annotations, i) }
        }
      }
    }
  }
}

seu <- AddMetaData(seu, new_annotations, col.name = "new_anno")

seu <- JoinLayers(seu)

Idents(seu) <- "new_anno"



celltypes <- c("Best4+ Enterocytes",
               "Stem",
               "Enterocyte Progenitors",
               "Immature Enterocytes",
               "Enterocytes",
               "Immature Goblet",
               "Goblet",
               "Enteroendocrine",
               "M cells",
               "Cycling TA",
               "Secretory TA",
               "TA")
subseu <- subset(seu, idents = celltypes)


subseu <-  NormalizeData(subseu)

# p <- DotPlot(subseu, c("APC","CTNNB1","OLFM4"),
#              split.by = "new_anno",
#              group.by = "Health",
#              dot.scale = 10,
#              cols = "Reds")  +
#   theme(axis.text.x = element_text(angle = 90,vjust = 0.5, hjust=0.5, size = 12)) +
#   geom_point(aes(size=pct.exp), shape = 21, colour="black", stroke=0.5) +
#   xlab("") + ylab("")+
#   scale_color_viridis_c(option="magma") +
#   guides(size=guide_legend(override.aes=list(shape=21, colour="black", fill="white"))) +
#   geom_hline(yintercept = c(12.5,24.5), linetype = "dotted", linewidth = 1.1)
#
#
# p$data$id <- factor(p$data$id,
#                     levels = c(paste("Healthy", celltypes, sep = "_"),
#                                paste("Inflamed", celltypes, sep = "_"),
#                                paste("Non-inflamed", celltypes, sep = "_")
#                     )
# )
#
# p
#
# ggsave(filename = "New_Supplementary_Figure8.png", plot = p, width = 8, height = 11, bg = "white")




p <- DotPlot(subseu, c("APC","CTNNB1","OLFM4","RUNX1"),
        split.by = "new_anno",
        group.by = "Health",
        dot.scale = 10,
        cols = "Reds")  +
  theme(axis.text.x = element_text(angle = 90,vjust = 0.5, hjust=0.5, size = 12)) +
  geom_point(aes(size=pct.exp), shape = 21, colour="black", stroke=0.5) +
  xlab("") + ylab("")+
  scale_color_viridis_c(option="magma") +
  guides(size=guide_legend(override.aes=list(shape=21, colour="black", fill="white"))) +
  geom_hline(yintercept = c(12.5,24.5), linetype = "dotted", linewidth = 1.1)


p$data$id <- factor(p$data$id,
                    levels = c(paste("Healthy", celltypes, sep = "_"),
                    paste("Inflamed", celltypes, sep = "_"),
                    paste("Non-inflamed", celltypes, sep = "_")
                    )
                    )

p

ggsave(filename = "Supplementary_Figure8.png", plot = p, width = 8, height = 11, bg = "white")




p <- DotPlot(subseu, c("RUNX2", "RUNX3"),
             split.by = "new_anno",
             group.by = "Health",
             dot.scale = 10,
             cols = "Reds")  +
  theme(axis.text.x = element_text(angle = 90,vjust = 0.5, hjust=0.5, size = 12)) +
  geom_point(aes(size=pct.exp), shape = 21, colour="black", stroke=0.5) +
  xlab("") + ylab("")+
  scale_color_viridis_c(option="magma") +
  guides(size=guide_legend(override.aes=list(shape=21, colour="black", fill="white"))) +
  geom_hline(yintercept = seq(from = 3.5, to = 33.5, by = 3), linetype = "dotted", linewidth = 1.1)


p$data$id <- factor(p$data$id,
                    levels = outer(c("Healthy", "Inflamed","Non-inflamed"), celltypes, paste,sep = "_")
                    )

p

ggsave(filename = "Supplementary_Figure9.png", plot = p, width = 7, height = 11, bg = "white")


#
#
# p <- DotPlot(subseu, c("RUNX1","RUNX2", "RUNX3"),
#              split.by = "new_anno",
#              group.by = "Health",
#              dot.scale = 10,
#              cols = "Reds")  +
#   theme(axis.text.x = element_text(angle = 90,vjust = 0.5, hjust=0.5, size = 12)) +
#   geom_point(aes(size=pct.exp), shape = 21, colour="black", stroke=0.5) +
#   xlab("") + ylab("")+
#   scale_color_viridis_c(option="magma") +
#   guides(size=guide_legend(override.aes=list(shape=21, colour="black", fill="white"))) +
#   geom_hline(yintercept = seq(from = 3.5, to = 33.5, by = 3), linetype = "dotted", linewidth = 1.1)
#
#
# p$data$id <- factor(p$data$id,
#                     levels = outer(c("Healthy", "Inflamed","Non-inflamed"), celltypes, paste,sep = "_")
# )
#
# plot_list <- list()
# for (i in c("RUNX1","RUNX2", "RUNX3")){
#
#   if (i == "RUNX1"){
#     p<- DotPlot(subseu, i,
#                 split.by = "new_anno",
#                 group.by = "Health",
#                 dot.scale = 10,
#                 cols = "Reds")  +
#       theme(axis.text.x = element_text(angle = 90,vjust = 0.5, hjust=0.5, size = 12)) +
#       geom_point(aes(size=pct.exp), shape = 21, colour="black", stroke=0.5) +
#       xlab("") + ylab("")+
#       scale_color_viridis_c(option="magma") +
#       guides(size=guide_legend(override.aes=list(shape=21, colour="black", fill="white"))) +
#       geom_hline(yintercept = seq(from = 3.5, to = 33.5, by = 3), linetype = "dotted", linewidth = 1.1)
#
#
#     p$data$id <- factor(p$data$id,
#                         levels = outer(c("Healthy", "Inflamed","Non-inflamed"), celltypes, paste,sep = "_")
#     )
#
#
#     plot_list[[i]] <- p
#
#   } else {
#
#     p<- DotPlot(subseu, i,
#                 split.by = "new_anno",
#                 group.by = "Health",
#                 dot.scale = 10,
#                 cols = "Reds")  +
#       theme(axis.text.x = element_text(angle = 90,vjust = 0.5, hjust=0.5, size = 12),
#             axis.text.y = element_blank()) +
#       geom_point(aes(size=pct.exp), shape = 21, colour="black", stroke=0.5) +
#       xlab("") + ylab("")+
#       scale_color_viridis_c(option="magma") +
#       guides(size=guide_legend(override.aes=list(shape=21, colour="black", fill="white"))) +
#       geom_hline(yintercept = seq(from = 3.5, to = 33.5, by = 3), linetype = "dotted", linewidth = 1.1)
#
#
#     p$data$id <- factor(p$data$id,
#                         levels = outer(c("Healthy", "Inflamed","Non-inflamed"), celltypes, paste,sep = "_")
#     )
#
#
#     plot_list[[i]] <- p
#   }
#
# }
#
# p <- cowplot::plot_grid(plotlist = plot_list, ncol = 3)
# ggsave(filename = "Alternative_ New_Supplementary_Figure9.png", plot = p, width = 30, height = 11, bg = "white")
