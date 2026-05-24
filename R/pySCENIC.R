library(Seurat)
library(tidyverse)
library(SCENIC)
library(SCopeLoomR)

# Preparing the Seurat object to be used for pySCENIC as a loom file
seu <- readRDS("Seurat_file.rds")

exprMat <- seu@assays$RNA@data
cellInfo <- seu@meta.data

loci1 <- which(rowSums(exprMat) > 1*.01*ncol(exprMat))

dim(exprMat)

exprMat_filter <- exprMat[loci1, ]

dim(exprMat_filter)

add_cell_annotation <- function(loom, cellAnnotation)
{
  cellAnnotation <- data.frame(cellAnnotation)
  if(any(c("nGene", "nUMI") %in% colnames(cellAnnotation)))
  {
    warning("Columns 'nGene' and 'nUMI' will not be added as annotations to the loom file.")
    cellAnnotation <- cellAnnotation[,colnames(cellAnnotation) != "nGene", drop=FALSE]
    cellAnnotation <- cellAnnotation[,colnames(cellAnnotation) != "nUMI", drop=FALSE]
  }
  
  if(ncol(cellAnnotation)<=0) stop("The cell annotation contains no columns")
  if(!all(get_cell_ids(loom) %in% rownames(cellAnnotation))) stop("Cell IDs are missing in the annotation")
  
  cellAnnotation <- cellAnnotation[get_cell_ids(loom),,drop=FALSE]
  # Add annotation
  for(cn in colnames(cellAnnotation))
  {
    add_col_attr(loom=loom, key=cn, value=cellAnnotation[,cn])
  }
  
  invisible(loom)
}

loom <- build_loom("pyscenic.loom", dgem=exprMat_filter)
loom <- add_cell_annotation(loom, cellInfo)
close_loom(loom)


# Running pySCNEIC on the terminal
# Terminals Commands:
# pyscenic grn pyscenic.loom /
#  allTFs_hg38.txt /
#  -o adj.csv 
#  --num_workers 10
#
# pyscenic ctx adj.csv *feather / 
#  --annotations_fname 	motifs-v10nr_clust-nr.hgnc-m0.001-o0.0.tbl /
#  --expression_mtx_fname pyscenic.loom /
#  --output reg.csv /
#  --mask_dropouts /
#  --num_workers 10

# pyscenic aucell pyscenic.loom /
#  reg.csv /
#  --output pyscenic_results.loom /
#  --num_workers 10

# Feathers used
# hg38_10kbp_up_10kbp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather
# hg38_500bp_up_100bp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather
#

# After finishing pyscenic, run the Binary_AUC.py 


seu$sample_type <- str_extract(seu$orig.ident, "[23]D$")

seu$condition <- str_remove(seu$orig.ident, "_[23]D$")

seu$disease <- str_remove(seu$condition, "[12]$")

Idents(seu) <- "condition"



#After doing the previous steps you can restart from here
scenic <- readRDS('pyscenic_AUC.rds')
regulons <- scenic$regulons
regulonAUC <-  scenic$regulonAUC
AUCmat <- AUCell::getAUC(regulonAUC)
rownames(AUCmat) <- gsub("[(+)]", "", rownames(AUCmat))



## binary auc scores
Binarymat <- read.csv('Binary_AUC.csv',  sep = ',', row.names = 1)
Binarymat <- t(Binarymat)
rownames(Binarymat) <- gsub("[...]", "", rownames(Binarymat))

seu[['AUC']] <- CreateAssayObject(data = AUCmat)
seu[['AUCBinary']] <- CreateAssayObject(data = Binarymat)

DefaultAssay(seu) <- 'AUC'
DefaultAssay(seu) <- 'AUCBinary'



p <- DotPlot(seu, 
             features = c("TCF7L2", "TCF12", "HES1", "TCF4"), 
             assay = "AUC" ,
             group.by = "disease", 
             scale = T, 
             dot.scale = 15, 
             dot.min = 0.05, 
             scale.min = 0)  +
  ggtitle("Regulon Enrichment") +
  theme(axis.text.x = element_text(angle = 90,
                                   vjust = 0.5, 
                                   hjust=0.5, 
                                   size = 12),
        plot.title = element_text(vjust = 0.5, 
                                  size = 15,
                                  hjust = 0.5)
        ) +
  geom_point(aes(size=pct.exp), 
             shape = 21, 
             colour="black", stroke=0.5) +
  xlab("") +
  ylab("")+
  scale_color_viridis_c(option="magma") +
  guides(size=guide_legend(override.aes=list(shape=21, 
                                             colour="black", 
                                             fill="white"))) 

p


cowplot::ggsave2(filename = "Suppl_Fig_10.png", 
                 dpi = 300, plot=p, 
                 device = "png", bg = "white",
                 height = 5.77, width = 5.5)




