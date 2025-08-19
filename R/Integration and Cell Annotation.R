
# Opening Libraries
library(BiocManager)
library(Seurat)
library(Signac)
library(EnsDb.Hsapiens.v86)
library(tidyverse)
library(hdf5r)
library(future)
library(GenomicRanges)
library(SingleCellExperiment)
library(data.table)


# Loading H5 files

files <- c("string_of_paths_to_h5files")

# Reading H5 files

h5_files <- lapply(files, Read10X_h5)


# Extracting RNA and ATAC data

list_counts <- lapply(h5_files, function(x){
  rna <- x["Gene Expression"]
  atac <- x["Peaks"]

  counts <- list(rna, atac)
  names(counts) <- c("RNA", "ATAC")
  return(counts)
})

# Hardcoded labeling, requires to read the files in this specific order
names(list_counts) <- c("Healthy1_2D","UC1_2D","Healthy1_3D",
                        "UC1_3D","Healthy2_2D","UC2_2D","Healthy2_3D","UC2_3D")

# Obtaining annotations

annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
genome(annotations) <- "hg38"
seqlevelsStyle(annotations) <- 'UCSC'


# Creating Seurat Objects

seurat_list <- list()

atac_fragment_paths <- c("Path_to_atac_fragment.tsv_files") # must match the order of the h5 files


for (i in seq_along(list_counts)){

  # Getting RNA and ATAC information
  rna_counts <- list_counts[[i]][[1]][[1]]
  atac_counts <- list_counts[[i]][[2]][[1]]

  # Getting the name to keep track of file
  working_name <- names(list_counts)[i]

  # Creating Seurat Object
  seu <- CreateSeuratObject(rna_counts, project = working_name)

  # Obtaining mitochondria %
  seu[["percent.mt"]] <- PercentageFeatureSet(seu, pattern = "^MT-")

  # Obtaining GRanges information from the ATAC matrix
  grange.counts <- StringToGRanges(rownames(atac_counts), sep = c(":", "-"))
  grange.use <- seqnames(grange.counts) %in% standardChromosomes(grange.counts)

  # Filtering ATAC matrix to keep only the ones found in the GRanges
  atac_counts <- atac_counts[as.vector(grange.use), ]

  # Obtianing the fragment file
  frag.file <- atac_fragment_paths[i]

  # Creating the ATAC Assay
  chrom_assay <- CreateChromatinAssay(
    counts = atac_counts,
    sep = c(":", "-"),
    fragments = frag.file,
    min.cells = 10,
    annotation = annotations
  )
  # Saving the ATAC Assay into Seurat Object
  seu[["ATAC"]] <- chrom_assay

  # Saving the Seurat object before starting the new one
  seurat_list[[working_name]] <- seu

}


# Preprocessing

# Subsetting, Cell Cycle Score, SCTransform, and reductions

seurat_list <- lapply(seurat_list, function(x){

  # Subsetting Seurat with the following thresholds
  seu <- x
  seu <- subset(x = seu,
                subset = nCount_ATAC < 1e5 &
                  nCount_ATAC > 1e3 &
                  nCount_RNA < 50000 &
                  nCount_RNA > 1000 &
                  percent.mt < 30 &
                  nFeature_RNA > 200
  )

  # Performing SCTransform to obtain Cell Cycle Score

  seu <- SCTransform(seu,
                     assay = "RNA",
                     new.assay.name = 'SCT',
                     vars.to.regress = c('percent.mt', 'nFeature_RNA', 'nCount_RNA')
  )

  # Obtaining Cell Cycle Score

  seu <- CellCycleScoring(seu,
                          s.features = cc.genes.updated.2019$s.genes,
                          g2m.features = cc.genes.updated.2019$g2m.genes,
                          set.ident = TRUE,
                          assay = "SCT"
  )

  # Re-doing SCTransform to regress out cell cycle information
  # Performing PCA and UMAP

  seu <- SCTransform(seu,
                     assay = "RNA",
                     new.assay.name = 'SCT',
                     vars.to.regress = c('percent.mt', 'nFeature_RNA',
                                         'nCount_RNA', 'S.Score', 'G2M.Score')
  ) %>%
    RunPCA() %>%
    RunUMAP(dims = 1:50)

  return(seu)


})


# Processing the ATAC Assay for TSSEnrichment and TFIDF


seurat_list <- lapply(seurat_list, function(x){

  seu <- x

  # Ensuring the ATAC Assay is the Active Assay
  DefaultAssay(seu) <- "ATAC"

  # Running TFIDF and obtaining TSS Enrichment
  seu <- RunTFIDF(seu) %>%
    FindTopFeatures(min.cutoff = 'q0') %>% # Using q0 cutoff for Top Features
    RunSVD() %>%
    RunUMAP(reduction = 'lsi',
            dims = 2:50,
            reduction.name = "umap.atac",
            reduction.key = "atacUMAP_") %>% #
    NucleosomeSignal() %>%
    TSSEnrichment()

  return(seu)

})



# Making WNN Graphs

seurat_list <- lapply(seurat_list, function(x){
  seu <- x
  # Getting multimodal Neighbors
  seu <- FindMultiModalNeighbors(seu,
                                 reduction.list = list("pca", "lsi"),
                                 dims.list = list(1:50, 2:50)
  ) %>% # Performing UMAP and Clustering
    RunUMAP(nn.name = "weighted.nn",
            reduction.name = "wnn.umap",
            reduction.key = "wnnUMAP_",
            n.neighbors = 50,
            min.dist = 0.001,
            spread = 1) %>%
    FindClusters(graph.name = "wsnn",
                 resolution = 1.4,
                 algorithm = 3)

  return(seu)

})


# Combining 3D and 2D of each sample for Integration

Healthy1 <- merge(seurat_list[["Healthy1_2D"]],
                  y = seurat_list[["Healthy1_3D"]],
                  add.cell.ids = c("2D", '3D'),
                  project = "Healthy1")

Healthy2 <- merge(seurat_list[["Healthy2_2D"]],
                  y = seurat_list[["Healthy2_3D"]],
                  add.cell.ids = c("2D", "3D"),
                  project = "Healthy2")

UC1 <- merge(seurat_list[["UC1_2D"]],
             y = seurat_list[["UC1_3D"]],
             add.cell.ids = c("2D", "3D"),
             project = "UC1")

UC2 <- merge(seurat_list[["UC2_2D"]],
             y = seurat_list[["UC2_3D"]],
             add.cell.ids = c("2D", "3D"),
             project = "UC2")

# Ensuring the SCT is the active Assay
DefaultAssay(Healthy1) <- "SCT"
DefaultAssay(Healthy2) <- "SCT"
DefaultAssay(UC1) <- "SCT"
DefaultAssay(UC2) <- "SCT"

# Placing everything in a list
s_standard <- list(Healthy1, Healthy2, UC1, UC2)


# Performing Integration

for (i in 1:length(s_standard)) {
  s_standard[[i]] <- FindVariableFeatures(s_standard[[i]],
                                          selection.method = "vst",
                                          nfeatures = 3000,
                                          verbose = FALSE)
}

# Getting Anchors
s.anchors_standard <- FindIntegrationAnchors(object.list = s_standard,
                                             dims = 1:30)

# Integrating Datasets
s.integrated_standard <- IntegrateData(anchorset = s.anchors_standard,
                                       dims = 1:30)

# Making the new "integrated" Assay the Default Assay
DefaultAssay(s.integrated_standard) <- "integrated"


# Scaling Data and performing PCA and UMAP reductions and clustering

s.integrated_standard <- ScaleData(s.integrated_standard) %>%
  RunPCA() %>%
  FindNeighbors(reduction = "pca", dims = 1:30) %>%
  RunUMAP(dims = 1:30) %>%
  FindClusters(resolution = 0.35, graph.name = "integrated_snn")

# Saving clustering information in the metadata
s.integrated_standard[["res_35"]] <- Idents(s.integrated_standard)


# Relabeling samples names for plotting

Idents(s.integrated_standard) <- "orig.ident"

s.integrated_standard <- RenameIdents(s.integrated_standard,
                                      "Healthy1_2D" = "Healthy 1 (2D)",
                                      "Healthy1_3D" = "Healthy 1 (3D)",
                                      "Healthy2_2D" = "Healthy 2 (2D)",
                                      "Healthy2_3D" = "Healthy 2 (3D)",
                                      "UC1_2D" = "Active UC (2D)",
                                      "UC1_3D" = "Active UC (3D)",
                                      "UC2_2D" = "Inactive UC (2D)",
                                      "UC2_3D" = "Inactive UC (3D)"
)


# Saving new names in the metadata

s.integrated_standard[["condition"]] <- Idents(s.integrated_standard)

# Plotting
p <- DimPlot(s.integrated_standard, group.by = "condition", label = FALSE) + theme(plot.title = element_blank())

ggsave("Figure_2C.png", device = "png", bg= "white", dpi = 600)




#Doing cell annotation

s.integrated_standard <- RenameIdents(s.integrated_standard,   "0" = "Transitional_2D",
                                      "1" = "Late_G1_3D", #  Differentiated_3D
                                      "2" = "Differentiated_2D", # Differentiated_2D
                                      "3" = "Early_G1_3D", #  Intermediate_3D
                                      "4" = "Intermediate_2D",
                                      "5" = "BMI+Stem_2D", #  Progenitor1_2D
                                      "6" = "LGR5+Stem_3D", #  Progenitor_3D
                                      "7" = "Undifferentiated_2D", #  Progenitor2_2D
                                      "8" = "M-like_2D", # Mcell
                                      "9" = "UC-specific_3D", # UC_3D
                                      "10" = "Hes+Stem_2D" # ProgenitorHes_2D
                                      )

# Saving cell annotation
s.integrated_standard[["cell_anno"]] <- Idents(s.integrated_standard)



p <- DimPlot(s.integrated_standard, group.by = "cell_anno",label = TRUE, repel = TRUE, label.box = TRUE) +
  theme_void() +
  theme(legend.position =  "none",
        title = element_blank())

ggsave("Figure_2D.png", plot = p ,device = "png", bg= "white", dpi = 600)



# Saving Integrated Object
write_rds(s.integrated_standard, "Seurat_file.rds")







