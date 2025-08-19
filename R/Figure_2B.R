
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

# Subsetting

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

  return(seu)


})


# Merging Files based on plating condition (2D or 3D)

# Making a list with 3 out of the 4 Seurat Objects to use in merge() function
other_2d <- list(seurat_list[[2]], seurat_list[[5]], seurat_list[[6]])
other_3d <- list(seurat_list[[4]], seurat_list[[7]], seurat_list[[8]])

# Merging the Seurat Objects
merge2D <- merge(x = seurat_list[[1]], y = other_2d, project = "2D")
merge3D <- merge(x = seurat_list[[3]], y = other_3d, project = "3D")

# Redoing SCTransform based on
merge2D <- merge2D %>%
  SCTransform() %>%
  FindVariableFeatures() %>%
  RunPCA(assay = "SCT")

merge3D <- merge3D %>%
  SCTransform() %>%
  FindVariableFeatures() %>%
  RunPCA(assay = "SCT")

p <- DimPlot(merge2D, reduction = "pca", group.by = "orig.ident")

p

p <- DimPlot(merge3D, reduction = "pca", group.by = "orig.ident")

ggsave(filename = "Figure_2B.png", plot = p, device = "png", bg = "white", dpi = 600)
