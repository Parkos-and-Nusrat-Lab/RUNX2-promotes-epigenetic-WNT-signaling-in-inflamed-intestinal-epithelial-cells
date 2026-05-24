# Opening libraries
library(ArchR)
library(pheatmap)
library(magick)
library(BSgenome.Hsapiens.UCSC.hg38)
library(ggseqlogo)

# Adding Genome to ArchR
addArchRGenome("hg38")

# Setting Working directory
setwd('/path/to/directory/')

# Setting seed
set.seed(1)

# Preparing Arrow Files


inputFiles3D <- c('3D_healthy1' = "/path/to/file/atac_fragments.tsv.gz",
                  '3D_UC1' = "/path/to/file/atac_fragments.tsv.gz",
                  '3D_healthy2' = "/path/to/file/atac_fragments.tsv.gz",
                  '3D_UC2' = "/path/to/file/atac_fragments.tsv.gz")

inputFiles2D <- c('2D_healthy1' = "/path/to/file/atac_fragments.tsv.gz",
                  '2D_UC1' = "/path/to/file/atac_fragments.tsv.gz",
                  '2D_healthy2' = "/path/to/file/atac_fragments.tsv.gz",
                  '2D_UC2' = "/path/to/file/atac_fragments.tsv.gz")

# Making Arrow files for 2D samples
ArrowFiles2D <- createArrowFiles(
  inputFiles = inputFiles2D,
  sampleNames = names(inputFiles2D),
  minTSS = 3, 
  minFrags = 1000,
  addTileMat = TRUE,
  addGeneScoreMat = TRUE, force = TRUE
)

# Making Arrow files for 3D samples
ArrowFiles3D <- createArrowFiles(
  inputFiles = inputFiles3D,
  sampleNames = names(inputFiles3D),
  minTSS = 3, 
  minFrags = 1000,
  addTileMat = TRUE,
  addGeneScoreMat = TRUE, force = TRUE
)


# Performing doubles scores

doubScores <- addDoubletScores(
  input = ArrowFiles2D,
  k = 10,
  knnMethod = "UMAP", 
  LSIMethod = 1,
)
doubScores <- addDoubletScores(
  input = ArrowFiles3D,
  k = 10,
  knnMethod = "UMAP", 
  LSIMethod = 1,
)

projCol <- ArchRProject(
  ArrowFiles = c(ArrowFiles2D,ArrowFiles3D),
  outputDirectory = "projCol",
  copyArrows = TRUE
)

# Extracting the names to save as metadata
bioNames <- gsub("_healthy2|_healthy1|_UC1|_UC2","",projCol$Sample)
bioNames2 <- gsub("3D_|2D_","",projCol$Sample)

projCol$bioNames <- bioNames
projCol$bioNames2 <- bioNames2


# SAving ArchR Project
saveArchRProject(ArchRProj = projCol, outputDirectory = "projCol", load = FALSE)

# Performing QC filtering TSS Enrichment
idxPass <- which(projCol$TSSEnrichment >= 5)
cellsPass <- projCol$cellNames[idxPass]
projCol2 <- filterDoublets(projCol[cellsPass, ]) #pass to second project
###################################LSI#############################

# Performing LSI
projCol2 <- addIterativeLSI(
  ArchRProj = projCol2,
  useMatrix = "TileMatrix",
  name = "IterativeLSI",
  iterations = 10,
  clusterParams = list(
    resolution = 0.8,
    sampleCells = 10000,
    n.start = 10
  ),
  varFeatures = 50000,
  dimsToUse = 1:30, force = TRUE
)

#Harmony batch correction
projCol2 <- addHarmony(ArchRProj = projCol2,
                       reducedDims = "IterativeLSI",
                       force = TRUE, 
                       name = "Harmony",
                       groupBy = "bioNames2")

#Adding clusters
projCol2 <- addClusters(input = projCol2,
                        reducedDims = "IterativeLSI",
                        force = TRUE, 
                        method = "Seurat",
                        name = "Clusters",
                        resolution = 0.8)

#UMAP
projCol2 <- addUMAP(ArchRProj = projCol2,
                    reducedDims = "IterativeLSI",
                    force = TRUE,  
                    name = "UMAP", 
                    nNeighbors = 30, 
                    minDist = 0.5, 
                    metric = "cosine")

#UMAP harmony
projCol2 <- addUMAP(ArchRProj = projCol2, 
                    reducedDims = "Harmony", 
                    force = TRUE, 
                    name = "UMAPHarmony", 
                    nNeighbors = 30, 
                    minDist = 0.5,
                    metric = "cosine")

#Adding clusters
projCol2 <- addClusters(input = projCol2,
                        reducedDims = "Harmony",
                        force = TRUE, 
                        method = "Seurat",
                        name = "HarmonyClusters",
                        resolution = 0.8)


# Plotting
p1 <- plotEmbedding(ArchRProj = projCol2, 
                    colorBy = "cellColData", 
                    name = "HarmonyClusters", 
                    embedding = "UMAPHarmony")
p2 <- plotEmbedding(ArchRProj = projCol2, 
                    colorBy = "cellColData", 
                    name = "bioNames2", 
                    embedding = "UMAPHarmony")
ggAlignPlots(p1, p2, type = "h")
plotPDF(p1,p2, name = "Plot-HarmonyUMAP-Sample-Clusters.pdf", 
        ArchRProj = projCol2, addDOC = TRUE, width = 10, height = 5)

#######Looking at markers of differentiation######
markersGS <- getMarkerFeatures(
  ArchRProj = projCol2,
  useMatrix = "GeneScoreMatrix",
  groupBy = "HarmonyClusters",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  testMethod = "wilcoxon"
)
markerGenes <- c(
  "LGR5", "BMI1",
  "MUC2", "REG4",
  "TNFAIP2",
  "ITPR2",
  "CHGA","LYZ",
  "MKI67",
  "HES1", 
  "KRT8", 
  "CDH17", "GPA33", "KRT20", 
  "THBS2", "THBS1"
)

heatmapGS <- markerHeatmap(
  seMarker = markersGS,
  cutOff = "FDR <= 0.01 & Log2FC >= 1.25",
  labelMarkers = markerGenes,
  transpose = TRUE
)
ComplexHeatmap::draw(heatmapGS, heatmap_legend_side = "bot", annotation_legend_side = "bot")
#heatmapGS@row_order <- c(10,11,12,3,4,5,1,2,6,8,7,9)
#ComplexHeatmap::draw(heatmapGS, heatmap_legend_side = "bot", annotation_legend_side = "bot")
p <- plotEmbedding(
  ArchRProj = projCol2,
  colorBy = "GeneScoreMatrix",
  name = markerGenes,
  embedding = "UMAPHarmony",
  quantCut = c(0.01, 0.95),
  imputeWeights = NULL
)
p2 <- lapply(p, function(x){
  x + guides(color = FALSE, fill = FALSE) +
    theme_ArchR(baseSize = 6.5) +
    theme(plot.margin = unit(c(0, 0, 0, 0), "cm")) +
    theme(
      axis.text.x=element_blank(),
      axis.ticks.x=element_blank(),
      axis.text.y=element_blank(),
      axis.ticks.y=element_blank()
    )
})
do.call(cowplot::plot_grid, c(list(ncol = 3),p2))
plotPDF(plotList = p,
        name = "Plot-UMAP-Marker-Genes-WO-Imputation.pdf",
        ArchRProj = projCol2,
        addDOC = FALSE, width = 5, height = 5)
projCol2 <- addImputeWeights(projCol2)
p <- plotEmbedding(
  ArchRProj = projCol2,
  colorBy = "GeneScoreMatrix",
  name = markerGenes,
  embedding = "UMAPHarmony",
  imputeWeights = getImputeWeights(projCol2)
)
p2 <- lapply(p, function(x){
  x + guides(color = FALSE, fill = FALSE) +
    theme_ArchR(baseSize = 6.5) +
    theme(plot.margin = unit(c(0, 0, 0, 0), "cm")) +
    theme(
      axis.text.x=element_blank(),
      axis.ticks.x=element_blank(),
      axis.text.y=element_blank(),
      axis.ticks.y=element_blank()
    )
})
do.call(cowplot::plot_grid, c(list(ncol = 3),p2))
plotPDF(plotList = p,
        name = "Plot-UMAP-Marker-Genes-W-Imputation.pdf",
        ArchRProj = projCol2,
        addDOC = FALSE, width = 5, height = 5)


######################################New Module scoring######################
features <- list(
  ISCScore = c("LGR5", 'LRIG1',"OLFM4"),
  DiffScore = c("KRT20", "KRT8", "MUC2", "VIL1", "GPA33")
)
projCol2 <- addModuleScore(projCol2,
                               useMatrix = "GeneScoreMatrix",
                               name = "Module",
                               features = features)
###################################PREP RNA DATASET###########################
scRNA_I_test <- read_rds("Path/to/scRNAseq/file/withcluster.rds")
Idents(scRNA_I_test) <- "seurat_clusters"
new.cluster.ids <- c('0' = "Transitional_2D",
                     '1' = "Differentiated_3D",
                     '2'="Differentiated_2D",
                     '3'="Intermediate_3D",
                     '4'= "Intermediate_2D",
                     '5'="Progenitor_2D",
                     '6'="Progenitor_3D",
                     '7'="Undifferentiated_2D",
                     '8'="M_like",
                     '9'="UC_specific",
                     '10'="HES1_progenitor")
scRNA_I_test <- RenameIdents(scRNA_I_test,new.cluster.ids)
scRNA_I_test$seurat_clusters <- Idents(scRNA_I_test)
######################CONSTRAINED INTEGRATION#################################
clustDiff_2D <- paste0(paste0("C",6)) 
clustDiff_3D <- paste0(paste0("C",3))
clustHes <- paste0(paste0("C",13))
clustIntermediate_2D <- paste0(c(paste0("C", 7:9),paste0("C", 14))) 
clustIntermediate_3D <- paste0((paste0("C", 1:2)))
clustMcell <- paste0(paste0("C",10)) 
clustProg_2D <- paste0(paste0("C",11))
clustUndiff_2D <- paste0(paste0("C",15))
clustProg_3D <- paste0(paste0("C",12))
clustTransitional <- paste0(paste0("C",5))
clustUC <- paste0(paste0("C",4))

rnaDiff_2D <- CellsByIdentities(scRNA_I_test, idents = 'Differentiated_2D') %>% 
  unlist() %>% as.character()
rnaDiff_3D <- CellsByIdentities(scRNA_I_test, idents = 'Differentiated_3D') %>% 
  unlist() %>% as.character()
rnaHes <- CellsByIdentities(scRNA_I_test, idents = 'HES1_progenitor') %>% 
  unlist() %>% as.character()
rnaIntermediate_2D <- CellsByIdentities(scRNA_I_test, idents = 'Intermediate_2D') %>% 
  unlist() %>% as.character()
rnaIntermediate_3D <- CellsByIdentities(scRNA_I_test, idents = 'Intermediate_3D') %>% 
  unlist() %>% as.character()
rnaMcell <- CellsByIdentities(scRNA_I_test, idents = 'M_like') %>% 
  unlist() %>% as.character()
rnaProg_2D <- CellsByIdentities(scRNA_I_test, idents = 'Progenitor_2D') %>% 
  unlist() %>% as.character()
rnaUndiff_2D <- CellsByIdentities(scRNA_I_test, idents = 'Undifferentiated_2D') %>% 
  unlist() %>% as.character()
rnaProg_3D <- CellsByIdentities(scRNA_I_test, idents = 'Progenitor_3D') %>% 
  unlist() %>% as.character()
rnaTransitional <- CellsByIdentities(scRNA_I_test, idents = 'Transitional_2D') %>% 
  unlist() %>% as.character()
rnaUC <- CellsByIdentities(scRNA_I_test, idents = 'UC_specific') %>% 
  unlist() %>% as.character()

groupList <- SimpleList(
  Differentiated_2D = SimpleList(ATAC = projCol2$cellNames[projCol2$HarmonyClusters %in% clustDiff_2D], RNA = rnaDiff_2D),
  Differentiated_3D = SimpleList(ATAC = projCol2$cellNames[projCol2$HarmonyClusters %in% clustDiff_3D],RNA = rnaDiff_3D),
  HES1_progenitor = SimpleList(ATAC = projCol2$cellNames[projCol2$HarmonyClusters %in% clustHes], RNA = rnaHes),
  Intermediate_2D = SimpleList(ATAC = projCol2$cellNames[projCol2$HarmonyClusters %in% clustIntermediate_2D],RNA = rnaIntermediate_2D),
  Intermediate_3D = SimpleList(ATAC = projCol2$cellNames[projCol2$HarmonyClusters %in% clustIntermediate_3D], RNA = rnaIntermediate_3D),
  M_like = SimpleList(ATAC = projCol2$cellNames[projCol2$HarmonyClusters %in% clustMcell],RNA = rnaMcell),
  Progenitor_2D = SimpleList(ATAC = projCol2$cellNames[projCol2$HarmonyClusters %in% clustProg_2D], RNA = rnaProg_2D),
  Undifferentiated_2D = SimpleList(ATAC = projCol2$cellNames[projCol2$HarmonyClusters %in% clustUndiff_2D],RNA = rnaUndiff_2D),
  Progenitor_3D = SimpleList(ATAC = projCol2$cellNames[projCol2$HarmonyClusters %in% clustProg_3D],RNA = rnaProg_3D),
  Transitional_2D = SimpleList(ATAC = projCol2$cellNames[projCol2$HarmonyClusters %in% clustTransitional],RNA = rnaTransitional),
  UC_specific = SimpleList(ATAC = projCol2$cellNames[projCol2$HarmonyClusters %in% clustUC],RNA = rnaUC)
)

#########CREATE a nested list of ATAC and RNA clusters for constraint
projCol2 <- addGeneIntegrationMatrix(
  ArchRProj = projCol2,
  useMatrix = "GeneScoreMatrix",
  matrixName = "GeneIntegrationMatrix",
  reducedDims = "Harmony",
  seRNA = scRNA_I_test,
  addToArrow = TRUE,
  force= TRUE,
  groupList = groupList,
  groupRNA = "seurat_clusters",
  nameCell = "predictedCell",
  nameGroup = "predictedGroup",
  nameScore = "predictedScore"
)

projCol2 <- addImputeWeights(projCol2)
p1 <- plotEmbedding(
  ArchRProj = projCol2,
  colorBy = "GeneIntegrationMatrix",
  name = markerGenes,
  continuousSet = "horizonExtra",
  embedding = "UMAPHarmony",
  imputeWeights = getImputeWeights(projCol2)
)
p2 <- plotEmbedding(
  ArchRProj = projCol2,
  colorBy = "GeneScoreMatrix",
  continuousSet = "horizonExtra",
  name = markerGenes,
  embedding = "UMAPHarmony",
  imputeWeights = getImputeWeights(projCol2)
)
p1c <- lapply(p1, function(x){
  x + guides(color = FALSE, fill = FALSE) +
    theme_ArchR(baseSize = 6.5) +
    theme(plot.margin = unit(c(0, 0, 0, 0), "cm")) +
    theme(
      axis.text.x=element_blank(),
      axis.ticks.x=element_blank(),
      axis.text.y=element_blank(),
      axis.ticks.y=element_blank()
    )
})

p2c <- lapply(p2, function(x){
  x + guides(color = FALSE, fill = FALSE) +
    theme_ArchR(baseSize = 6.5) +
    theme(plot.margin = unit(c(0, 0, 0, 0), "cm")) +
    theme(
      axis.text.x=element_blank(),
      axis.ticks.x=element_blank(),
      axis.text.y=element_blank(),
      axis.ticks.y=element_blank()
    )
})
do.call(cowplot::plot_grid, c(list(ncol = 3), p1c))
do.call(cowplot::plot_grid, c(list(ncol = 3), p2c))
plotPDF(plotList = p1,
        name = "Plot-UMAP-Marker-Genes-RNA-W-Imputation.pdf",
        ArchRProj = projCol2,
        addDOC = FALSE, width = 5, height = 5)
################################Add peaks to integrated dataset##############
projCol2 <- addGroupCoverages(ArchRProj = projCol2,
                                  maxFragments = 25 * 10^6,
                                  minCells = 40,
                                  maxCells = 2000,
                                  sampleRatio = 0.9,
                                  minReplicates = 2,
                                  maxReplicates = 4,
                                  force = TRUE,
                                  groupBy = "predictedGroup")
projCol2 <- addReproduciblePeakSet(
  ArchRProj = projCol2,
  groupBy = "predictedGroup",
  peakMethod = "Macs2",
  reproducibility = "1",
  genomesize = 2.7e9,
  maxPeaks = 1500000,
  minCells = 150,
  shift = -37,
  extsize = 73,
  extendSummits = 250,
  force = TRUE,
  verbose = TRUE,
  peaksPerCell = 5000
)
projCol2 <- addPeakMatrix(projCol2)
markerPeaks <- getMarkerFeatures(
  ArchRProj = projCol2,
  useMatrix = "PeakMatrix",
  groupBy = "predictedGroup",
  #bias = c("TSSEnrichment", "log10(nFrags)"),
  testMethod = "wilcoxon"
)
heatmapPeaks <- markerHeatmap(
  seMarker = markerPeaks,
  cutOff = "FDR <= 0.1 & Log2FC >= 0.5",
  transpose = FALSE
)
draw(heatmapPeaks, heatmap_legend_side = "bot", annotation_legend_side = "bot")
###################################Motifs##################################
projCol2 <- addMotifAnnotations(ArchRProj = projCol2, force = TRUE, motifSet = "cisbp", name = "Motif")
enrichMotifs <- peakAnnoEnrichment(
  seMarker = markerPeaks,
  ArchRProj = projCol2,
  peakAnnotation = "Motif",
  cutOff = "FDR <= 0.1 & Log2FC >= 0.25"
)
doMotifs <- peakAnnoEnrichment(
  seMarker = markerPeaks,
  ArchRProj = projCol2,
  peakAnnotation = "Motif",
  cutOff = "FDR <= 0.1 & Log2FC <= -0.5"
)

heatmapUp <- plotEnrichHeatmap(enrichMotifs, n = 30,
                               transpose = TRUE,rastr = F,labelRows = TRUE,
                               returnMatrix = T)

heatmapEM <- plotEnrichHeatmap(doMotifs, n = 50, transpose = FALSE)
ComplexHeatmap::draw(heatmapEM, heatmap_legend_side = "bot", annotation_legend_side = "bot")
plotPDF(heatmapEM, name = "Motifs-Enriched-Marker-Heatmap", width = 4, height = 6, ArchRProj = projCol2, addDOC = FALSE)
pdf(ComplexHeatmap::draw(heatmapUp, heatmap_legend_side = "bot", annotation_legend_side = "bot"))
Heatmap(heatmapUp)

####################Find motifs in genomic regions############################
pSet <- getPeakSet(ArchRProj = projCol2)
pSet$name <- paste(seqnames(pSet), start(pSet), end(pSet), sep = "_")
matches <- getMatches(ArchRProj = projCol2, name = "Motif")
rownames(matches) <- paste(seqnames(matches), start(matches), end(matches), sep = "_")
matches <- matches[pSet$name]
gr <- GRanges(seqnames = c("chr8"), ranges = IRanges(start = c(12166000), end = c(12169991))) 
queryHits <- queryHits(findOverlaps(query = pSet, subject = gr, type = "within"))
colnames(matches)[which(assay(matches[queryHits,]))]

#####################################MOTIF LETTERS########################
pwm <- getPeakAnnotation(projCol2, "Motif")$motifs[["RUNX2"]]
PWMatrixToProbMatrix <- function(x){
  if (class(x) != "PWMatrix") stop("x must be a TFBSTools::PWMatrix object")
  m <- (exp(as(x, "matrix"))) * TFBSTools::bg(x)/sum(TFBSTools::bg(x))
  m <- t(t(m)/colSums(m))
  m
}

ppm <- PWMatrixToProbMatrix(pwm)
colSums(ppm) %>% range

ggseqlogo(ppm, method = "bits")
ggseqlogo(ppm, method = "prob")


projCol2 <- saveArchRProject(ArchRProj = projCol2, outputDirectory = "Save-projCol2", load = TRUE)

