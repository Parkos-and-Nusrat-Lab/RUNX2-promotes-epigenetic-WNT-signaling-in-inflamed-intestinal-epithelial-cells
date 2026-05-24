#' @title GSEA Analysis
#'
#' @description
#' Perform GSEA  analysis of DEG results and generates plots with results with ClusterProfiler package
#' For more information please look clusterProfiler package
#'
#' @param x A dataframe with DEG results.
#' @param species string indicating the specie to use. ClusterProfiler can use human or mouse
#' @param category Which category to use. You can only choose one at a time
#' @param subcategory If the category can be separated in further categories. Indicate which one to use here
#' @param OrgDb the OrgDb file that will be use to obtain the information for GSEA analysis
#' @param log_FoldChange The log2 Fold Change threshold to subset the DEG dataset
#' @param pvalue_Cutoff The pvalue threshold to subset the DEG dataset
#' @param qvalue_Cutoff The pvalue cutoff to use during GSEA analysis
#' 
#'
#' @references 
#' Yu G, Wang L, Han Y, He Q (2012). “clusterProfiler: an R package for comparing biological themes among gene clusters.” 
#' OMICS: A Journal of Integrative Biology, 16(5), 284-287. doi:10.1089/omi.2011.0118.
#'
#' @returns 


GSEA_function <- function(x,
                          comparing, 
                          species = "Homo sapiens", 
                          category = "H", 
                          subcategory = NULL, 
                          OrgDb = org.Hs.eg.db, 
                          log_FoldChange = 1, 
                          pvalue_Cutoff = 1,
                          qvalue_Cutoff = 1) {
  
  
  name_change <- function(i){
    y <- strsplit(i, '/') [[1]]
    y <- entrez_id$SYMBOL[match(y, entrez_id$ENTREZID)]
    paste(y, collapse = '/')
  }
  
  
  comparison <- comparing
  
  hallmark_data <- msigdbr::msigdbr(species = species, 
                                    category = category, 
                                    subcategory = subcategory) %>% 
    dplyr::select(gs_name, entrez_gene)
  
  datframe = x
  
  entrez_id <- clusterProfiler::bitr(rownames(datframe), fromType = "SYMBOL" , toType = "ENTREZID", OrgDb = OrgDb)
  
  datframe$entrezid <- entrez_id$ENTREZID[match(rownames(datframe), entrez_id$SYMBOL)]
  
  
  #Getting the GSEA for the msigdbr category choosen (default is the Hallmark category)
  
  gene_list <- datframe %>% 
    dplyr::arrange(-log2FoldChange) %>% 
    dplyr::filter(abs(log2FoldChange) > log_FoldChange & p_val_adj < 0.05 &!is.na(.$entrezid)) %>%
    dplyr::select(entrezid, log2FoldChange) %>% 
    tibble::deframe()
  
  in_hallmark <- names(gene_list) %in% hallmark_data$entrez_gene
  
  
  gsea_msigdb <- clusterProfiler::GSEA(geneList = gene_list,
                                       TERM2GENE = hallmark_data,
                                       pAdjustMethod = "BH",
                                       pvalueCutoff = pvalue_Cutoff
                                       )
  

  gsea_msigdb@result$core_enrichment <- sapply(gsea_msigdb@result$core_enrichment, name_change)
      
      
  return(gsea_msigdb)
  
}
