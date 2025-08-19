#' @title GSEA Analysis
#'
#' @description
#' Perform GSEA  analysis of DEG results and generates plots with results with ClusterProfiler package
#' For more information please look clusterProfiler package
#'
#' @param x A dataframe with DEG results.
#' @param comparing A string indicating the comparison performed. This string will be use to label the plots
#' @param species string indicating the specie to use. ClusterProfiler can use human or mouse
#' @param category Which category to use. You can only choose one at a time
#' @param subcategory If the category can be separated in further categories. Indicate which one to use here
#' @param OrgDb the OrgDb file that will be use to obtain the information for GSEA analysis
#' @param top The top n pathways to plot
#' @param log_FoldChange The log2 Fold Change threshold to subset the DEG dataset
#' @param pvalue_Cutoff The pvalue threshold to subset the DEG dataset
#' @param qvalue_Cutoff The pvalue cutoff to use during GSEA analysis
#' @param saveData Logical. Indicate if you want the plots and GSEA results to be saved as png and csv, respectively in the working directory
#' @param is.seurat Logical. Indicate if the DEG dataset comes from `FindMarkers()` function from Seurat
#' @param return_object Logical. Return the GSEA results
#'
#' @references 
#' Yu G, Wang L, Han Y, He Q (2012). “clusterProfiler: an R package for comparing biological themes among gene clusters.” 
#' OMICS: A Journal of Integrative Biology, 16(5), 284-287. doi:10.1089/omi.2011.0118.
#'
#' @returns 


GSEA_analysis <- function(x, 
                          comparing,
                          species = "Homo sapiens", 
                          category = "H", 
                          subcategory = NULL,
                          OrgDb = org.Hs.eg.db, 
                          top = 20, 
                          log_FoldChange = 1, 
                          pvalue_Cutoff = 1, 
                          qvalue_Cutoff = 1, 
                          saveData = FALSE, 
                          is.seurat = FALSE, 
                          return_object= FALSE) {
  
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
  
  if(is.seurat == T) {
    datframe <- dplyr::rename(datframe, log2FoldChange = avg_log2FC)
  }
  
  entrez_id <- clusterProfiler::bitr(rownames(datframe), 
                                     fromType = "SYMBOL" , 
                                     toType = "ENTREZID",
                                     OrgDb = OrgDb)
  
  datframe$entrezid <- entrez_id$ENTREZID[match(rownames(datframe), entrez_id$SYMBOL)]
  
  #Getting the GSEA for the msigdbr category choosen (default is the Hallmark category)
  
  gene_list <- datframe %>% 
    dplyr::arrange(-log2FoldChange) %>% 
    dplyr::filter(abs(log2FoldChange) > log_FoldChange & !is.na(.$entrezid)) %>% 
    dplyr::select(entrezid, log2FoldChange) %>% 
    tibble::deframe()
  
  in_hallmark <- names(gene_list) %in% hallmark_data$entrez_gene
  
  if(sum(in_hallmark) > 0 ) {
    gsea_msigdb <- clusterProfiler::GSEA(geneList = gene_list, 
                                         TERM2GENE = hallmark_data,
                                         pAdjustMethod = "BH", 
                                         pvalueCutoff = pvalue_Cutoff)
  }
  #Changing ENTREZID name to Gene name
  if(exists("gsea_msigdb")){
    if(length(gsea_msigdb@result$NES)!= 0){
      gsea_msigdb@result$core_enrichment <- sapply(gsea_msigdb@result$core_enrichment, name_change)
    }
  }
  print("GSEA analysis for the chosen category done!")
  
  
  
  #GRAPHS
  
  print(paste("Printing the Graphs for top",top,"upregulated and downregulated GSEA results", sep=" "))
  list_plots <- list()
  # MSigDB Hallmark (upregulated)
  if(exists("gsea_msigdb")){
    if(length(gsea_msigdb@result$NES)!= 0){
      p1 <- gsea_msigdb %>% 
        dplyr::arrange(-NES) %>% 
        enrichplot::dotplot(showCategory = top, x = "NES") +
        ggplot2::ggtitle(paste(comparison, "Upregulated", sep = " ")) + 
        ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5))
      
      list_plots[[length(list_plots)+1]] <- p1
    }
  }
  
  # MSigDB Hallmark (downregulated)
  if(exists("gsea_msigdb")){
    if(length(gsea_msigdb@result$NES)!= 0){
      
      p2<- gsea_msigdb %>% 
        dplyr::arrange(NES)  %>% 
        enrichplot::dotplot(showCategory = top, x = "NES") + 
        ggplot2::ggtitle(paste(comparison, "Downregulated", sep = " ")) + 
        ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5))
      
      list_plots[[length(list_plots)+1]] <- p2
    }
  }
  
  
  if(length(list_plots) >= 1){
    print(cowplot::plot_grid(plotlist = list_plots))
  } else { print( "No graphs were created")}
  
  
  
  if(saveData == TRUE){
    
    print("Saving dataframes, objects, and graphs for GSEA results and graphs")
    print("Files stored in a folder inside working directory")
    
    dir.create(paste(getwd(), comparing, sep = "/"))
    
    for_saving <- stringr::str_replace_all(comparing, " ", "_")
    
    #Saving Objects
    if(exists("gsea_msigdb")){
      if(length(gsea_msigdb@result$NES)!= 0){
        if (is.null(subcategory)){
          save(gsea_msigdb, 
               file = paste(getwd(), 
                            comparing,
                            paste("Category", category,for_saving, "GSEA_MSIGBDR.RData", sep = "_"),
                            sep= "/"))
        } else {
          save(gsea_msigdb, 
               file = paste(getwd(), 
                            comparing, 
                            paste("Category", category, subcategory ,for_saving, "GSEA_MSIGBDR.RData", sep = "_"),
                            sep= "/"))
        }
      } else {print("No NES identified")}
    } else {print("There is no gsea object to save")}
    
    #Saving the GSEA results
    if(exists("gsea_msigdb")){
      if(length(gsea_msigdb@result$NES)!= 0){
        gsea_msigdb_results <- gsea_msigdb@result %>% dplyr::arrange(-NES)
        write.csv(gsea_msigdb_results, 
                  file =  paste(getwd(),
                                comparing,
                                paste0("Results ", 
                                      category, " ",
                                      "comparing ", comparison, ".csv"), 
                                sep= "/"))
      } else {print("There is no gsea data to make a spreadsheet")}
    } else {print("There is no gsea data to make a spreadsheet")}
    
    
    #Graphs
    if(exists("p1")){
      ggsave(plot = p1, 
             path = paste(getwd(), comparing, sep = "/"), 
             filename = paste0("Top ",
                              top," ", 
                              "upregulated pathways Category ", 
                              category," ",
                              "comparing ", 
                              comparison, 
                              ".png"), 
             height =8, width = 10, units = "in")
      
    } else {
      print("There is no graph with upregulated gsea genes to save")
      }
    
    if(exists("p2")){
      ggsave (plot = p2, 
              path = paste(getwd(), comparing, sep = "/"), 
              filename = paste0("Top ",
                                top," ", 
                                "downregulated pathways Category ", 
                                category," ",
                                "comparing ", 
                                comparison, 
                                ".png"), 
              height =8, width = 10, units = "in")
    } else {
      print("There is no graph with downregulated gsea genes to save")
      }
  }
  
  print("Complete!")
  
  
  if (return_object){
    to_return <- gsea_msigdb
    return(to_return)
    
  }
}



#Created by: Rodolfo Ismael Cabrera Silva & Zachary Wilson
#Required packages:
#tidyverse, ClusterProfiler, org.Hs.eg.db, enrichplot, msigdbr, cowplot, ggplot2
