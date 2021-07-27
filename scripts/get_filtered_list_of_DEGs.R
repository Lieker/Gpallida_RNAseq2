source("scripts/get_unfiltered_res_dds.R")

get_list_of_DEGs <- function(counts_csv_file = "input/counts.csv",
                             xp_design_csv_file = "input/xp_design.csv",
                             th = 8,
                             trtm = c("water","solA"),
                             log2FC_threshold = 0,
                             padj_threshold = 0.05) {
  res <- get_unfiltered_res_dds(counts_csv_file,
                                xp_design_csv_file,
                                th,
                                trtm) %>% as.data.frame()
  res_filtered <- dplyr::filter(res, padj < padj_threshold) %>% 
    dplyr::filter(., log2FoldChange > log2FC_threshold | log2FoldChange < -log2FC_threshold)
  
  
  return(res_filtered)
}


