source("scripts/get_unfiltered_res_dds.R")

get_shrunk_filtered_results <- function(d = dds,
                                        counts_csv_file = "input/counts.csv",
                                        xp_design_csv_file = "input/xp_design.csv",
                                        th = 8,
                                        trtm = c("water", "solA"),
                                        log2FC_threshold = 0,
                                        padj_threshold = 0.05) {
  res <- get_unfiltered_res_dds(counts_csv_file,
                                xp_design_csv_file,
                                th,
                                trtm)
  target <- xp_design %>% dplyr::filter(xp_design$tp %in% th) 
  target <- target %>% dplyr::filter(target$treatment %in% trtm)
  target <- unique(as.character(target$group))
  res_shr <- lfcShrink(dds = d,
                       res = res,
                       type = "ashr",
                       contrast = c("group", target[2], target[1])) %>% as.data.frame()
  res_shr$padj[is.na(res_shr$padj)] <- 0.99
  res_shr_filtered <- res_shr %>% as.data.frame() %>% dplyr::filter(., padj < padj_threshold) %>% 
    dplyr::filter(., log2FoldChange > log2FC_threshold | log2FoldChange < -log2FC_threshold)
  return(res_shr_filtered)
}
