library(DESeq2)
source("scripts/get_DESeq_dds.R")

get_unfiltered_res_dds <- function(counts_csv_file = "input/counts.csv",
                                   xp_design_csv_file = "input/xp_design.csv",
                                   th = 8,
                                   trtm = c("water", "solA")
                                   ) {
  dds <- get_DESeq_dds(counts_csv_file,
                       xp_design_csv_file)
  xp_design$treatment <- factor(xp_design$treatment, levels = c("water", "solA", "PRD"), ordered = TRUE)
  xp_design <- xp_design[order(xp_design$treatment),]
  target <- xp_design %>% dplyr::filter(xp_design$tp %in% th) 
  target <- target %>% dplyr::filter(target$treatment %in% trtm)
  target <- unique(as.character(target$group))
  t1 <- target[1]
  t2 <- target[2]
  res <- results(dds, contrast = c("group", t2, t1))
  return(res)
}
