source("scripts/get_DESeq_dds.R")
source("scripts/get_unfiltered_res_dds.R")

make_volcanoplot <- function(counts_csv_file = "input/counts.csv",
                             xp_design_csv_file = "input/xp_design.csv",
                             trtm = c("water","solA"),
                             th = 8,
                             log2FC_threshold = 0,
                             FCcutoff_volcano = 0,
                             padj_threshold = 0,
                             xsize = 2,
                             ttl = ""
) {
  dds <- get_DESeq_dds(counts_csv_file,
                       xp_design_csv_file)
  
  res <- get_unfiltered_res_dds(counts_csv_file,
                                xp_design_csv_file,
                                th,
                                trtm)
  
  target <- xp_design %>% dplyr::filter(xp_design$tp %in% th) 
  target <- target %>% dplyr::filter(target$treatment %in% trtm)
  target <- unique(as.character(target$group))
  
  shrunk <- lfcShrink(dds = dds,
                      res = res,
                      type = "ashr",
                      contrast = c("group", target[2], target[1]))
  v <- EnhancedVolcano(toptable = shrunk,
                       x = "log2FoldChange",
                       y = "padj",
                       title = ttl,
                       subtitle = "",
                       pCutoff = padj_threshold,
                       FCcutoff = FCcutoff_volcano,
                       lab = rownames(shrunk),
                       xlim = c(-(xsize), xsize))
  return(v)
}
