library(DESeq2)
library(tidyr)
library(tibble)

get_DESeq_dds <- function(counts_csv_file = "input/counts.csv",
                          xp_design_csv_file = "input/xp_design.csv",
                          trtm = c("water","solA"),
                          ref_treatment = "water",
                          treatment2 = "solA",
                          method = "treatment",
                          tp = 8
                          ) {
  counts <- read.csv(file = counts_csv_file, 
                     header = TRUE, 
                     stringsAsFactors = FALSE, 
                     fileEncoding = "UTF-8-BOM") %>% column_to_rownames("Geneid")
  counts <- as.data.frame(t(counts)) %>% rownames_to_column("sample")
  
  xp_design <- read.csv(file = xp_design_csv_file, header = TRUE, check.names = FALSE, fileEncoding = "UTF-8-BOM") 
  xp_design$tp <- as.factor(xp_design$tp)
  
  if(method == "treatment"){
    xp_design <- xp_design[which(xp_design$treatment %in% trtm),] %>%
    dplyr::filter(treatment == ref_treatment | treatment == treatment2)
    xp_design <- xp_design[which(xp_design$tp %in% tp),]
    counts <- counts %>% dplyr::filter(sample %in% xp_design$sample) %>% column_to_rownames("sample") %>% t()
    dds <- DESeqDataSetFromMatrix(countData = counts, colData = xp_design, design = ~ treatment)
  }
  
  dds <- DESeq(dds)
  return(dds)
}
