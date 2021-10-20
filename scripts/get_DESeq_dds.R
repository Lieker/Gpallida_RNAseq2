library(DESeq2)
library(tidyr)
library(tibble)

get_DESeq_dds <- function(counts_csv_file = "input/counts.csv",
                          xp_design_csv_file = "input/xp_design.csv"
                          ) {
  counts <- read.csv(file = counts_csv_file, 
                     header = TRUE, 
                     stringsAsFactors = FALSE, 
                     fileEncoding = "UTF-8-BOM") %>% column_to_rownames("Geneid")
  counts <- as.data.frame(t(counts)) %>% rownames_to_column("sample")
  
  xp_design <<- read.csv(file = xp_design_csv_file, header = TRUE, check.names = FALSE, fileEncoding = "UTF-8-BOM") 
  xp_design$tp <- as.factor(xp_design$tp)
  
  counts <- counts %>% column_to_rownames("sample") %>% t()
  dds <- DESeqDataSetFromMatrix(countData = counts, colData = xp_design, design = ~ group)
  
  keep <- rowSums(counts(dds)) >= 10
  dds <- dds[keep,]
  dds <- DESeq(dds)
  return(dds)
}

get_DESeq_dds2 <- function(counts_csv_file = "input/countsW.csv",
                          xp_design_csv_file = "input/xp_designW.csv"
) {
  counts <- read.csv(file = counts_csv_file, 
                     header = TRUE, 
                     stringsAsFactors = FALSE, 
                     fileEncoding = "UTF-8-BOM") %>% column_to_rownames("Geneid")
  counts <- as.data.frame(t(counts)) %>% rownames_to_column("sample")
  
  xp_design <- read.csv(file = xp_design_csv_file, header = TRUE, check.names = FALSE, fileEncoding = "UTF-8-BOM") 
  
  counts <- counts %>% column_to_rownames("sample") %>% t()
  dds <- DESeqDataSetFromMatrix(countData = counts, colData = xp_design, design = ~ group)
  
  keep <- rowSums(counts(dds)) >= 10
  dds <- dds[keep,]
  dds <- DESeq(dds)
  return(dds)
}
