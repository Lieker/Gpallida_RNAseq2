library(DESeq2)
source("scripts/get_DESeq_dds.R")

get_unfiltered_res_dds <- function(counts_csv_file = "input/counts.csv",
                                   xp_design_csv_file = "input/xp_design.csv",
                                   trtm = c("water","solA"),
                                   ref_treatment = "water",
                                   treatment2 = "solA",
                                   method = "treatment",
                                   tp = 8) {
  dds <- get_DESeq_dds(counts_csv_file,
                       xp_design_csv_file,
                       trtm,
                       ref_treatment,
                       treatment2,
                       method,
                       tp)
  res <- results(dds)
  return(res)
}
