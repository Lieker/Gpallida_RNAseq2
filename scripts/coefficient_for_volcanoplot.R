library(DESeq2)
source("scripts/get_DESeq_dds.R")

coefficient <- function(counts_csv_file = "input/counts.csv",
                        xp_design_csv_file = "input/xp_design.csv",
                        trtm = c("water","solA"),
                        ref_treatment = "water",
                        treatment2 = "solA",
                        method = "treatment", #this parameter chooses which formula design will be chosen: ~treatment or ~solA
                        tp = 8) {
  dds <- get_DESeq_dds(counts_csv_file,
                       xp_design_csv_file,
                       trtm,
                       ref_treatment,
                       treatment2,
                       method,
                       tp)
  if(method == "treatment"){
        r <- resultsNames(dds)[2]
    dds$condition <- relevel(dds$treatment, ref = ref_treatment)
  }
  
  return(r)
}
