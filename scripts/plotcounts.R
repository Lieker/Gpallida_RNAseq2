plotcounts <- function(gene, d = dds) {
  data <- plotCounts(d, gene,
                     intgroup = "group", 
                     returnData = TRUE) %>% rownames_to_column("sample") %>% separate(., 
                                                                                      col = group, 
                                                                                      into = c("treatment", "time"), 
                                                                                      sep = "_", )
  data$treatment <- factor(data$treatment, levels = c("water", "solA", "PRD"), ordered = TRUE)
  data <- data[order(data$treatment), ]
  data$time <- factor(data$time, levels = c(8,24,48), ordered = TRUE)
  data <- data[order(data$time), ]
  
  g <- ggplot(data, aes(x=time, y=count)) + 
    geom_jitter(aes(shape = time, color = treatment), 
                position = position_jitterdodge(jitter.width = 0.1, dodge.width = 0.8), size = 2) +
    theme_classic() +
    ggtitle(paste0(gene," expression levels")) +
    scale_y_log10(breaks=c(1,2,5,10,20,50,100,200,500,1000,2000,5000,10000,20000,50000)) +
    scale_color_brewer(palette="Set1")
  return(g)
}
