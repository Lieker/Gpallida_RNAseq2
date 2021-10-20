plotcounts2 <- function(gene, d = dds) {
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
  colnames(data)[2] <- "counts"
  
  data.summary <- data %>% group_by(treatment, time) %>% summarise(sd = sd(counts, na.rm = TRUE),
                                                                   counts = mean(counts))
  
  g <- ggplot(data, aes(x=time, y=counts, colour = treatment)) + 
    geom_jitter(position = position_jitter(0.2),
                alpha = 0.8) +
    geom_line(data = data.summary, 
              aes(colour = treatment, 
                  group = treatment),
              size = 1.2) +
    scale_color_brewer(palette="Set1") +
    geom_point(data = data.summary, size = 3) +
    theme_classic() +
    ggtitle(paste0(gene," expression levels")) +
    scale_y_log10(breaks=c(1,2,5,10,20,50,100,200,500,1000,2000,5000,10000,20000,50000))
  
  return(g)
}

plotcounts3 <- function(gene, d = dds, s = c("L2","dauerentry","dauer","dauerexit","L3")) {
  data <- plotCounts(d, gene,
                     intgroup = "group", 
                     returnData = TRUE) %>% rownames_to_column("sample")
  data$group <- factor(data$group, levels = c("L1", "L2", "dauerentry","dauer","dauerexit","L3","L4male","L4","L4soma","YA"), ordered = TRUE)
  data <- data[order(data$group), ]
  colnames(data)[2] <- "counts"
  
  data.summary <- data %>% group_by(group) %>% summarise(sd = sd(counts, na.rm = TRUE),
                                                                   counts = mean(counts))
  
  sel <- s
  data <- data[data$group %in% sel,]
  data.summary <- data.summary[data.summary$group %in% sel,]
  
  g <- ggplot(data, aes(x=group, y=counts)) + 
    geom_jitter(position = position_jitter(0.2),
                alpha = 0.8) +
    geom_line(data = data.summary,
              aes(group = 1),
              colour = "red3") +
    geom_point(data = data.summary, size = 3, colour = "red3") +
    theme_classic() +
    ggtitle(paste0(gene," expression levels")) +
    scale_y_log10(breaks=c(1,2,5,10,20,50,100,200,500,1000,2000,5000,10000,20000,50000)) +
    theme(axis.text.x = element_text(angle = 45, hjust=1))
  
  return(g)
}
