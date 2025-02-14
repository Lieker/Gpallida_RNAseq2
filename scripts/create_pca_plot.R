plot_pca <- function(count_csv_file = "input/counts.csv",
                     xp_design_csv_file = "input/xp_design.csv",
                     pc_x_axis = "PC1", 
                     pc_y_axis = "PC2",
                     tp = c(8,24,48),
                     pca_colour = "treatment",
                     pca_shape = "tp") {
  
  source("scripts/produce_scaled_counts_matrix.R")
  source("scripts/filter_counts_based_on_timepoint.R")
  source("scripts/extract_variance.R")
  
  produce_score_df <- function(counts = filtered_counts, 
                               xp_design_csv_file) {
    
    pca_results <- mypca(x = counts, center = TRUE, scale = TRUE)
    scores <- pca_results$score
    
    # keep only components of interest (numbers here to filter the score df)
    first_pc <- as.integer(gsub(pattern = "PC",replacement = "", x = pc_x_axis))
    second_pc <- as.integer(gsub(pattern = "PC",replacement = "", x = pc_y_axis))
    
    scores_filtered <- scores[,c(first_pc, second_pc)]
    
    # add xp_design extra info for future plot
    xp_design = read.csv(file = xp_design_csv_file, header = TRUE, check.names = FALSE, fileEncoding = "UTF-8-BOM") 
    xp_design$tp <- as.factor(xp_design$tp)
    
    scores_with_xp_design_info <- scores_filtered %>% 
      as.data.frame() %>% 
      rownames_to_column("sample") %>% 
      left_join(x = ., y = xp_design, by = "sample") %>% 
      column_to_rownames("sample")  
    
    return(scores_with_xp_design_info)
  }
  
  
  
  # scale raw counts the DESeq2 way
  scaled_counts = produce_scaled_counts_matrix(count_csv_file = count_csv_file,
                                               xp_design_csv_file = xp_design_csv_file)
  
  # based on a timepoint, filter the corresponding scaled_counts matrix
  filtered_counts <- filter_counts_based_on_timepoint(counts = scaled_counts, 
                                                      xp_design_csv_file = xp_design_csv_file, 
                                                      timep = tp) 
  
  # compute PCA and return scores as a dataframe with additional XP info
  # also returns explained variance per component 
   score_df <- produce_score_df(counts = filtered_counts, 
                               xp_design_csv_file = xp_design_csv_file)
  
  
  explained_variance_per_component <- extract_explained_variance_per_component(counts = filtered_counts) 
  
  # PCA plot
  variance_pc_x <- explained_variance_per_component[as.integer(gsub(pattern = "PC",replacement = "", x = pc_x_axis))]
  variance_pc_y <- explained_variance_per_component[as.integer(gsub(pattern = "PC",replacement = "", x = pc_y_axis))]
  
  
  pca_plot <- ggplot(data = score_df) +
    geom_point(aes_string(x = pc_x_axis, 
                          y = pc_y_axis, 
                          col = pca_colour,
                          shape = pca_shape), 
               size = 3) +
    xlab(paste0(pc_x_axis," ", variance_pc_x, "%")) +
    ylab(paste0(pc_y_axis," ", variance_pc_y, "%")) +
    scale_color_manual(values=c('#0099FF','#00b348', '#ff6700'))

  return(pca_plot)
}


