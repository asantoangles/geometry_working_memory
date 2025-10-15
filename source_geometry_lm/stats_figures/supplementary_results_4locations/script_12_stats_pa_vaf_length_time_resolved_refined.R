##############################
## connected scatter plots  ##
##############################

segments_start <- seq(from = 1, to = 3601, by = 300)
segments_end <- seq(from = 300, to = 3900, by = 300)
delay_segments = cbind(segments_start, segments_end)

for (metric_i in metrics) {
  
  path = paste(path_results, "/group_results/", metric_i, sep = "")
  
  for (fun in functions_data) {
    for (perf in performance) {
      
      ################################
      ## principal angle - between  ##
      ################################
      
      p_values = c()
      t_values = c()
      
      # stim
      time_window = paste(delay_segments[time_i,1], "to", delay_segments[time_i,2], sep = "")
      filename = paste(path, "/", fun, "/", perf,  "/design_matrix/between_within_stim_resolved_refined_length3_800to1100.txt", sep = "")
      design_matrix_l3 = read.table(filename, sep = ",", header = TRUE)
      colnames(design_matrix_l3) = c("response", "Subspaces", "subject", "session")
      design_matrix_l3$length = 3
      filename = paste(path, "/", fun, "/", perf,  "/design_matrix/between_within_stim_resolved_refined_length4_400to700.txt", sep = "")
      design_matrix_l4 = read.table(filename, sep = ",", header = TRUE)
      colnames(design_matrix_l4) = c("response", "Subspaces", "subject", "session")
      design_matrix_l4$length = 4
      design_matrix = rbind(design_matrix_l3, design_matrix_l4)
      design_matrix$Subspaces = factor(ifelse(design_matrix$Subspaces == 1, "Between ranks", "Within rank"))
      design_matrix$length = factor(design_matrix$length)
      design_matrix = subset(design_matrix, Subspaces == "Between ranks")
      model = summary(lm(response ~ length, data = design_matrix))
      p_values = c(p_values, model$coefficients[2,4])
      t_values = c(t_values, model$coefficients[2,3])
      
      
      time_window = paste(delay_segments[time_i,1], "to", delay_segments[time_i,2], sep = "")
      filename = paste(path, "/", fun, "/", perf,  "/design_matrix/between_within_stim_resolved_refined_length3_1200to1500.txt", sep = "")
      design_matrix_l3 = read.table(filename, sep = ",", header = TRUE)
      colnames(design_matrix_l3) = c("response", "Subspaces", "subject", "session")
      design_matrix_l3$length = 3
      filename = paste(path, "/", fun, "/", perf,  "/design_matrix/between_within_stim_resolved_refined_length4_800to1100.txt", sep = "")
      design_matrix_l4 = read.table(filename, sep = ",", header = TRUE)
      colnames(design_matrix_l4) = c("response", "Subspaces", "subject", "session")
      design_matrix_l4$length = 4
      design_matrix = rbind(design_matrix_l3, design_matrix_l4)
      design_matrix$Subspaces = factor(ifelse(design_matrix$Subspaces == 1, "Between ranks", "Within rank"))
      design_matrix$length = factor(design_matrix$length)
      design_matrix = subset(design_matrix, Subspaces == "Between ranks")
      model = summary(lm(response ~ length, data = design_matrix))
      p_values = c(p_values, model$coefficients[2,4])
      t_values = c(t_values, model$coefficients[2,3])
      
      # delay
      
      for (time_i in 1:nrow(delay_segments)) {
        time_window = paste(delay_segments[time_i,1], "to", delay_segments[time_i,2], sep = "")
        filename = paste(path, "/", fun, "/", perf,  "/design_matrix/between_within_delay_length3_", time_window, ".txt", sep = "")
        design_matrix_l3 = read.table(filename, sep = ",", header = TRUE)
        colnames(design_matrix_l3) = c("response", "Subspaces", "subject", "session")
        design_matrix_l3$length = 3
        filename = paste(path, "/", fun, "/", perf,  "/design_matrix/between_within_delay_length4_", time_window, ".txt", sep = "")
        design_matrix_l4 = read.table(filename, sep = ",", header = TRUE)
        colnames(design_matrix_l4) = c("response", "Subspaces", "subject", "session")
        design_matrix_l4$length = 4
        design_matrix = rbind(design_matrix_l3, design_matrix_l4)
        design_matrix$Subspaces = factor(ifelse(design_matrix$Subspaces == 1, "Between ranks", "Within rank"))
        design_matrix$length = factor(design_matrix$length)
        design_matrix = subset(design_matrix, Subspaces == "Between ranks")
        model = summary(lm(response ~ length, data = design_matrix))
        p_values = c(p_values, model$coefficients[2,4])
        t_values = c(t_values, model$coefficients[2,3])
      }
      
    }
  }
}
