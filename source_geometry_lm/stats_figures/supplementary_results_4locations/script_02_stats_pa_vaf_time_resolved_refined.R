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
      
      # length 3
      length_seq = 3
      p_values_l3 = c()
      t_values_l3 = c()
      between_avg_l3 = c()
      between_sd_l3 = c()
      within_avg_l3 = c()
      within_sd_l3 = c()
      
      # stim
      filename = paste(path, "/", fun, "/", perf,  "/design_matrix/between_within_stim_resolved_refined_length", length_seq, "_800to1100.txt", sep = "")
      design_matrix = read.table(filename, sep = ",", header = TRUE)
      colnames(design_matrix) = c("response", "Subspaces", "subject", "session")
      design_matrix$Subspaces = factor(ifelse(design_matrix$Subspaces == 1, "Between ranks", "Within rank"))
      model = summary(lm(response ~ Subspaces, data = design_matrix))
      p_values_l3 = c(p_values_l3, model$coefficients[2,4])
      t_values_l3 = c(t_values_l3, model$coefficients[2,3])

      filename = paste(path, "/", fun, "/", perf,  "/design_matrix/between_within_stim_resolved_refined_length", length_seq, "_1200to1500.txt", sep = "")
      design_matrix = read.table(filename, sep = ",", header = TRUE)
      colnames(design_matrix) = c("response", "Subspaces", "subject", "session")
      design_matrix$Subspaces = factor(ifelse(design_matrix$Subspaces == 1, "Between ranks", "Within rank"))
      model = summary(lm(response ~ Subspaces, data = design_matrix))
      p_values_l3 = c(p_values_l3, model$coefficients[2,4])
      t_values_l3 = c(t_values_l3, model$coefficients[2,3])

      # delay
      
      for (time_i in 1:nrow(delay_segments)) {
        time_window = paste(delay_segments[time_i,1], "to", delay_segments[time_i,2], sep = "")
        filename = paste(path, "/", fun, "/", perf,  "/design_matrix/between_within_delay_length", length_seq, "_", time_window, ".txt", sep = "")
        design_matrix = read.table(filename, sep = ",", header = TRUE)
        colnames(design_matrix) = c("response", "Subspaces", "subject", "session")
        design_matrix$Subspaces = factor(ifelse(design_matrix$Subspaces == 1, "Between ranks", "Within rank"))
        model = summary(lm(response ~ Subspaces, data = design_matrix))
        p_values_l3 = c(p_values_l3, model$coefficients[2,4])
        t_values_l3 = c(t_values_l3, model$coefficients[2,3])
        data_l3 = subset(design_matrix, Subspaces == "Between ranks")
        between_avg_l3 = c(between_avg_l3, mean(data_l3$response))
        between_sd_l3 = c(between_sd_l3, sd(data_l3$response))
        data_l3 = subset(design_matrix, Subspaces == "Within rank")
        within_avg_l3 = c(within_avg_l3, mean(data_l3$response))
        within_sd_l3 = c(within_sd_l3, sd(data_l3$response))
      }
      
      # length 4
      length_seq = 4
      p_values_l4 = c()
      t_values_l4 = c()
      between_avg_l4 = c()
      between_sd_l4 = c()
      within_avg_l4 = c()
      within_sd_l4 = c()

      
      # stim
      filename = paste(path, "/", fun, "/", perf,  "/design_matrix/between_within_stim_resolved_refined_length", length_seq, "_400to700.txt", sep = "")
      design_matrix = read.table(filename, sep = ",", header = TRUE)
      colnames(design_matrix) = c("response", "Subspaces", "subject", "session")
      design_matrix$Subspaces = factor(ifelse(design_matrix$Subspaces == 1, "Between ranks", "Within rank"))
      model = summary(lm(response ~ Subspaces, data = design_matrix))
      p_values_l4 = c(p_values_l4, model$coefficients[2,4])
      t_values_l4 = c(t_values_l4, model$coefficients[2,3])
      
      filename = paste(path, "/", fun, "/", perf,  "/design_matrix/between_within_stim_resolved_refined_length", length_seq, "_800to1100.txt", sep = "")
      design_matrix = read.table(filename, sep = ",", header = TRUE)
      colnames(design_matrix) = c("response", "Subspaces", "subject", "session")
      design_matrix$Subspaces = factor(ifelse(design_matrix$Subspaces == 1, "Between ranks", "Within rank"))
      model = summary(lm(response ~ Subspaces, data = design_matrix))
      p_values_l4 = c(p_values_l4, model$coefficients[2,4])
      t_values_l4 = c(t_values_l4, model$coefficients[2,3])
      
      filename = paste(path, "/", fun, "/", perf,  "/design_matrix/between_within_stim_resolved_refined_length", length_seq, "_1200to1500.txt", sep = "")
      design_matrix = read.table(filename, sep = ",", header = TRUE)
      colnames(design_matrix) = c("response", "Subspaces", "subject", "session")
      design_matrix$Subspaces = factor(ifelse(design_matrix$Subspaces == 1, "Between ranks", "Within rank"))
      model = summary(lm(response ~ Subspaces, data = design_matrix))
      p_values_l4 = c(p_values_l4, model$coefficients[2,4])
      t_values_l4 = c(t_values_l4, model$coefficients[2,3])
      
      for (time_i in 1:nrow(delay_segments)) {
        time_window = paste(delay_segments[time_i,1], "to", delay_segments[time_i,2], sep = "")
        filename = paste(path, "/", fun, "/", perf,  "/design_matrix/between_within_delay_length", length_seq, "_", time_window, ".txt", sep = "")
        design_matrix = read.table(filename, sep = ",", header = TRUE)
        colnames(design_matrix) = c("response", "Subspaces", "subject", "session")
        design_matrix$Subspaces = factor(ifelse(design_matrix$Subspaces == 1, "Between ranks", "Within rank"))
        model = summary(lm(response ~ Subspaces, data = design_matrix))
        p_values_l4 = c(p_values_l4, model$coefficients[2,4])
        t_values_l4 = c(t_values_l4, model$coefficients[2,3])
        data_l4 = subset(design_matrix, Subspaces == "Between ranks")
        between_avg_l4 = c(between_avg_l4, mean(data_l4$response))
        between_sd_l4 = c(between_sd_l4, sd(data_l4$response))
        data_l4 = subset(design_matrix, Subspaces == "Within rank")
        within_avg_l4 = c(within_avg_l4, mean(data_l4$response))
        within_sd_l4 = c(within_sd_l4, sd(data_l4$response))
      }
      
      
    }
  }
}

p_values = c(p_values_l3, p_values_l4)
t_values = c(t_values_l3, t_values_l4)
