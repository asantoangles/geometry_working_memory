##############################
## connected scatter plots  ##
##############################

segments_start <- seq(from = 1, to = 3601, by = 300)
segments_end <- seq(from = 300, to = 3900, by = 300)
delay_segments = cbind(segments_start, segments_end)

for (metric_i in metrics) {
  
  path = paste(path_results, "/group_results/", metric_i, sep = "")
  path_random = paste(path_results_random, "/group_results/", metric_i, sep = "")
  
  for (fun in functions_data) {
    for (perf in performance) {
      
      ##########################################
      ## principal angle - between vs within  ##
      ##########################################
      
      # length 3
      length_seq = 3
      p_values_l3 = c()
      t_values_l3 = c()

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
      }
      
      # length 4
      length_seq = 4
      p_values_l4 = c()
      t_values_l4 = c()
      
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
      }
      
      p_values_l3_empirical = p_values_l3
      p_values_l4_empirical = p_values_l4
      t_values_l3_empirical = t_values_l3
      t_values_l4_empirical = t_values_l4
      
      # set to one the p-value when PA (VAF) within is greater (smaller) than between
      if (metric_i == "vaf") {
        p_values_l3_empirical = ifelse(t_values_l3_empirical > 0, p_values_l3_empirical, 1)
        p_values_l4_empirical = ifelse(t_values_l4_empirical > 0, p_values_l4_empirical, 1)
      } else {
        p_values_l3_empirical = ifelse(t_values_l3_empirical < 0, p_values_l3_empirical, 1)
        p_values_l4_empirical = ifelse(t_values_l4_empirical < 0, p_values_l4_empirical, 1)
      }
      
      ##########################################
      ## principal angle - between vs random  ##
      ##########################################
      
      # length 3
      length_seq = 3
      p_values_l3 = c()
      t_values_l3 = c()

      # stim
      filename = paste(path, "/", fun, "/", perf,  "/design_matrix/between_within_stim_resolved_refined_length", length_seq, "_800to1100.txt", sep = "")
      design_matrix = read.table(filename, sep = ",", header = TRUE)
      colnames(design_matrix) = c("response", "Subspaces", "subject", "session")
      design_matrix$Subspaces = factor(ifelse(design_matrix$Subspaces == 1, "Between ranks", "Within rank"))
      design_matrix = subset(design_matrix, Subspaces == "Between ranks")
      filename = paste(path_random, "/", fun, "/", perf,  "/design_matrix/between_stim_resolved_refined_length", length_seq, "_800to1100.txt", sep = "")
      design_matrix_random = read.table(filename, sep = ",", header = TRUE)
      colnames(design_matrix_random) = c("response", "Subspaces", "subject", "session")
      design_matrix_random$Subspaces = factor(ifelse(design_matrix_random$Subspaces == 1, "Between ranks", "Within rank"))
      design_matrix_random = subset(design_matrix_random, Subspaces == "Between ranks")
      design_matrix_random$Subspaces = "Surrogate"
      design_matrix = rbind(design_matrix, design_matrix_random)
      model = summary(lm(response ~ Subspaces, data = design_matrix))
      p_values_l3 = c(p_values_l3, model$coefficients[2,4])
      t_values_l3 = c(t_values_l3, model$coefficients[2,3])
      
      filename = paste(path, "/", fun, "/", perf,  "/design_matrix/between_within_stim_resolved_refined_length", length_seq, "_1200to1500.txt", sep = "")
      design_matrix = read.table(filename, sep = ",", header = TRUE)
      colnames(design_matrix) = c("response", "Subspaces", "subject", "session")
      design_matrix$Subspaces = factor(ifelse(design_matrix$Subspaces == 1, "Between ranks", "Within rank"))
      design_matrix = subset(design_matrix, Subspaces == "Between ranks")
      filename = paste(path_random, "/", fun, "/", perf,  "/design_matrix/between_stim_resolved_refined_length", length_seq, "_1200to1500.txt", sep = "")
      design_matrix_random = read.table(filename, sep = ",", header = TRUE)
      colnames(design_matrix_random) = c("response", "Subspaces", "subject", "session")
      design_matrix_random$Subspaces = factor(ifelse(design_matrix_random$Subspaces == 1, "Between ranks", "Within rank"))
      design_matrix_random = subset(design_matrix_random, Subspaces == "Between ranks")
      design_matrix_random$Subspaces = "Surrogate"
      design_matrix = rbind(design_matrix, design_matrix_random)
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
        design_matrix = subset(design_matrix, Subspaces == "Between ranks")
        filename = paste(path_random, "/", fun, "/", perf,  "/design_matrix/between_delay_length", length_seq, "_", time_window, ".txt", sep = "")
        design_matrix_random = read.table(filename, sep = ",", header = TRUE)
        colnames(design_matrix_random) = c("response", "Subspaces", "subject", "session")
        design_matrix_random$Subspaces = factor(ifelse(design_matrix_random$Subspaces == 1, "Between ranks", "Within rank"))
        design_matrix_random = subset(design_matrix_random, Subspaces == "Between ranks")
        design_matrix_random$Subspaces = "Surrogate"
        design_matrix = rbind(design_matrix, design_matrix_random)
        model = summary(lm(response ~ Subspaces, data = design_matrix))
        p_values_l3 = c(p_values_l3, model$coefficients[2,4])
        t_values_l3 = c(t_values_l3, model$coefficients[2,3])
      }
      
      # length 4
      length_seq = 4
      p_values_l4 = c()
      t_values_l4 = c()

      
      # stim
      filename = paste(path, "/", fun, "/", perf,  "/design_matrix/between_within_stim_resolved_refined_length", length_seq, "_400to700.txt", sep = "")
      design_matrix = read.table(filename, sep = ",", header = TRUE)
      colnames(design_matrix) = c("response", "Subspaces", "subject", "session")
      design_matrix$Subspaces = factor(ifelse(design_matrix$Subspaces == 1, "Between ranks", "Within rank"))
      design_matrix = subset(design_matrix, Subspaces == "Between ranks")
      filename = paste(path_random, "/", fun, "/", perf,  "/design_matrix/between_stim_resolved_refined_length", length_seq, "_400to700.txt", sep = "")
      design_matrix_random = read.table(filename, sep = ",", header = TRUE)
      colnames(design_matrix_random) = c("response", "Subspaces", "subject", "session")
      design_matrix_random$Subspaces = factor(ifelse(design_matrix_random$Subspaces == 1, "Between ranks", "Within rank"))
      design_matrix_random = subset(design_matrix_random, Subspaces == "Between ranks")
      design_matrix_random$Subspaces = "Surrogate"
      design_matrix = rbind(design_matrix, design_matrix_random)
      model = summary(lm(response ~ Subspaces, data = design_matrix))
      p_values_l4 = c(p_values_l4, model$coefficients[2,4])
      t_values_l4 = c(t_values_l4, model$coefficients[2,3])
      
      filename = paste(path, "/", fun, "/", perf,  "/design_matrix/between_within_stim_resolved_refined_length", length_seq, "_800to1100.txt", sep = "")
      design_matrix = read.table(filename, sep = ",", header = TRUE)
      colnames(design_matrix) = c("response", "Subspaces", "subject", "session")
      design_matrix$Subspaces = factor(ifelse(design_matrix$Subspaces == 1, "Between ranks", "Within rank"))
      design_matrix = subset(design_matrix, Subspaces == "Between ranks")
      filename = paste(path_random, "/", fun, "/", perf,  "/design_matrix/between_stim_resolved_refined_length", length_seq, "_800to1100.txt", sep = "")
      design_matrix_random = read.table(filename, sep = ",", header = TRUE)
      colnames(design_matrix_random) = c("response", "Subspaces", "subject", "session")
      design_matrix_random$Subspaces = factor(ifelse(design_matrix_random$Subspaces == 1, "Between ranks", "Within rank"))
      design_matrix_random = subset(design_matrix_random, Subspaces == "Between ranks")
      design_matrix_random$Subspaces = "Surrogate"
      design_matrix = rbind(design_matrix, design_matrix_random)
      model = summary(lm(response ~ Subspaces, data = design_matrix))
      p_values_l4 = c(p_values_l4, model$coefficients[2,4])
      t_values_l4 = c(t_values_l4, model$coefficients[2,3])
      
      filename = paste(path, "/", fun, "/", perf,  "/design_matrix/between_within_stim_resolved_refined_length", length_seq, "_1200to1500.txt", sep = "")
      design_matrix = read.table(filename, sep = ",", header = TRUE)
      colnames(design_matrix) = c("response", "Subspaces", "subject", "session")
      design_matrix$Subspaces = factor(ifelse(design_matrix$Subspaces == 1, "Between ranks", "Within rank"))
      design_matrix = subset(design_matrix, Subspaces == "Between ranks")
      filename = paste(path_random, "/", fun, "/", perf,  "/design_matrix/between_stim_resolved_refined_length", length_seq, "_1200to1500.txt", sep = "")
      design_matrix_random = read.table(filename, sep = ",", header = TRUE)
      colnames(design_matrix_random) = c("response", "Subspaces", "subject", "session")
      design_matrix_random$Subspaces = factor(ifelse(design_matrix_random$Subspaces == 1, "Between ranks", "Within rank"))
      design_matrix_random = subset(design_matrix_random, Subspaces == "Between ranks")
      design_matrix_random$Subspaces = "Surrogate"
      design_matrix = rbind(design_matrix, design_matrix_random)
      model = summary(lm(response ~ Subspaces, data = design_matrix))
      p_values_l4 = c(p_values_l4, model$coefficients[2,4])
      t_values_l4 = c(t_values_l4, model$coefficients[2,3])
      
            
      for (time_i in 1:nrow(delay_segments)) {
        time_window = paste(delay_segments[time_i,1], "to", delay_segments[time_i,2], sep = "")
        filename = paste(path, "/", fun, "/", perf,  "/design_matrix/between_within_delay_length", length_seq, "_", time_window, ".txt", sep = "")
        design_matrix = read.table(filename, sep = ",", header = TRUE)
        colnames(design_matrix) = c("response", "Subspaces", "subject", "session")
        design_matrix$Subspaces = factor(ifelse(design_matrix$Subspaces == 1, "Between ranks", "Within rank"))
        design_matrix = subset(design_matrix, Subspaces == "Between ranks")
        filename = paste(path_random, "/", fun, "/", perf,  "/design_matrix/between_delay_length", length_seq, "_", time_window, ".txt", sep = "")
        design_matrix_random = read.table(filename, sep = ",", header = TRUE)
        colnames(design_matrix_random) = c("response", "Subspaces", "subject", "session")
        design_matrix_random$Subspaces = factor(ifelse(design_matrix_random$Subspaces == 1, "Between ranks", "Within rank"))
        design_matrix_random = subset(design_matrix_random, Subspaces == "Between ranks")
        design_matrix_random$Subspaces = "Surrogate"
        design_matrix = rbind(design_matrix, design_matrix_random)
        model = summary(lm(response ~ Subspaces, data = design_matrix))
        p_values_l4 = c(p_values_l4, model$coefficients[2,4])
        t_values_l4 = c(t_values_l4, model$coefficients[2,3])
      }
      
      p_values_l3_random = p_values_l3
      p_values_l4_random = p_values_l4
      t_values_l3_random = t_values_l3
      t_values_l4_random = t_values_l4
      
      # set to one the p-value when PA (VAF) within is greater (smaller) than between
      if (metric_i == "vaf") {
        p_values_l3_random = ifelse(t_values_l3_random > 0, p_values_l3_random, 1)
        p_values_l4_random = ifelse(t_values_l4_random > 0, p_values_l4_random, 1)
      } else {
        p_values_l3_random = ifelse(t_values_l3_random < 0, p_values_l3_random, 1)
        p_values_l4_random = ifelse(t_values_l4_random < 0, p_values_l4_random, 1)
      }
      
      
    }
  }
}

p_values = c(p_values_l3_empirical, p_values_l4_empirical)
t_values = c(t_values_l3_empirical, t_values_l4_empirical)

p_values_random = c(p_values_l3_random, p_values_l4_random)
t_values_random = c(t_values_l3_random, t_values_l4_random)

