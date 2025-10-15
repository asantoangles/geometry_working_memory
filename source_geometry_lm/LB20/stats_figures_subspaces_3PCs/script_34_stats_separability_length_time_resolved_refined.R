
#########################################
## separability - length contrast      ##
#########################################

segments_start <- seq(from = 1, to = 3601, by = 300)
segments_end <- seq(from = 300, to = 3900, by = 300)
delay_segments = cbind(segments_start, segments_end)

for (separability in separability_metrics) {
  
  path = paste(path_results, "/group_results/", separability, sep = "")
  
  for (fun in functions_data) {
    
    p_values = c()
    t_values = c()
    
    # stim
    
    time_window = paste(delay_segments[time_i,1], "to", delay_segments[time_i,2], sep = "")
    filename = paste(path, "/", fun, "/correct_trials/design_matrix/design_matrix_stim_resolved_refined_length3_800to1100.txt", sep = "")
    design_matrix_l3 = read.table(filename, sep = ",", header = TRUE)
    design_matrix_l3$length = 3
    filename = paste(path, "/", fun, "/correct_trials/design_matrix/design_matrix_stim_resolved_refined_length4_400to700.txt", sep = "")
    design_matrix_l4 = read.table(filename, sep = ",", header = TRUE)
    design_matrix_l4$length = 4
    design_matrix = rbind(design_matrix_l3, design_matrix_l4)
    design_matrix$length = factor(design_matrix$length)
    model = summary(lm(response ~ length, data = design_matrix))
    p_values = c(p_values, model$coefficients[2,4])
    t_values = c(t_values, model$coefficients[2,3])
    
    time_window = paste(delay_segments[time_i,1], "to", delay_segments[time_i,2], sep = "")
    filename = paste(path, "/", fun, "/correct_trials/design_matrix/design_matrix_stim_resolved_refined_length3_1200to1500.txt", sep = "")
    design_matrix_l3 = read.table(filename, sep = ",", header = TRUE)
    design_matrix_l3$length = 3
    filename = paste(path, "/", fun, "/correct_trials/design_matrix/design_matrix_stim_resolved_refined_length4_800to1100.txt", sep = "")
    design_matrix_l4 = read.table(filename, sep = ",", header = TRUE)
    design_matrix_l4$length = 4
    design_matrix = rbind(design_matrix_l3, design_matrix_l4)
    design_matrix$length = factor(design_matrix$length)
    model = summary(lm(response ~ length, data = design_matrix))
    p_values = c(p_values, model$coefficients[2,4])
    t_values = c(t_values, model$coefficients[2,3])
    
    
    # delay
    
    for (time_i in 1:nrow(delay_segments)) {
      time_window = paste(delay_segments[time_i,1], "to", delay_segments[time_i,2], sep = "")
      filename = paste(path, "/", fun, "/correct_trials/design_matrix/design_matrix_delay_length3_", time_window, ".txt", sep = "")
      design_matrix_l3 = read.table(filename, sep = ",", header = TRUE)
      design_matrix_l3$length = 3
      filename = paste(path, "/", fun, "/correct_trials/design_matrix/design_matrix_delay_length4_", time_window, ".txt", sep = "")
      design_matrix_l4 = read.table(filename, sep = ",", header = TRUE)
      design_matrix_l4$length = 4
      design_matrix = rbind(design_matrix_l3, design_matrix_l4)
      design_matrix$length = factor(design_matrix$length)
      model = summary(lm(response ~ length, data = design_matrix))
      p_values = c(p_values, model$coefficients[2,4])
      t_values = c(t_values, model$coefficients[2,3])
    }
  }
}

