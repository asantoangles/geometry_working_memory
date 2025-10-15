#################################
## separability - performance  ##
#################################

segments_start <- seq(from = 1, to = 3601, by = 300)
segments_end <- seq(from = 300, to = 3900, by = 300)
delay_segments = cbind(segments_start, segments_end)

for (separability_metric in separability_metrics) {
  
  path = paste(path_results, "/group_results/", separability_metric, sep = "")
  path_random = paste(path_results_random, "/group_results/", separability_metric, sep = "")
  
  for (fun in functions_data) {
    
    ############################
    ### correct vs incorrect ###
    ############################

    ## length 3
    length_seq = 3
    p_values_l3 = c()
    t_values_l3 = c()

    # stim
    
    filename = paste(path, "/", fun, "/correct_trials/design_matrix/design_matrix_stim_resolved_refined_length", length_seq, "_400to700.txt", sep = "")
    design_matrix_correct = read.table(filename, sep = ",", header = TRUE)
    design_matrix_correct$performance = "correct"
    filename = paste(path, "/", fun, "/incorrect_trials/design_matrix/design_matrix_stim_resolved_refined_length", length_seq, "_400to700.txt", sep = "")
    design_matrix_incorrect = read.table(filename, sep = ",", header = TRUE)
    design_matrix_incorrect$performance = "incorrect"
    design_matrix = rbind(design_matrix_correct, design_matrix_incorrect)
    design_matrix$performance = factor(design_matrix$performance, levels = c("incorrect", "correct"))
    model = summary(lm(response ~ performance, data = design_matrix))
    p_values_l3 = c(p_values_l3, model$coefficients[2,4])
    t_values_l3 = c(t_values_l3, model$coefficients[2,3])
    
    filename = paste(path, "/", fun, "/correct_trials/design_matrix/design_matrix_stim_resolved_refined_length", length_seq, "_800to1100.txt", sep = "")
    design_matrix_correct = read.table(filename, sep = ",", header = TRUE)
    design_matrix_correct$performance = "correct"
    filename = paste(path, "/", fun, "/incorrect_trials/design_matrix/design_matrix_stim_resolved_refined_length", length_seq, "_800to1100.txt", sep = "")
    design_matrix_incorrect = read.table(filename, sep = ",", header = TRUE)
    design_matrix_incorrect$performance = "incorrect"
    design_matrix = rbind(design_matrix_correct, design_matrix_incorrect)
    design_matrix$performance = factor(design_matrix$performance, levels = c("incorrect", "correct"))
    model = summary(lm(response ~ performance, data = design_matrix))
    p_values_l3 = c(p_values_l3, model$coefficients[2,4])
    t_values_l3 = c(t_values_l3, model$coefficients[2,3])
    
    filename = paste(path, "/", fun, "/correct_trials/design_matrix/design_matrix_stim_resolved_refined_length", length_seq, "_1200to1500.txt", sep = "")
    design_matrix_correct = read.table(filename, sep = ",", header = TRUE)
    design_matrix_correct$performance = "correct"
    filename = paste(path, "/", fun, "/incorrect_trials/design_matrix/design_matrix_stim_resolved_refined_length", length_seq, "_1200to1500.txt", sep = "")
    design_matrix_incorrect = read.table(filename, sep = ",", header = TRUE)
    design_matrix_incorrect$performance = "incorrect"
    design_matrix = rbind(design_matrix_correct, design_matrix_incorrect)
    design_matrix$performance = factor(design_matrix$performance, levels = c("incorrect", "correct"))
    model = summary(lm(response ~ performance, data = design_matrix))
    p_values_l3 = c(p_values_l3, model$coefficients[2,4])
    t_values_l3 = c(t_values_l3, model$coefficients[2,3])
    
    
    # delay
    
    for (time_i in 1:nrow(delay_segments)) {
      time_window = paste(delay_segments[time_i,1], "to", delay_segments[time_i,2], sep = "")
      
      filename = paste(path, "/", fun, "/correct_trials/design_matrix/design_matrix_delay_length", length_seq, "_", time_window, ".txt", sep = "")
      design_matrix_correct = read.table(filename, sep = ",", header = TRUE)
      design_matrix_correct$performance = "correct"

      filename = paste(path, "/", fun, "/incorrect_trials/design_matrix/design_matrix_delay_length", length_seq, "_", time_window, ".txt", sep = "")
      design_matrix_incorrect = read.table(filename, sep = ",", header = TRUE)
      design_matrix_incorrect$performance = "incorrect"
      
      design_matrix = rbind(design_matrix_correct, design_matrix_incorrect)
      
      design_matrix$performance = factor(design_matrix$performance, levels = c("incorrect", "correct"))
      
      model = summary(lm(response ~ performance, data = design_matrix))
      p_values_l3 = c(p_values_l3, model$coefficients[2,4])
      t_values_l3 = c(t_values_l3, model$coefficients[2,3])
      
    }
    
    ## length 4
    length_seq = 4
    p_values_l4 = c()
    t_values_l4 = c()
    
    # stim
    
    filename = paste(path, "/", fun, "/correct_trials/design_matrix/design_matrix_stim_resolved_refined_length", length_seq, "_1to300.txt", sep = "")
    design_matrix_correct = read.table(filename, sep = ",", header = TRUE)
    design_matrix_correct$performance = "correct"
    filename = paste(path, "/", fun, "/incorrect_trials/design_matrix/design_matrix_stim_resolved_refined_length", length_seq, "_1to300.txt", sep = "")
    design_matrix_incorrect = read.table(filename, sep = ",", header = TRUE)
    design_matrix_incorrect$performance = "incorrect"
    design_matrix = rbind(design_matrix_correct, design_matrix_incorrect)
    design_matrix$performance = factor(design_matrix$performance, levels = c("incorrect", "correct"))
    model = summary(lm(response ~ performance, data = design_matrix))
    p_values_l4 = c(p_values_l4, model$coefficients[2,4])
    t_values_l4 = c(t_values_l4, model$coefficients[2,3])
    
    filename = paste(path, "/", fun, "/correct_trials/design_matrix/design_matrix_stim_resolved_refined_length", length_seq, "_400to700.txt", sep = "")
    design_matrix_correct = read.table(filename, sep = ",", header = TRUE)
    design_matrix_correct$performance = "correct"
    filename = paste(path, "/", fun, "/incorrect_trials/design_matrix/design_matrix_stim_resolved_refined_length", length_seq, "_400to700.txt", sep = "")
    design_matrix_incorrect = read.table(filename, sep = ",", header = TRUE)
    design_matrix_incorrect$performance = "incorrect"
    design_matrix = rbind(design_matrix_correct, design_matrix_incorrect)
    design_matrix$performance = factor(design_matrix$performance, levels = c("incorrect", "correct"))
    model = summary(lm(response ~ performance, data = design_matrix))
    p_values_l4 = c(p_values_l4, model$coefficients[2,4])
    t_values_l4 = c(t_values_l4, model$coefficients[2,3])
    
    filename = paste(path, "/", fun, "/correct_trials/design_matrix/design_matrix_stim_resolved_refined_length", length_seq, "_800to1100.txt", sep = "")
    design_matrix_correct = read.table(filename, sep = ",", header = TRUE)
    design_matrix_correct$performance = "correct"
    filename = paste(path, "/", fun, "/incorrect_trials/design_matrix/design_matrix_stim_resolved_refined_length", length_seq, "_800to1100.txt", sep = "")
    design_matrix_incorrect = read.table(filename, sep = ",", header = TRUE)
    design_matrix_incorrect$performance = "incorrect"
    design_matrix = rbind(design_matrix_correct, design_matrix_incorrect)
    design_matrix$performance = factor(design_matrix$performance, levels = c("incorrect", "correct"))
    model = summary(lm(response ~ performance, data = design_matrix))
    p_values_l4 = c(p_values_l4, model$coefficients[2,4])
    t_values_l4 = c(t_values_l4, model$coefficients[2,3])
    
    filename = paste(path, "/", fun, "/correct_trials/design_matrix/design_matrix_stim_resolved_refined_length", length_seq, "_1200to1500.txt", sep = "")
    design_matrix_correct = read.table(filename, sep = ",", header = TRUE)
    design_matrix_correct$performance = "correct"
    filename = paste(path, "/", fun, "/incorrect_trials/design_matrix/design_matrix_stim_resolved_refined_length", length_seq, "_1200to1500.txt", sep = "")
    design_matrix_incorrect = read.table(filename, sep = ",", header = TRUE)
    design_matrix_incorrect$performance = "incorrect"
    design_matrix = rbind(design_matrix_correct, design_matrix_incorrect)
    design_matrix$performance = factor(design_matrix$performance, levels = c("incorrect", "correct"))
    model = summary(lm(response ~ performance, data = design_matrix))
    p_values_l4 = c(p_values_l4, model$coefficients[2,4])
    t_values_l4 = c(t_values_l4, model$coefficients[2,3])
    
    # delay
    
    for (time_i in 1:nrow(delay_segments)) {
      time_window = paste(delay_segments[time_i,1], "to", delay_segments[time_i,2], sep = "")
      
      filename = paste(path, "/", fun, "/correct_trials/design_matrix/design_matrix_delay_length", length_seq, "_", time_window, ".txt", sep = "")
      design_matrix_correct = read.table(filename, sep = ",", header = TRUE)
      design_matrix_correct$performance = "correct"
      
      filename = paste(path, "/", fun, "/incorrect_trials/design_matrix/design_matrix_delay_length", length_seq, "_", time_window, ".txt", sep = "")
      design_matrix_incorrect = read.table(filename, sep = ",", header = TRUE)
      design_matrix_incorrect$performance = "incorrect"
      
      design_matrix = rbind(design_matrix_correct, design_matrix_incorrect)
      
      design_matrix$performance = factor(design_matrix$performance, levels = c("incorrect", "correct"))
      
      model = summary(lm(response ~ performance, data = design_matrix))
      p_values_l4 = c(p_values_l4, model$coefficients[2,4])
      t_values_l4 = c(t_values_l4, model$coefficients[2,3])
      
    }
    
    p_values_l3_empirical = p_values_l3
    t_values_l3_empirical = t_values_l3
    p_values_l4_empirical = p_values_l4
    t_values_l4_empirical = t_values_l4
    
    
    
    
    
    

    
  }
}

p_values = c(p_values_l3_empirical, p_values_l4_empirical)
t_values = c(t_values_l3_empirical, t_values_l4_empirical)
