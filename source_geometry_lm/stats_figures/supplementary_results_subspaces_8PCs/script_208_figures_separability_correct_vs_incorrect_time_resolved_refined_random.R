##############################
## connected scatter plots  ##
##############################

segments_start <- seq(from = 1, to = 3601, by = 300)
segments_end <- seq(from = 300, to = 4000, by = 300)
delay_segments = cbind(segments_start, segments_end)

for (separability in separability_metrics) {
  
  path = paste(path_results, "/group_results/", separability, sep = "")

  cat(separability); cat(' ')
  
  if (separability == "distance") {
    p_values = p_values_all[1:33]
  }
  if (separability == "volume") {
    p_values = p_values_all[34:66]
  }
  
  for (fun in functions_data) {
    
    ## length 3
    
    seq_length = 3
    
    correct_avg = c()
    correct_sd = c()
    incorrect_avg = c()
    incorrect_sd = c()

    
    # stim
    
    filename = paste(path, "/", fun, "/correct_trials/design_matrix/design_matrix_stim_resolved_refined_length", seq_length, "_400to700.txt", sep = "")
    matrix_design = read.table(filename, sep = ",", header = TRUE)
    data = matrix_design
    colnames(data) = c("response", "condition", "subject", "session")
    correct_avg = c(correct_avg, mean(data$response))
    correct_sd = c(correct_sd, sd(data$response))
    filename = paste(path, "/", fun, "/incorrect_trials/design_matrix/design_matrix_stim_resolved_refined_length", seq_length, "_400to700.txt", sep = "")
    matrix_design = read.table(filename, sep = ",", header = TRUE)
    data = matrix_design
    colnames(data) = c("response", "condition", "subject", "session")
    incorrect_avg = c(incorrect_avg, mean(data$response))
    incorrect_sd = c(incorrect_sd, sd(data$response))
    
    filename = paste(path, "/", fun, "/correct_trials/design_matrix/design_matrix_stim_resolved_refined_length", seq_length, "_800to1100.txt", sep = "")
    matrix_design = read.table(filename, sep = ",", header = TRUE)
    data = matrix_design
    colnames(data) = c("response", "condition", "subject", "session")
    correct_avg = c(correct_avg, mean(data$response))
    correct_sd = c(correct_sd, sd(data$response))
    filename = paste(path, "/", fun, "/incorrect_trials/design_matrix/design_matrix_stim_resolved_refined_length", seq_length, "_800to1100.txt", sep = "")
    matrix_design = read.table(filename, sep = ",", header = TRUE)
    data = matrix_design
    colnames(data) = c("response", "condition", "subject", "session")
    incorrect_avg = c(incorrect_avg, mean(data$response))
    incorrect_sd = c(incorrect_sd, sd(data$response))
    
    filename = paste(path, "/", fun, "/correct_trials/design_matrix/design_matrix_stim_resolved_refined_length", seq_length, "_1200to1500.txt", sep = "")
    matrix_design = read.table(filename, sep = ",", header = TRUE)
    data = matrix_design
    colnames(data) = c("response", "condition", "subject", "session")
    correct_avg = c(correct_avg, mean(data$response))
    correct_sd = c(correct_sd, sd(data$response))
    filename = paste(path, "/", fun, "/incorrect_trials/design_matrix/design_matrix_stim_resolved_refined_length", seq_length, "_1200to1500.txt", sep = "")
    matrix_design = read.table(filename, sep = ",", header = TRUE)
    data = matrix_design
    colnames(data) = c("response", "condition", "subject", "session")
    incorrect_avg = c(incorrect_avg, mean(data$response))
    incorrect_sd = c(incorrect_sd, sd(data$response))
    
    # delay
    
    for (time_i in 1:nrow(delay_segments)) {
      time_window = paste(delay_segments[time_i,1], "to", delay_segments[time_i,2], sep = "")
      
      # correct_trials
      filename = paste(path, "/", fun, "/correct_trials/design_matrix/design_matrix_delay_length", seq_length, "_", time_window, ".txt", sep = "")
      matrix_design = read.table(filename, sep = ",", header = TRUE)
      data = matrix_design
      colnames(data) = c("response", "condition", "subject", "session")
      correct_avg = c(correct_avg, mean(data$response))
      correct_sd = c(correct_sd, sd(data$response))
      
      # incorrect_trials
      filename = paste(path, "/", fun, "/incorrect_trials/design_matrix/design_matrix_delay_length", seq_length, "_", time_window, ".txt", sep = "")
      matrix_design = read.table(filename, sep = ",", header = TRUE)
      data = matrix_design
      colnames(data) = c("response", "condition", "subject", "session")
      incorrect_avg = c(incorrect_avg, mean(data$response))
      incorrect_sd = c(incorrect_sd, sd(data$response))
      
    }
    

    # data
    time = c(-750, -450, -150, delay_segments[,2] - 150) / 1000
    correct_se = correct_sd / sqrt(15)
    incorrect_se = incorrect_sd / sqrt(15)
    avg_column = c(correct_avg, incorrect_avg)
    sd_column = c(correct_sd, incorrect_sd)
    se_column = c(correct_se, incorrect_se)
    performance = c(rep("correct trials", length(correct_avg)), 
                    rep("incorrect trials", length(correct_avg)))
    data = as.data.frame(cbind(avg_column, sd_column, se_column, time, performance))
    colnames(data) = c("avg", "sd", "se", "time", "performance")
    data$performance = factor(data$performance)
    data$avg = as.numeric(data$avg)
    data$sd = as.numeric(data$sd)
    data$se = as.numeric(data$se)
    data$time = as.numeric(data$time)
    
    # data p-values
    
    # Find segments with consecutive values below 0.05
    p_values_tmp = p_values[1:16]
    below_threshold <- p_values_tmp < 0.05
    start_points <- which(diff(c(FALSE, below_threshold)) == 1)
    end_points <- which(diff(c(below_threshold, FALSE)) == -1)
    if (separability == "distance") {
      y_pvalue = max(data$avg + data$se) + (max(data$avg + data$se)/200)
    } else {
      y_pvalue = max(data$avg + data$se) + (max(data$avg + data$se)/50)
    }
    
    # Create segments data frame
    segments <- data.frame(
      x = data$time[start_points],
      xend = data$time[end_points],
      y = rep(y_pvalue, length(start_points)),  # You can adjust the y-value as needed
      yend = rep(y_pvalue, length(start_points))  # Same y-value for end points
    )
    segments$x = segments$x - 0.05
    segments$xend = segments$xend + 0.05
    

    
    
    
    
    
    plot = 
      ggplot() +
      geom_vline(xintercept = 0, linetype = "dashed", color = "black", size = 0.3) +
      geom_segment(data = segments, aes(x = x, xend = xend, y = y, yend = yend), color = "black", size = 1) +
      geom_line(data = subset(data, performance == "incorrect trials"), 
                aes(x = time, y = avg, color = "incorrect trials"), size = 1) + 
      geom_line(data = subset(data, performance == "correct trials"), 
                aes(x = time, y = avg, color = "correct trials"), size = 1) + 
      geom_ribbon(data = subset(data, performance == "incorrect trials"), 
                  aes(x = time, ymin = avg - se, ymax = avg + se), alpha = 0.2, fill = "#00BFC4") + 
      geom_ribbon(data = subset(data, performance == "correct trials"), 
                  aes(x = time, ymin = avg - se, ymax = avg + se), alpha = 0.2, fill = "#F8766D") + 
      scale_color_manual(
        values = c("correct trials" = "#F8766D", "incorrect trials" = "#00BFC4"),
        labels = c("correct trials" = "Correct", "incorrect trials" = "Incorrect")) +
      labs(color = " ") +
      xlab("Time (s)") + ylab(separability) +
      scale_x_continuous(breaks = seq(0, 4, by = 1), limits = c(min(time) - 0.05, 4)) +
      theme(plot.title = element_text(hjust = 0.5), 
            legend.position = "none", 
            panel.background = element_blank(),
            plot.background = element_blank(), 
            axis.line = element_line(color = "black"), 
            axis.title = element_text(color = "black", size = 12), 
            axis.text.x = element_text(color = "black", size = 12), 
            axis.text.y = element_text(color = "black", size = 12), 
            panel.grid.major = element_blank(), 
            panel.grid.minor = element_blank())
    
    ggsave(paste(path_figures, "/connected_scatterplot_separability_", separability, "_length", seq_length, "_correct_vs_incorrect.png", sep = ""), plot, width = 3, height = 3, units = "in")
    
    
    plot = 
      ggplot() +
      geom_vline(xintercept = 0, linetype = "dashed", color = "black", size = 0.3) +
      geom_segment(data = segments, aes(x = x, xend = xend, y = y, yend = yend), color = "black", size = 1) +
      geom_line(data = subset(data, performance == "incorrect trials"), 
                aes(x = time, y = avg, color = "incorrect trials"), size = 1) + 
      geom_line(data = subset(data, performance == "correct trials"), 
                aes(x = time, y = avg, color = "correct trials"), size = 1) + 
      geom_ribbon(data = subset(data, performance == "incorrect trials"), 
                  aes(x = time, ymin = avg - se, ymax = avg + se), alpha = 0.2, fill = "#00BFC4") + 
      geom_ribbon(data = subset(data, performance == "correct trials"), 
                  aes(x = time, ymin = avg - se, ymax = avg + se), alpha = 0.2, fill = "#F8766D") + 
      scale_color_manual(
        values = c("correct trials" = "#F8766D", "incorrect trials" = "#00BFC4"),
        labels = c("correct trials" = "Correct", "incorrect trials" = "Incorrect")) +
      labs(color = " ") +
      xlab("Time (s)") + ylab(separability) +
      scale_x_continuous(breaks = seq(0, 4, by = 1), limits = c(min(time) - 0.05, 4)) +
      theme(plot.title = element_text(hjust = 0.5), 
            legend.position = "bottom", 
            panel.background = element_blank(),
            plot.background = element_blank(), 
            axis.line = element_line(color = "black"), 
            axis.title = element_text(color = "black", size = 12), 
            axis.text.x = element_text(color = "black", size = 12), 
            axis.text.y = element_text(color = "black", size = 12), 
            panel.grid.major = element_blank(), 
            panel.grid.minor = element_blank())
    
    ggsave(paste(path_figures, "/connected_scatterplot_separability_", separability, "_length", seq_length, "_correct_vs_incorrect_legend.png", sep = ""), plot, width = 3.5, height = 3.5, units = "in")
    
    
    
    
    
    
    
    
    ## length 4
    
    seq_length = 4
    
    correct_avg = c()
    correct_sd = c()
    incorrect_avg = c()
    incorrect_sd = c()

    
    # stim
    
    filename = paste(path, "/", fun, "/correct_trials/design_matrix/design_matrix_stim_resolved_refined_length", seq_length, "_1to300.txt", sep = "")
    matrix_design = read.table(filename, sep = ",", header = TRUE)
    data = matrix_design
    colnames(data) = c("response", "condition", "subject", "session")
    correct_avg = c(correct_avg, mean(data$response))
    correct_sd = c(correct_sd, sd(data$response))
    filename = paste(path, "/", fun, "/incorrect_trials/design_matrix/design_matrix_stim_resolved_refined_length", seq_length, "_1to300.txt", sep = "")
    matrix_design = read.table(filename, sep = ",", header = TRUE)
    data = matrix_design
    colnames(data) = c("response", "condition", "subject", "session")
    incorrect_avg = c(incorrect_avg, mean(data$response))
    incorrect_sd = c(incorrect_sd, sd(data$response))
    
    filename = paste(path, "/", fun, "/correct_trials/design_matrix/design_matrix_stim_resolved_refined_length", seq_length, "_400to700.txt", sep = "")
    matrix_design = read.table(filename, sep = ",", header = TRUE)
    data = matrix_design
    colnames(data) = c("response", "condition", "subject", "session")
    correct_avg = c(correct_avg, mean(data$response))
    correct_sd = c(correct_sd, sd(data$response))
    filename = paste(path, "/", fun, "/incorrect_trials/design_matrix/design_matrix_stim_resolved_refined_length", seq_length, "_400to700.txt", sep = "")
    matrix_design = read.table(filename, sep = ",", header = TRUE)
    data = matrix_design
    colnames(data) = c("response", "condition", "subject", "session")
    incorrect_avg = c(incorrect_avg, mean(data$response))
    incorrect_sd = c(incorrect_sd, sd(data$response))
    
    filename = paste(path, "/", fun, "/correct_trials/design_matrix/design_matrix_stim_resolved_refined_length", seq_length, "_800to1100.txt", sep = "")
    matrix_design = read.table(filename, sep = ",", header = TRUE)
    data = matrix_design
    colnames(data) = c("response", "condition", "subject", "session")
    correct_avg = c(correct_avg, mean(data$response))
    correct_sd = c(correct_sd, sd(data$response))
    filename = paste(path, "/", fun, "/incorrect_trials/design_matrix/design_matrix_stim_resolved_refined_length", seq_length, "_800to1100.txt", sep = "")
    matrix_design = read.table(filename, sep = ",", header = TRUE)
    data = matrix_design
    colnames(data) = c("response", "condition", "subject", "session")
    incorrect_avg = c(incorrect_avg, mean(data$response))
    incorrect_sd = c(incorrect_sd, sd(data$response))
    
    filename = paste(path, "/", fun, "/correct_trials/design_matrix/design_matrix_stim_resolved_refined_length", seq_length, "_1200to1500.txt", sep = "")
    matrix_design = read.table(filename, sep = ",", header = TRUE)
    data = matrix_design
    colnames(data) = c("response", "condition", "subject", "session")
    correct_avg = c(correct_avg, mean(data$response))
    correct_sd = c(correct_sd, sd(data$response))
    filename = paste(path, "/", fun, "/incorrect_trials/design_matrix/design_matrix_stim_resolved_refined_length", seq_length, "_1200to1500.txt", sep = "")
    matrix_design = read.table(filename, sep = ",", header = TRUE)
    data = matrix_design
    colnames(data) = c("response", "condition", "subject", "session")
    incorrect_avg = c(incorrect_avg, mean(data$response))
    incorrect_sd = c(incorrect_sd, sd(data$response))
    
    # delay
    
    for (time_i in 1:nrow(delay_segments)) {
      time_window = paste(delay_segments[time_i,1], "to", delay_segments[time_i,2], sep = "")
      
      # correct_trials
      filename = paste(path, "/", fun, "/correct_trials/design_matrix/design_matrix_delay_length", seq_length, "_", time_window, ".txt", sep = "")
      matrix_design = read.table(filename, sep = ",", header = TRUE)
      data = matrix_design
      colnames(data) = c("response", "condition", "subject", "session")
      correct_avg = c(correct_avg, mean(data$response))
      correct_sd = c(correct_sd, sd(data$response))
      
      # incorrect_trials
      filename = paste(path, "/", fun, "/incorrect_trials/design_matrix/design_matrix_delay_length", seq_length, "_", time_window, ".txt", sep = "")
      matrix_design = read.table(filename, sep = ",", header = TRUE)
      data = matrix_design
      colnames(data) = c("response", "condition", "subject", "session")
      incorrect_avg = c(incorrect_avg, mean(data$response))
      incorrect_sd = c(incorrect_sd, sd(data$response))
      
    }
    

    # data
    time = c(-1050, -750, -450, -150, delay_segments[,2] - 150) / 1000
    correct_se = correct_sd / sqrt(15)
    incorrect_se = incorrect_sd / sqrt(15)
    avg_column = c(correct_avg, incorrect_avg)
    sd_column = c(correct_sd, incorrect_sd)
    se_column = c(correct_se, incorrect_se)
    performance = c(rep("correct trials", length(correct_avg)), 
                    rep("incorrect trials", length(correct_avg)))
    data = as.data.frame(cbind(avg_column, sd_column, se_column, time, performance))
    colnames(data) = c("avg", "sd", "se", "time", "performance")
    data$performance = factor(data$performance)
    data$avg = as.numeric(data$avg)
    data$sd = as.numeric(data$sd)
    data$se = as.numeric(data$se)
    data$time = as.numeric(data$time)
    
    # data p-values
    
    # Find segments with consecutive values below 0.05
    p_values_tmp = p_values[17:33]
    below_threshold <- p_values_tmp < 0.05
    start_points <- which(diff(c(FALSE, below_threshold)) == 1)
    end_points <- which(diff(c(below_threshold, FALSE)) == -1)
    if (separability == "distance") {
      y_pvalue = max(data$avg + data$se) + (max(data$avg + data$se)/200)
    } else {
      y_pvalue = max(data$avg + data$se) + (max(data$avg + data$se)/50)
    }
    
    # Create segments data frame
    segments <- data.frame(
      x = data$time[start_points],
      xend = data$time[end_points],
      y = rep(y_pvalue, length(start_points)),  # You can adjust the y-value as needed
      yend = rep(y_pvalue, length(start_points))  # Same y-value for end points
    )
    segments$x = segments$x - 0.05
    segments$xend = segments$xend + 0.05
    

    
    
    
    plot = 
      ggplot() +
      geom_vline(xintercept = 0, linetype = "dashed", color = "black", size = 0.3) +
      geom_segment(data = segments, aes(x = x, xend = xend, y = y, yend = yend), color = "black", size = 1) +
      geom_line(data = subset(data, performance == "incorrect trials"), 
                aes(x = time, y = avg, color = "incorrect trials"), size = 1) + 
      geom_line(data = subset(data, performance == "correct trials"), 
                aes(x = time, y = avg, color = "correct trials"), size = 1) + 
      geom_ribbon(data = subset(data, performance == "incorrect trials"), 
                  aes(x = time, ymin = avg - se, ymax = avg + se), alpha = 0.2, fill = "#00BFC4") + 
      geom_ribbon(data = subset(data, performance == "correct trials"), 
                  aes(x = time, ymin = avg - se, ymax = avg + se), alpha = 0.2, fill = "#F8766D") + 
      scale_color_manual(
        values = c("correct trials" = "#F8766D", "incorrect trials" = "#00BFC4"),
        labels = c("correct trials" = "Correct", "incorrect trials" = "Incorrect")) +
      labs(color = " ") +
      xlab("Time (s)") + ylab(separability) +
      scale_x_continuous(breaks = seq(0, 4, by = 1), limits = c(min(time) - 0.05, 4)) +
      theme(plot.title = element_text(hjust = 0.5), 
            legend.position = "none", 
            panel.background = element_blank(),
            plot.background = element_blank(), 
            axis.line = element_line(color = "black"), 
            axis.title = element_text(color = "black", size = 12), 
            axis.text.x = element_text(color = "black", size = 12), 
            axis.text.y = element_text(color = "black", size = 12), 
            panel.grid.major = element_blank(), 
            panel.grid.minor = element_blank())
    
    
    ggsave(paste(path_figures, "/connected_scatterplot_separability_", separability, "_length", seq_length, "_correct_vs_incorrect.png", sep = ""), plot, width = 3, height = 3, units = "in")
    
    
    plot = 
      ggplot() +
      geom_vline(xintercept = 0, linetype = "dashed", color = "black", size = 0.3) +
      geom_segment(data = segments, aes(x = x, xend = xend, y = y, yend = yend), color = "black", size = 1) +
      geom_line(data = subset(data, performance == "incorrect trials"), 
                aes(x = time, y = avg, color = "incorrect trials"), size = 1) + 
      geom_line(data = subset(data, performance == "correct trials"), 
                aes(x = time, y = avg, color = "correct trials"), size = 1) + 
      geom_ribbon(data = subset(data, performance == "incorrect trials"), 
                  aes(x = time, ymin = avg - se, ymax = avg + se), alpha = 0.2, fill = "#00BFC4") + 
      geom_ribbon(data = subset(data, performance == "correct trials"), 
                  aes(x = time, ymin = avg - se, ymax = avg + se), alpha = 0.2, fill = "#F8766D") + 
      scale_color_manual(
        values = c("correct trials" = "#F8766D", "incorrect trials" = "#00BFC4"),
        labels = c("correct trials" = "Correct", "incorrect trials" = "Incorrect")) +
      labs(color = " ") +
      xlab("Time (s)") + ylab(separability) +
      scale_x_continuous(breaks = seq(0, 4, by = 1), limits = c(min(time) - 0.05, 4)) +
      theme(plot.title = element_text(hjust = 0.5), 
            legend.position = "bottom", 
            panel.background = element_blank(),
            plot.background = element_blank(), 
            axis.line = element_line(color = "black"), 
            axis.title = element_text(color = "black", size = 12), 
            axis.text.x = element_text(color = "black", size = 12), 
            axis.text.y = element_text(color = "black", size = 12), 
            panel.grid.major = element_blank(), 
            panel.grid.minor = element_blank())
    
    ggsave(paste(path_figures, "/connected_scatterplot_separability_", separability, "_length", seq_length, "_correct_vs_incorrect_legend.png", sep = ""), plot, width = 3.5, height = 3.5, units = "in")
    
    
    
    
  }
}
