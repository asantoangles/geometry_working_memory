##############################
## connected scatter plots  ##
##############################

segments_start <- seq(from = 1, to = 3601, by = 300)
segments_end <- seq(from = 300, to = 4000, by = 300)
delay_segments = cbind(segments_start, segments_end)

for (metric_i in metrics) {
  
  path = paste(path_results, "/group_results/", metric_i, sep = "")
  
  if (metric_i == "principal_angle") {
    ylab_metric = "Principal angle (°)"
    p_values = p_values_all[1:31]
    y_axis_step = 20
  } else if (metric_i == "vaf") {
    ylab_metric = "VAF"
    p_values = p_values_all[32:62]
    y_axis_step = 0.2
  }
  
  for (fun in functions_data) {

    ### length 3
    
    length_seq = 3
    correct_avg = c()
    correct_sd = c()
    incorrect_avg = c()
    incorrect_sd = c()
    
    correct_subject_data = data.frame(
      subject = numeric(0),
      time = numeric(0),
      value = numeric(0)
    )
    
    incorrect_subject_data = data.frame(
      subject = numeric(0),
      time = numeric(0),
      value = numeric(0)
    )
    
    # stim
    
    # correct trials
    filename = paste(path, "/", fun, "/correct_trials/design_matrix/stim_length", length_seq, "_400to700.txt", sep = "")
    design_matrix = read.table(filename, sep = ",", header = TRUE)
    colnames(design_matrix) = c("response", "Subspaces", "subject", "session")
    design_matrix$Subspaces = factor(ifelse(design_matrix$Subspaces == 1, "Between ranks", "Within rank"))
    data = subset(design_matrix, Subspaces == "Between ranks")
    correct_avg = c(correct_avg, mean(data$response))
    correct_sd = c(correct_sd, sd(data$response))
    data_correct = data
    # incorrect trials
    filename = paste(path, "/", fun, "/incorrect_trials/design_matrix/stim_length", length_seq, "_400to700.txt", sep = "")
    design_matrix = read.table(filename, sep = ",", header = TRUE)
    colnames(design_matrix) = c("response", "Subspaces", "subject", "session")
    design_matrix$Subspaces = factor(ifelse(design_matrix$Subspaces == 1, "Between ranks", "Within rank"))
    data = subset(design_matrix, Subspaces == "Between ranks")
    incorrect_avg = c(incorrect_avg, mean(data$response))
    incorrect_sd = c(incorrect_sd, sd(data$response))
    data_incorrect = data
    
    current_row = nrow(correct_subject_data)
    t = (current_row / n_subjects) + 1
    correct_subject_data[(current_row + 1):(current_row + n_subjects), ] = data.frame(
      subject = 1:n_subjects,
      time = rep(t, n_subjects),
      value = tapply(data_correct$response, data_correct$subject, mean) )
    
    current_row = nrow(incorrect_subject_data)
    t = (current_row / n_subjects) + 1
    incorrect_subject_data[(current_row + 1):(current_row + n_subjects), ] = data.frame(
      subject = 1:n_subjects,
      time = rep(t, n_subjects),
      value = tapply(data_incorrect$response, data_incorrect$subject, mean) )
    
    # correct trials
    filename = paste(path, "/", fun, "/correct_trials/design_matrix/stim_length", length_seq, "_800to1100.txt", sep = "")
    design_matrix = read.table(filename, sep = ",", header = TRUE)
    colnames(design_matrix) = c("response", "Subspaces", "subject", "session")
    design_matrix$Subspaces = factor(ifelse(design_matrix$Subspaces == 1, "Between ranks", "Within rank"))
    data = subset(design_matrix, Subspaces == "Between ranks")
    correct_avg = c(correct_avg, mean(data$response))
    correct_sd = c(correct_sd, sd(data$response))
    data_correct = data
    # incorrect trials
    filename = paste(path, "/", fun, "/incorrect_trials/design_matrix/stim_length", length_seq, "_800to1100.txt", sep = "")
    design_matrix = read.table(filename, sep = ",", header = TRUE)
    colnames(design_matrix) = c("response", "Subspaces", "subject", "session")
    design_matrix$Subspaces = factor(ifelse(design_matrix$Subspaces == 1, "Between ranks", "Within rank"))
    data = subset(design_matrix, Subspaces == "Between ranks")
    incorrect_avg = c(incorrect_avg, mean(data$response))
    incorrect_sd = c(incorrect_sd, sd(data$response))
    data_incorrect = data
    
    current_row = nrow(correct_subject_data)
    t = (current_row / n_subjects) + 1
    correct_subject_data[(current_row + 1):(current_row + n_subjects), ] = data.frame(
      subject = 1:n_subjects,
      time = rep(t, n_subjects),
      value = tapply(data_correct$response, data_correct$subject, mean) )
    
    current_row = nrow(incorrect_subject_data)
    t = (current_row / n_subjects) + 1
    incorrect_subject_data[(current_row + 1):(current_row + n_subjects), ] = data.frame(
      subject = 1:n_subjects,
      time = rep(t, n_subjects),
      value = tapply(data_incorrect$response, data_incorrect$subject, mean) )
    
    # delay
    
    for (time_i in 1:nrow(delay_segments)) {
      time_window = paste(delay_segments[time_i,1], "to", delay_segments[time_i,2], sep = "")
      # correct trials
      filename = paste(path, "/", fun, "/correct_trials/design_matrix/delay_length", length_seq, "_", time_window, ".txt", sep = "")
      design_matrix = read.table(filename, sep = ",", header = TRUE)
      colnames(design_matrix) = c("response", "Subspaces", "subject", "session")
      design_matrix$Subspaces = factor(ifelse(design_matrix$Subspaces == 1, "Between ranks", "Within rank"))
      data = subset(design_matrix, Subspaces == "Between ranks")
      correct_avg = c(correct_avg, mean(data$response))
      correct_sd = c(correct_sd, sd(data$response))
      data_correct = data
      # incorrect trials
      filename = paste(path, "/", fun, "/incorrect_trials/design_matrix/delay_length", length_seq, "_", time_window, ".txt", sep = "")
      design_matrix = read.table(filename, sep = ",", header = TRUE)
      colnames(design_matrix) = c("response", "Subspaces", "subject", "session")
      design_matrix$Subspaces = factor(ifelse(design_matrix$Subspaces == 1, "Between ranks", "Within rank"))
      data = subset(design_matrix, Subspaces == "Between ranks")
      incorrect_avg = c(incorrect_avg, mean(data$response))
      incorrect_sd = c(incorrect_sd, sd(data$response))
      data_incorrect = data
      
      current_row = nrow(correct_subject_data)
      t = (current_row / n_subjects) + 1
      correct_subject_data[(current_row + 1):(current_row + n_subjects), ] = data.frame(
        subject = 1:n_subjects,
        time = rep(t, n_subjects),
        value = tapply(data_correct$response, data_correct$subject, mean) )
      
      current_row = nrow(incorrect_subject_data)
      t = (current_row / n_subjects) + 1
      incorrect_subject_data[(current_row + 1):(current_row + n_subjects), ] = data.frame(
        subject = 1:n_subjects,
        time = rep(t, n_subjects),
        value = tapply(data_incorrect$response, data_incorrect$subject, mean) )
      
    }
    

    # data
    time = c(-450, -150, delay_segments[,2] - 150) / 1000
    
    correct_subject_data$time = time[correct_subject_data$time]
    incorrect_subject_data$time = time[incorrect_subject_data$time]
    
    correct_se = correct_sd / sqrt(15)
    incorrect_se = incorrect_sd / sqrt(15)
    avg_column = c(correct_avg, incorrect_avg)
    sd_column = c(correct_sd, incorrect_sd)
    se_column = c(correct_se, incorrect_se)
    performance = c(rep("correct_trials", length(correct_avg)), rep("incorrect_trials", length(correct_avg)))
    data = as.data.frame(cbind(avg_column, sd_column, se_column, time, performance))
    colnames(data) = c("avg", "sd", "se", "time", "performance")
    data$time = as.numeric(data$time)
    data$avg = as.numeric(data$avg)
    data$sd = as.numeric(data$sd)
    data$se = as.numeric(data$se)
    data$performance = factor(data$performance)
    
    ## data p-values
    
    # Find segments with consecutive values below 0.05
    p_values_tmp = p_values[1:15]
    below_threshold <- p_values_tmp < 0.025
    start_points <- which(diff(c(FALSE, below_threshold)) == 1)
    end_points <- which(diff(c(below_threshold, FALSE)) == -1)
    y_pvalue = max(c(correct_subject_data$value, incorrect_subject_data$value))
    
    # Create segments data frame
    segments <- data.frame(
      x = data$time[start_points],
      xend = data$time[end_points],
      y = rep(y_pvalue, length(start_points)),  # You can adjust the y-value as needed
      yend = rep(y_pvalue, length(start_points))  # Same y-value for end points
    )
    segments$x = segments$x - 0.05
    segments$xend = segments$xend + 0.05
    
    # Create segments data frame
    segments <- data.frame(
      x = data$time[start_points],
      xend = data$time[end_points],
      y = rep(y_pvalue, length(start_points)),  # You can adjust the y-value as needed
      yend = rep(y_pvalue, length(start_points))  # Same y-value for end points
    )
    segments$x = segments$x - 0.05
    segments$xend = segments$xend + 0.05

    ylim_metric = c(min(c(correct_subject_data$value, incorrect_subject_data$value)),
                    max(c(correct_subject_data$value, incorrect_subject_data$value)))
    
    plot = 
      ggplot() + 
      geom_vline(xintercept = 0, linetype = "dashed", color = "black", size = 0.3) +
      geom_segment(data = segments, aes(x = x, xend = xend, y = y, yend = yend), color = "black", size = 1) +
      geom_point(data = correct_subject_data, aes(x = time, y = value), color = "red", alpha = 0.3, size = 0.3) +
      geom_point(data = incorrect_subject_data, aes(x = time, y = value), color = "blue", alpha = 0.3, size = 0.3) +
      geom_line(data = correct_subject_data, aes(x = time, y = value, group = subject), color = "red", alpha = 0.2, size = 0.3) +
      geom_line(data = incorrect_subject_data, aes(x = time, y = value, group = subject), color = "blue", alpha = 0.2, size = 0.3) +
      geom_line(data = subset(data, performance == "correct_trials"), aes(x = time, y = avg, color = "correct trials"), size = 1.5) + 
      geom_line(data = subset(data, performance == "incorrect_trials"), aes(x = time, y = avg, color = "incorrect trials"), size = 1.5) + 
      geom_ribbon(data = subset(data, performance == "correct_trials"), 
                  aes(x = time, ymin = avg - se, ymax = avg + se, fill = "correct trials"), alpha = 0.2) + 
      geom_ribbon(data = subset(data, performance == "incorrect_trials"), 
                  aes(x = time, ymin = avg - se, ymax = avg + se, fill = "incorrect trials"), alpha = 0.2) + 
      scale_color_manual(values = c("correct trials" = "red", "incorrect trials" = "blue")) + 
      scale_fill_manual(values = c("correct trials" = "red", "incorrect trials" = "blue")) + 
      labs(color = " ", fill = " ") +  # Set the legend titles
      xlab(NULL) + ylab(NULL) +
      scale_x_continuous(breaks = seq(0, 4, by = 1), limits = c(min(time) - 0.05, 4)) +
      scale_y_continuous(
        limits = ylim_metric,
        breaks = seq(
          floor(ylim_metric[1] / y_axis_step) * y_axis_step,
          ceiling(ylim_metric[2] / y_axis_step) * y_axis_step,
          by = y_axis_step
        )
      ) + 
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
    
    
    ggsave(paste(path_figures, "/connected_scatterplot_", metric_i, "_length", length_seq, "_correct_vs_incorrect.png", sep = ""), plot, width = 3, height = 3, units = "in")
    
    plot = 
      ggplot() + 
      geom_vline(xintercept = 0, linetype = "dashed", color = "black", size = 0.3) +
      geom_segment(data = segments, aes(x = x, xend = xend, y = y, yend = yend), color = "black", size = 1) +
      geom_point(data = correct_subject_data, aes(x = time, y = value), color = "red", alpha = 0.3, size = 0.3) +
      geom_point(data = incorrect_subject_data, aes(x = time, y = value), color = "blue", alpha = 0.3, size = 0.3) +
      geom_line(data = correct_subject_data, aes(x = time, y = value, group = subject), color = "red", alpha = 0.2, size = 0.3) +
      geom_line(data = incorrect_subject_data, aes(x = time, y = value, group = subject), color = "blue", alpha = 0.2, size = 0.3) +
      geom_line(data = subset(data, performance == "correct_trials"), aes(x = time, y = avg, color = "correct trials"), size = 1.5) + 
      geom_line(data = subset(data, performance == "incorrect_trials"), aes(x = time, y = avg, color = "incorrect trials"), size = 1.5) + 
      geom_ribbon(data = subset(data, performance == "correct_trials"), 
                  aes(x = time, ymin = avg - se, ymax = avg + se, fill = "correct trials"), alpha = 0.2) + 
      geom_ribbon(data = subset(data, performance == "incorrect_trials"), 
                  aes(x = time, ymin = avg - se, ymax = avg + se, fill = "incorrect trials"), alpha = 0.2) + 
      scale_color_manual(values = c("correct trials" = "red", "incorrect trials" = "blue")) + 
      scale_fill_manual(values = c("correct trials" = "red", "incorrect trials" = "blue")) + 
      labs(color = " ", fill = " ") +  # Set the legend titles
      xlab("Time (s)") + ylab(ylab_metric) +
      scale_x_continuous(breaks = seq(0, 4, by = 1), limits = c(min(time) - 0.05, 4)) +
      scale_y_continuous(
        limits = ylim_metric,
        breaks = seq(
          floor(ylim_metric[1] / y_axis_step) * y_axis_step,
          ceiling(ylim_metric[2] / y_axis_step) * y_axis_step,
          by = y_axis_step
        )
      ) + 
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
    
    
    ggsave(paste(path_figures, "/connected_scatterplot_", metric_i, "_length", length_seq, "_correct_vs_incorrect_legend.png", sep = ""), plot, width = 3, height = 3, units = "in")
    
    
    
    
    
    
    
    
    ### length 4
    
    length_seq = 4
    correct_avg = c()
    correct_sd = c()
    incorrect_avg = c()
    incorrect_sd = c()
    
    correct_subject_data = data.frame(
      subject = numeric(0),
      time = numeric(0),
      value = numeric(0)
    )
    
    incorrect_subject_data = data.frame(
      subject = numeric(0),
      time = numeric(0),
      value = numeric(0)
    )
    
    # stim
    
    # correct trials
    filename = paste(path, "/", fun, "/correct_trials/design_matrix/stim_length", length_seq, "_400to700.txt", sep = "")
    design_matrix = read.table(filename, sep = ",", header = TRUE)
    colnames(design_matrix) = c("response", "Subspaces", "subject", "session")
    design_matrix$Subspaces = factor(ifelse(design_matrix$Subspaces == 1, "Between ranks", "Within rank"))
    data = subset(design_matrix, Subspaces == "Between ranks")
    correct_avg = c(correct_avg, mean(data$response))
    correct_sd = c(correct_sd, sd(data$response))
    data_correct = data
    # incorrect trials
    filename = paste(path, "/", fun, "/incorrect_trials/design_matrix/stim_length", length_seq, "_400to700.txt", sep = "")
    design_matrix = read.table(filename, sep = ",", header = TRUE)
    colnames(design_matrix) = c("response", "Subspaces", "subject", "session")
    design_matrix$Subspaces = factor(ifelse(design_matrix$Subspaces == 1, "Between ranks", "Within rank"))
    data = subset(design_matrix, Subspaces == "Between ranks")
    incorrect_avg = c(incorrect_avg, mean(data$response))
    incorrect_sd = c(incorrect_sd, sd(data$response))
    data_incorrect = data
    
    current_row = nrow(correct_subject_data)
    t = (current_row / n_subjects) + 1
    correct_subject_data[(current_row + 1):(current_row + n_subjects), ] = data.frame(
      subject = 1:n_subjects,
      time = rep(t, n_subjects),
      value = tapply(data_correct$response, data_correct$subject, mean) )
    
    current_row = nrow(incorrect_subject_data)
    t = (current_row / n_subjects) + 1
    incorrect_subject_data[(current_row + 1):(current_row + n_subjects), ] = data.frame(
      subject = 1:n_subjects,
      time = rep(t, n_subjects),
      value = tapply(data_incorrect$response, data_incorrect$subject, mean) )
    
    # correct trials
    filename = paste(path, "/", fun, "/correct_trials/design_matrix/stim_length", length_seq, "_800to1100.txt", sep = "")
    design_matrix = read.table(filename, sep = ",", header = TRUE)
    colnames(design_matrix) = c("response", "Subspaces", "subject", "session")
    design_matrix$Subspaces = factor(ifelse(design_matrix$Subspaces == 1, "Between ranks", "Within rank"))
    data = subset(design_matrix, Subspaces == "Between ranks")
    correct_avg = c(correct_avg, mean(data$response))
    correct_sd = c(correct_sd, sd(data$response))
    data_correct = data
    # incorrect trials
    filename = paste(path, "/", fun, "/incorrect_trials/design_matrix/stim_length", length_seq, "_800to1100.txt", sep = "")
    design_matrix = read.table(filename, sep = ",", header = TRUE)
    colnames(design_matrix) = c("response", "Subspaces", "subject", "session")
    design_matrix$Subspaces = factor(ifelse(design_matrix$Subspaces == 1, "Between ranks", "Within rank"))
    data = subset(design_matrix, Subspaces == "Between ranks")
    incorrect_avg = c(incorrect_avg, mean(data$response))
    incorrect_sd = c(incorrect_sd, sd(data$response))
    data_incorrect = data
    
    current_row = nrow(correct_subject_data)
    t = (current_row / n_subjects) + 1
    correct_subject_data[(current_row + 1):(current_row + n_subjects), ] = data.frame(
      subject = 1:n_subjects,
      time = rep(t, n_subjects),
      value = tapply(data_correct$response, data_correct$subject, mean) )
    
    current_row = nrow(incorrect_subject_data)
    t = (current_row / n_subjects) + 1
    incorrect_subject_data[(current_row + 1):(current_row + n_subjects), ] = data.frame(
      subject = 1:n_subjects,
      time = rep(t, n_subjects),
      value = tapply(data_incorrect$response, data_incorrect$subject, mean) )
    
    
    # correct trials
    filename = paste(path, "/", fun, "/correct_trials/design_matrix/stim_length", length_seq, "_1200to1500.txt", sep = "")
    design_matrix = read.table(filename, sep = ",", header = TRUE)
    colnames(design_matrix) = c("response", "Subspaces", "subject", "session")
    design_matrix$Subspaces = factor(ifelse(design_matrix$Subspaces == 1, "Between ranks", "Within rank"))
    data = subset(design_matrix, Subspaces == "Between ranks")
    correct_avg = c(correct_avg, mean(data$response))
    correct_sd = c(correct_sd, sd(data$response))
    data_correct = data
    # incorrect trials
    filename = paste(path, "/", fun, "/incorrect_trials/design_matrix/stim_length", length_seq, "_1200to1500.txt", sep = "")
    design_matrix = read.table(filename, sep = ",", header = TRUE)
    colnames(design_matrix) = c("response", "Subspaces", "subject", "session")
    design_matrix$Subspaces = factor(ifelse(design_matrix$Subspaces == 1, "Between ranks", "Within rank"))
    data = subset(design_matrix, Subspaces == "Between ranks")
    incorrect_avg = c(incorrect_avg, mean(data$response))
    incorrect_sd = c(incorrect_sd, sd(data$response))
    data_incorrect = data
    
    current_row = nrow(correct_subject_data)
    t = (current_row / n_subjects) + 1
    correct_subject_data[(current_row + 1):(current_row + n_subjects), ] = data.frame(
      subject = 1:n_subjects,
      time = rep(t, n_subjects),
      value = tapply(data_correct$response, data_correct$subject, mean) )
    
    current_row = nrow(incorrect_subject_data)
    t = (current_row / n_subjects) + 1
    incorrect_subject_data[(current_row + 1):(current_row + n_subjects), ] = data.frame(
      subject = 1:n_subjects,
      time = rep(t, n_subjects),
      value = tapply(data_incorrect$response, data_incorrect$subject, mean) )
    
    # delay
    
    
    
    for (time_i in 1:nrow(delay_segments)) {
      time_window = paste(delay_segments[time_i,1], "to", delay_segments[time_i,2], sep = "")
      # correct trials
      filename = paste(path, "/", fun, "/correct_trials/design_matrix/delay_length", length_seq, "_", time_window, ".txt", sep = "")
      design_matrix = read.table(filename, sep = ",", header = TRUE)
      colnames(design_matrix) = c("response", "Subspaces", "subject", "session")
      design_matrix$Subspaces = factor(ifelse(design_matrix$Subspaces == 1, "Between ranks", "Within rank"))
      data = subset(design_matrix, Subspaces == "Between ranks")
      correct_avg = c(correct_avg, mean(data$response))
      correct_sd = c(correct_sd, sd(data$response))
      data_correct = data
      # correct trials
      filename = paste(path, "/", fun, "/incorrect_trials/design_matrix/delay_length", length_seq, "_", time_window, ".txt", sep = "")
      design_matrix = read.table(filename, sep = ",", header = TRUE)
      colnames(design_matrix) = c("response", "Subspaces", "subject", "session")
      design_matrix$Subspaces = factor(ifelse(design_matrix$Subspaces == 1, "Between ranks", "Within rank"))
      data = subset(design_matrix, Subspaces == "Between ranks")
      incorrect_avg = c(incorrect_avg, mean(data$response))
      incorrect_sd = c(incorrect_sd, sd(data$response))
      data_incorrect = data
      
      current_row = nrow(correct_subject_data)
      t = (current_row / n_subjects) + 1
      correct_subject_data[(current_row + 1):(current_row + n_subjects), ] = data.frame(
        subject = 1:n_subjects,
        time = rep(t, n_subjects),
        value = tapply(data_correct$response, data_correct$subject, mean) )
      
      current_row = nrow(incorrect_subject_data)
      t = (current_row / n_subjects) + 1
      incorrect_subject_data[(current_row + 1):(current_row + n_subjects), ] = data.frame(
        subject = 1:n_subjects,
        time = rep(t, n_subjects),
        value = tapply(data_incorrect$response, data_incorrect$subject, mean) )
      
    }
    
    
    # data
    time = c(-750, -450, -150, delay_segments[,2] - 150) / 1000
    
    correct_subject_data$time = time[correct_subject_data$time]
    incorrect_subject_data$time = time[incorrect_subject_data$time]

    correct_se = correct_sd / sqrt(15)
    incorrect_se = incorrect_sd / sqrt(15)
    avg_column = c(correct_avg, incorrect_avg)
    sd_column = c(correct_sd, incorrect_sd)
    se_column = c(correct_se, incorrect_se)
    performance = c(rep("correct_trials", length(correct_avg)), rep("incorrect_trials", length(correct_avg)))
    data = as.data.frame(cbind(avg_column, sd_column, se_column, time, performance))
    colnames(data) = c("avg", "sd", "se", "time", "performance")
    data$time = as.numeric(data$time)
    data$avg = as.numeric(data$avg)
    data$sd = as.numeric(data$sd)
    data$se = as.numeric(data$se)
    data$performance = factor(data$performance)
    
    ## data p-values
    
    # Find segments with consecutive values below 0.05
    p_values_tmp = p_values[16:31]
    below_threshold <- p_values_tmp < 0.025
    start_points <- which(diff(c(FALSE, below_threshold)) == 1)
    end_points <- which(diff(c(below_threshold, FALSE)) == -1)
    y_pvalue = max(c(correct_subject_data$value, incorrect_subject_data$value))
    
    # Create segments data frame
    segments <- data.frame(
      x = data$time[start_points],
      xend = data$time[end_points],
      y = rep(y_pvalue, length(start_points)),  # You can adjust the y-value as needed
      yend = rep(y_pvalue, length(start_points))  # Same y-value for end points
    )
    segments$x = segments$x - 0.05
    segments$xend = segments$xend + 0.05
    
    ylim_metric = c(min(c(correct_subject_data$value, incorrect_subject_data$value)),
                    max(c(correct_subject_data$value, incorrect_subject_data$value)))
    
    plot = 
      ggplot() + 
      geom_vline(xintercept = 0, linetype = "dashed", color = "black", size = 0.3) +
      geom_segment(data = segments, aes(x = x, xend = xend, y = y, yend = yend), color = "black", size = 1) +
      geom_point(data = correct_subject_data, aes(x = time, y = value), color = "red", alpha = 0.3, size = 0.3) +
      geom_point(data = incorrect_subject_data, aes(x = time, y = value), color = "blue", alpha = 0.3, size = 0.3) +
      geom_line(data = correct_subject_data, aes(x = time, y = value, group = subject), color = "red", alpha = 0.2, size = 0.3) +
      geom_line(data = incorrect_subject_data, aes(x = time, y = value, group = subject), color = "blue", alpha = 0.2, size = 0.3) +
      geom_line(data = subset(data, performance == "correct_trials"), aes(x = time, y = avg, color = "correct trials"), size = 1.5) + 
      geom_line(data = subset(data, performance == "incorrect_trials"), aes(x = time, y = avg, color = "incorrect trials"), size = 1.5) + 
      geom_ribbon(data = subset(data, performance == "correct_trials"), 
                  aes(x = time, ymin = avg - se, ymax = avg + se, fill = "correct trials"), alpha = 0.2) + 
      geom_ribbon(data = subset(data, performance == "incorrect_trials"), 
                  aes(x = time, ymin = avg - se, ymax = avg + se, fill = "incorrect trials"), alpha = 0.2) + 
      scale_color_manual(values = c("correct trials" = "red", "incorrect trials" = "blue")) + 
      scale_fill_manual(values = c("correct trials" = "red", "incorrect trials" = "blue")) + 
      labs(color = " ", fill = " ") +  # Set the legend titles
      xlab(NULL) + ylab(NULL) +
      scale_x_continuous(breaks = seq(0, 4, by = 1), limits = c(min(time) - 0.05, 4)) +
      scale_y_continuous(
        limits = ylim_metric,
        breaks = seq(
          floor(ylim_metric[1] / y_axis_step) * y_axis_step,
          ceiling(ylim_metric[2] / y_axis_step) * y_axis_step,
          by = y_axis_step
        )
      ) + 
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
    
    
    ggsave(paste(path_figures, "/connected_scatterplot_", metric_i, "_length", length_seq, "_correct_vs_incorrect.png", sep = ""), plot, width = 3, height = 3, units = "in")
    
    plot = 
      ggplot() + 
      geom_vline(xintercept = 0, linetype = "dashed", color = "black", size = 0.3) +
      geom_segment(data = segments, aes(x = x, xend = xend, y = y, yend = yend), color = "black", size = 1) +
      geom_point(data = correct_subject_data, aes(x = time, y = value), color = "red", alpha = 0.3, size = 0.3) +
      geom_point(data = incorrect_subject_data, aes(x = time, y = value), color = "blue", alpha = 0.3, size = 0.3) +
      geom_line(data = correct_subject_data, aes(x = time, y = value, group = subject), color = "red", alpha = 0.2, size = 0.3) +
      geom_line(data = incorrect_subject_data, aes(x = time, y = value, group = subject), color = "blue", alpha = 0.2, size = 0.3) +
      geom_line(data = subset(data, performance == "correct_trials"), aes(x = time, y = avg, color = "correct trials"), size = 1.5) + 
      geom_line(data = subset(data, performance == "incorrect_trials"), aes(x = time, y = avg, color = "incorrect trials"), size = 1.5) + 
      geom_ribbon(data = subset(data, performance == "correct_trials"), 
                  aes(x = time, ymin = avg - se, ymax = avg + se, fill = "correct trials"), alpha = 0.2) + 
      geom_ribbon(data = subset(data, performance == "incorrect_trials"), 
                  aes(x = time, ymin = avg - se, ymax = avg + se, fill = "incorrect trials"), alpha = 0.2) + 
      scale_color_manual(values = c("correct trials" = "red", "incorrect trials" = "blue")) + 
      scale_fill_manual(values = c("correct trials" = "red", "incorrect trials" = "blue")) + 
      labs(color = " ", fill = " ") +  # Set the legend titles
      xlab("Time (s)") + ylab(ylab_metric) +
      scale_x_continuous(breaks = seq(0, 4, by = 1), limits = c(min(time) - 0.05, 4)) +
      scale_y_continuous(
        limits = ylim_metric,
        breaks = seq(
          floor(ylim_metric[1] / y_axis_step) * y_axis_step,
          ceiling(ylim_metric[2] / y_axis_step) * y_axis_step,
          by = y_axis_step
        )
      ) + 
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
    
    
    ggsave(paste(path_figures, "/connected_scatterplot_", metric_i, "_length", length_seq, "_correct_vs_incorrect_legend.png", sep = ""), plot, width = 3, height = 3, units = "in")
    
    
      
            

  }
}
