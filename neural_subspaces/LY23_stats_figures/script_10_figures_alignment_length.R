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
    p_values = p_values_all[1:15]
    y_axis_step = 20
  } else if (metric_i == "principal_angle_min") {
    ylab_metric = "Principal angle min (°)"
    p_values = p_values_all[16:30]
    y_axis_step = 20
  } else if (metric_i == "vaf") {
    ylab_metric = "VAF"
    p_values = p_values_all[31:45]
    y_axis_step = 0.05
  }
  
  for (fun in functions_data) {

    length3_avg = c()
    length3_sd = c()
    length4_avg = c()
    length4_sd = c()
    
    length3_subject_data = data.frame(
      subject = numeric(0),
      time = numeric(0),
      value = numeric(0)
    )
    
    length4_subject_data = data.frame(
      subject = numeric(0),
      time = numeric(0),
      value = numeric(0)
    )
    
    # stim
    
    filename = paste(path, "/", fun, "/correct_trials/design_matrix/stim_length3_400to700.txt", sep = "")
    design_matrix = read.table(filename, sep = ",", header = TRUE)
    colnames(design_matrix) = c("response", "Subspaces", "subject", "session")
    design_matrix$Subspaces = factor(ifelse(design_matrix$Subspaces == 1, "Between ranks", "Within rank"))
    data = subset(design_matrix, Subspaces == "Between ranks")
    length3_avg = c(length3_avg, mean(data$response))
    length3_sd = c(length3_sd, sd(data$response))
    data_length3 = data
    filename = paste(path, "/", fun, "/correct_trials/design_matrix/stim_length4_400to700.txt", sep = "")
    design_matrix = read.table(filename, sep = ",", header = TRUE)
    colnames(design_matrix) = c("response", "Subspaces", "subject", "session")
    design_matrix$Subspaces = factor(ifelse(design_matrix$Subspaces == 1, "Between ranks", "Within rank"))
    data = subset(design_matrix, Subspaces == "Between ranks")
    length4_avg = c(length4_avg, mean(data$response))
    length4_sd = c(length4_sd, sd(data$response))
    data_length4 = data
    
    current_row = nrow(length3_subject_data)
    t = (current_row / n_subjects) + 1
    length3_subject_data[(current_row + 1):(current_row + n_subjects), ] = data.frame(
      subject = 1:n_subjects,
      time = rep(t, n_subjects),
      value = tapply(data_length3$response, data_length3$subject, mean) )
    
    current_row = nrow(length4_subject_data)
    t = (current_row / n_subjects) + 1
    length4_subject_data[(current_row + 1):(current_row + n_subjects), ] = data.frame(
      subject = 1:n_subjects,
      time = rep(t, n_subjects),
      value = tapply(data_length4$response, data_length4$subject, mean) )
    
    filename = paste(path, "/", fun, "/correct_trials/design_matrix/stim_length3_800to1100.txt", sep = "")
    design_matrix = read.table(filename, sep = ",", header = TRUE)
    colnames(design_matrix) = c("response", "Subspaces", "subject", "session")
    design_matrix$Subspaces = factor(ifelse(design_matrix$Subspaces == 1, "Between ranks", "Within rank"))
    data = subset(design_matrix, Subspaces == "Between ranks")
    length3_avg = c(length3_avg, mean(data$response))
    length3_sd = c(length3_sd, sd(data$response))
    data_length3 = data
    filename = paste(path, "/", fun, "/correct_trials/design_matrix/stim_length4_800to1100.txt", sep = "")
    design_matrix = read.table(filename, sep = ",", header = TRUE)
    colnames(design_matrix) = c("response", "Subspaces", "subject", "session")
    design_matrix$Subspaces = factor(ifelse(design_matrix$Subspaces == 1, "Between ranks", "Within rank"))
    data = subset(design_matrix, Subspaces == "Between ranks")
    length4_avg = c(length4_avg, mean(data$response))
    length4_sd = c(length4_sd, sd(data$response))
    data_length4 = data
    
    current_row = nrow(length3_subject_data)
    t = (current_row / n_subjects) + 1
    length3_subject_data[(current_row + 1):(current_row + n_subjects), ] = data.frame(
      subject = 1:n_subjects,
      time = rep(t, n_subjects),
      value = tapply(data_length3$response, data_length3$subject, mean) )
    
    current_row = nrow(length4_subject_data)
    t = (current_row / n_subjects) + 1
    length4_subject_data[(current_row + 1):(current_row + n_subjects), ] = data.frame(
      subject = 1:n_subjects,
      time = rep(t, n_subjects),
      value = tapply(data_length4$response, data_length4$subject, mean) )
    
    # delay
    
    for (time_i in 1:nrow(delay_segments)) {
      time_window = paste(delay_segments[time_i,1], "to", delay_segments[time_i,2], sep = "")
      
      # length-3
      filename = paste(path, "/", fun, "/correct_trials/design_matrix/delay_length3_", time_window, ".txt", sep = "")
      design_matrix = read.table(filename, sep = ",", header = TRUE)
      colnames(design_matrix) = c("response", "Subspaces", "subject", "session")
      design_matrix$Subspaces = factor(ifelse(design_matrix$Subspaces == 1, "Between ranks", "Within rank"))
      data = subset(design_matrix, Subspaces == "Between ranks")
      data_length3 = data
      length3_avg = c(length3_avg, mean(data$response))
      length3_sd = c(length3_sd, sd(data$response))
      
      # length-4
      filename = paste(path, "/", fun, "/correct_trials/design_matrix/delay_length4_", time_window, ".txt", sep = "")
      design_matrix = read.table(filename, sep = ",", header = TRUE)
      colnames(design_matrix) = c("response", "Subspaces", "subject", "session")
      design_matrix$Subspaces = factor(ifelse(design_matrix$Subspaces == 1, "Between ranks", "Within rank"))
      data = subset(design_matrix, Subspaces == "Between ranks")
      data_length4 = data
      length4_avg = c(length4_avg, mean(data$response))
      length4_sd = c(length4_sd, sd(data$response))
      
      current_row = nrow(length3_subject_data)
      t = (current_row / n_subjects) + 1
      length3_subject_data[(current_row + 1):(current_row + n_subjects), ] = data.frame(
        subject = 1:n_subjects,
        time = rep(t, n_subjects),
        value = tapply(data_length3$response, data_length3$subject, mean) )
      
      current_row = nrow(length4_subject_data)
      t = (current_row / n_subjects) + 1
      length4_subject_data[(current_row + 1):(current_row + n_subjects), ] = data.frame(
        subject = 1:n_subjects,
        time = rep(t, n_subjects),
        value = tapply(data_length4$response, data_length4$subject, mean) )
    }
    
    # data
    time = c(-450, -150, delay_segments[,2] - 150) / 1000
    
    length3_subject_data$time = time[length3_subject_data$time]
    length4_subject_data$time = time[length4_subject_data$time]
    
    length3_se = length3_sd / sqrt(15)
    length4_se = length4_sd / sqrt(15)
    avg_column = c(length3_avg, length4_avg)
    sd_column = c(length3_sd, length4_sd)
    se_column = c(length3_se, length4_se)
    length = c(rep(3, length(length3_avg)), rep(4, length(length3_avg)))
    data = as.data.frame(cbind(avg_column, sd_column, se_column, time, length))
    colnames(data) = c("avg", "sd", "se", "time", "length")
    data$length = ifelse(data$length == 3, "length-3", "length-4")
    data$length = factor(data$length)
    
    # data p-values
    
    # Find segments with consecutive values below 0.05
    below_threshold <- p_values < 0.025
    start_points <- which(diff(c(FALSE, below_threshold)) == 1)
    end_points <- which(diff(c(below_threshold, FALSE)) == -1)
    y_pvalue = max(c(length3_subject_data$value, length4_subject_data$value))
    
    # Create segments data frame
    segments <- data.frame(
      x = data$time[start_points],
      xend = data$time[end_points],
      y = rep(y_pvalue, length(start_points)),  # You can adjust the y-value as needed
      yend = rep(y_pvalue, length(start_points))  # Same y-value for end points
    )
    segments$x = segments$x - 0.05
    segments$xend = segments$xend + 0.05
    
    ylim_metric = c(min(c(length3_subject_data$value, length4_subject_data$value)),
                    max(c(length3_subject_data$value, length4_subject_data$value)))
    
    plot = 
      ggplot() + 
      geom_vline(xintercept = 0, linetype = "dashed", color = "black", size = 0.3) +
      geom_segment(data = segments, aes(x = x, xend = xend, y = y, yend = yend), color = "black", size = 1) +
      geom_point(data = length3_subject_data, aes(x = time, y = value), color = "red", alpha = 0.3, size = 0.3) +
      geom_point(data = length4_subject_data, aes(x = time, y = value), color = "blue", alpha = 0.3, size = 0.3) +
      geom_line(data = length3_subject_data, aes(x = time, y = value, group = subject), color = "red", alpha = 0.2, size = 0.3) +
      geom_line(data = length4_subject_data, aes(x = time, y = value, group = subject), color = "blue", alpha = 0.2, size = 0.3) +
      geom_line(data = subset(data, length == "length-3"), aes(x = time, y = avg, color = "length-3"), size = 1.5) + 
      geom_line(data = subset(data, length == "length-4"), aes(x = time, y = avg, color = "length-4"), size = 1.5) + 
      geom_ribbon(data = subset(data, length == "length-3"), 
                  aes(x = time, ymin = avg - se, ymax = avg + se), alpha = 0.2, fill = "red") + 
      geom_ribbon(data = subset(data, length == "length-4"), 
                  aes(x = time, ymin = avg - se, ymax = avg + se), alpha = 0.2, fill = "blue") + 
      scale_color_manual(values = c("length-3" = "red", "length-4" = "blue")) + 
      labs(color = " ") +  # Set the legend title
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
      

    ggsave(paste(path_figures, "/connected_scatterplot_", metric_i, "_length3vs4_legend.png", sep = ""), plot, width = 3, height = 3, units = "in")
    
    
    plot = 
      ggplot() + 
      geom_vline(xintercept = 0, linetype = "dashed", color = "black", size = 0.3) +
      geom_segment(data = segments, aes(x = x, xend = xend, y = y, yend = yend), color = "black", size = 1) +
      geom_point(data = length3_subject_data, aes(x = time, y = value), color = "red", alpha = 0.3, size = 0.3) +
      geom_point(data = length4_subject_data, aes(x = time, y = value), color = "blue", alpha = 0.3, size = 0.3) +
      geom_line(data = length3_subject_data, aes(x = time, y = value, group = subject), color = "red", alpha = 0.2, size = 0.3) +
      geom_line(data = length4_subject_data, aes(x = time, y = value, group = subject), color = "blue", alpha = 0.2, size = 0.3) +
      geom_line(data = subset(data, length == "length-3"), aes(x = time, y = avg, color = "length-3"), size = 1.5) + 
      geom_line(data = subset(data, length == "length-4"), aes(x = time, y = avg, color = "length-4"), size = 1.5) + 
      geom_ribbon(data = subset(data, length == "length-3"), 
                  aes(x = time, ymin = avg - se, ymax = avg + se), alpha = 0.2, fill = "red") + 
      geom_ribbon(data = subset(data, length == "length-4"), 
                  aes(x = time, ymin = avg - se, ymax = avg + se), alpha = 0.2, fill = "blue") + 
      scale_color_manual(values = c("length-3" = "red", "length-4" = "blue")) + 
      labs(color = " ") +  # Set the legend title
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
    
    
    ggsave(paste(path_figures, "/connected_scatterplot_", metric_i, "_length3vs4.png", sep = ""), plot, width = 3, height = 3, units = "in")
    
    
    
    
    
    
    
    
    
    
    
      
            

  }
}
