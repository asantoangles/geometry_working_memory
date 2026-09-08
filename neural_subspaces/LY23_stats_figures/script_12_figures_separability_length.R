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
    p_values = p_values_all[1:15]
    separability_name = "Distance"  
    y_axis_step = 0.1
    }
  if (separability == "volume") {
    p_values = p_values_all[16:30]
    separability_name = "Volume"  
    y_axis_step = 0.5
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
    matrix_design = read.table(filename, sep = ",", header = TRUE)
    data = matrix_design
    colnames(data) = c("response", "condition", "subject", "session")
    length3_avg = c(length3_avg, mean(data$response))
    length3_sd = c(length3_sd, sd(data$response))
    data_length3 = data
    filename = paste(path, "/", fun, "/correct_trials/design_matrix/stim_length4_400to700.txt", sep = "")
    matrix_design = read.table(filename, sep = ",", header = TRUE)
    data = matrix_design
    colnames(data) = c("response", "condition", "subject", "session")
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
    matrix_design = read.table(filename, sep = ",", header = TRUE)
    data = matrix_design
    colnames(data) = c("response", "condition", "subject", "session")
    length3_avg = c(length3_avg, mean(data$response))
    length3_sd = c(length3_sd, sd(data$response))
    data_length3 = data
    filename = paste(path, "/", fun, "/correct_trials/design_matrix/stim_length4_800to1100.txt", sep = "")
    matrix_design = read.table(filename, sep = ",", header = TRUE)
    data = matrix_design
    colnames(data) = c("response", "condition", "subject", "session")
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
      
      # length 3
      filename = paste(path, "/", fun, "/correct_trials/design_matrix/delay_length3_", time_window, ".txt", sep = "")
      matrix_design = read.table(filename, sep = ",", header = TRUE)
      data = matrix_design
      colnames(data) = c("response", "condition", "subject", "session")
      length3_avg = c(length3_avg, mean(data$response))
      length3_sd = c(length3_sd, sd(data$response))
      data_length3 = data

      # length 4
      filename = paste(path, "/", fun, "/correct_trials/design_matrix/delay_length4_", time_window, ".txt", sep = "")
      matrix_design = read.table(filename, sep = ",", header = TRUE)
      data = matrix_design
      colnames(data) = c("response", "condition", "subject", "session")
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
    length = c(rep("length-3", length(length3_avg)), rep("length-4", length(length3_avg)))
    data = as.data.frame(cbind(avg_column, sd_column, se_column, time, length))
    colnames(data) = c("avg", "sd", "se", "time", "length")
    data$length = factor(data$length)
    data$avg = as.numeric(data$avg)
    data$sd = as.numeric(data$sd)
    data$se = as.numeric(data$se)
    data$time = as.numeric(data$time)
    
    # data p-values
    
    # Find segments with consecutive values below 0.05
    p_values_tmp = p_values[1:15]
    below_threshold <- p_values_tmp < 0.025
    start_points <- which(diff(c(FALSE, below_threshold)) == 1)
    end_points <- which(diff(c(below_threshold, FALSE)) == -1)

    y_pvalue = max(c(length3_subject_data$value, length3_subject_data$value))
    
    # Create segments data frame
    segments <- data.frame(
      x = data$time[start_points],
      xend = data$time[end_points],
      y = rep(y_pvalue, length(start_points)),  # You can adjust the y-value as needed
      yend = rep(y_pvalue, length(start_points))  # Same y-value for end points
    )
    segments$x = segments$x - 0.05
    segments$xend = segments$xend + 0.05
    
    ylim_metric = c(min(c(length3_subject_data$value, length4_subject_data$value, min(data$avg - data$se))),
                    max(c(length3_subject_data$value, length4_subject_data$value, max(data$avg + data$se))))
    
    plot = 
      ggplot() +
      geom_vline(xintercept = 0, linetype = "dashed", color = "black", size = 0.3) +
      geom_segment(data = segments, aes(x = x, xend = xend, y = y, yend = yend), color = "black", size = 1) +
      geom_point(data = length3_subject_data, aes(x = time, y = value), color = "#F8766D", alpha = 0.3, size = 0.3) +
      geom_point(data = length4_subject_data, aes(x = time, y = value), color = "#00BFC4", alpha = 0.3, size = 0.3) +
      geom_line(data = length3_subject_data, aes(x = time, y = value, group = subject), color = "#F8766D", alpha = 0.2, size = 0.3) +
      geom_line(data = length4_subject_data, aes(x = time, y = value, group = subject), color = "#00BFC4", alpha = 0.2, size = 0.3) +
      geom_line(data = subset(data, length == "length-3"), aes(x = time, y = avg, color = "length-3"), size = 1.5) + 
      geom_line(data = subset(data, length == "length-4"), aes(x = time, y = avg, color = "length-4"), size = 1.5) + 
      geom_ribbon(data = subset(data, length == "length-3"), 
                  aes(x = time, ymin = avg - se, ymax = avg + se), alpha = 0.2, fill = "#F8766D") + 
      geom_ribbon(data = subset(data, length == "length-4"), 
                  aes(x = time, ymin = avg - se, ymax = avg + se), alpha = 0.2, fill = "#00BFC4") + 
      scale_color_manual(values = c("length-3" = "#F8766D", "length-4" = "#00BFC4")) + 
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
    
    
    ggsave(paste(path_figures, "/connected_scatterplot_separability_", separability, "_length3vs4.png", sep = ""), plot, width = 3, height = 3, units = "in")
    
    
    plot = 
      ggplot() +
      geom_vline(xintercept = 0, linetype = "dashed", color = "black", size = 0.3) +
      geom_segment(data = segments, aes(x = x, xend = xend, y = y, yend = yend), color = "black", size = 1) +
      geom_point(data = length3_subject_data, aes(x = time, y = value), color = "#F8766D", alpha = 0.3, size = 0.3) +
      geom_point(data = length4_subject_data, aes(x = time, y = value), color = "#00BFC4", alpha = 0.3, size = 0.3) +
      geom_line(data = length3_subject_data, aes(x = time, y = value, group = subject), color = "#F8766D", alpha = 0.2, size = 0.3) +
      geom_line(data = length4_subject_data, aes(x = time, y = value, group = subject), color = "#00BFC4", alpha = 0.2, size = 0.3) +
      geom_line(data = subset(data, length == "length-3"), aes(x = time, y = avg, color = "length-3"), size = 1.5) + 
      geom_line(data = subset(data, length == "length-4"), aes(x = time, y = avg, color = "length-4"), size = 1.5) + 
      geom_ribbon(data = subset(data, length == "length-3"), 
                  aes(x = time, ymin = avg - se, ymax = avg + se), alpha = 0.2, fill = "#F8766D") + 
      geom_ribbon(data = subset(data, length == "length-4"), 
                  aes(x = time, ymin = avg - se, ymax = avg + se), alpha = 0.2, fill = "#00BFC4") + 
      scale_color_manual(values = c("length-3" = "#F8766D", "length-4" = "#00BFC4")) + 
      labs(color = " ") +  # Set the legend title
      xlab("Time (s)") + ylab(separability_name) +
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
    
    
    ggsave(paste(path_figures, "/connected_scatterplot_separability_", separability, "_length3vs4_legend.png", sep = ""), plot, width = 3, height = 3, units = "in")
    
  }
}
