##############################
## connected scatter plots  ##
##############################

segments_start <- seq(from = 1, to = 3601, by = 300)
segments_end <- seq(from = 300, to = 4000, by = 300)
delay_segments = cbind(segments_start, segments_end)

for (metric_i in metrics) {
  
  path = paste(path_results, "/group_results/", metric_i, sep = "")
  
  if (metric_i == "principal_angle") {
    ylim_metric = c(10,90)
    ylab_metric = "Principal angle (°)"
    p_values = p_values_all[1:15]
  } else if (metric_i == "vaf") {
    ylim_metric = c(0,1)
    ylab_metric = "VAF"
    p_values = p_values_all[16:30]
  }
  
  for (fun in functions_data) {

    length3_avg = c()
    length3_sd = c()
    length4_avg = c()
    length4_sd = c()
    
    # stim
    
    filename = paste(path, "/", fun, "/correct_trials/design_matrix/between_within_stim_resolved_refined_length3_800to1100.txt", sep = "")
    design_matrix = read.table(filename, sep = ",", header = TRUE)
    colnames(design_matrix) = c("response", "Subspaces", "subject", "session")
    design_matrix$Subspaces = factor(ifelse(design_matrix$Subspaces == 1, "Between ranks", "Within rank"))
    data = subset(design_matrix, Subspaces == "Between ranks")
    length3_avg = c(length3_avg, mean(data$response))
    length3_sd = c(length3_sd, sd(data$response))
    filename = paste(path, "/", fun, "/correct_trials/design_matrix/between_within_stim_resolved_refined_length4_400to700.txt", sep = "")
    design_matrix = read.table(filename, sep = ",", header = TRUE)
    colnames(design_matrix) = c("response", "Subspaces", "subject", "session")
    design_matrix$Subspaces = factor(ifelse(design_matrix$Subspaces == 1, "Between ranks", "Within rank"))
    data = subset(design_matrix, Subspaces == "Between ranks")
    length4_avg = c(length4_avg, mean(data$response))
    length4_sd = c(length4_sd, sd(data$response))
    
    filename = paste(path, "/", fun, "/correct_trials/design_matrix/between_within_stim_resolved_refined_length3_1200to1500.txt", sep = "")
    design_matrix = read.table(filename, sep = ",", header = TRUE)
    colnames(design_matrix) = c("response", "Subspaces", "subject", "session")
    design_matrix$Subspaces = factor(ifelse(design_matrix$Subspaces == 1, "Between ranks", "Within rank"))
    data = subset(design_matrix, Subspaces == "Between ranks")
    length3_avg = c(length3_avg, mean(data$response))
    length3_sd = c(length3_sd, sd(data$response))
    filename = paste(path, "/", fun, "/correct_trials/design_matrix/between_within_stim_resolved_refined_length4_800to1100.txt", sep = "")
    design_matrix = read.table(filename, sep = ",", header = TRUE)
    colnames(design_matrix) = c("response", "Subspaces", "subject", "session")
    design_matrix$Subspaces = factor(ifelse(design_matrix$Subspaces == 1, "Between ranks", "Within rank"))
    data = subset(design_matrix, Subspaces == "Between ranks")
    length4_avg = c(length4_avg, mean(data$response))
    length4_sd = c(length4_sd, sd(data$response))
    
    
    # delay
    
    for (time_i in 1:nrow(delay_segments)) {
      time_window = paste(delay_segments[time_i,1], "to", delay_segments[time_i,2], sep = "")
      
      # length-3
      filename = paste(path, "/", fun, "/correct_trials/design_matrix/between_within_delay_length3_", time_window, ".txt", sep = "")
      design_matrix = read.table(filename, sep = ",", header = TRUE)
      colnames(design_matrix) = c("response", "Subspaces", "subject", "session")
      design_matrix$Subspaces = factor(ifelse(design_matrix$Subspaces == 1, "Between ranks", "Within rank"))
      data = subset(design_matrix, Subspaces == "Between ranks")
      length3_avg = c(length3_avg, mean(data$response))
      length3_sd = c(length3_sd, sd(data$response))
      
      # length-4
      filename = paste(path, "/", fun, "/correct_trials/design_matrix/between_within_delay_length4_", time_window, ".txt", sep = "")
      design_matrix = read.table(filename, sep = ",", header = TRUE)
      colnames(design_matrix) = c("response", "Subspaces", "subject", "session")
      design_matrix$Subspaces = factor(ifelse(design_matrix$Subspaces == 1, "Between ranks", "Within rank"))
      data = subset(design_matrix, Subspaces == "Between ranks")
      length4_avg = c(length4_avg, mean(data$response))
      length4_sd = c(length4_sd, sd(data$response))
    }
    
    # data
    time = c(-450, -150, delay_segments[,2] - 150) / 1000
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
    below_threshold <- p_values < 0.05
    start_points <- which(diff(c(FALSE, below_threshold)) == 1)
    end_points <- which(diff(c(below_threshold, FALSE)) == -1)
    y_pvalue = max(data$avg + data$se) + (max(data$avg + data$se)/40)
    
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
