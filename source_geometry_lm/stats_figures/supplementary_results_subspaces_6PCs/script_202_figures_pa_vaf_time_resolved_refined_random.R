##############################
## connected scatter plots  ##
##############################

segments_start <- seq(from = 1, to = 3601, by = 300)
segments_end <- seq(from = 300, to = 4000, by = 300)
delay_segments = cbind(segments_start, segments_end)

for (metric_i in metrics) {
  
  path = paste(path_results, "/group_results/", metric_i, sep = "")
  path_random = paste(path_results_random, "/group_results/", metric_i, sep = "")
  
  # between vs within
  if (metric_i == "principal_angle") {
    ylab_metric = "Principal angle (°)"
    p_values = p_values_all[1:31]
    p_values_random = p_values_all_random[1:31]
  } else if (metric_i == "principal_angle_min") {
    ylab_metric = "Principal angle min (°)"
    p_values = p_values_all[1:31]
    p_values_random = p_values_all_random[1:31]
  } else if (metric_i == "vaf") {
    ylab_metric = "VAF"
    p_values = p_values_all[32:62]
    p_values_random = p_values_all_random[32:62]
  }
  
  
  
  for (fun in functions_data) {
    for (perf in performance) {
      
      ### length 3
      
      length_seq = 3
      between_avg = c()
      between_sd = c()
      within_avg = c()
      within_sd = c()
      between_random_avg = c()
      between_random_sd = c()
      
      # stim
      
      time_window = paste(delay_segments[time_i,1], "to", delay_segments[time_i,2], sep = "")
      filename = paste(path, "/", fun, "/", perf,  "/design_matrix/between_within_stim_resolved_refined_length", length_seq, "_800to1100.txt", sep = "")
      design_matrix = read.table(filename, sep = ",", header = TRUE)
      colnames(design_matrix) = c("response", "Subspaces", "subject", "session")
      design_matrix$Subspaces = factor(ifelse(design_matrix$Subspaces == 1, "Between ranks", "Within rank"))
      data = subset(design_matrix, Subspaces == "Between ranks")
      between_avg = c(between_avg, mean(data$response))
      between_sd = c(between_sd, sd(data$response))
      data = subset(design_matrix, Subspaces == "Within rank")
      within_avg = c(within_avg, mean(data$response))
      within_sd = c(within_sd, sd(data$response))
      
      time_window = paste(delay_segments[time_i,1], "to", delay_segments[time_i,2], sep = "")
      filename = paste(path, "/", fun, "/", perf,  "/design_matrix/between_within_stim_resolved_refined_length", length_seq, "_1200to1500.txt", sep = "")
      design_matrix = read.table(filename, sep = ",", header = TRUE)
      colnames(design_matrix) = c("response", "Subspaces", "subject", "session")
      design_matrix$Subspaces = factor(ifelse(design_matrix$Subspaces == 1, "Between ranks", "Within rank"))
      data = subset(design_matrix, Subspaces == "Between ranks")
      between_avg = c(between_avg, mean(data$response))
      between_sd = c(between_sd, sd(data$response))
      data = subset(design_matrix, Subspaces == "Within rank")
      within_avg = c(within_avg, mean(data$response))
      within_sd = c(within_sd, sd(data$response))
      
      
      # delay
      
      for (time_i in 1:nrow(delay_segments)) {
        time_window = paste(delay_segments[time_i,1], "to", delay_segments[time_i,2], sep = "")
        filename = paste(path, "/", fun, "/", perf,  "/design_matrix/between_within_delay_length", length_seq, "_", time_window, ".txt", sep = "")
        design_matrix = read.table(filename, sep = ",", header = TRUE)
        colnames(design_matrix) = c("response", "Subspaces", "subject", "session")
        design_matrix$Subspaces = factor(ifelse(design_matrix$Subspaces == 1, "Between ranks", "Within rank"))
        data = subset(design_matrix, Subspaces == "Between ranks")
        between_avg = c(between_avg, mean(data$response))
        between_sd = c(between_sd, sd(data$response))
        data = subset(design_matrix, Subspaces == "Within rank")
        within_avg = c(within_avg, mean(data$response))
        within_sd = c(within_sd, sd(data$response))
      }
      
      
      
      # stim random
      
      time_window = paste(delay_segments[time_i,1], "to", delay_segments[time_i,2], sep = "")
      filename = paste(path_random, "/", fun, "/", perf,  "/design_matrix/between_stim_resolved_refined_length", length_seq, "_800to1100.txt", sep = "")
      design_matrix = read.table(filename, sep = ",", header = TRUE)
      colnames(design_matrix) = c("response", "Subspaces", "subject", "session")
      design_matrix$Subspaces = factor(ifelse(design_matrix$Subspaces == 1, "Between ranks", "Within rank"))
      data = subset(design_matrix, Subspaces == "Between ranks")
      between_random_avg = c(between_random_avg, mean(data$response))
      between_random_sd = c(between_random_sd, sd(data$response))

      time_window = paste(delay_segments[time_i,1], "to", delay_segments[time_i,2], sep = "")
      filename = paste(path_random, "/", fun, "/", perf,  "/design_matrix/between_stim_resolved_refined_length", length_seq, "_1200to1500.txt", sep = "")
      design_matrix = read.table(filename, sep = ",", header = TRUE)
      colnames(design_matrix) = c("response", "Subspaces", "subject", "session")
      design_matrix$Subspaces = factor(ifelse(design_matrix$Subspaces == 1, "Between ranks", "Within rank"))
      data = subset(design_matrix, Subspaces == "Between ranks")
      between_random_avg = c(between_random_avg, mean(data$response))
      between_random_sd = c(between_random_sd, sd(data$response))
      
      
      # delay random
      
      for (time_i in 1:nrow(delay_segments)) {
        time_window = paste(delay_segments[time_i,1], "to", delay_segments[time_i,2], sep = "")
        filename = paste(path_random, "/", fun, "/", perf,  "/design_matrix/between_delay_length", length_seq, "_", time_window, ".txt", sep = "")
        design_matrix = read.table(filename, sep = ",", header = TRUE)
        colnames(design_matrix) = c("response", "Subspaces", "subject", "session")
        design_matrix$Subspaces = factor(ifelse(design_matrix$Subspaces == 1, "Between ranks", "Within rank"))
        data = subset(design_matrix, Subspaces == "Between ranks")
        between_random_avg = c(between_random_avg, mean(data$response))
        between_random_sd = c(between_random_sd, sd(data$response))
      }
      
      
      
      # set in data frame
      time = c(-450, -150, delay_segments[,2] - 150) / 1000

      between_values = between_avg
      between_values_sd = between_sd
      between_values_se = between_values_sd / sqrt(15)
      
      within_values = within_avg
      within_values_sd = within_sd
      within_values_se = within_values_sd / sqrt(15)
      
      between_random = between_random_avg
      between_random_sd = between_random_sd
      between_random_se = between_random_sd / sqrt(15)
      
      data = as.data.frame(cbind(between_values, between_values_sd, between_values_se, within_values, within_values_sd, within_values_se, between_random, between_random_sd, between_random_se, time))

      # Using complete.cases()
      data <- data[complete.cases(data),]
      
      ## data p-values between vs within
      
      # Find segments with consecutive values below 0.05
      p_values_tmp = p_values[1:15]
      below_threshold <- p_values_tmp < 0.05
      start_points <- which(diff(c(FALSE, below_threshold)) == 1)
      end_points <- which(diff(c(below_threshold, FALSE)) == -1)
      y_pvalue = max(c(data$between_values + data$between_values_se, data$within_values + data$within_values_se, data$between_random + data$between_random_se)) + (max(c(data$between_values + data$between_values_se, data$within_values + data$within_values_se, data$between_random + data$between_random_se))/100)
      
      # Create segments data frame
      segments <- data.frame(
        x = data$time[start_points],
        xend = data$time[end_points],
        y = rep(y_pvalue, length(start_points)),  # You can adjust the y-value as needed
        yend = rep(y_pvalue, length(start_points))  # Same y-value for end points
      )
      segments$x = segments$x - 0.05
      segments$xend = segments$xend + 0.05
      
      ## data p-values between vs random
      
      # Find segments with consecutive values below 0.05
      p_values_tmp = p_values_random[1:15]
      below_threshold <- p_values_tmp < 0.05
      start_points <- which(diff(c(FALSE, below_threshold)) == 1)
      end_points <- which(diff(c(below_threshold, FALSE)) == -1)
      y_pvalue = y_pvalue + (y_pvalue/50)
      
      # Create segments data frame
      segments_random <- data.frame(
        x = data$time[start_points],
        xend = data$time[end_points],
        y = rep(y_pvalue, length(start_points)),  # You can adjust the y-value as needed
        yend = rep(y_pvalue, length(start_points))  # Same y-value for end points
      )
      segments_random$x = segments_random$x - 0.05
      segments_random$xend = segments_random$xend + 0.05
      
      ylim_metric = c(min(c(between_values - between_values_se, within_values - within_values_se, between_random - between_random_se)) - min(c(between_values - between_values_se, within_values - within_values_se, between_random - between_random_se))/100,
                      max(c(segments$y, segments_random$y)) + max(c(segments$y, segments_random$y))/20)

      # Plot between and within with standard errors
      
      plot = 
        ggplot() +
        geom_vline(xintercept = 0, linetype = "dashed", color = "black", size = 0.3) +
        geom_segment(data = segments, aes(x = x, xend = xend, y = y, yend = yend), color = "black", size = 1) +
        geom_segment(data = segments_random, aes(x = x, xend = xend, y = y, yend = yend), color = "grey", size = 1) +
        geom_line(data = data, aes(x = time, y = between_values, color = "Between"), linetype = "solid", size = 1.5) +
        geom_line(data = data, aes(x = time, y = within_values, color = "Within"), linetype = "solid", size = 1.5) +
        geom_line(data = data, aes(x = time, y = between_random, color = "Surrogate"), linetype = "solid", size = 1.5) +
        geom_ribbon(data = data,
                    aes(x = time, ymin = between_values - between_values_se, ymax = between_values + between_values_se),
                    alpha = 0.2, fill = "red") +
        geom_ribbon(data = data,
                    aes(x = time, ymin = within_values - within_values_se, ymax = within_values + within_values_se),
                    alpha = 0.2, fill = "blue") +
        geom_ribbon(data = data,
                    aes(x = time, ymin = between_random - between_random_se, ymax = between_random + between_random_se),
                    alpha = 0.2, fill = "green") +
        ylab(ylab_metric) +
        xlab("Time (s)") +
        scale_color_manual(values = c("Between" = "red", "Within" = "blue", "Surrogate" = "green"), breaks = c("Between", "Within", "Surrogate")) + 
        scale_x_continuous(breaks = seq(0, 4, by = 1), limits = c(min(time) - 0.05, 4)) +
        scale_y_continuous(limits = ylim_metric) + 
        theme(plot.title = element_text(hjust = 0.5), 
              legend.position = "none", 
              panel.background = element_blank(),
              plot.background = element_blank(),
              axis.line = element_line(color = "black"),
              axis.title = element_text(color = "black", size = 12),
              axis.text.x = element_text(color = "black", size = 12),
              axis.text.y = element_text(color = "black", size = 12),
              panel.grid.major = element_blank(),
              panel.grid.minor = element_blank()) +
        labs(color = "")
      
      ggsave(paste(path_figures, "/connected_scatterplot_", metric_i, "_between_vs_within_length", length_seq, ".png", sep = ""), plot, width = 3, height = 3, units = "in")
      
      plot = 
        ggplot() +
        geom_vline(xintercept = 0, linetype = "dashed", color = "black", size = 0.3) +
        geom_segment(data = segments, aes(x = x, xend = xend, y = y, yend = yend), color = "black", size = 1) +
        geom_segment(data = segments_random, aes(x = x, xend = xend, y = y, yend = yend), color = "grey", size = 1) +
        geom_line(data = data, aes(x = time, y = between_values, color = "Between"), linetype = "solid", size = 1.5) +
        geom_line(data = data, aes(x = time, y = within_values, color = "Within"), linetype = "solid", size = 1.5) +
        geom_line(data = data, aes(x = time, y = between_random, color = "Surrogate"), linetype = "solid", size = 1.5) +
        geom_ribbon(data = data,
                    aes(x = time, ymin = between_values - between_values_se, ymax = between_values + between_values_se),
                    alpha = 0.2, fill = "red") +
        geom_ribbon(data = data,
                    aes(x = time, ymin = within_values - within_values_se, ymax = within_values + within_values_se),
                    alpha = 0.2, fill = "blue") +
        geom_ribbon(data = data,
                    aes(x = time, ymin = between_random - between_random_se, ymax = between_random + between_random_se),
                    alpha = 0.2, fill = "green") +
        ylab(ylab_metric) +
        xlab("Time (s)") +
        scale_color_manual(values = c("Between" = "red", "Within" = "blue", "Surrogate" = "green"), breaks = c("Between", "Within", "Surrogate")) + 
        scale_x_continuous(breaks = seq(0, 4, by = 1), limits = c(min(time) - 0.05, 4)) +
        scale_y_continuous(limits = ylim_metric) + 
        theme(plot.title = element_text(hjust = 0.5), 
              legend.position = "bottom", 
              panel.background = element_blank(),
              plot.background = element_blank(),
              axis.line = element_line(color = "black"),
              axis.title = element_text(color = "black", size = 12),
              axis.text.x = element_text(color = "black", size = 12),
              axis.text.y = element_text(color = "black", size = 12),
              panel.grid.major = element_blank(),
              panel.grid.minor = element_blank()) +
        guides(linetype = "none", color = guide_legend(title = ""))  # <- this line is key
      
      ggsave(paste(path_figures, "/connected_scatterplot_", metric_i, "_between_vs_within_length", length_seq, "_legend.png", sep = ""), plot, width = 3, height = 3, units = "in")
      
      
      
      
      
      
      
      
      
      
      
      
      ### length 4
      length_seq = 4
      between_avg = c()
      between_sd = c()
      within_avg = c()
      within_sd = c()
      between_random_avg = c()
      between_random_sd = c()

      # stim
      
      time_window = paste(delay_segments[time_i,1], "to", delay_segments[time_i,2], sep = "")
      filename = paste(path, "/", fun, "/", perf,  "/design_matrix/between_within_stim_resolved_refined_length", length_seq, "_400to700.txt", sep = "")
      design_matrix = read.table(filename, sep = ",", header = TRUE)
      colnames(design_matrix) = c("response", "Subspaces", "subject", "session")
      design_matrix$Subspaces = factor(ifelse(design_matrix$Subspaces == 1, "Between ranks", "Within rank"))
      data = subset(design_matrix, Subspaces == "Between ranks")
      between_avg = c(between_avg, mean(data$response))
      between_sd = c(between_sd, sd(data$response))
      data = subset(design_matrix, Subspaces == "Within rank")
      within_avg = c(within_avg, mean(data$response))
      within_sd = c(within_sd, sd(data$response))
      
      time_window = paste(delay_segments[time_i,1], "to", delay_segments[time_i,2], sep = "")
      filename = paste(path, "/", fun, "/", perf,  "/design_matrix/between_within_stim_resolved_refined_length", length_seq, "_800to1100.txt", sep = "")
      design_matrix = read.table(filename, sep = ",", header = TRUE)
      colnames(design_matrix) = c("response", "Subspaces", "subject", "session")
      design_matrix$Subspaces = factor(ifelse(design_matrix$Subspaces == 1, "Between ranks", "Within rank"))
      data = subset(design_matrix, Subspaces == "Between ranks")
      between_avg = c(between_avg, mean(data$response))
      between_sd = c(between_sd, sd(data$response))
      data = subset(design_matrix, Subspaces == "Within rank")
      within_avg = c(within_avg, mean(data$response))
      within_sd = c(within_sd, sd(data$response))
      
      time_window = paste(delay_segments[time_i,1], "to", delay_segments[time_i,2], sep = "")
      filename = paste(path, "/", fun, "/", perf,  "/design_matrix/between_within_stim_resolved_refined_length", length_seq, "_1200to1500.txt", sep = "")
      design_matrix = read.table(filename, sep = ",", header = TRUE)
      colnames(design_matrix) = c("response", "Subspaces", "subject", "session")
      design_matrix$Subspaces = factor(ifelse(design_matrix$Subspaces == 1, "Between ranks", "Within rank"))
      data = subset(design_matrix, Subspaces == "Between ranks")
      between_avg = c(between_avg, mean(data$response))
      between_sd = c(between_sd, sd(data$response))
      data = subset(design_matrix, Subspaces == "Within rank")
      within_avg = c(within_avg, mean(data$response))
      within_sd = c(within_sd, sd(data$response))
      
      # delay

      for (time_i in 1:nrow(delay_segments)) {
        time_window = paste(delay_segments[time_i,1], "to", delay_segments[time_i,2], sep = "")
        filename = paste(path, "/", fun, "/", perf,  "/design_matrix/between_within_delay_length", length_seq, "_", time_window, ".txt", sep = "")
        design_matrix = read.table(filename, sep = ",", header = TRUE)
        colnames(design_matrix) = c("response", "Subspaces", "subject", "session")
        design_matrix$Subspaces = factor(ifelse(design_matrix$Subspaces == 1, "Between ranks", "Within rank"))
        data = subset(design_matrix, Subspaces == "Between ranks")
        between_avg = c(between_avg, mean(data$response))
        between_sd = c(between_sd, sd(data$response))
        data = subset(design_matrix, Subspaces == "Within rank")
        within_avg = c(within_avg, mean(data$response))
        within_sd = c(within_sd, sd(data$response))
      }
      
      
      # stim random
      
      time_window = paste(delay_segments[time_i,1], "to", delay_segments[time_i,2], sep = "")
      filename = paste(path_random, "/", fun, "/", perf,  "/design_matrix/between_stim_resolved_refined_length", length_seq, "_400to700.txt", sep = "")
      design_matrix = read.table(filename, sep = ",", header = TRUE)
      colnames(design_matrix) = c("response", "Subspaces", "subject", "session")
      design_matrix$Subspaces = factor(ifelse(design_matrix$Subspaces == 1, "Between ranks", "Within rank"))
      data = subset(design_matrix, Subspaces == "Between ranks")
      between_random_avg = c(between_random_avg, mean(data$response))
      between_random_sd = c(between_random_sd, sd(data$response))

      time_window = paste(delay_segments[time_i,1], "to", delay_segments[time_i,2], sep = "")
      filename = paste(path_random, "/", fun, "/", perf,  "/design_matrix/between_stim_resolved_refined_length", length_seq, "_800to1100.txt", sep = "")
      design_matrix = read.table(filename, sep = ",", header = TRUE)
      colnames(design_matrix) = c("response", "Subspaces", "subject", "session")
      design_matrix$Subspaces = factor(ifelse(design_matrix$Subspaces == 1, "Between ranks", "Within rank"))
      data = subset(design_matrix, Subspaces == "Between ranks")
      between_random_avg = c(between_random_avg, mean(data$response))
      between_random_sd = c(between_random_sd, sd(data$response))
      
      time_window = paste(delay_segments[time_i,1], "to", delay_segments[time_i,2], sep = "")
      filename = paste(path_random, "/", fun, "/", perf,  "/design_matrix/between_stim_resolved_refined_length", length_seq, "_1200to1500.txt", sep = "")
      design_matrix = read.table(filename, sep = ",", header = TRUE)
      colnames(design_matrix) = c("response", "Subspaces", "subject", "session")
      design_matrix$Subspaces = factor(ifelse(design_matrix$Subspaces == 1, "Between ranks", "Within rank"))
      data = subset(design_matrix, Subspaces == "Between ranks")
      between_random_avg = c(between_random_avg, mean(data$response))
      between_random_sd = c(between_random_sd, sd(data$response))
      
      # delay random
      
      for (time_i in 1:nrow(delay_segments)) {
        time_window = paste(delay_segments[time_i,1], "to", delay_segments[time_i,2], sep = "")
        filename = paste(path_random, "/", fun, "/", perf,  "/design_matrix/between_delay_length", length_seq, "_", time_window, ".txt", sep = "")
        design_matrix = read.table(filename, sep = ",", header = TRUE)
        colnames(design_matrix) = c("response", "Subspaces", "subject", "session")
        design_matrix$Subspaces = factor(ifelse(design_matrix$Subspaces == 1, "Between ranks", "Within rank"))
        data = subset(design_matrix, Subspaces == "Between ranks")
        between_random_avg = c(between_random_avg, mean(data$response))
        between_random_sd = c(between_random_sd, sd(data$response))
      }
      
      
      
      # set in data frame
      time = c(-750, -450, -150, delay_segments[,2] - 150) / 1000

      between_values = between_avg
      between_values_sd = between_sd
      between_values_se = between_values_sd / sqrt(15)
      
      within_values = within_avg
      within_values_sd = within_sd
      within_values_se = within_values_sd / sqrt(15)
      
      between_random = between_random_avg
      between_random_sd = between_random_sd
      between_random_se = between_random_sd / sqrt(15)
      
      data = as.data.frame(cbind(between_values, between_values_sd, between_values_se, within_values, within_values_sd, within_values_se, between_random, between_random_sd, between_random_se, time))
      
      # Using complete.cases()
      data <- data[complete.cases(data),]
      
      ## data p-values between vs within
      
      # Find segments with consecutive values below 0.05
      p_values_tmp = p_values[16:31]
      below_threshold <- p_values_tmp < 0.05
      start_points <- which(diff(c(FALSE, below_threshold)) == 1)
      end_points <- which(diff(c(below_threshold, FALSE)) == -1)
      y_pvalue = max(c(data$between_values + data$between_values_se, data$within_values + data$within_values_se, data$between_random + data$between_random_se)) + (max(c(data$between_values + data$between_values_se, data$within_values + data$within_values_se, data$between_random + data$between_random_se))/100)
      
      # Create segments data frame
      segments <- data.frame(
        x = data$time[start_points],
        xend = data$time[end_points],
        y = rep(y_pvalue, length(start_points)),  # You can adjust the y-value as needed
        yend = rep(y_pvalue, length(start_points))  # Same y-value for end points
      )
      segments$x = segments$x - 0.05
      segments$xend = segments$xend + 0.05
      
      ## data p-values between vs random
      
      # Find segments with consecutive values below 0.05
      p_values_tmp = p_values_random[16:31]
      below_threshold <- p_values_tmp < 0.05
      start_points <- which(diff(c(FALSE, below_threshold)) == 1)
      end_points <- which(diff(c(below_threshold, FALSE)) == -1)
      y_pvalue = y_pvalue + (y_pvalue/50)
      
      # Create segments data frame
      segments_random <- data.frame(
        x = data$time[start_points],
        xend = data$time[end_points],
        y = rep(y_pvalue, length(start_points)),  # You can adjust the y-value as needed
        yend = rep(y_pvalue, length(start_points))  # Same y-value for end points
      )
      segments_random$x = segments_random$x - 0.05
      segments_random$xend = segments_random$xend + 0.05
      
      ylim_metric = c(min(c(between_values - between_values_se, within_values - within_values_se, between_random - between_random_se)) - min(c(between_values - between_values_se, within_values - within_values_se, between_random - between_random_se))/100,
                      max(c(segments$y, segments_random$y)) + max(c(segments$y, segments_random$y))/20)
      
      # Plot between and within with standard errors
      
      plot = 
        ggplot() +
        geom_vline(xintercept = 0, linetype = "dashed", color = "black", size = 0.3) +
        geom_segment(data = segments, aes(x = x, xend = xend, y = y, yend = yend), color = "black", size = 1) +
        geom_segment(data = segments_random, aes(x = x, xend = xend, y = y, yend = yend), color = "grey", size = 1) +
        geom_line(data = data, aes(x = time, y = between_values, color = "Between"), linetype = "solid", size = 1.5) +
        geom_line(data = data, aes(x = time, y = within_values, color = "Within"), linetype = "solid", size = 1.5) +
        geom_line(data = data, aes(x = time, y = between_random, color = "Surrogate"), linetype = "solid", size = 1.5) +
        geom_ribbon(data = data,
                    aes(x = time, ymin = between_values - between_values_se, ymax = between_values + between_values_se),
                    alpha = 0.2, fill = "red") +
        geom_ribbon(data = data,
                    aes(x = time, ymin = within_values - within_values_se, ymax = within_values + within_values_se),
                    alpha = 0.2, fill = "blue") +
        geom_ribbon(data = data,
                    aes(x = time, ymin = between_random - between_random_se, ymax = between_random + between_random_se),
                    alpha = 0.2, fill = "green") +
        ylab(ylab_metric) +
        xlab("Time (s)") +
        scale_color_manual(values = c("Between" = "red", "Within" = "blue", "Surrogate" = "green"), breaks = c("Between", "Within", "Surrogate")) + 
        scale_x_continuous(breaks = seq(0, 4, by = 1), limits = c(min(time) - 0.05, 4)) +
        scale_y_continuous(limits = ylim_metric) + 
        theme(plot.title = element_text(hjust = 0.5), 
              legend.position = "none", 
              panel.background = element_blank(),
              plot.background = element_blank(),
              axis.line = element_line(color = "black"),
              axis.title = element_text(color = "black", size = 12),
              axis.text.x = element_text(color = "black", size = 12),
              axis.text.y = element_text(color = "black", size = 12),
              panel.grid.major = element_blank(),
              panel.grid.minor = element_blank()) +
        labs(color = "")
      
      ggsave(paste(path_figures, "/connected_scatterplot_", metric_i, "_between_vs_within_length", length_seq, ".png", sep = ""), plot, width = 3, height = 3, units = "in")
      
      plot = 
        ggplot() +
        geom_vline(xintercept = 0, linetype = "dashed", color = "black", size = 0.3) +
        geom_segment(data = segments, aes(x = x, xend = xend, y = y, yend = yend), color = "black", size = 1) +
        geom_segment(data = segments_random, aes(x = x, xend = xend, y = y, yend = yend), color = "grey", size = 1) +
        geom_line(data = data, aes(x = time, y = between_values, color = "Between"), linetype = "solid", size = 1.5) +
        geom_line(data = data, aes(x = time, y = within_values, color = "Within"), linetype = "solid", size = 1.5) +
        geom_line(data = data, aes(x = time, y = between_random, color = "Surrogate"), linetype = "solid", size = 1.5) +
        geom_ribbon(data = data,
                    aes(x = time, ymin = between_values - between_values_se, ymax = between_values + between_values_se),
                    alpha = 0.2, fill = "red") +
        geom_ribbon(data = data,
                    aes(x = time, ymin = within_values - within_values_se, ymax = within_values + within_values_se),
                    alpha = 0.2, fill = "blue") +
        geom_ribbon(data = data,
                    aes(x = time, ymin = between_random - between_random_se, ymax = between_random + between_random_se),
                    alpha = 0.2, fill = "green") +
        ylab(ylab_metric) +
        xlab("Time (s)") +
        scale_color_manual(values = c("Between" = "red", "Within" = "blue", "Surrogate" = "green"), breaks = c("Between", "Within", "Surrogate")) + 
        scale_x_continuous(breaks = seq(0, 4, by = 1), limits = c(min(time) - 0.05, 4)) +
        scale_y_continuous(limits = ylim_metric) + 
        theme(plot.title = element_text(hjust = 0.5), 
              legend.position = "bottom", 
              panel.background = element_blank(),
              plot.background = element_blank(),
              axis.line = element_line(color = "black"),
              axis.title = element_text(color = "black", size = 12),
              axis.text.x = element_text(color = "black", size = 12),
              axis.text.y = element_text(color = "black", size = 12),
              panel.grid.major = element_blank(),
              panel.grid.minor = element_blank()) +
        guides(linetype = "none", color = guide_legend(title = ""))  # <- this line is key
      
      ggsave(paste(path_figures, "/connected_scatterplot_", metric_i, "_between_vs_within_length", length_seq, "_legend.png", sep = ""), plot, width = 3, height = 3, units = "in")
      
      
      
      
      
      
            

          
      
            
    }
  }
}
