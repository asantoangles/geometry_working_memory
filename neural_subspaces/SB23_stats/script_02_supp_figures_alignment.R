##############################
## connected scatter plots  ##
##############################

segments_start <- seq(from = 1, to = 3601, by = 300)
segments_end <- seq(from = 300, to = 4000, by = 300)
delay_segments = cbind(segments_start, segments_end)

for (metric_i in metrics) {
  
  # between vs within
  if (metric_i == "principal_angle") {
    ylab_metric = "Principal angle (°)"
    y_axis_step = 10
  } else if (metric_i == "vaf") {
    ylab_metric = "VAF"
    y_axis_step = 0.1
  }
  
  
  
  for (fun in functions_data) {
    for (perf in performance) {
      
      ### length 3
      
      length_seq = 3
      
      between_supp <- list()
      between_se_supp <- list()
      
      for (supp_i in supplementaries) {
        
        between_avg = c()
        between_sd = c()
        
        path = paste(path_results, "/supp", supp_i, "/group_results/", metric_i, sep = "")
        
        # stim
        
        time_window = paste(delay_segments[time_i,1], "to", delay_segments[time_i,2], sep = "")
        filename = paste(path, "/", fun, "/", perf,  "/design_matrix/stim_length", length_seq, "_400to700.txt", sep = "")
        design_matrix = read.table(filename, sep = ",", header = TRUE)
        colnames(design_matrix) = c("response", "Subspaces", "subject", "session")
        design_matrix$Subspaces = factor(ifelse(design_matrix$Subspaces == 1, "Between ranks", "Within rank"))
        data = subset(design_matrix, Subspaces == "Between ranks")
        between_avg = c(between_avg, mean(data$response))
        between_sd = c(between_sd, sd(data$response))
        
        time_window = paste(delay_segments[time_i,1], "to", delay_segments[time_i,2], sep = "")
        filename = paste(path, "/", fun, "/", perf,  "/design_matrix/stim_length", length_seq, "_800to1100.txt", sep = "")
        design_matrix = read.table(filename, sep = ",", header = TRUE)
        colnames(design_matrix) = c("response", "Subspaces", "subject", "session")
        design_matrix$Subspaces = factor(ifelse(design_matrix$Subspaces == 1, "Between ranks", "Within rank"))
        data = subset(design_matrix, Subspaces == "Between ranks")
        between_avg = c(between_avg, mean(data$response))
        between_sd = c(between_sd, sd(data$response))
        
        
        # delay
        
        for (time_i in 1:nrow(delay_segments)) {
          time_window = paste(delay_segments[time_i,1], "to", delay_segments[time_i,2], sep = "")
          filename = paste(path, "/", fun, "/", perf,  "/design_matrix/delay_length", length_seq, "_", time_window, ".txt", sep = "")
          design_matrix = read.table(filename, sep = ",", header = TRUE)
          colnames(design_matrix) = c("response", "Subspaces", "subject", "session")
          design_matrix$Subspaces = factor(ifelse(design_matrix$Subspaces == 1, "Between ranks", "Within rank"))
          data = subset(design_matrix, Subspaces == "Between ranks")
          between_avg = c(between_avg, mean(data$response))
          between_sd = c(between_sd, sd(data$response))
        }
        
        between_supp[[supp_i]] <- between_avg
        between_se_supp[[supp_i]] <- between_sd / sqrt(15)     
        
      }
      
      # set in data frame
      time = c(-450, -150, delay_segments[,2] - 150) / 1000
      data = as.data.frame(cbind(between_supp[[1]], between_supp[[2]], between_supp[[3]], between_supp[[4]], between_supp[[5]],
                                 between_se_supp[[1]], between_se_supp[[2]], between_se_supp[[3]], between_se_supp[[4]], between_se_supp[[5]],
                                 time))
      colnames(data) = c("between_supp1", "between_supp2", "between_supp3", "between_supp4", "between_supp5", 
                         "between_se_supp1", "between_se_supp2", "between_se_supp3", "between_se_supp4", "between_se_supp5", 
                         "time")
      
      # Plot between and within with standard errors
      
      plot =
        ggplot() +
        geom_vline(xintercept = 0, linetype = "dashed", color = "black", size = 0.3) +
        
        # Lines
        geom_line(data = data, aes(x = time, y = between_supp1, color = "supp1"), size = 1) +
        geom_line(data = data, aes(x = time, y = between_supp2, color = "supp2"), size = 1) +
        geom_line(data = data, aes(x = time, y = between_supp3, color = "supp3"), size = 1) +
        geom_line(data = data, aes(x = time, y = between_supp4, color = "supp4"), size = 1) +
        geom_line(data = data, aes(x = time, y = between_supp5, color = "supp5"), size = 1) +
        
        # Ribbons
        geom_ribbon(data = data, aes(x = time,
                                     ymin = between_supp1 - between_se_supp1,
                                     ymax = between_supp1 + between_se_supp1,
                                     fill = "supp1"), alpha = 0.2) +
        geom_ribbon(data = data, aes(x = time,
                                     ymin = between_supp2 - between_se_supp2,
                                     ymax = between_supp2 + between_se_supp2,
                                     fill = "supp2"), alpha = 0.2) +
        geom_ribbon(data = data, aes(x = time,
                                     ymin = between_supp3 - between_se_supp3,
                                     ymax = between_supp3 + between_se_supp3,
                                     fill = "supp3"), alpha = 0.2) +
        geom_ribbon(data = data, aes(x = time,
                                     ymin = between_supp4 - between_se_supp4,
                                     ymax = between_supp4 + between_se_supp4,
                                     fill = "supp4"), alpha = 0.2) +
        geom_ribbon(data = data, aes(x = time,
                                     ymin = between_supp5 - between_se_supp5,
                                     ymax = between_supp5 + between_se_supp5,
                                     fill = "supp5"), alpha = 0.2) +
        
        # Color scales (matched)
        scale_color_manual(values = c(
          supp1 = "#67000d",
          supp2 = "#a50f15",
          supp3 = "#cb181d",
          supp4 = "#ef3b2c",
          supp5 = "#fb6a4a"
        )) +
        
        scale_fill_manual(values = c(
          supp1 = "#67000d",
          supp2 = "#a50f15",
          supp3 = "#cb181d",
          supp4 = "#ef3b2c",
          supp5 = "#fb6a4a"
        )) +   
        
        
        
        ylab(NULL) + xlab(NULL) +
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
      
      ggsave(paste(path_figures, "/connected_scatterplot_", metric_i, "_length", length_seq, ".png", sep = ""), plot, width = 3, height = 3, units = "in")
      
      plot =
        ggplot() +
        geom_vline(xintercept = 0, linetype = "dashed", color = "black", size = 0.3) +
        
        # Lines
        geom_line(data = data, aes(x = time, y = between_supp1, color = "supp1"), size = 1) +
        geom_line(data = data, aes(x = time, y = between_supp2, color = "supp2"), size = 1) +
        geom_line(data = data, aes(x = time, y = between_supp3, color = "supp3"), size = 1) +
        geom_line(data = data, aes(x = time, y = between_supp4, color = "supp4"), size = 1) +
        geom_line(data = data, aes(x = time, y = between_supp5, color = "supp5"), size = 1) +
        
        # Ribbons
        geom_ribbon(data = data, aes(x = time,
                                     ymin = between_supp1 - between_se_supp1,
                                     ymax = between_supp1 + between_se_supp1,
                                     fill = "supp1"), alpha = 0.2) +
        geom_ribbon(data = data, aes(x = time,
                                     ymin = between_supp2 - between_se_supp2,
                                     ymax = between_supp2 + between_se_supp2,
                                     fill = "supp2"), alpha = 0.2) +
        geom_ribbon(data = data, aes(x = time,
                                     ymin = between_supp3 - between_se_supp3,
                                     ymax = between_supp3 + between_se_supp3,
                                     fill = "supp3"), alpha = 0.2) +
        geom_ribbon(data = data, aes(x = time,
                                     ymin = between_supp4 - between_se_supp4,
                                     ymax = between_supp4 + between_se_supp4,
                                     fill = "supp4"), alpha = 0.2) +
        geom_ribbon(data = data, aes(x = time,
                                     ymin = between_supp5 - between_se_supp5,
                                     ymax = between_supp5 + between_se_supp5,
                                     fill = "supp5"), alpha = 0.2) +
        
        # Color scales with labels
        scale_color_manual(
          values = c(
            supp1 = "#67000d",
            supp2 = "#a50f15",
            supp3 = "#cb181d",
            supp4 = "#ef3b2c",
            supp5 = "#fb6a4a"
          ),
          labels = c("20%", "40%", "60%", "80%", "100%"),
          name = "Condition"
        ) +
        
        scale_fill_manual(
          values = c(
            supp1 = "#67000d",
            supp2 = "#a50f15",
            supp3 = "#cb181d",
            supp4 = "#ef3b2c",
            supp5 = "#fb6a4a"
          ),
          labels = c("20%", "40%", "60%", "80%", "100%"),
          name = "Condition"
        ) +
        
        ylab(NULL) + xlab(NULL) +
        
        theme(
          plot.title = element_text(hjust = 0.5),
          legend.position = "right",   # ← show legend
          panel.background = element_blank(),
          plot.background = element_blank(),
          axis.line = element_line(color = "black"),
          axis.title = element_text(color = "black", size = 12),
          axis.text.x = element_text(color = "black", size = 12),
          axis.text.y = element_text(color = "black", size = 12),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank()
        )
      
      ggsave(paste(path_figures, "/connected_scatterplot_", metric_i, "_length", length_seq, "_legend.png", sep = ""), plot, width = 3, height = 3, units = "in")
      
      
      
      
      
      
      
      
      
      
      
      
      ### length 4
      length_seq = 4
      
      between_supp <- list()
      between_se_supp <- list()
      
      for (supp_i in supplementaries) {
        
        between_avg = c()
        between_sd = c()
        
        path = paste(path_results, "/supp", supp_i, "/group_results/", metric_i, sep = "")
        
        # stim
        
        time_window = paste(delay_segments[time_i,1], "to", delay_segments[time_i,2], sep = "")
        filename = paste(path, "/", fun, "/", perf,  "/design_matrix/stim_length", length_seq, "_400to700.txt", sep = "")
        design_matrix = read.table(filename, sep = ",", header = TRUE)
        colnames(design_matrix) = c("response", "Subspaces", "subject", "session")
        design_matrix$Subspaces = factor(ifelse(design_matrix$Subspaces == 1, "Between ranks", "Within rank"))
        data = subset(design_matrix, Subspaces == "Between ranks")
        between_avg = c(between_avg, mean(data$response))
        between_sd = c(between_sd, sd(data$response))
        
        time_window = paste(delay_segments[time_i,1], "to", delay_segments[time_i,2], sep = "")
        filename = paste(path, "/", fun, "/", perf,  "/design_matrix/stim_length", length_seq, "_800to1100.txt", sep = "")
        design_matrix = read.table(filename, sep = ",", header = TRUE)
        colnames(design_matrix) = c("response", "Subspaces", "subject", "session")
        design_matrix$Subspaces = factor(ifelse(design_matrix$Subspaces == 1, "Between ranks", "Within rank"))
        data = subset(design_matrix, Subspaces == "Between ranks")
        between_avg = c(between_avg, mean(data$response))
        between_sd = c(between_sd, sd(data$response))
        
        time_window = paste(delay_segments[time_i,1], "to", delay_segments[time_i,2], sep = "")
        filename = paste(path, "/", fun, "/", perf,  "/design_matrix/stim_length", length_seq, "_1200to1500.txt", sep = "")
        design_matrix = read.table(filename, sep = ",", header = TRUE)
        colnames(design_matrix) = c("response", "Subspaces", "subject", "session")
        design_matrix$Subspaces = factor(ifelse(design_matrix$Subspaces == 1, "Between ranks", "Within rank"))
        data = subset(design_matrix, Subspaces == "Between ranks")
        between_avg = c(between_avg, mean(data$response))
        between_sd = c(between_sd, sd(data$response))
        
        # delay
        
        for (time_i in 1:nrow(delay_segments)) {
          time_window = paste(delay_segments[time_i,1], "to", delay_segments[time_i,2], sep = "")
          filename = paste(path, "/", fun, "/", perf,  "/design_matrix/delay_length", length_seq, "_", time_window, ".txt", sep = "")
          design_matrix = read.table(filename, sep = ",", header = TRUE)
          colnames(design_matrix) = c("response", "Subspaces", "subject", "session")
          design_matrix$Subspaces = factor(ifelse(design_matrix$Subspaces == 1, "Between ranks", "Within rank"))
          data = subset(design_matrix, Subspaces == "Between ranks")
          between_avg = c(between_avg, mean(data$response))
          between_sd = c(between_sd, sd(data$response))
        }
        
        
        between_supp[[supp_i]] <- between_avg
        between_se_supp[[supp_i]] <- between_sd / sqrt(15)     
        
      }
      
      
      
      # set in data frame
      time = c(-750, -450, -150, delay_segments[,2] - 150) / 1000
      data = as.data.frame(cbind(between_supp[[1]], between_supp[[2]], between_supp[[3]], between_supp[[4]], between_supp[[5]],
                                 between_se_supp[[1]], between_se_supp[[2]], between_se_supp[[3]], between_se_supp[[4]], between_se_supp[[5]],
                                 time))
      colnames(data) = c("between_supp1", "between_supp2", "between_supp3", "between_supp4", "between_supp5", 
                         "between_se_supp1", "between_se_supp2", "between_se_supp3", "between_se_supp4", "between_se_supp5", 
                         "time")
      
      # Plot between and within with standard errors
      
      plot =
        ggplot() +
        geom_vline(xintercept = 0, linetype = "dashed", color = "black", size = 0.3) +
        
        # Lines
        geom_line(data = data, aes(x = time, y = between_supp1, color = "supp1"), size = 1) +
        geom_line(data = data, aes(x = time, y = between_supp2, color = "supp2"), size = 1) +
        geom_line(data = data, aes(x = time, y = between_supp3, color = "supp3"), size = 1) +
        geom_line(data = data, aes(x = time, y = between_supp4, color = "supp4"), size = 1) +
        geom_line(data = data, aes(x = time, y = between_supp5, color = "supp5"), size = 1) +
        
        # Ribbons
        geom_ribbon(data = data, aes(x = time,
                                     ymin = between_supp1 - between_se_supp1,
                                     ymax = between_supp1 + between_se_supp1,
                                     fill = "supp1"), alpha = 0.2) +
        geom_ribbon(data = data, aes(x = time,
                                     ymin = between_supp2 - between_se_supp2,
                                     ymax = between_supp2 + between_se_supp2,
                                     fill = "supp2"), alpha = 0.2) +
        geom_ribbon(data = data, aes(x = time,
                                     ymin = between_supp3 - between_se_supp3,
                                     ymax = between_supp3 + between_se_supp3,
                                     fill = "supp3"), alpha = 0.2) +
        geom_ribbon(data = data, aes(x = time,
                                     ymin = between_supp4 - between_se_supp4,
                                     ymax = between_supp4 + between_se_supp4,
                                     fill = "supp4"), alpha = 0.2) +
        geom_ribbon(data = data, aes(x = time,
                                     ymin = between_supp5 - between_se_supp5,
                                     ymax = between_supp5 + between_se_supp5,
                                     fill = "supp5"), alpha = 0.2) +
        
        # Color scales (matched)
        scale_color_manual(values = c(
          supp1 = "#67000d",
          supp2 = "#a50f15",
          supp3 = "#cb181d",
          supp4 = "#ef3b2c",
          supp5 = "#fb6a4a"
        )) +
        
        scale_fill_manual(values = c(
          supp1 = "#67000d",
          supp2 = "#a50f15",
          supp3 = "#cb181d",
          supp4 = "#ef3b2c",
          supp5 = "#fb6a4a"
        )) +   
        
        
        
        ylab(NULL) + xlab(NULL) +
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
      
      ggsave(paste(path_figures, "/connected_scatterplot_", metric_i, "_length", length_seq, ".png", sep = ""), plot, width = 3, height = 3, units = "in")
      
      plot =
        ggplot() +
        geom_vline(xintercept = 0, linetype = "dashed", color = "black", size = 0.3) +
        
        # Lines
        geom_line(data = data, aes(x = time, y = between_supp1, color = "supp1"), size = 1) +
        geom_line(data = data, aes(x = time, y = between_supp2, color = "supp2"), size = 1) +
        geom_line(data = data, aes(x = time, y = between_supp3, color = "supp3"), size = 1) +
        geom_line(data = data, aes(x = time, y = between_supp4, color = "supp4"), size = 1) +
        geom_line(data = data, aes(x = time, y = between_supp5, color = "supp5"), size = 1) +
        
        # Ribbons
        geom_ribbon(data = data, aes(x = time,
                                     ymin = between_supp1 - between_se_supp1,
                                     ymax = between_supp1 + between_se_supp1,
                                     fill = "supp1"), alpha = 0.2) +
        geom_ribbon(data = data, aes(x = time,
                                     ymin = between_supp2 - between_se_supp2,
                                     ymax = between_supp2 + between_se_supp2,
                                     fill = "supp2"), alpha = 0.2) +
        geom_ribbon(data = data, aes(x = time,
                                     ymin = between_supp3 - between_se_supp3,
                                     ymax = between_supp3 + between_se_supp3,
                                     fill = "supp3"), alpha = 0.2) +
        geom_ribbon(data = data, aes(x = time,
                                     ymin = between_supp4 - between_se_supp4,
                                     ymax = between_supp4 + between_se_supp4,
                                     fill = "supp4"), alpha = 0.2) +
        geom_ribbon(data = data, aes(x = time,
                                     ymin = between_supp5 - between_se_supp5,
                                     ymax = between_supp5 + between_se_supp5,
                                     fill = "supp5"), alpha = 0.2) +
        
        # Color scales with labels
        scale_color_manual(
          values = c(
            supp1 = "#67000d",
            supp2 = "#a50f15",
            supp3 = "#cb181d",
            supp4 = "#ef3b2c",
            supp5 = "#fb6a4a"
          ),
          labels = c("20%", "40%", "60%", "80%", "100%"),
          name = "Condition"
        ) +
        
        scale_fill_manual(
          values = c(
            supp1 = "#67000d",
            supp2 = "#a50f15",
            supp3 = "#cb181d",
            supp4 = "#ef3b2c",
            supp5 = "#fb6a4a"
          ),
          labels = c("20%", "40%", "60%", "80%", "100%"),
          name = "Condition"
        ) +
        
        ylab(NULL) + xlab(NULL) +
        
        theme(
          plot.title = element_text(hjust = 0.5),
          legend.position = "right",   # ← show legend
          panel.background = element_blank(),
          plot.background = element_blank(),
          axis.line = element_line(color = "black"),
          axis.title = element_text(color = "black", size = 12),
          axis.text.x = element_text(color = "black", size = 12),
          axis.text.y = element_text(color = "black", size = 12),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank()
        )
      
      ggsave(paste(path_figures, "/connected_scatterplot_", metric_i, "_length", length_seq, "_legend.png", sep = ""), plot, width = 3, height = 3, units = "in")
      
      
      
      
      
      
      
    }
    
    
  }
}

