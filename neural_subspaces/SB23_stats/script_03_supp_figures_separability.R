##############################
## connected scatter plots  ##
##############################

segments_start <- seq(from = 1, to = 3601, by = 300)
segments_end <- seq(from = 300, to = 4000, by = 300)
delay_segments = cbind(segments_start, segments_end)

for (separability in separability_metrics) {
  
  cat(separability); cat(' ')
  
  if (separability == "distance") {
    y_axis_step = 0.05
  }
  if (separability == "volume") {
    y_axis_step = 0.5
  }

  for (fun in functions_data) {
    
    ## length 3
    
    seq_length = 3
    
    separability_supp <- list()
    separability_se_supp <- list()
    
    for (supp_i in supplementaries) {
      
      separability_avg = c()
      separability_sd = c()
      
      path = paste(path_results, "/supp", supp_i, "/group_results/", separability, sep = "")
      
      # stim
      
      filename = paste(path, "/", fun, "/correct_trials/design_matrix/stim_length", seq_length, "_1to300.txt", sep = "")
      matrix_design = read.table(filename, sep = ",", header = TRUE)
      data = matrix_design
      colnames(data) = c("response", "condition", "subject", "session")
      separability_avg = c(separability_avg, mean(data$response))
      separability_sd = c(separability_sd, sd(data$response))

      filename = paste(path, "/", fun, "/correct_trials/design_matrix/stim_length", seq_length, "_400to700.txt", sep = "")
      matrix_design = read.table(filename, sep = ",", header = TRUE)
      data = matrix_design
      colnames(data) = c("response", "condition", "subject", "session")
      separability_avg = c(separability_avg, mean(data$response))
      separability_sd = c(separability_sd, sd(data$response))
      
      filename = paste(path, "/", fun, "/correct_trials/design_matrix/stim_length", seq_length, "_800to1100.txt", sep = "")
      matrix_design = read.table(filename, sep = ",", header = TRUE)
      data = matrix_design
      colnames(data) = c("response", "condition", "subject", "session")
      separability_avg = c(separability_avg, mean(data$response))
      separability_sd = c(separability_sd, sd(data$response))
      
      # delay
      
      for (time_i in 1:nrow(delay_segments)) {
        time_window = paste(delay_segments[time_i,1], "to", delay_segments[time_i,2], sep = "")
        
        filename = paste(path, "/", fun, "/correct_trials/design_matrix/delay_length", seq_length, "_", time_window, ".txt", sep = "")
        matrix_design = read.table(filename, sep = ",", header = TRUE)
        data = matrix_design
        colnames(data) = c("response", "condition", "subject", "session")
        separability_avg = c(separability_avg, mean(data$response))
        separability_sd = c(separability_sd, sd(data$response))
        

      }
      
      separability_supp[[supp_i]] <- separability_avg
      separability_se_supp[[supp_i]] <- separability_sd / sqrt(15)     
      
      
    }
    
    # set in data frame
    time = c(-750, -450, -150, delay_segments[,2] - 150) / 1000
    data = as.data.frame(cbind(separability_supp[[1]], separability_supp[[2]], separability_supp[[3]], separability_supp[[4]], separability_supp[[5]],
                               separability_se_supp[[1]], separability_se_supp[[2]], separability_se_supp[[3]], separability_se_supp[[4]], separability_se_supp[[5]],
                               time))
    colnames(data) = c("separability_supp1", "separability_supp2", "separability_supp3", "separability_supp4", "separability_supp5", 
                       "separability_se_supp1", "separability_se_supp2", "separability_se_supp3", "separability_se_supp4", "separability_se_supp5", 
                       "time")
    
    # Plot between and within with standard errors
    
    plot =
      ggplot() +
      geom_vline(xintercept = 0, linetype = "dashed", color = "black", size = 0.3) +
      
      # Lines
      geom_line(data = data, aes(x = time, y = separability_supp1, color = "supp1"), size = 1) +
      geom_line(data = data, aes(x = time, y = separability_supp2, color = "supp2"), size = 1) +
      geom_line(data = data, aes(x = time, y = separability_supp3, color = "supp3"), size = 1) +
      geom_line(data = data, aes(x = time, y = separability_supp4, color = "supp4"), size = 1) +
      geom_line(data = data, aes(x = time, y = separability_supp5, color = "supp5"), size = 1) +
      
      # Ribbons
      geom_ribbon(data = data, aes(x = time,
                                   ymin = separability_supp1 - separability_se_supp1,
                                   ymax = separability_supp1 + separability_se_supp1,
                                   fill = "supp1"), alpha = 0.2) +
      geom_ribbon(data = data, aes(x = time,
                                   ymin = separability_supp2 - separability_se_supp2,
                                   ymax = separability_supp2 + separability_se_supp2,
                                   fill = "supp2"), alpha = 0.2) +
      geom_ribbon(data = data, aes(x = time,
                                   ymin = separability_supp3 - separability_se_supp3,
                                   ymax = separability_supp3 + separability_se_supp3,
                                   fill = "supp3"), alpha = 0.2) +
      geom_ribbon(data = data, aes(x = time,
                                   ymin = separability_supp4 - separability_se_supp4,
                                   ymax = separability_supp4 + separability_se_supp4,
                                   fill = "supp4"), alpha = 0.2) +
      geom_ribbon(data = data, aes(x = time,
                                   ymin = separability_supp5 - separability_se_supp5,
                                   ymax = separability_supp5 + separability_se_supp5,
                                   fill = "supp5"), alpha = 0.2) +
      
      # Color scales (matched)
      scale_color_manual(values = c(
        supp1 = "#00441b",
        supp2 = "#1b7837",
        supp3 = "#5aae61",
        supp4 = "#a6dba0",
        supp5 = "#d9f0d3"
      )) +
      
      scale_fill_manual(values = c(
        supp1 = "#00441b",
        supp2 = "#1b7837",
        supp3 = "#5aae61",
        supp4 = "#a6dba0",
        supp5 = "#d9f0d3"
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
    
    ggsave(paste(path_figures, "/connected_scatterplot_", separability, "_length", seq_length, ".png", sep = ""), plot, width = 3, height = 3, units = "in")
    
    plot =
      ggplot() +
      geom_vline(xintercept = 0, linetype = "dashed", color = "black", size = 0.3) +
      
      # Lines
      geom_line(data = data, aes(x = time, y = separability_supp1, color = "supp1"), size = 1) +
      geom_line(data = data, aes(x = time, y = separability_supp2, color = "supp2"), size = 1) +
      geom_line(data = data, aes(x = time, y = separability_supp3, color = "supp3"), size = 1) +
      geom_line(data = data, aes(x = time, y = separability_supp4, color = "supp4"), size = 1) +
      geom_line(data = data, aes(x = time, y = separability_supp5, color = "supp5"), size = 1) +
      
      # Ribbons
      geom_ribbon(data = data, aes(x = time,
                                   ymin = separability_supp1 - separability_se_supp1,
                                   ymax = separability_supp1 + separability_se_supp1,
                                   fill = "supp1"), alpha = 0.2) +
      geom_ribbon(data = data, aes(x = time,
                                   ymin = separability_supp2 - separability_se_supp2,
                                   ymax = separability_supp2 + separability_se_supp2,
                                   fill = "supp2"), alpha = 0.2) +
      geom_ribbon(data = data, aes(x = time,
                                   ymin = separability_supp3 - separability_se_supp3,
                                   ymax = separability_supp3 + separability_se_supp3,
                                   fill = "supp3"), alpha = 0.2) +
      geom_ribbon(data = data, aes(x = time,
                                   ymin = separability_supp4 - separability_se_supp4,
                                   ymax = separability_supp4 + separability_se_supp4,
                                   fill = "supp4"), alpha = 0.2) +
      geom_ribbon(data = data, aes(x = time,
                                   ymin = separability_supp5 - separability_se_supp5,
                                   ymax = separability_supp5 + separability_se_supp5,
                                   fill = "supp5"), alpha = 0.2) +
      
      # Color scales with labels
      scale_color_manual(
        values = c(
          supp1 = "#00441b",
          supp2 = "#1b7837",
          supp3 = "#5aae61",
          supp4 = "#a6dba0",
          supp5 = "#d9f0d3"
        ),
        labels = c("20%", "40%", "60%", "80%", "100%"),
        name = "Condition"
      ) +
      
      scale_fill_manual(
        values = c(
          supp1 = "#00441b",
          supp2 = "#1b7837",
          supp3 = "#5aae61",
          supp4 = "#a6dba0",
          supp5 = "#d9f0d3"
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
    
    ggsave(paste(path_figures, "/connected_scatterplot_", separability, "_length", seq_length, "_legend.png", sep = ""), plot, width = 3, height = 3, units = "in")
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    ## length 4
    
    seq_length = 4
    
    separability_supp <- list()
    separability_se_supp <- list()
    
    for (supp_i in supplementaries) {
      
      separability_avg = c()
      separability_sd = c()
      
      path = paste(path_results, "/supp", supp_i, "/group_results/", separability, sep = "")
      
      # stim
      
      filename = paste(path, "/", fun, "/correct_trials/design_matrix/stim_length", seq_length, "_1to300.txt", sep = "")
      matrix_design = read.table(filename, sep = ",", header = TRUE)
      data = matrix_design
      colnames(data) = c("response", "condition", "subject", "session")
      separability_avg = c(separability_avg, mean(data$response))
      separability_sd = c(separability_sd, sd(data$response))
      
      filename = paste(path, "/", fun, "/correct_trials/design_matrix/stim_length", seq_length, "_400to700.txt", sep = "")
      matrix_design = read.table(filename, sep = ",", header = TRUE)
      data = matrix_design
      colnames(data) = c("response", "condition", "subject", "session")
      separability_avg = c(separability_avg, mean(data$response))
      separability_sd = c(separability_sd, sd(data$response))
      
      filename = paste(path, "/", fun, "/correct_trials/design_matrix/stim_length", seq_length, "_800to1100.txt", sep = "")
      matrix_design = read.table(filename, sep = ",", header = TRUE)
      data = matrix_design
      colnames(data) = c("response", "condition", "subject", "session")
      separability_avg = c(separability_avg, mean(data$response))
      separability_sd = c(separability_sd, sd(data$response))
      
      filename = paste(path, "/", fun, "/correct_trials/design_matrix/stim_length", seq_length, "_1200to1500.txt", sep = "")
      matrix_design = read.table(filename, sep = ",", header = TRUE)
      data = matrix_design
      colnames(data) = c("response", "condition", "subject", "session")
      separability_avg = c(separability_avg, mean(data$response))
      separability_sd = c(separability_sd, sd(data$response))
      
      # delay
      
      for (time_i in 1:nrow(delay_segments)) {
        time_window = paste(delay_segments[time_i,1], "to", delay_segments[time_i,2], sep = "")
        
        filename = paste(path, "/", fun, "/correct_trials/design_matrix/delay_length", seq_length, "_", time_window, ".txt", sep = "")
        matrix_design = read.table(filename, sep = ",", header = TRUE)
        data = matrix_design
        colnames(data) = c("response", "condition", "subject", "session")
        separability_avg = c(separability_avg, mean(data$response))
        separability_sd = c(separability_sd, sd(data$response))
        
        
      }
      
      separability_supp[[supp_i]] <- separability_avg
      separability_se_supp[[supp_i]] <- separability_sd / sqrt(15)     
      
      
    }
    
    # set in data frame
    time = c(-1150, -750, -450, -150, delay_segments[,2] - 150) / 1000
    data = as.data.frame(cbind(separability_supp[[1]], separability_supp[[2]], separability_supp[[3]], separability_supp[[4]], separability_supp[[5]],
                               separability_se_supp[[1]], separability_se_supp[[2]], separability_se_supp[[3]], separability_se_supp[[4]], separability_se_supp[[5]],
                               time))
    colnames(data) = c("separability_supp1", "separability_supp2", "separability_supp3", "separability_supp4", "separability_supp5", 
                       "separability_se_supp1", "separability_se_supp2", "separability_se_supp3", "separability_se_supp4", "separability_se_supp5", 
                       "time")
    
    # Plot between and within with standard errors
    
    plot =
      ggplot() +
      geom_vline(xintercept = 0, linetype = "dashed", color = "black", size = 0.3) +
      
      # Lines
      geom_line(data = data, aes(x = time, y = separability_supp1, color = "supp1"), size = 1) +
      geom_line(data = data, aes(x = time, y = separability_supp2, color = "supp2"), size = 1) +
      geom_line(data = data, aes(x = time, y = separability_supp3, color = "supp3"), size = 1) +
      geom_line(data = data, aes(x = time, y = separability_supp4, color = "supp4"), size = 1) +
      geom_line(data = data, aes(x = time, y = separability_supp5, color = "supp5"), size = 1) +
      
      # Ribbons
      geom_ribbon(data = data, aes(x = time,
                                   ymin = separability_supp1 - separability_se_supp1,
                                   ymax = separability_supp1 + separability_se_supp1,
                                   fill = "supp1"), alpha = 0.2) +
      geom_ribbon(data = data, aes(x = time,
                                   ymin = separability_supp2 - separability_se_supp2,
                                   ymax = separability_supp2 + separability_se_supp2,
                                   fill = "supp2"), alpha = 0.2) +
      geom_ribbon(data = data, aes(x = time,
                                   ymin = separability_supp3 - separability_se_supp3,
                                   ymax = separability_supp3 + separability_se_supp3,
                                   fill = "supp3"), alpha = 0.2) +
      geom_ribbon(data = data, aes(x = time,
                                   ymin = separability_supp4 - separability_se_supp4,
                                   ymax = separability_supp4 + separability_se_supp4,
                                   fill = "supp4"), alpha = 0.2) +
      geom_ribbon(data = data, aes(x = time,
                                   ymin = separability_supp5 - separability_se_supp5,
                                   ymax = separability_supp5 + separability_se_supp5,
                                   fill = "supp5"), alpha = 0.2) +
      
      # Color scales (matched)
      scale_color_manual(values = c(
        supp1 = "#00441b",
        supp2 = "#1b7837",
        supp3 = "#5aae61",
        supp4 = "#a6dba0",
        supp5 = "#d9f0d3"
        )) +
      
      scale_fill_manual(values = c(
        supp1 = "#00441b",
        supp2 = "#1b7837",
        supp3 = "#5aae61",
        supp4 = "#a6dba0",
        supp5 = "#d9f0d3"
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
    
    ggsave(paste(path_figures, "/connected_scatterplot_", separability, "_length", seq_length, ".png", sep = ""), plot, width = 3, height = 3, units = "in")
    
    plot =
      ggplot() +
      geom_vline(xintercept = 0, linetype = "dashed", color = "black", size = 0.3) +
      
      # Lines
      geom_line(data = data, aes(x = time, y = separability_supp1, color = "supp1"), size = 1) +
      geom_line(data = data, aes(x = time, y = separability_supp2, color = "supp2"), size = 1) +
      geom_line(data = data, aes(x = time, y = separability_supp3, color = "supp3"), size = 1) +
      geom_line(data = data, aes(x = time, y = separability_supp4, color = "supp4"), size = 1) +
      geom_line(data = data, aes(x = time, y = separability_supp5, color = "supp5"), size = 1) +
      
      # Ribbons
      geom_ribbon(data = data, aes(x = time,
                                   ymin = separability_supp1 - separability_se_supp1,
                                   ymax = separability_supp1 + separability_se_supp1,
                                   fill = "supp1"), alpha = 0.2) +
      geom_ribbon(data = data, aes(x = time,
                                   ymin = separability_supp2 - separability_se_supp2,
                                   ymax = separability_supp2 + separability_se_supp2,
                                   fill = "supp2"), alpha = 0.2) +
      geom_ribbon(data = data, aes(x = time,
                                   ymin = separability_supp3 - separability_se_supp3,
                                   ymax = separability_supp3 + separability_se_supp3,
                                   fill = "supp3"), alpha = 0.2) +
      geom_ribbon(data = data, aes(x = time,
                                   ymin = separability_supp4 - separability_se_supp4,
                                   ymax = separability_supp4 + separability_se_supp4,
                                   fill = "supp4"), alpha = 0.2) +
      geom_ribbon(data = data, aes(x = time,
                                   ymin = separability_supp5 - separability_se_supp5,
                                   ymax = separability_supp5 + separability_se_supp5,
                                   fill = "supp5"), alpha = 0.2) +
      
      # Color scales with labels
      scale_color_manual(
        values = c(
          supp1 = "#00441b",
          supp2 = "#1b7837",
          supp3 = "#5aae61",
          supp4 = "#a6dba0",
          supp5 = "#d9f0d3"
        ),
        labels = c("20%", "40%", "60%", "80%", "100%"),
        name = "Condition"
      ) +
      
      scale_fill_manual(
        values = c(
          supp1 = "#00441b",
          supp2 = "#1b7837",
          supp3 = "#5aae61",
          supp4 = "#a6dba0",
          supp5 = "#d9f0d3"
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
    
    ggsave(paste(path_figures, "/connected_scatterplot_", separability, "_length", seq_length, "_legend.png", sep = ""), plot, width = 3, height = 3, units = "in")
    
    
    
    
    
    
  }
}
