# libraries
library(ggplot2)
library(lmerTest)
library(gridExtra)
library(cowplot)
library(GGally)
library(zoo)

setwd("/path_to_local/github/scripts/source_geometry_lm/RB23_stats")

####################
##### figures ######
####################

folder_figures = "RB23_stats"
path_figures = paste("/path_to_local/github/results/source_geometry_lm/", folder_figures, sep = "")

if (!dir.exists(path_figures)) {
  dir.create(path_figures, recursive = TRUE)
}

### alignment

folder = "RB23"
path_results = paste("/path_to_local/github/results/source_geometry_lm/", folder, sep = "")

performance = c("correct_trials")
functions_data = c("mean")
metrics = c("principal_angle", "vaf")
supplementaries = seq(1,5)

source('script_02_supp_figures_alignment.R')







### separability

folder = "RB23"
path_results = paste("/path_to_local/github/results/source_geometry_lm/", folder, sep = "")

functions_data = c("mean")
separability_metrics = c("distance", "volume")

source('script_03_supp_figures_separability.R')




