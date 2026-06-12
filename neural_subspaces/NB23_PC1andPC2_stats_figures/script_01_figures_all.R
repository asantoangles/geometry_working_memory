# libraries
library(ggplot2)
library(lmerTest)
library(gridExtra)
library(cowplot)
library(GGally)
library(zoo)

####################
#####  stats  ######
####################

main_folder = "NB23_PC1andPC2_stats_figures"
setwd(paste("/path_to_local/github/scripts/source_geometry_lm/", main_folder, sep = ""))
source('script_02_stats_all.R')

####################
##### figures ######
####################

path_figures = paste("/path_to_local/github/results/source_geometry_lm/", main_folder, sep = "")

if (!dir.exists(path_figures)) {
  dir.create(path_figures, recursive = TRUE)
}

last_idx = 0

### principal angle / vaf (between vs within) in stim and delay for correct trials time-resolved

folder = "NB23_PC1andPC2"
path_results = paste("/path_to_local/github/results/source_geometry_lm/", folder, sep = "")

folder_random = "NB23_PC1andPC2_random"
path_results_random = paste("/path_to_local/github/results/source_geometry_lm/", folder_random, sep = "")

performance = c("correct_trials")
functions_data = c("mean")
metrics = c("principal_angle", "vaf")

first_idx = last_idx + 1
last_idx = length(p_values_003_pa)
p_values_fdr_003_pa = p_values_fdr_plots[first_idx:last_idx]

first_idx = last_idx + 1
last_idx = last_idx + length(p_values_003_vaf) 
p_values_fdr_003_vaf = p_values_fdr_plots[first_idx:last_idx]

p_values_all = c(p_values_fdr_003_pa, p_values_fdr_003_vaf)

first_idx = last_idx + 1
last_idx = last_idx + length(p_values_003_pa_random) 
p_values_fdr_003_pa_random = p_values_fdr_plots[first_idx:last_idx]

first_idx = last_idx + 1
last_idx = last_idx + length(p_values_003_vaf_random) 
p_values_fdr_003_vaf_random = p_values_fdr_plots[first_idx:last_idx]

p_values_all_random = c(p_values_fdr_003_pa_random, p_values_fdr_003_vaf_random)

source('script_08_figures_alignment.R')









### principal angle / vaf (correct vs incorrect trials) in stim and delay - time-resolved

folder = "NB23_PC1andPC2_controlled_resampling"
path_results = paste("/path_to_local/github/results/source_geometry_lm/", folder, sep = "")

functions_data = c("mean")
metrics = c("principal_angle", "vaf")

first_idx = last_idx + 1
last_idx = last_idx + length(p_values_004_pa) 
p_values_fdr_004_pa = p_values_fdr_plots[first_idx:last_idx]

first_idx = last_idx + 1
last_idx = last_idx + length(p_values_004_vaf) 
p_values_fdr_004_vaf = p_values_fdr_plots[first_idx:last_idx]

p_values_all = c(p_values_fdr_004_pa, p_values_fdr_004_vaf)

source('script_09_figures_alignment_performance.R')







### principal angle / vaf (length 3 vs 4) in stim and delay - time-resolved

folder = "NB23_PC1andPC2"
path_results = paste("/path_to_local/github/results/source_geometry_lm/", folder, sep = "")

functions_data = c("mean")
metrics = c("principal_angle", "vaf")

first_idx = last_idx + 1
last_idx = last_idx + length(p_values_005_pa) 
p_values_fdr_005_pa = p_values_fdr_plots[first_idx:last_idx]

first_idx = last_idx + 1
last_idx = last_idx + length(p_values_005_vaf) 
p_values_fdr_005_vaf = p_values_fdr_plots[first_idx:last_idx]

p_values_all = c(p_values_fdr_005_pa, p_values_fdr_005_vaf)

source('script_10_figures_alignment_length.R')










### separability (distance, distance_noclose_pairs, volume) - correct vs incorrect trials - time-resolved (bootstrap)

folder = "NB23_PC1andPC2_controlled_resampling"
path_results = paste("/path_to_local/github/results/source_geometry_lm/", folder, sep = "")

functions_data = c("mean")
separability_metrics = c("distance", "volume")

first_idx = last_idx + 1
last_idx = last_idx + length(p_values_006_distance) 
p_values_fdr_006_distance = p_values_fdr_plots[first_idx:last_idx]

first_idx = last_idx + 1
last_idx = last_idx + length(p_values_006_volume) 
p_values_fdr_006_volume = p_values_fdr_plots[first_idx:last_idx]

p_values_all = c(p_values_fdr_006_distance, p_values_fdr_006_volume)

source('script_11_figures_separability_performance.R')







### separability (distance) - length 3 vs length 4 in correct trials - time-resolved (bootstrap)


folder = "NB23_PC1andPC2"
path_results = paste("/path_to_local/github/results/source_geometry_lm/", folder, sep = "")

functions_data = c("mean")
separability_metrics = c("distance", "volume")

first_idx = last_idx + 1
last_idx = last_idx + length(p_values_007_distance) 
p_values_fdr_007_distance = p_values_fdr_plots[first_idx:last_idx]

first_idx = last_idx + 1
last_idx = last_idx + length(p_values_007_volume) 
p_values_fdr_007_volume = p_values_fdr_plots[first_idx:last_idx]

p_values_all = c(p_values_fdr_007_distance, p_values_fdr_007_volume)

source('script_12_figures_separability_length.R')





