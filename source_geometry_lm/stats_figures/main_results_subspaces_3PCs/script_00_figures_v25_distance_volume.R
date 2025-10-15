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

setwd("/path_to_local/scripts/source_geometry_lm/LB20/main_results_subspaces_3PCs")
source('script_00_stats_v25_distance_volume.R')

####################
##### figures ######
####################

folder_figures = "LB20"
path_figures = paste("/path_to_local/results/source_geometry_lm/", folder_figures, "/main_results_subspaces_3PCs", sep = "")






### principal angle / vaf (between vs within) in stim and delay for correct trials time-resolved

folder = "LB20"
path_results = paste("/path_to_local/results/source_geometry_lm/", folder, sep = "")

folder_random = "LB20_random"
path_results_random = paste("/path_to_local/results/source_geometry_lm/", folder_random, sep = "")

performance = c("correct_trials")
functions_data = c("mean")
metrics = c("principal_angle", "vaf")

p_values_fdr_02 = p_values_fdr_plots[1:31]
p_values_fdr_04 = p_values_fdr_plots[32:62]
p_values_all = c(p_values_fdr_02, p_values_fdr_04)

p_values_fdr_02_random = p_values_fdr_plots[63:93]
p_values_fdr_04_random = p_values_fdr_plots[94:124]
p_values_all_random = c(p_values_fdr_02_random, p_values_fdr_04_random)

source('script_202_figures_pa_vaf_time_resolved_refined_random.R')









### principal angle / vaf (correct vs incorrect trials) in stim and delay - time-resolved

folder = "LB20_controlled_resampling"
path_results = paste("/path_to_local/results/source_geometry_lm/", folder, sep = "")

functions_data = c("mean")
metrics = c("principal_angle", "vaf")

p_values_fdr_06 = p_values_fdr_plots[125:155]
p_values_fdr_08 = p_values_fdr_plots[156:186]
p_values_all = c(p_values_fdr_06, p_values_fdr_08)

source('script_204_figures_pa_vaf_performance_time_resolved_refined.R')







### principal angle / vaf (length 3 vs 4) in stim and delay - time-resolved

folder = "LB20"
path_results = paste("/path_to_local/results/source_geometry_lm/", folder, sep = "")

functions_data = c("mean")
metrics = c("principal_angle", "vaf")

p_values_fdr_12 = p_values_fdr_plots[187:201]
p_values_fdr_16 = p_values_fdr_plots[202:216]
p_values_all = c(p_values_fdr_12, p_values_fdr_16)

source('script_206_figures_pa_vaf_length3vs4_time_resolved_refined.R')










### separability (distance, distance_noclose_pairs, volume) - correct vs incorrect trials - time-resolved (bootstrap)

folder = "LB20_controlled_resampling"
path_results = paste("/path_to_local/results/source_geometry_lm/", folder, sep = "")

functions_data = c("mean")
separability_metrics = c("distance", "volume")

p_values_fdr_24a = p_values_fdr_plots[217:249]
p_values_fdr_24c = p_values_fdr_plots[250:282]
p_values_all = c(p_values_fdr_24a, p_values_fdr_24c)

source('script_208_figures_separability_correct_vs_incorrect_time_resolved_refined_random.R')







### separability (distance) - length 3 vs length 4 in correct trials - time-resolved (bootstrap)


folder = "LB20"
path_results = paste("/path_to_local/results/source_geometry_lm/", folder, sep = "")

functions_data = c("mean")
separability_metrics = c("distance", "volume")

p_values_fdr_34a = p_values_fdr_plots[283:297]
p_values_fdr_34c = p_values_fdr_plots[298:312]
p_values_all = c(p_values_fdr_34a, p_values_fdr_34c)

source('script_210_figures_separability_length3vs4_time_resolved_refined.R')





