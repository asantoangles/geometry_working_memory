
# libraries
library(ggplot2)
library(lmerTest)
library(gridExtra)
library(cowplot)
library(GGally)
library(zoo)


setwd("/path_to_local/scripts/source_geometry_lm/LB20/stats_figures_subspaces_3PCs")

##################
##### stats ######
##################



### principal angle (between vs within) in stim and delay for correct trials time-resolved

folder_source = "LB20"
path_results = paste("/path_to_local/results/source_geometry_lm/", folder_source, sep = "")

folder_source_random = "LB20_random"
path_results_random = paste("/path_to_local/results/source_geometry_lm/", folder_source_random, sep = "")

performance = c("correct_trials")
functions_data = c("mean")
metrics = c("principal_angle")

source('script_02_stats_pa_vaf_time_resolved_refined_random.R')

p_values_02 = p_values
t_values_02 = t_values

p_values_02_random = p_values_random
t_values_02_random = t_values_random






### vaf (between vs within) in stim and delay for correct trials time-resolved

folder_source = "LB20"
path_results = paste("/path_to_local/results/source_geometry_lm/", folder_source, sep = "")

folder_source_random = "LB20_random"
path_results_random = paste("/path_to_local/results/source_geometry_lm/", folder_source_random, sep = "")

performance = c("correct_trials")
functions_data = c("mean")
metrics = c("vaf")

source('script_02_stats_pa_vaf_time_resolved_refined_random.R')

p_values_04 = p_values
t_values_04 = t_values

p_values_04_random = p_values_random
t_values_04_random = t_values_random














### principal angle (correct vs incorrect) in stim and delay time-resolved

folder = "LB20_controlled_resampling"
path_results = paste("/path_to_local/results/source_geometry_lm/", folder, sep = "")

functions_data = c("mean")
metrics = c("principal_angle")

source('script_06_stats_pa_vaf_performance_time_resolved_refined.R')

p_values_06 = p_values
t_values_06 = t_values







### vaf (correct vs incorrect) in stim and delay time-resolved

folder = "LB20_controlled_resampling"
path_results = paste("/path_to_local/results/source_geometry_lm/", folder, sep = "")

functions_data = c("mean")
metrics = c("vaf")

source('script_06_stats_pa_vaf_performance_time_resolved_refined.R')

p_values_08 = p_values
t_values_08 = t_values






















### principal angle (between) (length 3 vs length 4) in correct trials in stim and delay time-resolved (bootstrap)

folder_source = "LB20"
path_results = paste("/path_to_local/results/source_geometry_lm/", folder_source, sep = "")

performance = c("correct_trials")
functions_data = c("mean")
metrics = c("principal_angle")

source('script_12_stats_pa_vaf_length_time_resolved_refined.R')

p_values_12 = p_values
t_values_12 = t_values













### vaf (length 3 vs length 4) in correct trials in stim and delay time-resolved (bootstrap)

folder_source = "LB20"
path_results = paste("/path_to_local/results/source_geometry_lm/", folder_source, sep = "")

performance = c("correct_trials")
functions_data = c("mean")
metrics = c("vaf")

source('script_12_stats_pa_vaf_length_time_resolved_refined.R')

p_values_16 = p_values
t_values_16 = t_values

















### separability (distance) - correct vs incorrect trials - time-resolved (bootstrap)

folder = "LB20_controlled_resampling"
path_results = paste("/path_to_local/results/source_geometry_lm/", folder, sep = "")

folder_source_random = "LB20_random"
path_results_random = paste("/path_to_local/results/source_geometry_lm/", folder_source_random, sep = "")

separability_metrics = c("distance")
source('script_24_stats_separability_time_resolved_refined_random.R')

p_values_24a = p_values
t_values_24a = t_values

# no surrogate here
# p_values_24a_random = p_values_random
# t_values_24a_random = t_values_random

separability_metrics = c("volume")
source('script_24_stats_separability_time_resolved_refined_random.R')

p_values_24c = p_values
t_values_24c = t_values

# no surrogate here
# p_values_24c_random = p_values_random
# t_values_24c_random = t_values_random














### separability (distance) - length 3 vs length 4 in correct trials - time-resolved (bootstrap)

folder = "LB20"
path_results = paste("/path_to_local/results/source_geometry_lm/", folder, sep = "")
functions_data = c("mean")

separability_metrics = c("distance")
source('script_34_stats_separability_length_time_resolved_refined.R')
p_values_34a = p_values
t_values_34a = t_values

separability_metrics = c("volume")
source('script_34_stats_separability_length_time_resolved_refined.R')
p_values_34c = p_values
t_values_34c = t_values






### FDR correction

p_values_uncorrected = c(p_values_02, p_values_04, p_values_02_random, p_values_04_random,
                         p_values_06, p_values_08, 
                         p_values_12, p_values_16, 
                         p_values_24a, p_values_24c,
                         p_values_34a, p_values_34c)

t_values = c(t_values_02, t_values_04, t_values_02_random, t_values_04_random, 
             t_values_06, t_values_08, 
             t_values_12, t_values_16, 
             t_values_24a, t_values_24c, 
             t_values_34a, t_values_34c)


list_stats = c(rep('pa_bt_wt', length(p_values_02)), rep('vaf_bt_wt', length(p_values_04)),
               rep('pa_bt_su', length(p_values_02_random)), rep('vaf_bt_su', length(p_values_04_random)),
               rep('pa_perf', length(p_values_06)), rep('vaf_perf', length(p_values_08)),
               rep('pa_length', length(p_values_12)), rep('vaf_perf', length(p_values_16)),
               rep('distance_perf', length(p_values_24a)), rep('volume_perf', length(p_values_24c)),
               rep('distance_length', length(p_values_34a)), rep('volume_length', length(p_values_34c))
               )

p_values_fdr <- p.adjust(p_values_uncorrected, method = "fdr")
p_values_uncorrected = ifelse(p_values_uncorrected > 0.025, 1, p_values_uncorrected)
p_values_fdr = ifelse(p_values_fdr > 0.025, 1, p_values_fdr)
p_values = cbind(list_stats, round(p_values_uncorrected,5), round(p_values_fdr,5), t_values)

colnames(p_values) = c("stats", "uncorr", "fdr", "t_value")


# Save matrix to a .txt file
write.table(p_values, file = "p_values_distance_volume_LB20_noavg4lme.csv", sep = ",", row.names = FALSE, col.names = TRUE)


p_values_fdr_plots <- p.adjust(p_values_uncorrected, method = "fdr")

