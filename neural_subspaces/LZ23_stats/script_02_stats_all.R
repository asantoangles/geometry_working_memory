
# libraries
library(ggplot2)
library(lmerTest)
library(gridExtra)
library(cowplot)
library(GGally)
library(zoo)


setwd("/path_to_local/github/scripts/source_geometry_lm/LZ23_stats")

##################
##### stats ######
##################



### principal angle (between vs within) in stim and delay for correct trials time-resolved

folder_source = "LZ23"
path_results = paste("/path_to_local/github/results/source_geometry_lm/", folder_source, sep = "")

folder_source_random = "LZ23_random"
path_results_random = paste("/path_to_local/github/results/source_geometry_lm/", folder_source_random, sep = "")

performance = c("correct_trials")
functions_data = c("mean")
metrics = c("principal_angle")

source('script_03_stats_alignment_random.R')

p_values_003_pa = p_values
t_values_003_pa = t_values

p_values_003_pa_random = p_values_random
t_values_003_pa_random = t_values_random









### principal angle min (between vs within) in stim and delay for correct trials time-resolved

folder_source = "LZ23"
path_results = paste("/path_to_local/github/results/source_geometry_lm/", folder_source, sep = "")

folder_source_random = "LZ23_random"
path_results_random = paste("/path_to_local/github/results/source_geometry_lm/", folder_source_random, sep = "")

performance = c("correct_trials")
functions_data = c("mean")
metrics = c("principal_angle_min")

source('script_03_stats_alignment_random.R')

p_values_003_pa_min = p_values
t_values_003_pa_min = t_values

p_values_003_pa_min_random = p_values_random
t_values_003_pa_min_random = t_values_random








### vaf (between vs within) in stim and delay for correct trials time-resolved

folder_source = "LZ23"
path_results = paste("/path_to_local/github/results/source_geometry_lm/", folder_source, sep = "")

folder_source_random = "LZ23_random"
path_results_random = paste("/path_to_local/github/results/source_geometry_lm/", folder_source_random, sep = "")

performance = c("correct_trials")
functions_data = c("mean")
metrics = c("vaf")

source('script_03_stats_alignment_random.R')

p_values_003_vaf = p_values
t_values_003_vaf = t_values

p_values_003_vaf_random = p_values_random
t_values_003_vaf_random = t_values_random














### principal angle (correct vs incorrect) in stim and delay time-resolved

folder = "LZ23_controlled_resampling"
path_results = paste("/path_to_local/github/results/source_geometry_lm/", folder, sep = "")

functions_data = c("mean")
metrics = c("principal_angle")

source('script_04_stats_alignment_performance.R')

p_values_004_pa = p_values
t_values_004_pa = t_values








### principal angle min (correct vs incorrect) in stim and delay time-resolved

folder = "LZ23_controlled_resampling"
path_results = paste("/path_to_local/github/results/source_geometry_lm/", folder, sep = "")

functions_data = c("mean")
metrics = c("principal_angle_min")

source('script_04_stats_alignment_performance.R')

p_values_004_pa_min = p_values
t_values_004_pa_min = t_values








### vaf (correct vs incorrect) in stim and delay time-resolved

folder = "LZ23_controlled_resampling"
path_results = paste("/path_to_local/github/results/source_geometry_lm/", folder, sep = "")

functions_data = c("mean")
metrics = c("vaf")

source('script_04_stats_alignment_performance.R')

p_values_004_vaf = p_values
t_values_004_vaf = t_values






















### principal angle (between) (length 3 vs length 4) in correct trials in stim and delay time-resolved (bootstrap)

folder_source = "LZ23"
path_results = paste("/path_to_local/github/results/source_geometry_lm/", folder_source, sep = "")

performance = c("correct_trials")
functions_data = c("mean")
metrics = c("principal_angle")

source('script_05_stats_alignment_length.R')

p_values_005_pa = p_values
t_values_005_pa = t_values








### principal angle min (between) (length 3 vs length 4) in correct trials in stim and delay time-resolved (bootstrap)

folder_source = "LZ23"
path_results = paste("/path_to_local/github/results/source_geometry_lm/", folder_source, sep = "")

performance = c("correct_trials")
functions_data = c("mean")
metrics = c("principal_angle_min")

source('script_05_stats_alignment_length.R')

p_values_005_pa_min = p_values
t_values_005_pa_min = t_values











### vaf (length 3 vs length 4) in correct trials in stim and delay time-resolved (bootstrap)

folder_source = "LZ23"
path_results = paste("/path_to_local/github/results/source_geometry_lm/", folder_source, sep = "")

performance = c("correct_trials")
functions_data = c("mean")
metrics = c("vaf")

source('script_05_stats_alignment_length.R')

p_values_005_vaf = p_values
t_values_005_vaf = t_values

















### separability (distance) - correct vs incorrect trials - time-resolved (bootstrap)

folder = "LZ23_controlled_resampling"
path_results = paste("/path_to_local/github/results/source_geometry_lm/", folder, sep = "")

folder_source_random = "LZ23_random"
path_results_random = paste("/path_to_local/github/results/source_geometry_lm/", folder_source_random, sep = "")

separability_metrics = c("distance")
source('script_06_stats_separability.R')

p_values_006_distance = p_values
t_values_006_distance = t_values

# no surrogate here
# p_values_006_distance_random = p_values_random
# t_values_006_distance_random = t_values_random

separability_metrics = c("volume")
source('script_06_stats_separability.R')

p_values_006_volume = p_values
t_values_006_volume = t_values

# no surrogate here
# p_values_006_volume_random = p_values_random
# t_values_006_volume_random = t_values_random














### separability (distance) - length 3 vs length 4 in correct trials - time-resolved (bootstrap)

folder = "LZ23"
path_results = paste("/path_to_local/github/results/source_geometry_lm/", folder, sep = "")
functions_data = c("mean")

separability_metrics = c("distance")
source('script_07_stats_separability_length.R')
p_values_007_distance = p_values
t_values_007_distance = t_values

separability_metrics = c("volume")
source('script_07_stats_separability_length.R')
p_values_007_volume = p_values
t_values_007_volume = t_values






### FDR correction

p_values_uncorrected = c(p_values_003_pa, p_values_003_pa_min, p_values_003_vaf, 
                         p_values_003_pa_random, p_values_003_pa_min_random, p_values_003_vaf_random,
                         p_values_004_pa, p_values_004_pa_min, p_values_004_vaf, 
                         p_values_005_pa, p_values_005_pa_min, p_values_005_vaf, 
                         p_values_006_distance, p_values_006_volume,
                         p_values_007_distance, p_values_007_volume)

t_values = c(t_values_003_pa, t_values_003_pa_min, t_values_003_vaf, 
             t_values_003_pa_random, t_values_003_pa_min_random, t_values_003_vaf_random, 
             t_values_004_pa, t_values_004_pa_min, t_values_004_vaf, 
             t_values_005_pa, t_values_005_pa_min, t_values_005_vaf, 
             t_values_006_distance, t_values_006_volume, 
             t_values_007_distance, t_values_007_volume)


list_stats = c(rep('pa_bt_wt', length(p_values_003_pa)), rep('pa__min_bt_wt', length(p_values_003_pa_min)), rep('vaf_bt_wt', length(p_values_003_vaf)),
               rep('pa_bt_su', length(p_values_003_pa_random)), rep('pa_min_bt_su', length(p_values_003_pa_min_random)), rep('vaf_bt_su', length(p_values_003_vaf_random)),
               rep('pa_perf', length(p_values_004_pa)), rep('pa_min_perf', length(p_values_004_pa_min)), rep('vaf_perf', length(p_values_004_vaf)),
               rep('pa_length', length(p_values_005_pa)), rep('pa_min_length', length(p_values_005_pa_min)), rep('vaf_perf', length(p_values_005_vaf)),
               rep('distance_perf', length(p_values_006_distance)), rep('volume_perf', length(p_values_006_volume)),
               rep('distance_length', length(p_values_007_distance)), rep('volume_length', length(p_values_007_volume))
               )

p_values_fdr <- p.adjust(p_values_uncorrected, method = "fdr")
p_values_uncorrected = ifelse(p_values_uncorrected > 0.025, 1, p_values_uncorrected)
p_values_fdr = ifelse(p_values_fdr > 0.025, 1, p_values_fdr)
p_values = cbind(list_stats, round(p_values_uncorrected,5), round(p_values_fdr,5), t_values)

colnames(p_values) = c("stats", "uncorr", "fdr", "t_value")


# Save matrix to a .txt file
write.table(p_values, file = "p_values.csv", sep = ",", row.names = FALSE, col.names = TRUE)


p_values_fdr_plots <- p.adjust(p_values_uncorrected, method = "fdr")

