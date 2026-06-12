##### summary t-stats

### Figure 2 - LB23
### Figure S2 - LB23
### Figure S3 - LB23

## Principal Angle - Figure 2A

print('Principal Angle (between vs within) - length 3')
p_values = p_values_fdr_003_pa[1:15]
t_values = t_values_003_pa[1:15]
t_values = t_values * -1
summarize_significant_segments(t_values, p_values)


print('Principal Angle (between vs within) - length 4')
p_values = p_values_fdr_003_pa[16:31]
t_values = t_values_003_pa[16:31]
t_values = t_values * -1
summarize_significant_segments(t_values, p_values)


print('Principal Angle (between vs surrogate) - length 3')
p_values = p_values_fdr_003_pa_random[1:15]
t_values = t_values_003_pa_random[1:15]
t_values = t_values * -1
summarize_significant_segments(t_values, p_values)


print('Principal Angle (between vs surrogate) - length 4')
p_values = p_values_fdr_003_pa_random[16:31]
t_values = t_values_003_pa_random[16:31]
t_values = t_values * -1
summarize_significant_segments(t_values, p_values)



## VAF - Figure 2A


print('VAF (between vs within) - length 3')
p_values = p_values_fdr_003_vaf[1:15]
t_values = t_values_003_vaf[1:15]
t_values = t_values * -1
summarize_significant_segments(t_values, p_values)


print('VAF (between vs within) - length 4')
p_values = p_values_fdr_003_vaf[16:31]
t_values = t_values_003_vaf[16:31]
t_values = t_values * -1
summarize_significant_segments(t_values, p_values)


print('VAF (between vs surrogate) - length 3')
p_values = p_values_fdr_003_vaf_random[1:15]
t_values = t_values_003_vaf_random[1:15]
t_values = t_values * -1
summarize_significant_segments(t_values, p_values)


print('VAF (between vs surrogate) - length 4')
p_values = p_values_fdr_003_vaf_random[16:31]
t_values = t_values_003_vaf_random[16:31]
t_values = t_values * -1
summarize_significant_segments(t_values, p_values)



## Principal Angle - Figure 2B

print('Principal Angle (correct vs incorrect) - length 3')
p_values = p_values_fdr_004_pa[1:15]
t_values = t_values_004_pa[1:15]
t_values = t_values * -1
summarize_significant_segments(t_values, p_values)

print('Principal Angle (correct vs incorrect) - length 4')
p_values = p_values_fdr_004_pa[16:31]
t_values = t_values_004_pa[16:31]
t_values = t_values * -1
summarize_significant_segments(t_values, p_values)



print('VAF (correct vs incorrect) - length 3')
p_values = p_values_fdr_004_vaf[1:15]
t_values = t_values_004_vaf[1:15]
t_values = t_values * -1
summarize_significant_segments(t_values, p_values)

print('VAF (correct vs incorrect) - length 4')
p_values = p_values_fdr_004_vaf[16:31]
t_values = t_values_004_vaf[16:31]
t_values = t_values * -1
summarize_significant_segments(t_values, p_values)




## Figure S2


print('Volume (correct vs incorrect) - length 4')
p_values = p_values_fdr_006_volume[17:33]
t_values = t_values_006_volume[17:33]
t_values = t_values * -1
summarize_significant_segments(t_values, p_values)






## Figure S3


print('Principal Angle (length 3 vs length 4)')
p_values = p_values_fdr_005_pa
t_values = t_values_005_pa
t_values = t_values * -1
summarize_significant_segments(t_values, p_values)


print('VAF (length 3 vs length 4)')
p_values = p_values_fdr_005_vaf
t_values = t_values_005_vaf
t_values = t_values * -1
summarize_significant_segments(t_values, p_values)









### Figure S6 - LB24 / LB25 / LB21 / LB27

## Principal Angle

print('Principal Angle (between vs within) - length 3')
p_values = p_values_fdr_003_pa[1:15]
t_values = t_values_003_pa[1:15]
t_values = t_values * -1
t_values[t_values < 0] = 0
mask = ifelse(t_values > 0, 1, 0)
p_values = ifelse(mask == 0, 1, p_values)
summarize_significant_segments(t_values, p_values)


print('Principal Angle (between vs within) - length 4')
p_values = p_values_fdr_003_pa[16:31]
t_values = t_values_003_pa[16:31]
t_values = t_values * -1
t_values[t_values < 0] = 0
mask = ifelse(t_values > 0, 1, 0)
p_values = ifelse(mask == 0, 1, p_values)
summarize_significant_segments(t_values, p_values)


print('Principal Angle (between vs surrogate) - length 3')
p_values = p_values_fdr_003_pa_random[1:15]
t_values = t_values_003_pa_random[1:15]
t_values = t_values * -1
t_values[t_values < 0] = 0
mask = ifelse(t_values > 0, 1, 0)
p_values = ifelse(mask == 0, 1, p_values)
summarize_significant_segments(t_values, p_values)


print('Principal Angle (between vs surrogate) - length 4')
p_values = p_values_fdr_003_pa_random[16:31]
t_values = t_values_003_pa_random[16:31]
t_values = t_values * -1
t_values[t_values < 0] = 0
mask = ifelse(t_values > 0, 1, 0)
p_values = ifelse(mask == 0, 1, p_values)
summarize_significant_segments(t_values, p_values)



## VAF


print('VAF (between vs within) - length 3')
p_values = p_values_fdr_003_vaf[1:15]
t_values = t_values_003_vaf[1:15]
t_values[t_values < 0] = 0
mask = ifelse(t_values > 0, 1, 0)
p_values = ifelse(mask == 0, 1, p_values)
summarize_significant_segments(t_values, p_values)


print('VAF (between vs within) - length 4')
p_values = p_values_fdr_003_vaf[16:31]
t_values = t_values_003_vaf[16:31]
t_values[t_values < 0] = 0
mask = ifelse(t_values > 0, 1, 0)
p_values = ifelse(mask == 0, 1, p_values)
summarize_significant_segments(t_values, p_values)


print('VAF (between vs surrogate) - length 3')
p_values = p_values_fdr_003_vaf_random[1:15]
t_values = t_values_003_vaf_random[1:15]
t_values[t_values < 0] = 0
mask = ifelse(t_values > 0, 1, 0)
p_values = ifelse(mask == 0, 1, p_values)
summarize_significant_segments(t_values, p_values)


print('VAF (between vs surrogate) - length 4')
p_values = p_values_fdr_003_vaf_random[16:31]
t_values = t_values_003_vaf_random[16:31]
t_values[t_values < 0] = 0
mask = ifelse(t_values > 0, 1, 0)
p_values = ifelse(mask == 0, 1, p_values)
summarize_significant_segments(t_values, p_values)






### Figure S9 - LY23
### Figure S10 - LZ23

## Principal Angle

print('Principal Angle (between vs within) - length 3')
p_values = p_values_fdr_003_pa[1:15]
t_values = t_values_003_pa[1:15]
t_values = t_values * -1
t_values[t_values < 0] = 0
mask = ifelse(t_values > 0, 1, 0)
p_values = ifelse(mask == 0, 1, p_values)
summarize_significant_segments(t_values, p_values)


print('Principal Angle (between vs within) - length 4')
p_values = p_values_fdr_003_pa[16:31]
t_values = t_values_003_pa[16:31]
t_values = t_values * -1
t_values[t_values < 0] = 0
mask = ifelse(t_values > 0, 1, 0)
p_values = ifelse(mask == 0, 1, p_values)
summarize_significant_segments(t_values, p_values)


print('Principal Angle (between vs surrogate) - length 3')
p_values = p_values_fdr_003_pa_random[1:15]
t_values = t_values_003_pa_random[1:15]
t_values = t_values * -1
t_values[t_values < 0] = 0
mask = ifelse(t_values > 0, 1, 0)
p_values = ifelse(mask == 0, 1, p_values)
summarize_significant_segments(t_values, p_values)


print('Principal Angle (between vs surrogate) - length 4')
p_values = p_values_fdr_003_pa_random[16:31]
t_values = t_values_003_pa_random[16:31]
t_values = t_values * -1
t_values[t_values < 0] = 0
mask = ifelse(t_values > 0, 1, 0)
p_values = ifelse(mask == 0, 1, p_values)
summarize_significant_segments(t_values, p_values)




## Principal Angle min

print('Principal Angle (between vs within) - length 3')
p_values = p_values_fdr_003_pa_min[1:15]
t_values = t_values_003_pa_min[1:15]
t_values = t_values * -1
t_values[t_values < 0] = 0
mask = ifelse(t_values > 0, 1, 0)
p_values = ifelse(mask == 0, 1, p_values)
summarize_significant_segments(t_values, p_values)


print('Principal Angle (between vs within) - length 4')
p_values = p_values_fdr_003_pa_min[16:31]
t_values = t_values_003_pa_min[16:31]
t_values = t_values * -1
t_values[t_values < 0] = 0
mask = ifelse(t_values > 0, 1, 0)
p_values = ifelse(mask == 0, 1, p_values)
summarize_significant_segments(t_values, p_values)


print('Principal Angle (between vs surrogate) - length 3')
p_values = p_values_fdr_003_pa_min_random[1:15]
t_values = t_values_003_pa_min_random[1:15]
t_values = t_values * -1
t_values[t_values < 0] = 0
mask = ifelse(t_values > 0, 1, 0)
p_values = ifelse(mask == 0, 1, p_values)
summarize_significant_segments(t_values, p_values)


print('Principal Angle (between vs surrogate) - length 4')
p_values = p_values_fdr_003_pa_min_random[16:31]
t_values = t_values_003_pa_min_random[16:31]
t_values = t_values * -1
t_values[t_values < 0] = 0
mask = ifelse(t_values > 0, 1, 0)
p_values = ifelse(mask == 0, 1, p_values)
summarize_significant_segments(t_values, p_values)





## VAF


print('VAF (between vs within) - length 3')
p_values = p_values_fdr_003_vaf[1:15]
t_values = t_values_003_vaf[1:15]
t_values = t_values * -1
summarize_significant_segments(t_values, p_values)


print('VAF (between vs within) - length 4')
p_values = p_values_fdr_003_vaf[16:31]
t_values = t_values_003_vaf[16:31]
t_values = t_values * -1
summarize_significant_segments(t_values, p_values)


print('VAF (between vs surrogate) - length 3')
p_values = p_values_fdr_003_vaf_random[1:15]
t_values = t_values_003_vaf_random[1:15]
t_values = t_values * -1
summarize_significant_segments(t_values, p_values)


print('VAF (between vs surrogate) - length 4')
p_values = p_values_fdr_003_vaf_random[16:31]
t_values = t_values_003_vaf_random[16:31]
t_values = t_values * -1
summarize_significant_segments(t_values, p_values)



## Principal Angle (correct vs incorrect)

print('Principal Angle (correct vs incorrect) - length 3')
p_values = p_values_fdr_004_pa[1:15]
t_values = t_values_004_pa[1:15]
t_values = t_values * -1
summarize_significant_segments(t_values, p_values)

print('Principal Angle (correct vs incorrect) - length 4')
p_values = p_values_fdr_004_pa[16:31]
t_values = t_values_004_pa[16:31]
t_values = t_values * -1
summarize_significant_segments(t_values, p_values)

## Principal Angle min (correct vs incorrect)

print('Principal Angle (correct vs incorrect) - length 3')
p_values = p_values_fdr_004_pa_min[1:15]
t_values = t_values_004_pa_min[1:15]
t_values = t_values * -1
summarize_significant_segments(t_values, p_values)

print('Principal Angle (correct vs incorrect) - length 4')
p_values = p_values_fdr_004_pa_min[16:31]
t_values = t_values_004_pa_min[16:31]
t_values = t_values * -1
summarize_significant_segments(t_values, p_values)

## VAF (correct vs incorrect)


print('VAF (correct vs incorrect) - length 3')
p_values = p_values_fdr_004_vaf[1:15]
t_values = t_values_004_vaf[1:15]
t_values = t_values * -1
summarize_significant_segments(t_values, p_values)

print('VAF (correct vs incorrect) - length 4')
p_values = p_values_fdr_004_vaf[16:31]
t_values = t_values_004_vaf[16:31]
t_values = t_values * -1
summarize_significant_segments(t_values, p_values)


## Separability

print('Volume (correct vs incorrect) - length 3')
p_values = p_values_fdr_006_volume[1:16]
t_values = t_values_006_volume[1:16]
t_values = t_values * -1
summarize_significant_segments(t_values, p_values)

print('Volume (correct vs incorrect) - length 4')
p_values = p_values_fdr_006_volume[17:33]
t_values = t_values_006_volume[17:33]
t_values = t_values * -1
summarize_significant_segments(t_values, p_values)






### Figure S11 ABC - LC23
### Figure S11 DEF - LD23
### Figure S12 ABC - NB23_PC1
### Figure S12 DEF - NB23_PC1andPC2
### Figure S13 ABC - MB23
### Figure S14 ABC - JB23


## Principal Angle

print('Principal Angle (between vs within) - length 3')
p_values = p_values_fdr_003_pa[1:15]
t_values = t_values_003_pa[1:15]
t_values = t_values * -1
summarize_significant_segments(t_values, p_values)


print('Principal Angle (between vs within) - length 4')
p_values = p_values_fdr_003_pa[16:31]
t_values = t_values_003_pa[16:31]
t_values = t_values * -1
summarize_significant_segments(t_values, p_values)


print('Principal Angle (between vs surrogate) - length 3')
p_values = p_values_fdr_003_pa_random[1:15]
t_values = t_values_003_pa_random[1:15]
t_values = t_values * -1
summarize_significant_segments(t_values, p_values)


print('Principal Angle (between vs surrogate) - length 4')
p_values = p_values_fdr_003_pa_random[16:31]
t_values = t_values_003_pa_random[16:31]
t_values = t_values * -1
summarize_significant_segments(t_values, p_values)



## VAF


print('VAF (between vs within) - length 3')
p_values = p_values_fdr_003_vaf[1:15]
t_values = t_values_003_vaf[1:15]
t_values = t_values * -1
summarize_significant_segments(t_values, p_values)


print('VAF (between vs within) - length 4')
p_values = p_values_fdr_003_vaf[16:31]
t_values = t_values_003_vaf[16:31]
t_values = t_values * -1
summarize_significant_segments(t_values, p_values)


print('VAF (between vs surrogate) - length 3')
p_values = p_values_fdr_003_vaf_random[1:15]
t_values = t_values_003_vaf_random[1:15]
t_values = t_values * -1
summarize_significant_segments(t_values, p_values)


print('VAF (between vs surrogate) - length 4')
p_values = p_values_fdr_003_vaf_random[16:31]
t_values = t_values_003_vaf_random[16:31]
t_values = t_values * -1
summarize_significant_segments(t_values, p_values)



## Principal Angle - performance

print('Principal Angle (correct vs incorrect) - length 3')
p_values = p_values_fdr_004_pa[1:15]
t_values = t_values_004_pa[1:15]
t_values = t_values * -1
summarize_significant_segments(t_values, p_values)

print('Principal Angle (correct vs incorrect) - length 4')
p_values = p_values_fdr_004_pa[16:31]
t_values = t_values_004_pa[16:31]
t_values = t_values * -1
summarize_significant_segments(t_values, p_values)

## VAF - performance


print('VAF (correct vs incorrect) - length 3')
p_values = p_values_fdr_004_vaf[1:15]
t_values = t_values_004_vaf[1:15]
t_values = t_values * -1
summarize_significant_segments(t_values, p_values)

print('VAF (correct vs incorrect) - length 4')
p_values = p_values_fdr_004_vaf[16:31]
t_values = t_values_004_vaf[16:31]
t_values = t_values * -1
summarize_significant_segments(t_values, p_values)




## Separability

print('Volume (correct vs incorrect) - length 3')
p_values = p_values_fdr_006_volume[1:16]
t_values = t_values_006_volume[1:16]
t_values = t_values * -1
summarize_significant_segments(t_values, p_values)

print('Volume (correct vs incorrect) - length 4')
p_values = p_values_fdr_006_volume[17:33]
t_values = t_values_006_volume[17:33]
t_values = t_values * -1
summarize_significant_segments(t_values, p_values)









### Figure S15A - LH23
### Figure S15B - LG23
### Figure S15C - LI23


## Principal Angle - performance

print('Principal Angle (correct vs incorrect) - length 3')
p_values = p_values_fdr_004_pa[1:15]
t_values = t_values_004_pa[1:15]
t_values = t_values * -1
summarize_significant_segments(t_values, p_values)

print('Principal Angle (correct vs incorrect) - length 4')
p_values = p_values_fdr_004_pa[16:31]
t_values = t_values_004_pa[16:31]
t_values = t_values * -1
summarize_significant_segments(t_values, p_values)

## VAF - performance


print('VAF (correct vs incorrect) - length 3')
p_values = p_values_fdr_004_vaf[1:15]
t_values = t_values_004_vaf[1:15]
t_values = t_values * -1
summarize_significant_segments(t_values, p_values)

print('VAF (correct vs incorrect) - length 4')
p_values = p_values_fdr_004_vaf[16:31]
t_values = t_values_004_vaf[16:31]
t_values = t_values * -1
summarize_significant_segments(t_values, p_values)




## Separability

print('Distance (correct vs incorrect) - length 3')
p_values = p_values_fdr_006_distance[1:16]
t_values = t_values_006_distance[1:16]
t_values = t_values * -1
summarize_significant_segments(t_values, p_values)

print('Distance (correct vs incorrect) - length 4')
p_values = p_values_fdr_006_distance[17:33]
t_values = t_values_006_distance[17:33]
t_values = t_values * -1
summarize_significant_segments(t_values, p_values)


print('Volume (correct vs incorrect) - length 3')
p_values = p_values_fdr_006_volume[1:16]
t_values = t_values_006_volume[1:16]
t_values = t_values * -1
summarize_significant_segments(t_values, p_values)

print('Volume (correct vs incorrect) - length 4')
p_values = p_values_fdr_006_volume[17:33]
t_values = t_values_006_volume[17:33]
t_values = t_values * -1
summarize_significant_segments(t_values, p_values)

















### function

summarize_significant_segments = function(t_values, p_values, threshold = 0.025){
  
  n = length(t_values)
  
  # Base time vector
  base_time0 = seq(0, 3.6, by = 0.3)
  
  # Add preceding values if needed
  if(n > length(base_time0)){
    
    n_extra = n - length(base_time0)
    
    extra_time0 = seq(
      from = -0.3 * n_extra,
      to = -0.3,
      by = 0.3
    )
    
    time0 = c(extra_time0, base_time0)
    
  } else {
    time0 = base_time0[1:n]
  }
  
  time1 = time0 + 0.3
  
  # Significant timepoints
  sig = p_values < threshold
  
  # Identify contiguous segments
  r = rle(sig)
  
  ends = cumsum(r$lengths)
  starts = ends - r$lengths + 1
  
  sig_segments = which(r$values)
  
  results = data.frame(
    interval = character(),
    peak_t = numeric(),
    peak_latency = numeric(),
    t_range = character(),
    p_min = numeric(),
    stringsAsFactors = FALSE
  )
  
  for(i in sig_segments){
    
    idx = starts[i]:ends[i]
    
    t_seg = t_values[idx]
    p_seg = p_values[idx]
    
    # Interval boundaries
    interval_start = time0[min(idx)]
    interval_end   = time1[max(idx)]
    
    interval = paste0(
      "[",
      sprintf("%.2f", interval_start),
      " ",
      sprintf("%.2f", interval_end),
      "]"
    )
    
    # Peak t
    peak_idx_local = which.max(abs(t_seg))
    peak_idx_global = idx[peak_idx_local]
    
    peak_t = t_values[peak_idx_global]
    
    # Peak latency (midpoint)
    peak_latency = mean(c(
      time0[peak_idx_global],
      time1[peak_idx_global]
    ))
    
    # t-range
    t_range = paste0(
      "[",
      sprintf("%.2f", min(t_seg)),
      " ",
      sprintf("%.2f", max(t_seg)),
      "]"
    )
    
    # NEW: minimum p-value in interval
    p_min = min(p_seg, na.rm = TRUE)
    
    results[nrow(results) + 1, ] = list(
      interval,
      round(peak_t, 2),
      round(peak_latency, 2),
      t_range,
      signif(p_min, 3)
    )
  }
  
  # Excel-friendly output
  cat(
    paste(
      c("Significant interval (s)",
        "Peak t",
        "Peak latency (s)",
        "t-range",
        "p-min"),
      collapse = "\t"
    ),
    "\n"
  )
  
  apply(results, 1, function(x){
    cat(paste(x, collapse = "\t"), "\n")
  })
  
  invisible(results)
}