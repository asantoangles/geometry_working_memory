# geometry_working_memory
Code of manuscript: Santo-Angles A., Yang J., Zhou Y., Chu W.K.H., Lindsay G.W., Sreenivasan K.K. Neural Subspaces Encode Sequential Working Memory, but Neural Sequences Do Not. bioRxiv (2025). doi: https://doi.org/10.1101/2025.09.05.674385

## task_design

Psychotoolbox code of the sequential working memory task and functional localizer.

## behavioral

Behavioral analysis of sequential working memory task.

## meg_preprocessing

Preprocessing of MEG data, based on FLUX pipeline (neuosc.com/flux)

## source_reconstruction

Source reconstruction using beamforming and cortical parcellation.

## replay

Analysis of neural sequences, using Temporally delayed linear modelling (TDLM). 

- Trial averaged approach (Figure 4B): sequence_trial_avg_delay and sequence_trial_avg_stim folders.
- Trial concatenate approach (Figure 4C): sequence_concat_delay and sequence_concat_stim folders.

## source_geometry_lm

Geometric analysis of neural subspaces.

Main results:

- LB20: main analysis of subspaces of dimensionality 3 PCs.
- LB20_controlled_resampling: resampling approach to control for imbalanced number of trials of alignment metrics (PA and VAF) and separability metrics (volume and distance) in the comparison between correct and incorrect trials.
- LB20_random: null distribution of subspace alignment metrics.
- LB20_surrogate: null distribution of eigenvalue spectra.

Supplementary results:

- LB20_allPCs_accumulated: analysis of subspaces of dimensionality 6 and 8 PCs.
- LB20_controlled_resampling_allPCs_accumulated: resampling approach to control for imbalanced number of trials for dimensionality 6 and 8 PCs.
- LB20_random_comp1to6: null distribution of subspace alignment metrics for subspaces of 6 PCs.
- LB20_random_comp1to8: null distribution of subspace alignment metrics for subspaces of 8 PCs.
- MB20: analysis of subspaces of dimensionality 3 PCs with 4 locations.
- MB20_controlled_resampling: resampling approach to control for imbalanced number of trials for subspaces with 4 locations.
- MB20_random: null distribution of subspace alignment metrics for subspaces with 4 locations.

Stats and figures:

- main_results_subspaces_3PCs: statistical analysis and figures of subspaces of dimensionality 3 PCs. Figure 2, Figure S2, Figure S3.
- supplementary_results_subspaces_6PCs: Figure S5.
- supplementary_results_subspaces_8PCs: Figure S6.
- supplementary_results_4locations: Figure S7.




