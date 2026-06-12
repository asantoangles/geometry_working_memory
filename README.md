# geometry_working_memory

Code of manuscript: Santo-Angles A., Yang J., Zhou Y., Chu W.K.H., Lindsay G.W., Sreenivasan K.K. Neural Subspaces Encode Sequential Working Memory, but Neural Sequences Do Not. bioRxiv (2025). doi: https://doi.org/10.1101/2025.09.05.674385

This repository contains all scripts used in Santo-Angles et al. (2025).

To facilitate reproducibility, intermediate outputs from both the neural geometry and neural sequence analyses are available at the OSF repository https://osf.io/hvq8p.

For the neural geometry pipeline, we provide the neural activity matrices (X) used to compute the neural subspaces, as well as the resulting subspaces (principal component scores) derived from these matrices. For the neural sequences pipeline, we provide the decoded state spaces (classifier prediction outputs) that serve as inputs to the Temporal Delayed Linear Modelling (TDLM) analyses. The scripts in folders neural_subspaces and neural_sequences reproduce the figures reported in the manuscript.

To run the analyses, download the GitHub repository and place it in path_to_local/scripts, and download the OSF repository and place it in path_to_local/results.

## task_design

Psychotoolbox code of the sequential working memory task and functional localizer.

## behavioral

Behavioral analysis of sequential working memory task.

## meg_preprocessing

Preprocessing of MEG data, based on FLUX pipeline (neuosc.com/flux)

## source_reconstruction

Source reconstruction using beamforming and cortical parcellation.

## neural_sequences

Analysis of neural sequences, using Temporally Delayed Linear Modelling (TDLM). 

- Across-trials concatenated TDLM approach (Figure 4B): sequence_concat_delay and sequence_concat_stim folders.

- Trial-averaged TDLM approach (Figure 4C): sequence_trial_avg_delay and sequence_trial_avg_stim folders.

## neural_subspaces

Geometric analysis of neural subspaces.

Folders define distinct blocks of analysis. For example, LB23 contains empirical results for a given analysis block. Folders with the suffix _controlled_resampling contain resampling analyses used to compare correct and incorrect trials. Folders with the suffix _random contain surrogate datasets used for geometry control analyses. Finally, folders with the suffix _stats_figures contain statistical analyses and figure-generation scripts, which take outputs from the other folders as input.

- LB23: Figure 2, 3, S2, S3, S4, S5, S7, S8

- LB24: Figure S6A

- LB25: Figure S6B

- LB21: Figure S6C

- LB27: Figure S6D

- LY23: Figure S9

- LZ23: Figure S10

- LC23: Figure S11 ABC

- LD23: Figure S11 DEF

- NB23_PC1: Figure S12 ABC

- NB23_PC1andPC2: Figure S12 DEF

- MB23: Figure S13

- LB23: Figure S14

- LH23: Figure S15 A

- LG23: Figure S15 B

- LI23: Figure S15 C

- RB23: Figure S16 A

- SB23: Figure S16 B
