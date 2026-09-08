# geometry_working_memory

Code of manuscript: Santo-Angles A., Yang J., Zhou Y., Chu W.K.H., Lindsay G.W., Sreenivasan K.K. Neural Subspaces Encode Sequential Working Memory, but Neural Sequences Do Not. bioRxiv (2025). doi: https://doi.org/10.1101/2025.09.05.674385

To facilitate reproducibility, intermediate outputs from both the neural geometry and neural sequence analyses are available at the OSF repository https://osf.io/hvq8p.

For the neural geometry pipeline, we provide the neural activity matrices (X) used to compute the neural subspaces, as well as the resulting subspaces (principal component scores) derived from these matrices. For the neural sequences pipeline, we provide the decoded state spaces (classifier prediction outputs) that serve as inputs to the Temporal Delayed Linear Modelling (TDLM) analyses. The scripts in folders neural_subspaces and neural_sequences reproduce the figures reported in the manuscript.

To run the analyses, download the GitHub repository and place it in path_to_local/scripts, and download the OSF repository and place it in path_to_local/results.

The expected runtime of the entire pipeline is approximately 24-36 hours, depending on the computational environment and configuration. When starting from the intermediate outputs made available, the remaining pipeline takes approximately 2-4 hours to complete.

## software

The code was developed and tested using MATLAB R2022a on a Mac laptop and the NYU Abu Dhabi (NYUAD) HPC cluster. R version 4.2.1 was used on the Mac laptop. The software is made publicly available under the MIT License.

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

- LB23: Figure 2, 3, S2, S3, S5, S6

- LB24: Figure S4A

- LB25: Figure S4B

- LB21: Figure S4C

- LB27: Figure S4D

- LY23: Figure S7

- LZ23: Figure S8

- LC23: Figure S9 ABC

- LD23: Figure S9 DEF

- NB23_PC1: Figure S10 ABC

- NB23_PC1andPC2: Figure S10 DEF

- MB23: Figure S11 ABC

- JB23: Figure S11 DEF

- LH23: Figure S12 A

- LG23: Figure S12 B

- LI23: Figure S12 C

- RB23: Figure S13 A

- SB23: Figure S13 B


