function [z_k, z_k_expl, planes] = compute_plane_components_alltrials(X_all_trials, X, components)

% compute_plane_components_alltrials performs a dimensionality reduction
% and subspace decomposition of neural population activity across all
% trials, identifying low-dimensional planes that best capture structured
% variance within subsets of conditions.
%
% SYNTAX:
%   [z_k, z_k_expl, planes] = compute_plane_components_alltrials(X_all_trials, X, components)
%
% INPUTS:
%   X_all_trials : [N × D] matrix
%       Data matrix containing all trials or conditions, where N is the
%       number of observations and D is the number of recorded features
%       (e.g., MEG sensors, electrodes, or neural units). Used to estimate
%       the full covariance structure of the data for PCA.
%
%   X : [M × D] matrix
%       Subset of trials or conditions (e.g., averaged or condition-specific
%       activity) to be projected onto the low-dimensional space defined by
%       the principal components of X_all_trials.
%
%   components : vector or integer
%       Indices of the principal components to retain (e.g., 1:3 for the
%       first three components). These define the reduced basis used for
%       projection and subsequent plane estimation.
%
% OUTPUTS:
%   z_k : [M × length(components)] matrix
%       Low-dimensional representation of X obtained by projecting it onto
%       the selected principal components of X_all_trials.
%
%   z_k_expl : [D × 1] vector
%       Fraction of variance explained by each component, computed from the
%       squared singular values of the SVD of the covariance matrix.
%
%   planes : struct
%       Structure containing, for each identified subspace (rank r):
%           - planes.Y_r#: demeaned population activity matrix for subset r
%           - planes.plane_r#: two principal component vectors defining the
%             best-fitting plane for subset r
%           - planes.score_r#: projections of Y_r onto the best-fitting plane
%           - planes.explained_r#: variance explained by the two plane components
%
% DESCRIPTION:
%   The function first computes the covariance matrix of X_all_trials and
%   applies singular value decomposition (SVD) to extract the principal
%   components (equivalent to PCA). The data matrix X is then projected
%   onto a reduced basis defined by the selected components, producing
%   z_k — a low-dimensional representation of neural activity.
%
%   The projected data are divided into subsets (ranks), each corresponding
%   to a distinct condition group (ordinal position in a sequence or rank). 
%   For each subset, the mean across conditions
%   is subtracted (demeaning), and a second PCA is applied
%   to extract the two eigenvectors that best define the plane capturing
%   local structure within that subset.
%
%   The number of subsets (ranks) is determined adaptively:
%       - If fewer than 20 conditions: 4 conditions per subset (N/4).
%       - If 20 or more conditions: 8 conditions per subset (N/8).
%
%   Each plane thus represents a two-dimensional subspace that locally
%   approximates the neural manifold structure across groups of trials or
%   conditions.
%
%
% EXAMPLE:
%   % Perform dimensionality reduction and plane estimation
%   [z_k, z_k_expl, planes] = compute_plane_components_alltrials(X_all, X, 1:3);
%
% REFERENCES:
%   Santo-Angles A., Yang J., Zhou Y., Chu W.K.H., Lindsay G.W., Sreenivasan K.K. 
%   Neural Subspaces Encode Sequential Working Memory, but Neural Sequences Do Not. 
%   bioRxiv (2025). doi: https://doi.org/10.1101/2025.09.05.674385

% settings
if size(X,1) < 20
    ranks = size(X,1)/4;
else
    ranks = size(X,1)/8;
end

%% decomposition of X with SVD (PCA)

warning off

% Calculate the covariance matrix of the data
covariance_matrix = cov(X_all_trials);

% Perform Singular Value Decomposition (SVD) on the covariance matrix
[U, S, ~] = svd(covariance_matrix);

% Choose the first three columns of U as the reduced basis
U_reduced = U(:, components);

% Project the data onto the reduced basis
z_k = X * U_reduced;

% Calculate the singular values from SVD
singular_values = diag(S);

% Compute the total sum of squared singular values
total_variance = sum(singular_values.^2);

% Compute the variance explained by each component
z_k_expl = (singular_values.^2) / total_variance;

%% identify the best-fitting planes

% new population activity matrix Y_r in which z_k is the population
% vector for condition location l and rank r in reduced
% dimensionality space, and demeaned by the average across
% locations in each rank (that is, the mean of each column of Y_r
% is zero)

if size(X,1) < 20

    Y_r1 = z_k(1:4,:) - mean(z_k(1:4,:),1);
    Y_r2 = z_k(5:8,:) - mean(z_k(5:8,:),1);
    Y_r3 = z_k(9:12,:) - mean(z_k(9:12,:),1);
    if ranks == 4
        Y_r4 = z_k(13:16,:) - mean(z_k(13:16,:),1);
    end

else

    Y_r1 = z_k(1:8,:) - mean(z_k(1:8,:),1);
    Y_r2 = z_k(9:16,:) - mean(z_k(9:16,:),1);
    Y_r3 = z_k(17:24,:) - mean(z_k(17:24,:),1);
    if ranks == 4
        Y_r4 = z_k(25:32,:) - mean(z_k(25:32,:),1);
    end

end

% pca on T_r<rank> (the first two eigenvectors identify the best
% fitting plane
[coeff, score, ~, ~, explained, ~] = pca(Y_r1);
plane_r1 = coeff(:,1:2);
score_r1 = score(:,1:2);
explained_r1 = sum(explained(1:2));

[coeff, score, ~, ~, explained, ~] = pca(Y_r2);
plane_r2 = coeff(:,1:2);
score_r2 = score(:,1:2);
explained_r2 = sum(explained(1:2));

[coeff, score, ~, ~, explained, ~] = pca(Y_r3);
plane_r3 = coeff(:,1:2);
score_r3 = score(:,1:2);
explained_r3 = sum(explained(1:2));

if ranks == 4
    [coeff, score, ~, ~, explained, ~] = pca(Y_r4);
    plane_r4 = coeff(:,1:2);
    score_r4 = score(:,1:2);
    explained_r4 = sum(explained(1:2));
end

planes = [];
planes.plane_r1 = plane_r1;
planes.plane_r2 = plane_r2;
planes.plane_r3 = plane_r3;
planes.score_r1 = score_r1;
planes.score_r2 = score_r2;
planes.score_r3 = score_r3;
planes.explained_r1 = explained_r1;
planes.explained_r2 = explained_r2;
planes.explained_r3 = explained_r3;
planes.Y_r1 = Y_r1;
planes.Y_r2 = Y_r2;
planes.Y_r3 = Y_r3;

if ranks == 4
    planes.plane_r4 = plane_r4;
    planes.explained_r4 = explained_r4;
    planes.score_r4 = score_r4;
    planes.Y_r4 = Y_r4;
end

warning on


end