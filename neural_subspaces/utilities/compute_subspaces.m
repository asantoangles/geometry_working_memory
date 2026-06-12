function [z_k, z_k_expl, planes] = compute_subspaces(X_eigvecs, X, components, number_locations)

% COMPUTE_SUBSPACES  Project data into PCA subspaces and extract best-fitting planes.
%
%   [z_k, z_k_expl, planes] = compute_subspaces(X_eigvecs, X, components)
%
%   This function performs a two-stage dimensionality reduction and geometric
%   analysis of neural population activity. It first computes a PCA basis using
%   X_eigvecs and then projects a (possibly different) data matrix X into that
%   subspace. This allows, for example, projecting activity from incorrect
%   trials into principal components derived from correct trials.
%
%   The function additionally computes best-fitting planes for each item-rank
%   within the sequence by performing PCA separately on rank-specific,
%   location-demeaned population vectors.
%
%   INPUTS:
%       X_eigvecs   Data matrix used to compute PCA eigenvectors
%                   (conditions × channels).
%
%       X           Data matrix to be projected into the PCA space derived
%                   from X_eigvecs. Dimensions must match in the channel
%                   dimension.
%
%       components  Indices of principal components to retain (e.g., 1:3).
%
%   OUTPUTS:
%       z_k         Projection of X onto the selected PCA components.
%                   Size: (conditions × length(components)).
%
%       z_k_expl    Fraction of variance explained by each PCA component,
%                   computed from the eigenvalues of X_eigvecs.
%
%       planes       Struct containing, for each rank:
%                       • plane_r<k>     – 2D plane basis (first two PCs)
%                       • score_r<k>     – projected rank-specific scores
%                       • explained_r<k> – variance explained by the plane
%                       • Y_r<k>         – demeaned population vectors
%

% settings
ranks = size(X,1)/number_locations;

%% decomposition of X with SVD (PCA)

warning off

% Calculate the covariance matrix of the data
covariance_matrix = cov(X_eigvecs);

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