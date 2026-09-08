function [planes_split1, planes_split2, Y_r, score] = compute_subspaces_within(X_all_trials, X_split1, X_split2, components, number_locations)

% COMPUTE_SUBSPACES_WITHIN  Split-half PCA projection and rank-specific plane estimation.
%
%   [planes_split1, planes_split2, Y_r, score] = ...
%       compute_subspaces_within(X_all_trials, X_split1, X_split2, components)
%
%   This function performs a *within-condition split-half analysis* to assess
%   the reliability of neural subspaces across two halves of the data. It uses
%   PCA computed from **all trials** (X_all_trials) to define a shared set of
%   principal axes, and then projects two independent splits of the data
%   (X_split1 and X_split2) into that common PCA space.
%
%   After projection, the function extracts rank-specific, location-demeaned
%   population vectors and fits a 2-D plane to each rank via PCA. This yields
%   two independent sets of subspaces (one per split).
%
%   ===========================================================
%   INPUTS
%   ===========================================================
%
%   X_all_trials     Full dataset used to compute PCA eigenvectors.
%                    (conditions × channels)
%
%   X_split1         First split-half dataset projected into PCA space.
%   X_split2         Second split-half dataset projected into PCA space.
%                    Splits must match X_all_trials in channel dimension.
%
%   components       Indices of principal components to retain 
%                    (typically 1:3).
%
%   ===========================================================
%   OUTPUTS
%   ===========================================================
%
%   planes_split1    Struct with best-fitting 2-D planes for ranks 1–R
%                    estimated from split 1:
%                       • plane_r1, plane_r2, ...
%
%   planes_split2    Same as above, estimated from split 2.
%
%   Y_r              Struct containing demeaned, rank-specific population
%                    activity for both splits:
%                       • Y_r1_split1, Y_r2_split1, ...
%                       • Y_r1_split2, Y_r2_split2, ...
%
%   score            Struct containing PCA scores in reduced space
%                    for each rank and each split:
%                       • score_r1_split1, score_r2_split1, ...
%                       • score_r1_split2, score_r2_split2, ...
%
%

% settings
ranks = size(X_all_trials,1)/number_locations;

%% decomposition of X with SVD (PCA)

warning off

for split_i = 1:2

    % set inputs
    if split_i == 1

        X = X_split1;

    else

        X = X_split2;

    end

    % Calculate the covariance matrix of the data
    covariance_matrix = cov(X_all_trials);
    
    % Perform Singular Value Decomposition (SVD) on the covariance matrix
    [U, S, ~] = svd(covariance_matrix);
    
    % Choose the first three columns of U as the reduced basis
    U_reduced = U(:, components);
    
    % Project the data onto the reduced basis
    z_k = X * U_reduced;

    
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
        Y_r3 = z_k(17:24,:) - mean(z_k(1:24,:),1);
        if ranks == 4
            Y_r4 = z_k(25:32,:) - mean(z_k(25:32,:),1);
        end
    
    end
    
    % pca on T_r<rank> (the first two eigenvectors identify the best
    % fitting plane
    [coeff, score] = pca(Y_r1);
    plane_r1 = coeff(:,1:2);
    score_r1 = score(:,1:2);
    
    [coeff, score] = pca(Y_r2);
    plane_r2 = coeff(:,1:2);
    score_r2 = score(:,1:2);
    
    [coeff, score] = pca(Y_r3);
    plane_r3 = coeff(:,1:2);
    score_r3 = score(:,1:2);
    
    if ranks == 4
        [coeff, score] = pca(Y_r4);
        plane_r4 = coeff(:,1:2);
        score_r4 = score(:,1:2);
    end
    
    planes = [];
    planes.plane_r1 = plane_r1;
    planes.plane_r2 = plane_r2;
    planes.plane_r3 = plane_r3;
    if ranks == 4
        planes.plane_r4 = plane_r4;
    end
    
    % set outputss
    if split_i == 1

        z_k_split1 = z_k;
        planes_split1 = planes;

        Y_r1_split1 = Y_r1;
        Y_r2_split1 = Y_r2;
        Y_r3_split1 = Y_r3;
        if ranks == 4
            Y_r4_split1 = Y_r4;
        end

        score_r1_split1 = score_r1;
        score_r2_split1 = score_r2;
        score_r3_split1 = score_r3;
        if ranks == 4
            score_r4_split1 = score_r4;
        end

    else

        z_k_split2 = z_k;
        planes_split2 = planes;
        Y_r1_split2 = Y_r1;
        Y_r2_split2 = Y_r2;
        Y_r3_split2 = Y_r3;
        if ranks == 4
            Y_r4_split2 = Y_r4;
        end

        score_r1_split2 = score_r1;
        score_r2_split2 = score_r2;
        score_r3_split2 = score_r3;
        if ranks == 4
            score_r4_split2 = score_r4;
        end

    end

end

warning on

% merge outputs
Y_r = [];
Y_r.Y_r1_split1 = Y_r1_split1;
Y_r.Y_r2_split1 = Y_r2_split1;
Y_r.Y_r3_split1 = Y_r3_split1;
Y_r.Y_r1_split2 = Y_r1_split2;
Y_r.Y_r2_split2 = Y_r2_split2;
Y_r.Y_r3_split2 = Y_r3_split2;

if ranks == 4
    Y_r.Y_r4_split1 = Y_r4_split1;
    Y_r.Y_r4_split2 = Y_r4_split2;
end

score = [];
score.score_r1_split1 = score_r1_split1;
score.score_r2_split1 = score_r2_split1;
score.score_r3_split1 = score_r3_split1;
score.score_r1_split2 = score_r1_split2;
score.score_r2_split2 = score_r2_split2;
score.score_r3_split2 = score_r3_split2;

if ranks == 4
    score.score_r4_split1 = score_r4_split1;
    score.score_r4_split2 = score_r4_split2;
end




end