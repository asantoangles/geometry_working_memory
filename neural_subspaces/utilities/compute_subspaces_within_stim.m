function [planes_split1, planes_split2, Y_r, score] = compute_subspaces_within_stim(X_all_trials, X_split1, X_split2, components, stim_order)

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

    % coeff = pca(X);
    % % coeff are eigenvectors
    % % latent are eigenvalues (sorted)
    % 
    % Pk = coeff(:,components);
    % 
    % % projection of data for a given condition into the reduced
    % % dimensionality space
    % 
    % z_k = nan(size(X,1),length(components));
    % 
    % for cond_i = 1:size(X, 1)
    %     z_k(cond_i,:) = Pk' * X(cond_i,:)';
    % end
    
    %% identify the best-fitting planes
    
    % new population activity matrix Y_r in which z_k is the population
    % vector for condition location l and rank r in reduced
    % dimensionality space, and demeaned by the average across
    % locations in each rank (that is, the mean of each column of Y_r
    % is zero)
    
    if size(X_all_trials,1) < 20
    
        if stim_order == 1
    
            Y_r1 = z_k(1:4,:) - mean(z_k(1:4,:),1);
    
        elseif stim_order == 2
    
            Y_r1 = z_k(1:4,:) - mean(z_k(1:4,:),1);
            Y_r2 = z_k(5:8,:) - mean(z_k(5:8,:),1);
    
        elseif stim_order == 3
    
            Y_r1 = z_k(1:4,:) - mean(z_k(1:4,:),1);
            Y_r2 = z_k(5:8,:) - mean(z_k(5:8,:),1);
            Y_r3 = z_k(9:12,:) - mean(z_k(9:12,:),1);
    
        elseif stim_order == 4
    
            Y_r1 = z_k(1:4,:) - mean(z_k(1:4,:),1);
            Y_r2 = z_k(5:8,:) - mean(z_k(5:8,:),1);
            Y_r3 = z_k(9:12,:) - mean(z_k(9:12,:),1);
            Y_r4 = z_k(13:16,:) - mean(z_k(13:16,:),1);
    
        end
    
    
    
    else
    
        if stim_order == 1
    
            Y_r1 = z_k(1:8,:) - mean(z_k(1:8,:),1);
    
        elseif stim_order == 2
    
            Y_r1 = z_k(1:8,:) - mean(z_k(1:8,:),1);
            Y_r2 = z_k(9:16,:) - mean(z_k(9:16,:),1);
    
        elseif stim_order == 3
    
            Y_r1 = z_k(1:8,:) - mean(z_k(1:8,:),1);
            Y_r2 = z_k(9:16,:) - mean(z_k(9:16,:),1);
            Y_r3 = z_k(17:24,:) - mean(z_k(1:24,:),1);
    
        elseif stim_order == 4
    
            Y_r1 = z_k(1:8,:) - mean(z_k(1:8,:),1);
            Y_r2 = z_k(9:16,:) - mean(z_k(9:16,:),1);
            Y_r3 = z_k(17:24,:) - mean(z_k(1:24,:),1);
            Y_r4 = z_k(25:32,:) - mean(z_k(25:32,:),1);
    
        end
    
    end
    
    % pca on T_r<rank> (the first two eigenvectors identify the best
    % fitting plane
    [coeff, score, ~, ~, explained, ~] = pca(Y_r1);
    plane_r1 = coeff(:,1:2);
    score_r1 = score(:,1:2);
    explained_r1 = sum(explained(1:2));
    
    if stim_order > 1
    
        [coeff, score, ~, ~, explained, ~] = pca(Y_r2);
        plane_r2 = coeff(:,1:2);
        score_r2 = score(:,1:2);
        explained_r2 = sum(explained(1:2));
    
    end
    
    if stim_order > 2
    
        [coeff, score, ~, ~, explained, ~] = pca(Y_r3);
        plane_r3 = coeff(:,1:2);
        score_r3 = score(:,1:2);
        explained_r3 = sum(explained(1:2));
    end
    
    if stim_order > 3
    
        [coeff, score, ~, ~, explained, ~] = pca(Y_r4);
        plane_r4 = coeff(:,1:2);
        score_r4 = score(:,1:2);
        explained_r4 = sum(explained(1:2));
    
    end
    
    planes = [];
    planes.plane_r1 = plane_r1;
    planes.score_r1 = score_r1;
    planes.explained_r1 = explained_r1;
    planes.Y_r1 = Y_r1;
    
    if stim_order > 1
    
        planes.plane_r2 = plane_r2;
        planes.score_r2 = score_r2;
        planes.explained_r2 = explained_r2;
        planes.Y_r2 = Y_r2;
    
    end
    
    if stim_order > 2
    
        planes.plane_r3 = plane_r3;
        planes.score_r3 = score_r3;
        planes.explained_r3 = explained_r3;
        planes.Y_r3 = Y_r3;
    
    end
    
    if stim_order > 3
    
        planes.plane_r4 = plane_r4;
        planes.explained_r4 = explained_r4;
        planes.score_r4 = score_r4;
        planes.Y_r4 = Y_r4;
        
    end


    % set outputss
    if split_i == 1

        z_k_split1 = z_k;
        planes_split1 = planes;

        Y_r1_split1 = Y_r1;
        if stim_order > 1
            Y_r2_split1 = Y_r2;
        end
        if stim_order > 2
            Y_r3_split1 = Y_r3;
        end
        if stim_order > 3
            Y_r4_split1 = Y_r4;
        end

        score_r1_split1 = score_r1;
        if stim_order > 1
            score_r2_split1 = score_r2;
        end
        if stim_order > 2
            score_r3_split1 = score_r3;
        end
        if stim_order > 3
            score_r4_split1 = score_r4;
        end

    else

        z_k_split2 = z_k;
        planes_split2 = planes;

        Y_r1_split2 = Y_r1;
        if stim_order > 1
            Y_r2_split2 = Y_r2;
        end
        if stim_order > 2
            Y_r3_split2 = Y_r3;
        end
        if stim_order > 3
            Y_r4_split2 = Y_r4;
        end

        score_r1_split2 = score_r1;
        if stim_order > 1
            score_r2_split2 = score_r2;
        end
        if stim_order > 2
            score_r3_split2 = score_r3;
        end
        if stim_order > 3
            score_r4_split2 = score_r4;
        end

    end




end

warning on

% merge outputs
Y_r = [];
Y_r.Y_r1_split1 = Y_r1_split1;
Y_r.Y_r1_split2 = Y_r1_split2;

if stim_order > 1
    Y_r.Y_r2_split1 = Y_r2_split1;
    Y_r.Y_r2_split2 = Y_r2_split2;
end

if stim_order > 2
    Y_r.Y_r3_split1 = Y_r3_split1;
    Y_r.Y_r3_split2 = Y_r3_split2;
end

if stim_order > 3
    Y_r.Y_r4_split1 = Y_r4_split1;
    Y_r.Y_r4_split2 = Y_r4_split2;
end

score = [];
score.score_r1_split1 = score_r1_split1;
score.score_r1_split2 = score_r1_split2;

if stim_order > 1
    score.score_r2_split1 = score_r2_split1;
    score.score_r2_split2 = score_r2_split2;
end

if stim_order > 2
    score.score_r3_split1 = score_r3_split1;
    score.score_r3_split2 = score_r3_split2;
end

if stim_order > 3
    score.score_r4_split1 = score_r4_split1;
    score.score_r4_split2 = score_r4_split2;
end




end