function [angle, vaf] = compute_alignment_metrics_within(X_split1, planes_split1, planes_split2, Y_r, score, ranks)

% COMPUTE_ALIGNMENT_METRICS_WITHIN  Compare subspace alignment across data splits.
%
%   [angle, vaf] = compute_alignment_metrics_within( ...
%                       X_split1, planes_split1, planes_split2, Y_r, score, ranks)
%
%   This function quantifies how similar the low-dimensional neural
%   subspaces are when computed separately on two independent data splits
%   (e.g., split-half reliability analysis). It computes:
%
%       (1) PRINCIPAL ANGLE: the geometric angle between 2-D planes
%           estimated for each rank in split 1 vs. split 2, using the
%           normal vectors of those planes.
%
%       (2) VAF (Variance Accounted For): how much variance in the
%           split-1 neural population activity for a given rank is
%           explained by projecting it onto the subspace estimated from
%           split-2 (and vice versa).
%
%   INPUTS:
%       X_split1        Population activity matrix for split 1
%                       (baseline-subtracted).
%
%       planes_split1   Struct from compute_subspaces_within containing:
%                           • plane_r<k> – 2-D basis vectors for rank k (split 1)
%
%       planes_split2   Same as above, but for split 2.
%
%       Y_r             Struct containing demeaned rank-specific population 
%                       activity for both splits:
%                           • Y_r<k>_split1, Y_r<k>_split2
%
%       score           Struct containing PCA scores for both splits:
%                           • score_r<k>_split1, score_r<k>_split2
%
%       ranks           Number of ranks (3 or 4).
%
%   OUTPUTS:
%       angle           ranks × ranks matrix of pairwise angles (degrees)
%                       between split-1 and split-2 planes for each rank.
%                       Values are bounded to [0, 90] degrees.
%
%       vaf             ranks × ranks matrix of Variance Accounted For (VAF):
%
%                           vaf(i,j) = variance in rank i (split 1)
%                                      explained by subspace of rank j (split 2)
%
%                       Values are bounded to [0, 1].





%% angle between splits

angle = zeros(ranks,ranks);

for rank_i = 1:ranks

    eval(['plane_split1 = planes_split1.plane_r' num2str(rank_i) ';']);

    for rank_ii = 1:ranks

        eval(['plane_split2 = planes_split2.plane_r' num2str(rank_ii) ';']);

        cross_product1 = cross(plane_split1(:,1), plane_split1(:,2));
        cross_product2 = cross(plane_split2(:,1), plane_split2(:,2));

        % Calculate the dot product of the cross products
        dot_product = dot(cross_product1, cross_product2);
        
        % Calculate the magnitudes of the cross products
        magnitude1 = norm(cross_product1);
        magnitude2 = norm(cross_product2);
        
        % Calculate the cosine of the angle between the two planes
        angle(rank_i, rank_ii) = dot_product / (magnitude1 * magnitude2);

    end

end

% convert to degrees
angle = real(acosd(angle));

% reduce to range [0 90]
% if angle > 90, 
angle(angle > 90) = abs(180 - angle(angle > 90));


%% variance accounted for (VAF)

vaf = zeros(ranks,ranks);

for rank_i = 1:ranks

    eval(['plane_split1 = planes_split1.plane_r' num2str(rank_i) ';']);
    eval(['data_split1 = Y_r.Y_r' num2str(rank_i) '_split1;']);
    eval(['score_split1 = score.score_r' num2str(rank_i) '_split1;']);

    for rank_ii = 1:ranks

        eval(['plane_split2 = planes_split2.plane_r' num2str(rank_ii) ';']);
        eval(['data_split2 = Y_r.Y_r' num2str(rank_ii) '_split2;']);
        eval(['score_split2 = score.score_r' num2str(rank_ii) '_split2;']);

        vaf(rank_i, rank_ii) = (norm(plane_split2*plane_split2'*plane_split1*score_split1','fro')/norm(plane_split1*score_split1','fro')).^2;

    end

end

end