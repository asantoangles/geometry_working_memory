function [angle, vaf] = compute_angle_vaf_within(planes_split1, planes_split2, Y_r, score, ranks)

% compute_angle_vaf_within extends compute_angle_vaf.m
% to compute principal angles (PA) and variance accounted for (VAF) between
% subsets of subspaces, enabling the assessment of within-subspace alignment.
% see compute_angle_vaf.m for documentation.
%
% REFERENCES:
%   Santo-Angles A., Yang J., Zhou Y., Chu W.K.H., Lindsay G.W., Sreenivasan K.K. 
%   Neural Subspaces Encode Sequential Working Memory, but Neural Sequences Do Not. 
%   bioRxiv (2025). doi: https://doi.org/10.1101/2025.09.05.674385

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