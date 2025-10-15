function [angle_min, angle_max, vaf] = compute_angle_vaf_accumulated_within(planes_split1, planes_split2, Y_r, score, ranks)

% compute_angle_vaf_accumulated_within extends compute_angle_vaf_accumulated.m
% to compute principal angles (PA) and variance accounted for (VAF) between
% subsets of subspaces, enabling the assessment of within-subspace alignment.
% see compute_angle_vaf_accumulated.m for documentation.
%
% REFERENCES:
%   Santo-Angles A., Yang J., Zhou Y., Chu W.K.H., Lindsay G.W., Sreenivasan K.K. 
%   Neural Subspaces Encode Sequential Working Memory, but Neural Sequences Do Not. 
%   bioRxiv (2025). doi: https://doi.org/10.1101/2025.09.05.674385

%% min angle between splits

angle_min = zeros(ranks,ranks);
angle_max = zeros(ranks,ranks);

for rank_i = 1:ranks

    eval(['plane_split1 = planes_split1.plane_r' num2str(rank_i) ';']);

    for rank_ii = 1:ranks

        eval(['plane_split2 = planes_split2.plane_r' num2str(rank_ii) ';']);

        [U1,~] = qr(plane_split1,0);
        [U2,~] = qr(plane_split2,0);

        C = U1' * U2;
        sigma = svd(C);

        angle_min(rank_i, rank_ii) = min(acosd(min(max(sigma, -1), 1)));
        angle_max(rank_i, rank_ii) = max(acosd(min(max(sigma, -1), 1)));

    end

end


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