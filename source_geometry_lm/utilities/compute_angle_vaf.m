function [angle, vaf] = compute_angle_vaf(planes, ranks)

% compute_angle_vaf computes the geometric similarity between subspaces
% of different ranks using two complementary metrics: the *angle* between
% subspaces and the *variance accounted for (VAF)* ratio.
%
% SYNTAX:
%   [angle, vaf] = compute_angle_vaf(planes, ranks)
%
% INPUTS:
%   planes : struct
%       A structure containing, for each rank r, the following fields:
%           - planes.plane_r#: basis vectors (columns) defining the subspace 
%             of rank # (e.g., plane_r1, plane_r2, …)
%           - planes.score_r#: component scores associated with the subspace
%           - planes.Y_r#: original data projected onto that subspace
%
%   ranks : integer
%       The number of subspaces (ranks) to compare. The function assumes
%       that fields in 'planes' are named accordingly (e.g., plane_r1, plane_r2, …).
%
% OUTPUTS:
%   angle : [ranks × ranks] matrix
%       Pairwise angles (in degrees) between subspaces of different ranks.
%       Angles are restricted to the range [0, 90].
%
%   vaf : [ranks × ranks] matrix
%       Pairwise variance accounted for (VAF) ratios, following Xie (2022).
%       Each element vaf(i,j) quantifies how much variance from subspace i
%       is explained by subspace j. Computed as the squared ratio of Frobenius norms:
%           vaf(i,j) = ( || Pj * Pi * Si' ||_F / || Pi * Si' ||_F )^2
%       where Pi and Pj are the bases (planes) and Si are the component scores.
%
% DESCRIPTION:
%   The function first computes the principal angles between subspaces by
%   taking the dot product of the normal vectors (cross products) of the
%   plane basis vectors. This yields an angular similarity measure that
%   reflects the orientation difference between the subspaces.
%
%   The VAF computation estimates how much of the variance captured by one
%   subspace can be represented within another subspace, providing a measure
%   of their representational overlap or shared structure.
%
% REFERENCES:
%   Santo-Angles A., Yang J., Zhou Y., Chu W.K.H., Lindsay G.W., Sreenivasan K.K. 
%   Neural Subspaces Encode Sequential Working Memory, but Neural Sequences Do Not. 
%   bioRxiv (2025). doi: https://doi.org/10.1101/2025.09.05.674385

%% angle between rank subspaces

angle = zeros(ranks,ranks);

for rank_i = 1:ranks

    eval(['plane1 = planes.plane_r' num2str(rank_i) ';']);

    for rank_ii = 1:ranks

        eval(['plane2 = planes.plane_r' num2str(rank_ii) ';']);

        cross_product1 = cross(plane1(:,1), plane1(:,2));
        cross_product2 = cross(plane2(:,1), plane2(:,2));

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

    eval(['plane1 = planes.plane_r' num2str(rank_i) ';']);
    eval(['score1 = planes.score_r' num2str(rank_i) ';']);
    eval(['data1 = planes.Y_r' num2str(rank_i) ';']);

    for rank_ii = 1:ranks

        eval(['plane2 = planes.plane_r' num2str(rank_ii) ';']);
        eval(['score2 = planes.score_r' num2str(rank_ii) ';']);

        vaf(rank_i, rank_ii) = (norm(plane2*plane2'*plane1*score1','fro')/norm(plane1*score1','fro')).^2;

    end

end

end

