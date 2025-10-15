function [angle_min, angle_max, vaf] = compute_angle_vaf_accumulated(planes, ranks)

% compute_angle_vaf_accumulated computes the geometric similarity between
% subspaces of different ranks using principal angles (PAs) and variance
% accounted for (VAF).
%
% SYNTAX:
%   [angle_min, angle_max, vaf] = compute_angle_vaf_accumulated(planes, ranks)
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
%   angle_min : [ranks × ranks] matrix
%       Minimum principal angles (in degrees) between each pair of subspaces.
%       These represent the smallest angular deviation between the closest
%       directions shared by the two subspaces.
%
%   angle_max : [ranks × ranks] matrix
%       Maximum principal angles (in degrees) between each pair of subspaces.
%       These represent the largest angular deviation between the most
%       divergent directions of the two subspaces.
%
%   vaf : [ranks × ranks] matrix
%       Pairwise variance accounted for (VAF) ratios, following Xie (2021).
%       Each element vaf(i,j) quantifies how much variance from subspace i
%       is explained by subspace j, computed as the squared ratio of
%       Frobenius norms:
%           vaf(i,j) = ( || Pj * Pj' * Pi * Si' ||_F / || Pi * Si' ||_F )^2
%       where Pi and Pj are the basis matrices (planes) and Si are the
%       component scores.
%
% DESCRIPTION:
%   This function computes principal angles between subspaces defined by
%   basis matrices of different ranks. Unlike simple normal vector
%   comparisons (compute_angle_vaf.m), principal angles quantify the full geometric relationship
%   between two multidimensional subspaces.
%
%   For each pair of subspaces, the function:
%       1. Orthogonalizes the basis matrices using QR decomposition.
%       2. Computes singular values of their cross-projection matrix (U1' * U2).
%       3. Converts singular values to principal angles via arccosine.
%       4. Extracts the smallest (angle_min) and largest (angle_max) principal angles.
%
%   The VAF computation estimates how much of the variance captured by one
%   subspace can be represented within another, providing a complementary
%   measure of representational overlap or shared structure.
%
%   For subspaces of dimensionality 3, this principal angle computation
%   yields identical results to the compute_angle_vaf.m function.
%   However, it also generalizes correctly to higher-dimensional subspaces.
%   The computation of VAF remains unchanged between this and compute_angle_vaf.m function.
%
% REFERENCES:
%   Santo-Angles A., Yang J., Zhou Y., Chu W.K.H., Lindsay G.W., Sreenivasan K.K. 
%   Neural Subspaces Encode Sequential Working Memory, but Neural Sequences Do Not. 
%   bioRxiv (2025). doi: https://doi.org/10.1101/2025.09.05.674385

%% angle between rank subspaces

angle_min = zeros(ranks,ranks);
angle_max = zeros(ranks,ranks);

for rank_i = 1:ranks

    eval(['plane1 = planes.plane_r' num2str(rank_i) ';']);

    for rank_ii = 1:ranks

        eval(['plane2 = planes.plane_r' num2str(rank_ii) ';']);

        [U1,~] = qr(plane1,0);
        [U2,~] = qr(plane2,0);

        C = U1' * U2;
        sigma = svd(C);

        angle_min(rank_i, rank_ii) = min(acosd(min(max(sigma, -1), 1)));
        angle_max(rank_i, rank_ii) = max(acosd(min(max(sigma, -1), 1)));

    end

end

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

